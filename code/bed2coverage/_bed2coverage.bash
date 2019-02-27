#!/bin/bash

WIG2BIGWIG="wigToBigWig" # should be in $PATH


# PARSE OPTIONS
OPTIND=1         # Reset in case getopts has been used previously in the shell.
USAGE=
usage() {
  echo >&2 "usage: ${SCRIPTNAME} -apmcs SuRE-count-files"
  echo >&2 "OPTIONS:"
  echo >&2 "  -a: name of 'all' output bigwig file"
  echo >&2 "  -p: name of 'plus' output bigwig file"
  echo >&2 "  -m: name of 'minus' output bigwig file"
  echo >&2 "  -c: name of cDNA column to extract from input table"
  echo >&2 "  -s: filename with chomosome sizes"
  echo >&2 "  -h: print this message"
  echo >&2 ""
  exit 1;
}
while getopts "h?a:p:m:c:s:" opt; do
  case $opt in
    a)
      OUTALL=$OPTARG;
      ;;
    p)
      OUTPLUS=$OPTARG;
      ;;
    m)
      OUTMINUS=$OPTARG;
      ;;
    c)
      COLUMN=$OPTARG;
      ;;
    s)
      CHROMSIZES=$OPTARG;
      ;;
    h)
      usage;
      ;;
    \?)
      echo "option not recognized: "$opt
      usage
      ;;
  esac
done
shift $(( OPTIND - 1 ))

INPUT="$@"

tmpfileall=$(mktemp wig.all.XXXXXX)
tmpfileplus=$(mktemp wig.plus.XXXXXX)
tmpfileminus=$(mktemp wig.minus.XXXXXX)

# trick to enforce waiting for subprocesses created by tee
# from https://unix.stackexchange.com/questions/351780/wait-for-bash-subshells
# and https://unix.stackexchange.com/questions/29851/shell-script-mktemp-whats-the-best-method-to-create-temporary-named-pipe
tempdir=$(mktemp -d waitfifo.tempdir.XXXXXX)
trap 'rm -rf ${tempdir}' EXIT INT TERM HUP
waitfifo="${tempdir}/pipe"
mkfifo "${waitfifo}"

zcat $INPUT | awk -v col="${COLUMN}" ' 
                  BEGIN { OFS="\t"; FS="\t"; }
		  NR==1 { 
		    for (i=1; i<=NF; i++) {
		      headers[$i]=i
		    }
		    chr=headers["chrom"]
		    start=headers["start_hg19"]
		    end=headers["end_hg19"]
		    strand=headers["strand"]
		    # the column with iPCR counts is called "count"
		    if (col=="iPCR") { col="count" }
		    # the column names in the SuRE-counts files miss the "-T1"
		    # extension, so I need to remove that here as well.
		    col=headers[gensub("-T1","",1,col)]
		    next
		  }
		  NR==2 {
		    # rewrite the name of the chromosome to hg19
		    chr = gensub("(.*)_[pm]aternal","chr\\1","g",$headers["chrom"])
		  }
		  !/^BC/{
		    if( $start>0 && $end>0) {
		      print chr, $start, $end, $strand, $col}
		    next
		  } ' |\
  tee >( { 
         awk '
           BEGIN { OFS="\t"; FS="\t" }
           NR==1 { next }
	   { print $1, $2, $3, $5 ; }' |\
	 tee >( awk 'BEGIN{OFS="\t"; FS="\t"} $4>0{$4=1} 1' | ${BED2COVERAGE_EXE} > "flat.${tmpfileall}" ; : >${waitfifo} ) |\
         ${BED2COVERAGE_EXE} > $tmpfileall; : >${waitfifo}
       } ) \
      >( { awk '
           BEGIN { OFS="\t"; FS="\t" }
           NR==1 { next }
           { print $1, $2, $3, ($4=="-"?$5:0) };' |\
	 tee >( awk 'BEGIN{OFS="\t"; FS="\t"} $4>0{$4=1} 1' | ${BED2COVERAGE_EXE} > "flat.${tmpfileminus}" ; : >${waitfifo} ) |\
	 ${BED2COVERAGE_EXE} > $tmpfileminus; : >${waitfifo}
       } ) \
      >( { awk '
           BEGIN { OFS="\t"; FS="\t" }
           NR==1 { next }
           { print $1, $2, $3, ($4=="+"?$5:0) };' |\
	 tee >( awk 'BEGIN{OFS="\t"; FS="\t"} $4>0{$4=1} 1' | ${BED2COVERAGE_EXE} > "flat.${tmpfileplus}" ; : >${waitfifo} ) |\
	 ${BED2COVERAGE_EXE} > $tmpfileplus; : >${waitfifo}
       } ) > /dev/null

# wait for all subprocesses to finish
for (( i=0;i<6;i++ )); do read <${waitfifo}; done

echo "$tmpfileall ${CHROMSIZES} $OUTALL
$tmpfileminus ${CHROMSIZES} $OUTMINUS
$tmpfileplus ${CHROMSIZES} $OUTPLUS" |\
  parallel --colsep ' ' ${WIG2BIGWIG}

echo "flat.$tmpfileall ${CHROMSIZES} ${OUTALL/.cov./.covflat.}
flat.$tmpfileminus ${CHROMSIZES} ${OUTMINUS/.cov./.covflat.}
flat.$tmpfileplus ${CHROMSIZES} ${OUTPLUS/.cov./.covflat.}" |\
  parallel --colsep ' ' ${WIG2BIGWIG}

rm -f $tmpfileall $tmpfileminus $tmpfileplus
rm -f flat.$tmpfileall flat.$tmpfileminus flat.$tmpfileplus

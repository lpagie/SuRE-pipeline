# AUTHOR / DATE
#   Ludo Pagie; Nov 21, 2018; LP180525_SNPannot-altRefReads.py

# INTRO / BACKGROUND

#
# USAGE / INPUT / ARGUMENTS / OUTPUT
# USAGE:
#   required:
#     -bp: bam file with reads aligned on paternal genome, sorted on readID
#     -bm: bam file with reads aligned on maternal genome, sorted on readID
#     -basename: template used to create filenames for ouput bam files
#   optional:
#     -d: directory for output file
#     -l: name of logfile
#   fastq files are expected as last arguments on the commandline
# INPUT:
#   2 bam files, with same set of reads aligned to paternal genome ('bp') and
#     maternal genome ('mb')
# OUTPUT:
#   2 bam files, with reads classified as paternal in bam file
#     'basename_paternal.bam' and reads classified as maternal in
#     'basename_maternal.bam'. Reads classified differently (lowQ, equal,
#     ambiguous) are dicarded.
# USAGE:
#   python indel_LP181120.py --bp XYZ_paternal.bam --bm XYZ_maternal.bam --basename XYZ_classified

# VERSIONS:

import sys
import argparse
import os.path
import glob
import pysam

def parse_options():
# parse user options:

    parser = argparse.ArgumentParser(
        description="Annotate SNPs in reads, given alignments on paternal and maternal genome")
    parser.add_argument('--bp', required=True,
                        help=('bam file with aligned reads using paternal genome'))
    parser.add_argument('--bm', required=True,
                        help=('bam file with aligned reads using maternal genome'))
    parser.add_argument('--basename', required=True,
                        help=('basename for output bam files'))

    parser.add_argument('-d', '--dirout', required=False,default='.',
                        help=('output directory'))
    parser.add_argument('-l', '--log',required=False,default=None,
                        help=('log'))
    options = parser.parse_args()
#    if options.log is None:
#        options.log=options.vcfout+'.stats'
    return options


def import_bam(options):
    bams = {'paternal':pysam.AlignmentFile(options.bp, "rb"),
            'maternal':pysam.AlignmentFile(options.bm, "rb")}
    return(bams)

def export_bam(options, input_bams):
    bams = {'paternal':pysam.AlignmentFile(options.basename+'_paternal.bam', 'wb', template=input_bams['paternal']),
            'maternal':pysam.AlignmentFile(options.basename+'_maternal.bam', 'wb',input_bams['maternal'])}
    return(bams)

def classify_parent(rpf,rmf,rpr,rmr):
    # return:
        # - lowQ:  the best MAPq score below threshold; pair should be discarded
        # - equal: pair has equal MAPQ and aln score for forw and rev reads; discard pair
        # - ambiguous: forw and rev reads suggest different parent; pair should be discarded
        # - paternal/maternal: read is classified and should be annotated with INDELs and SNPs
        # - discord: both read pairs are dicordant; discard pair
    min_mapq = 40
    # concordant pair or not?
    if (not (rpf.is_proper_pair | rmf.is_proper_pair)):
        # neither read pair is concordant, so discard read
        return ('discord',{'forw':rpf,'rev':rpr})
    # if MAPQ < (?) 40 return lowQ
    mapq = [r.mapping_quality for r in (rpf,rmf,rpr,rmr)]
    if max(min(mapq[0],mapq[2]),min(mapq[1],mapq[3])) < min_mapq:
        return ('lowQ',{'forw':rpf,'rev':rpr})

    # if both pairs are not concordant, return 'noCP'
    if (not((rpf.get_tag('YT')=='CP') or (rmf.get_tag('YT')=='CP'))):
        return ('noCP',{'forw':rpf,'rev':rpr})

    # if exactly 1 of the readpairs is a concordant pair this is the parent
    if ((rpf.get_tag('YT')=='CP') ^ (rmf.get_tag('YT')=='CP')):
        if (rpf.get_tag('YT')=='CP'):
            return('paternal',{'forw':rpf,'rev':rpr})
        else:
            return('maternal',{'forw':rmf,'rev':rmr})

    # MAPQ: select pair with higher MAPQ, or continue
    if (mapq[0] != mapq[1]):
        if ((mapq[0]-mapq[1]) > 0):
            return('paternal',{'forw':rpf,'rev':rpr})
        else:
            return('maternal',{'forw':rmf,'rev':rmr})

    # AS: select pair if both aln scores are equal or higher, or continue
    try:
        AS = [r.get_tag('AS') for r in (rpf,rmf,rpr,rmr)]
    except:
        print(rpf)
        print(rmf)
        print(rpr)
        print(rmr)

    # if ASdiff are equal values -> ambiguous
    if (AS[0]==AS[1]) & (AS[2]==AS[3]):
        return('equal',{'forw':rpf,'rev':rpr})
    ASdiff = [AS[0]-AS[1], AS[2]-AS[3]]
    # if ASdiff has opposite sign -> ambiguous
    if (ASdiff[0]*ASdiff[1]) < 0:
        return('ambiguous',{'forw':rpf,'rev':rpr})
    # sum of signs -> paternal or maternal
    if (ASdiff[0]+ASdiff[1])>0:
        return('paternal',{'forw':rpf,'rev':rpr})
    else:
        return('maternal',{'forw':rmf,'rev':rmr})

    # at this point there is no way to classify the read parent
    # in fact we should never reach this point :-)
    print('we are not here!!!')


def classify_reads(ibam, obam):
    reads = zip(ibam['paternal'].fetch(until_eof=True),
                ibam['maternal'].fetch(until_eof=True))
    for rpf, rmf in reads:
        rpr, rmr = next(reads)
        read_class,read = classify_parent(rpf,rmf,rpr,rmr)
        if (read_class=='paternal' or read_class=='maternal'):
            obam[read_class].write(read['forw'])
            obam[read_class].write(read['rev'])

def main(options):
    # open input bam files
    input_bams = import_bam(options)
    # open output bam files
    output_bams = export_bam(options, input_bams)
    # loop over input reads and write to output
    classify_reads(input_bams, output_bams)



if __name__ == "__main__":
    sys.stderr.write("command: %s\n" % " ".join(sys.argv))
    options = parse_options()
    main(options)


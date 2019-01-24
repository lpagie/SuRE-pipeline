# AUTHOR / DATE
#   Ludo Pagie; November 30, 2018; bam_to_annotBed.py

# INTRO / BACKGROUND

# USAGE / INPUT / ARGUMENTS / OUTPUT
# USAGE:
#   required:
#     -bam: bam file, sorted on readID
#     -vcf: vcf file, corresponding to bam file (cell line + paternal/maternal)
#   optional:
#     -d: directory for output file
#     -l: name of logfile
# INPUT:
#   
# OUTPUT:
#   

# VERSIONS:

import sys
import argparse
import os.path
import glob
import gzip
import re
import pysam
import vcf
from bx.intervals.intersection import Interval, IntervalTree

CHRS = ['1','2','3','4','5','6','7','8','9',
       '10','11','12','13','14','15','16','17','18','19',
       '20','21','22','X']

def parse_options():
# parse user options:
    parser = argparse.ArgumentParser(
        description="Annotate SNPs in reads, given alignments on paternal and maternal genome")
    parser.add_argument('-b', '--bam', required=True,
                        help=('bam file'))
    parser.add_argument('-v', '--vcf', required=True,
                        help=('vcf file for corresponding (parental) genome'))
    parser.add_argument('-o', '--out', required=True,
                        help=('output bedpe-like file'))
    parser.add_argument('-p', '--patmat',required=True,choices=['paternal', 'maternal'],
                        help=('paternal (1st) or maternal (2nd) genotype'))
    parser.add_argument('-l', '--log',required=False,default=None,
                        help=('genomeID'))
    options = parser.parse_args()
    return options


def import_VCF(vcf_fname):
    # using try/except mostly as an exercise
    try:
        vcf_reader = vcf.Reader(open(vcf_fname, 'r'))
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        sys.exit()
    except IOError as io_error:
        print(io_error)
        sys.exit()
    tree = IntervalTree()
    _c='' 
    for record in vcf_reader:
        # check we are reading first entry and therefor set _c to chrom
        # or we read the same chrom as before.
        # otherwise abort
        if _c is '':
            _c = record.CHROM
        assert _c == record.CHROM
        # insert record in tree if it is an INDEL
        # if record.is_indel:
        tree.insert_interval(Interval(record.start, record.end, record))
    return(tree)


def open_bam(bam_fname):
    bam = pysam.AlignmentFile(bam_fname, "rb")
    # check whether bam is name sorted
    head = bam.head(100)
    c = 0
    for r1 in head:
        r2 = next (head)
        c += 1
        try:
            assert r1.query_name == r2.query_name
        except AssertionError as error:
            print(r1)
            print(r2)
            print(error)
            print("error opening bam file, read %d pairs" % c)
            sys.exit()
    return bam

def open_output(options):
    out = gzip.open(options.out, 'wb')
    header = "readID chrom start end length strand iend istart MAPQ1 MAPQ2 MD1 MD2 XS1 XS2 SEQ1 SEQ2 CIGAR1 CIGAR2 SNP_ABS_POS SNP_REL_POS SNP_ID SNP_SEQ SNP_VAR SNP_PARENT SNP_TYPE SNP_SUBTYPE"
    # replace spaces in header string with tabs
    header = header.replace(" ", "\t")
    # add a newline to the headerline
    header = header + "\n"
    out.write(header.encode('utf-8'))
    return(out)

def iupac(base1, base2):
    if (base1 == base2): return base1
    if ((base1=="A" and base2=="G") or (base2=="A" and base1=="G")): return "R"
    if ((base1=="A" and base2=="T") or (base2=="A" and base1=="T")): return "W"
    if ((base1=="A" and base2=="C") or (base2=="A" and base1=="C")): return "M"
    if ((base1=="G" and base2=="C") or (base2=="G" and base1=="C")): return "S"
    if ((base1=="G" and base2=="T") or (base2=="G" and base1=="T")): return "K"
    if ((base1=="C" and base2=="T") or (base2=="C" and base1=="T")): return "Y"

def _map(func, iterable):
    return [func(x) if x is not None else '.' for x in iterable]

def _stringify(iterable, sep=','):
    return sep.join(_map(str, iterable))

def annotate_snp(snp, r1, r2, patmat):
    def annotate_snp_in_read(r):
        # overlap with read 'r'
        try:
        # check if genomic position is covered by read and not deleted in read
            rel_pos = r.get_reference_positions().index(snp.start)
        except ValueError:
            snp_base = None
            snp_var = -2
            snp_patmat = 'pos deleted in read'
            return (snp_base, snp_var, snp_patmat)

        snp_base = r.query_sequence[rel_pos]
        try:
            # check if observed base is either reference or alternative (1, 2, etc)
            snp_var = snp.alleles.index(snp_base) # 0: reference, 1..: 1st (2nd, 3rd, etc) allele
            if snp_var == int(snp.samples[0].data.GT.split('|')[patmat=='maternal']):
                # base is what is expected according to parent 'patmat'
                snp_patmat = patmat
            elif snp_var == int(snp.samples[0].data.GT.split('|')[patmat!='maternal']):
                # base is on parental chromosome other than expected
                snp_patmat = 'maternal_unexpected' if patmat=='paternal' else 'paternal_unexpected'
            else:
                # base is known allele but neither of the two parental alleles
                snp_patmat = "non_paternal_allele"
        except ValueError:
            snp_var = -1 # base is not a known allele
            snp_patmat = "unknown_allele"
        return (snp_base, snp_var, snp_patmat)

    snp_ID = snp.ID
    snp_abs_pos = snp.start
    snp_rel_pos = snp.start - r1.reference_start # both coord systems are 0-based
    snp_type    = snp.var_type
    snp_subtype = snp.var_subtype
    annot = []

    if snp.POS>r1.reference_end and snp.POS<=r2.reference_start:
        # snp in between r1 and r2
        try: 
            snp_base = iupac(snp.alleles[int(snp.samples[0].data.GT.split('|')[patmat=='paternal'])], 
                             snp.alleles[int(snp.samples[0].data.GT.split('|')[patmat=='maternal'])])
        except NameError:
            print("error in annotate_snp_in_read with SNP %s" % snp)
            print(snp.alleles)
            sys.exit()

        snp_var = -2 # base is not covered by either of the two reads
        snp_patmat = "unread"
        annot = [(snp_base, snp_var, snp_patmat)]
    else:
        # snp overlaps with either one or both reads
        if snp.POS<r1.reference_end:
            annot.append(annotate_snp_in_read(r1))
        if snp.POS > r2.reference_start:
            annot.append(annotate_snp_in_read(r2))

    return (snp_rel_pos, snp_ID, annot, snp_type, snp_subtype, snp_abs_pos)


def annotate_indel(snp, r1, r2, patmat):
    def annotate_indel_in_read(r):
        # overlap with read 'r'
        snp_var = int(snp.samples[0].data.GT.split('|')[patmat=='maternal'])
        expect_seq = snp.alleles[snp_var]
        # if alignment contains INDELs the read- and genome-positions are not synchronous.
        # In that case the CIGAR string needs to be used to infer corresponding positions
        # if I:D:N:S:P in CIGAR ....
        if re.match(".*[IDNSP].*", r.cigarstring):
            poss =  range(snp.start, snp.start+len(expect_seq))
            aln_pairs = r.get_aligned_pairs(with_seq=True)
            bases = [x[2] for x in aln_pairs if x[1] in poss]
            obs_seq = ''.join(bases)
        else:
            rel_pos = snp.start - r.reference_start # both coord systems are 0-based
            obs_seq = r.query_sequence[rel_pos:min(rel_pos+len(expect_seq), len(r.query_sequence))]
        if obs_seq == expect_seq:
            # snp_var doesn't change
            snp_base = expect_seq
            snp_patmat = patmat
        else:
            snp_var = -1 # observed sequence differs from expected sequence, unclear whether observed is another allele or sequencing error or something else
            snp_patmat = "unexpected"
            snp_base = obs_seq

        return (snp_base, snp_var, snp_patmat)

    snp_ID      = snp.ID
    snp_abs_pos = snp.start
    snp_rel_pos = snp.start - r1.reference_start # both coord systems are 0-based
    snp_max_end = snp.start + max(len(allele) for allele in snp.alleles) 
    snp_type    = snp.var_type
    snp_subtype = snp.var_subtype
    annot = []

    if snp.start < r1.reference_start or snp_max_end > r2.reference_end:
        # snp overlaps fragment boundaries; discard snp completely
        return None
    elif snp_max_end > r1.reference_end and snp.start < r2.reference_start:
        # snp does not overlap completely with either read; check if homolgous alleles
        if re.match(r'^(.*)\|\1$',snp.samples[0].data.GT): # both parents have same allele
            snp_base = snp.alleles[int(snp.samples[0].data.GT.split('|')[patmat=='maternal'])]
        else:
            snp_base = ""
        snp_var = -2 # base is not fully covered by either of the two reads
        snp_patmat = "unread"
        annot = [(snp_base, snp_var, snp_patmat)]
        return (snp_rel_pos, snp_ID, annot, snp_type, snp_subtype, snp_abs_pos)
    else:
        if snp_max_end < r1.reference_end:
            # snp overlaps completely with read1
            annot.append(annotate_indel_in_read(r1))
        if snp.start > r2.reference_start:
            # snp overlaps completely with read2
            annot.append(annotate_indel_in_read(r2))
    return (snp_rel_pos, snp_ID, annot, snp_type, snp_subtype, snp_abs_pos)

def annotate_fragment(r1, r2, vcf, patmat):
    # find overlapping SNPs
    # sam coordinates are 0-based
    snps = vcf.find(r1.reference_start, r2.reference_end-1) # reference_end gives "one past the last aligned residue" but is also 0-based (pysam manual). Nevertheless I found I need the 'end-1' anyway, not to include SNPs beyond the end of the alignment
    snp_annot = []
    # check SNP base/sequence identity
    if len(snps) is 0:
        return
    for s in snps:
        if s.value.is_snp:
            snp_annot.append(annotate_snp(s.value, r1, r2, patmat))
        elif s.value.is_indel:
            annot = annotate_indel(s.value,r1,r2, patmat)
            # indel may overlap boundaries of reads/fragment in which case annot may be None; discard in that case
            if annot is not None:
                snp_annot.append(annot)
        else:
            print (s.value)
    return(snp_annot)

def _stringify_fragment(r1, r2, snp_annot):
    # write read and snp info into bedpe-like format:
    # - chr, start, end (1-based, fully closed)
    # - length
    # - strand
    # - barcode sequence
    # - read count
    # - internal start, end
    # - MAPQ1 (now also MAPQ2)
    # - MD1, MD2
    # - alternative alignment exists 1,2
    # - read1 seq, read2 seq
    # - CIGAR1, 2
    # - rel_snp_pos (0-based), snp_base, abs_snp_pos (1-based), snp_var, snp_ind_in_vcf
    # - inf_base, inf_snp_var, SNP_ID, patmat

    # BC = ''
    alt1 = str(r1.get_tag('XS') if 'XS' in [e[0] for e in r1.get_tags()] else '')
    alt2 = str(r2.get_tag('XS') if 'XS' in [e[0] for e in r2.get_tags()] else '')
    md1 = str(r1.get_tag('MD') if 'MD' in [e[0] for e in r1.get_tags()] else '')
    md2 = str(r2.get_tag('MD') if 'MD' in [e[0] for e in r2.get_tags()] else '')
    cigar1 = r1.cigarstring or ''
    cigar2 = r2.cigarstring or ''
    if None is not snp_annot:
        snp_abs_pos = _stringify([pos  for annot in snp_annot for pos in [annot[5]]*len(annot[2])])
        snp_rel_pos = _stringify([pos  for annot in snp_annot for pos in [annot[0]]*len(annot[2])])
        snp_ID      = _stringify([ID   for annot in snp_annot for ID  in [annot[1]]*len(annot[2])])
        snp_base    = _stringify([a[0] for annot in snp_annot for a in annot[2]])
        snp_var     = _stringify([a[1] for annot in snp_annot for a in annot[2]])
        snp_patmat  = _stringify([a[2] for annot in snp_annot for a in annot[2]])
        snp_type    = _stringify([stype for annot in snp_annot for stype in [annot[3]]*len(annot[2])])
        snp_subtype = _stringify([stype for annot in snp_annot for stype in [annot[4]]*len(annot[2])])
    else:
        snp_abs_pos    = ""
        snp_rel_pos    = ""
        snp_ID         = ""
        snp_base       = ""
        snp_var        = ""
        snp_patmat     = ""
        snp_type       = ""
        snp_subtype    = ""

    try:
        fragment = '\t'.join([r1.query_name, r1.reference_name, str(r1.reference_start+1), str(r2.reference_end), # reference_start is 0-based
                          str(r2.reference_end-r1.reference_start), ('-' if r1.is_read2 else '+'), 
                          str(r1.reference_end), str(r2.reference_start), 
                          str(r1.mapping_quality), str(r2.mapping_quality), 
                          md1, md2, alt1, alt2, r1.query_sequence, r2.query_sequence, cigar1, cigar2, 
                          snp_abs_pos, snp_rel_pos, snp_ID, snp_base, snp_var, snp_patmat, snp_type, snp_subtype])
    except TypeError:
        print(md1)
        print(md2)
        print(type(md1))
        for e in [r1.query_name, r1.reference_name, str(r1.reference_start+1), str(r2.reference_end), # reference_start is 0-based
                                             str(r2.reference_end-r1.reference_start), ('-' if r1.is_reverse else '+'), 
                                             str(r1.reference_end), str(r2.reference_start), 
                                             str(r1.mapping_quality), str(r2.mapping_quality), 
                                             md1, md2, alt1, alt2, cigar1, cigar2, 
                                             snp_abs_pos, snp_rel_pos, snp_ID, snp_base, snp_var, snp_patmat, snp_type, snp_subtype]:
            print(type(e))
        print("failing fragment with readID %s" % r1.query_name)
        return

    fragment += "\n"
    return fragment

def write_fragment(r1, r2, snp_annot, out):
    fragment = _stringify_fragment(r1, r2, snp_annot)
    out.write(fragment.encode('utf-8'))
    return

def main(options):
    # some checks on cell line, chromosome, paternal/maternal, bam is sorted on readname
    # get interface to reads in bamfiles
    print("opening bam file %s" % options.bam)
    reads = open_bam(options.bam)
    print("bam file opened")
    # import SNPs files, store in interval tree
    print("importing vcf file %s" % options.vcf)
    vcf = import_VCF(options.vcf) # import_SNPs returns an IntervalTree
    print("import of vcf done")
    # open output file
    print("opening bedpe file for output: %s" % options.out)
    out = open_output(options)
    print("output opened")
    # iterate over reads, annotate with SNPs, output to (compressed) tabular txt file
    print("start processing reads")
    for r1 in reads:
        r2 = next(reads)
        # checks?
        # assert r1.is_read1
        if r1.is_reverse:
            r_tmp = r2
            r2 = r1
            r1 = r_tmp
        
        # check strand 
        snp_annot = annotate_fragment(r1, r2, vcf, options.patmat)
        if snp_annot is not None:
            write_fragment(r1, r2, snp_annot, out)
    print("processing done")
    out.close()
    print("output closed")
    print("done")
    return


# run -n LP180726_SNPannot-altRefReads.py

if __name__ == "__main__":
    sys.stderr.write("command: %s\n" % " ".join(sys.argv))
    options = parse_options()
    main(options)



# 163, 2nd, plus
# 83,  1st, minus
# 
# 99, 1st, plus
# 147, 2nd, minus
# 

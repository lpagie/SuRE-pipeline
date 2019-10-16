# AUTHOR / DATE
#   Ludo Pagie; January 15, 2019; iPCR_addBC2bedpe.py

# INTRO / BACKGROUND
#   Take a bedpe-like and an info file as input and associate barcode sequences
#   from the info file with fragments in the bedpe-like file, through common
#   readIDs in both

# USAGE / INPUT / ARGUMENTS / OUTPUT
# USAGE:
#   required:
#     -bedpe: bedpe file
#     -info: info file, corresponding to bedpe file (cell line + paternal/maternal)
#     -out: bedpe filename for output
#   optional:
#     -d: directory for output file
#     -l: name of logfile
# INPUT:
#   -bedpe: tabular text file, containing info and annotation on SuRE fragments, among which readIDs
#   -info: tabular text file, containing info on trimming of fastq files, among which readID and barcode sequence
# OUTPUT:
#   -out: tabular text file, compressed; contains merge of both input files 

# VERSIONS:

import sys
import argparse
import os.path
import glob
import gzip
import pandas

import re
import pysam

def parse_options():
# parse user options:
    parser = argparse.ArgumentParser(
        description="Annotate SNPs in reads, given alignments on paternal and maternal genome")
    parser.add_argument('-b', '--bedpe', required=True, nargs='+',
                        help=('bedpe file'))
    parser.add_argument('-i', '--info', required=True,
                        help=('info file for corresponding sample'))
    parser.add_argument('-o', '--out', required=True, nargs='+',
                        help=('output bedpe-like file'))
    parser.add_argument('-l', '--log', required=True,
                        help=('log file'))
    parser.add_argument('-s', '--stats', required=True,
                        help=('stats file'))

    options = parser.parse_args()
    return options

def open_log(log_fname):
    return open(log_fname, "a+")

def open_stats(stats_fname, log):
    log.write("iPCR_addBC2bedpe.py: opening statsfile (%s) for writing\n" % stats_fname )
    return open(stats_fname, "a+")

def import_bedpe(bedpe_fname, log, stats):
    log.write("iPCR_addBC2bedpe.py: opening bedpe %s for input\n" % bedpe_fname)
    bedpe = pandas.read_table(bedpe_fname,header=0)
    bedpe['readID'] = bedpe['readID'].str.extract(r'^(\S+?)(?:/[12])?(?:\s.*)?$', expand=False) # match (1) non-space string (non-greedy) (2) possible "/1,/1" (3) possible space plus rest of string
    return bedpe
    # return pandas.read_table(bedpe_fname,header=0) # LP190326; Need to trim the readID before returning

def import_barcodes(info_fname, log, stats):
    log.write("iPCR_addBC2bedpe.py: opening info %s for input\n" % info_fname)
    # open info file, store in pandas data frame
    # info file contains '-1' in 2nd column if read was not parsed
    # pandas.read_table doesn't handle this correctly if it occurs on the first line
    # first find first line in info file which does not contain a '-1' in 2nd column
    nskip=0
    with gzip.open(info_fname,'rt') as f:
        for line in f:
            if line.split("\t")[1]=="-1":
                nskip=nskip+1
            else:
                break
    log.write("iPCR_addBC2bedpe.py: skipping from info input file = %s\n" % str(nskip))
    barcodes = pandas.read_table(info_fname, sep="\t", header=None, usecols=[0,4], names=['readID','BC'], skiprows=nskip)
    # trim readID to make format compatible with bedpe file
    # barcodes['readID'] = barcodes['readID'].str.extract(r'^(.*) .*$', expand=False) # LP190320; this match fails if there is no white spec in readID at all
    # barcodes['readID'] = barcodes['readID'].str.extract(r'^(\S+).*$', expand=False) # LP190320; found that I need to be able to trim "/1" and "/2" from readID
    barcodes['readID'] = barcodes['readID'].str.extract(r'^(\S+?)(?:/[12])?(?:\s.*)?$', expand=False) # match (1) non-space string (non-greedy) (2) possible "/1,/1" (3) possible space plus rest of string
    log.write("iPCR_addBC2bedpe.py: from info file after rewriting 'readID':\n %s\n" % barcodes.head())
    # extract a table of barcode lengths from pandas column
    BClen_tbl = [len_count for len_count in barcodes['BC'].str.len().value_counts().sort_index().items()]
    # extract a table of number of N's in barcode from pandas column
    BCN_tbl = [NNN_count for NNN_count in barcodes['BC'].str.count('N').value_counts().sort_index().items()]
    tot = len(barcodes.index) # total number of barcodes items
    idx = (~barcodes['BC'].isnull()) & (barcodes['BC'].str.len() == 20) & ~barcodes['BC'].str.contains('N', na=False)
    log.write("iPCR_addBC2bedpe.py: length index: %s\n" % str(len(idx)))
    log.write("iPCR_addBC2bedpe.py: index: %s\n" % str(idx.value_counts()))
    incl = idx.value_counts()[True]
    excl = idx.value_counts()[False]

    # write stats to stats file
    stats.write("while filtering for proper barcodes: included = %d, discarded = %d\n\n" % (incl, excl))
    stats.write("bedpeFragmentCount\t%d\n" % incl)
    stats.write("iPCR_BC_lengths\t")
    stats.write("\t".join([str(i) for i,v in BClen_tbl])+"\n")
    stats.write("iPCR_BC_lengths_counts\t")
    stats.write("\t".join([str(v) for i,v in BClen_tbl])+"\n")
    stats.write("iPCR_BC_NNN\t")
    stats.write("\t".join([str(i) for i,v in BCN_tbl])+"\n")
    stats.write("iPCR_BC_NNN_counts\t")
    stats.write("\t".join([str(v) for i,v in BCN_tbl])+"\n")
    stats.write("\n\n") 
    stats.flush()

    return barcodes[idx]


def addBC2bedpe(bedpe, barcodes, log, stats):
    bedpe_merged = pandas.merge(barcodes, bedpe)
    stats.write("%d readIDs (of %d, %.2f%%) found in info file\n" % (len(bedpe_merged.index), 
                                                                     len(bedpe.index), 
                                                                     len(bedpe_merged.index)/len(bedpe.index)*100.0 if len(bedpe.index)>0 else float('NaN')))
    return (bedpe_merged)

def write_merged_bedpe(bedpe, bedpe_fname, log, stats):
    bedpe.to_csv(bedpe_fname, index=None, sep="\t", compression='gzip', header=True)
    return

def main(options):
    # open log file
    log = open_log(options.log)
    # open stats file
    stats = open_stats(options.stats, log)
    # import infofile with barcodes
    barcodes = import_barcodes(options.info, log, stats)
    # iterate over input bedpe files and merge each with the barcodes
    for bedpe_in, bedpe_out in zip(options.bedpe, options.out):
        # import bedpe file; store in pandas dataframe
        bedpe = import_bedpe(bedpe_in, log, stats)
        # merge bedpe and barcodes
        bedpe = addBC2bedpe(bedpe, barcodes, log, stats)
        # write merged bedpe to output
        write_merged_bedpe(bedpe, bedpe_out, log, stats)

    log.write("iPCR_addBC2bedpe.py: processing done\n")
    stats.close()
    log.close()
    return


# run -n LP180726_SNPannot-altRefReads.py

if __name__ == "__main__":
    sys.stderr.write("command: %s\n" % " ".join(sys.argv))
    options = parse_options()
    main(options)


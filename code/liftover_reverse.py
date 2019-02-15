# AUTHOR / DATE
#   Ludo Pagie; April 25, 2018; LP180416_liftover_VCFs.py
#   Ludo Pagie; February 13, 2019; liftover_reverse.py

# INTRO / BACKGROUND
#   A script to add hg19 positions to a tabular text file, given genome
#   positions in one of the 4(*2) genomes of the SuRE42-45 project.
#   make VCF/SNP files specific for each of the individual genomes in the SuRE
#   project.
#   The python script uses module pyliftover to apply LiftOver to regions.

#
# USAGE / INPUT / ARGUMENTS / OUTPUT
# USAGE:
#   required:
#     -i: tabular text file for input
#     -o: tabular text file for output
#     -c: chainfile
#     -p: specify whether genome is 1st (paternal) or 2nd (maternal)
#   optional:
#     -l: name of logfile
# INPUT:
#   tabular text file, with at least columns 'start', 'end'
# OUTPUT:
#   same tabular text file with added columns start_hg19, end_hg19

# VERSIONS:

import sys
import argparse
from pyliftover import LiftOver
import os
import re
import gzip

def parse_options():
    # parse user options:

    parser = argparse.ArgumentParser(
        description="Liftover a VCF, given a chainfile")
    parser.add_argument('-i', '--inn', required=True,
                        help=('tabular text to be liftovered'))
    parser.add_argument('-o', '--out', required=True,
                        help=('output abular text file'))
    parser.add_argument('-c', '--chainfile',required=True,
                        help=('chainfile'))
    parser.add_argument('-p', '--parent',required=True,choices=['paternal', 'maternal'],
                        help=('paternal (1st) or maternal (2nd) genotype'))
#    parser.add_argument('-s', '--sample',required=True,default=None,
#                        help=('sample'))
    parser.add_argument('-g', '--genomeTag', required=True,
                        help=('genomeID tag, added to new position columns'))
    parser.add_argument('-l', '--log',required=False,default=None,
                        help=('log file'))
    options = parser.parse_args()
    return options


def main(options):
    # open chainfile for input
    lo      = LiftOver(options.chainfile)
    # open files for input and for output
    log_out = open(options.log, "wt")
    outfile = gzip.open(options.out, "wt")
    infile  = gzip.open(options.inn, "rt")
    # check header contains expected columns
    header  = infile.readline()
    columns = header.split()
    new_columns = ["%s_%s" % (p, "hg19") for p in ["start","end","SNP_ABS_POS"]]
    outfile.write("\t".join(columns+new_columns)+"\n")
    expected_cols = ['chrom','start','end','strand','SNP_ABS_POS']
    assert all(e in columns for e in expected_cols)
    # get indices for columns of interest
    col2indices = {k:columns.index(k) for k in expected_cols}
    # iterate over input data, convert each position in a single line, check
    # that conversion makes sense,write output
    log_out.write("finished initialization")
    for line in infile:
        records = line.rstrip().split("\t")
        pos_lo = {p:None for p in ['start','end','SNP_ABS_POS']}
        for p in ['start','end']:
            pos =records[col2indices[p]]
            new_coord = lo.convert_coordinate(records[col2indices['chrom']],int(pos)-1, records[col2indices['strand']])
            #print(new_coord)
            if new_coord is None or len(new_coord) != 1 or new_coord[0][2] != records[col2indices['strand']] or new_coord[0][0] != re.sub('_[pm]aternal','',records[col2indices['chrom']]):
                pos_lo[p] = None
                log_out.write("no, or multiple positions after liftover (%d)\n" % (new_coord is not None and len(new_coord) or 0))
                continue
            pos_lo[p] = str(new_coord[0][1])
        snp_pos = records[col2indices['SNP_ABS_POS']]
        pos_new = []
        if snp_pos == '':
            pos_lo.update(SNP_ABS_POS=None)
        else:
            for pos in snp_pos.split(','):
                new_coord = lo.convert_coordinate(records[col2indices['chrom']],int(pos)-1, records[col2indices['strand']])
                if new_coord is None or len(new_coord) != 1 or new_coord[0][2] != records[col2indices['strand']] or new_coord[0][0] != re.sub('_[pm]aternal','',records[col2indices['chrom']]):
                    pos_lo.update(SNP_ABS_POS=None)
                    log_out.write("no, or multiple positions after liftover (%d)\n" % (new_coord is not None and len(new_coord) or 0))
                    log_out.write(line)
                    pos_new.append('')
                else:
                    pos_new.append(str(new_coord[0][1]))
            pos_lo['SNP_ABS_POS'] = ','.join(pos_new)
        records += pos_lo.values()
        outfile.write("\t".join(['' if e is None else str(e) for e in records])+"\n")
        # outfile.flush()

    # close output files
    infile.close()
    outfile.close()
    log_out.close()


if __name__ == "__main__":
    sys.stderr.write("command: %s\n" % " ".join(sys.argv))
    options = parse_options()
    main(options)


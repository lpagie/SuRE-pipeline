#!/bin/bash

# AUTHOR / DATE
#   Ludo Pagie; Feb 19, 2019; sort_cnt_tbls.sh

# INTRO / BACKGROUND
#   A short bash script to sort SuRE count tables on genomic coordinate after
#   reverse liftover. The script (in all its simplicity) is necessary to deal
#   with nested quotes and escapes. It first replaces the tab delimiters by a
#   'normal character (@)', then sorts with the delimiter specification, and
#   the replaces the '@' delimiter back to tab.
#
# USAGE / INPUT / ARGUMENTS / OUTPUT
# USAGE:
#   required:
#     The script reads from STDIN; a tabular text file with start/end
#       coordinates in columns 19/20. All entries are assumed to be on one
#       chromosome.
#     The script writes output (gzipped) to STDOUT


cat | sed 's/\t/@/g' | awk '
NR>1 {
  print | "sort --parallel 5  -k19,19n -k20,20n -t@"; 
  next
}
{
  print
}' | sed 's/@/\t/g' | gzip -c 


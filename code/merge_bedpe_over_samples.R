# AUTHOR / DATE
#   Ludo Pagie; January 29, 2019; merge_bedpe_over_samples.R

# INTRO / BACKGROUND
#   Take a bedpe-like from all samples as input and merge the data. Next group
#   the data on barcode, and subsequently on the set of SNP_ID's; replace this
#   set by the most frequent element.
#   The exported data consists of unique barcodes (per chromosome) and per
#   barcode the most prominent annotation regarding position and SNPs.

# USAGE / INPUT / ARGUMENTS / OUTPUT
#   This script is exclusively for use by snakemake. The script gets its input
#   from the S4 snakemake object 'snakemake'. The required snakemake elements
#   are:
#   - input:bedpe, a vector with filenames of the bedpe files to be imported
#   - output: log, filename of the logfile
#             bedpe, filename for exported bedpe data
# USAGE:
#   required:
#   optional:
# INPUT:
#   via S4 object snakemake, see above
# OUTPUT:
#   via S4 object snakemake, see above

# VERSIONS:


library(argparse)
library(data.table)

setDTthreads(threads = 4)

# argument parsing
args <- commandArgs(trailingOnly=TRUE)
output.bedpe <- args[1]
output.log   <- args[2]
input        <- tail(args, -2)

# setup logging function
logout <- NULL
log <- function(msg="") {
  if (is.null(logout)) {
    logout <<- file(description = output.log, open = "wt")
    .Last <- function() {
      close(logout)
    }
  }
  write(x = msg, file = logout)
  flush(logout)
}


bedpe_fnames <- input
log(c("feature-merge-bedpe-R: merging the following bedpe files:",bedpe_fnames)) 
bedpe <- rbindlist(lapply(bedpe_fnames, fread, nrows=-100, nThread=4))
log(c("feature-merge-bedpe-R: done importing, start grouping"))
bedpe <- bedpe[,.(chrom, start=max(start), end=min(end), count=.N, strand, SNP_ABS_POS, SNP_SEQ, SNP_PARENT, SNP_VAR, SNP_TYPE, SNP_SUBTYPE),by=.(BC,SNP_ID)
                   ][order(BC,-count), .SD[1,], by=BC]
log(c("feature-merge-bedpe-R: done grouping, start export"))
fwrite(x=bedpe, file=output.bedpe, quote=FALSE, sep="\t")
log("feature-merge-bedpe-R: done export")

log("feature-merge-bedpe-R: Session Info:")
log(capture.output(print(sessionInfo())))
log("feature-merge-bedpe-R: finished\n\n")

quit(save='no')


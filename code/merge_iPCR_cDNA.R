# AUTHOR / DATE
#   Ludo Pagie; February 7, 2019; merge_iPCR_cDNA.R

# INTRO / BACKGROUND
#   Take a bedpe-like and merge it with the counts of a set of inputcDNA coutn
#   data files. The merge is done on barcode sequence in both data files. The
#   elements of the input bedpe (ie the SuREfragments) are all retained whereas
#   the cDNA counts are discarded if no corresponding barcode is present in the
#   SuRE fragments.

# USAGE / INPUT / ARGUMENTS / OUTPUT
# USAGE:
#   required:
#     Arguments should be given in the following order:
#     - output bedpe filename
#     - output log filename
#     - input bedpe filename
#     - remaining args are input cdna count filenames
#   optional:

# VERSIONS:


library(argparse)
library(data.table)
library(R.utils)
# set max number of threads to be used by any data.table function:
setDTthreads(threads=4)

# argument parsing
args <- commandArgs(trailingOnly=TRUE)
output.bedpe <- args[1]
output.log   <- args[2]
input.bedpe  <- args[3]
input.cdna   <- tail(args, -3)

# setup function for logging progress
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

# prev.time <- cur.time <- start.time <- Sys.time()
# cat("starting, time is", start.time)

# fns <- 
#   c(
#     "/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE42/SuRE-pipelineOutput/cDNA/SuRE42_B1_T1/SuRE42_B1_T1_trimmed_table.txt.gz",
#     "/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE42/SuRE-pipelineOutput/cDNA/SuRE42_B2_T1/SuRE42_B2_T1_trimmed_table.txt.gz",
#     "/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE42/SuRE-pipelineOutput/cDNA/SuRE42_B3_T1/SuRE42_B3_T1_trimmed_table.txt.gz",
#     "/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE42/SuRE-pipelineOutput/cDNA/SuRE42_HEPG2_B1_T1/SuRE42_HEPG2_B1_T1_trimmed_table.txt.gz",
#     "/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE42/SuRE-pipelineOutput/cDNA/SuRE42_HEPG2_B2_T1/SuRE42_HEPG2_B2_T1_trimmed_table.txt.gz"
# 
#     )
# dt <- lapply(fns, fread, header=F, nrow=-100000, nThread=1, col.names=c('count','BC'),key='BC')
# names(dt) <- sub('-T1_trimmed_table.txt.gz','',basename(fns))
# for(i in seq_along(dt)){colnames(dt[[i]])[1]=names(dt)[i]}

# cat("finished reading cDNA, time is", cur.time <- Sys.time())
# cat("time elapsed in this step is",cur.time-prev.time); prev.time <- cur.time
# cat("time elapsed since start is",cur.time-start.time); prev.time <- cur.time

# dt.merged <- Reduce(function(...) merge(..., all=T), yy)

# bedpe.fn <- '/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP181220_SuRE42-45-INDEL/SuRE42-1_pipelineOutput/iPCR/samples_merged/08_bedpe_merged/equal/5.bedpe.gz'
# import bedpe data
log("reading bedpe input")
bedpe <- fread(input.bedpe, nrow=-10,key='BC')

# cat("finished reading bedpe, time is", cur.time <- Sys.time())
# cat("time elapsed in this step is",cur.time-prev.time); prev.time <- cur.time
# cat("time elapsed since start is",cur.time-start.time); prev.time <- cur.time

# import cdna data and merge with bedpe data
log("reading cdna input")
for (f in input.cdna) {
  if (nrow(bedpe)>0)
    bedpe <- merge(bedpe, fread(f, header=F, nrow=-100000, nThread=1, col.names=c(sub('[_-]T1_trimmed_table.txt.gz','',basename(f)),'BC'),key='BC'), all.x=TRUE)
  else
    bedpe[[sub('[_-]T1_trimmed_table.txt.gz','',basename(f))]]=integer()
}
# set all counts not observed in cDNA data to 0
bedpe[is.na(bedpe)] <- 0

# cat("finished merging cDNA, time is", cur.time <- Sys.time())
# cat("time elapsed in this step is",cur.time-prev.time); prev.time <- cur.time
# cat("time elapsed since start is",cur.time-start.time); prev.time <- cur.time

# export bedpe data to file
log("writing output")
if(endsWith(output.bedpe, ".gz")) {
  fname <- sub(".gz$","",output.bedpe)
  fwrite(bedpe, fname, sep="\t", quote=FALSE)
  gzip(fname, destname=output.bedpe)
  unlink(fname)
} else {
  fwrite(bedpe, output.bedpe, sep="\t", quote=FALSE)
}

# cat("finished writing and compressing output, time is", cur.time <- Sys.time())
# cat("time elapsed in this step is",cur.time-prev.time); prev.time <- cur.time
# cat("time elapsed since start is",cur.time-start.time); prev.time <- cur.time
log("done")
quit(save='no')


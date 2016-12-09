#!/usr/bin/env Rscript

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

getScriptPath <- function() {
    argv <-commandArgs()
    dir <- na.omit(stringr::str_match(argv, "^--file=(.*)$")[,2])[1]
    if (!is.na(dir) && !is.null(dir))
        return(dir)
}
tryCatch(setwd(file.path(dirname(getScriptPath()), "..")),
         error=function(...) tsmsg("WARNING: Could not determine script path. Ensure that you are already in the correct directory."))

library(stringr)
library(magrittr)
library(GenomicRanges)
library(rtracklayer)
library(SummarizedExperiment)
library(dplyr)
library(csaw)
library(Matrix)
library(assertthat)

library(doParallel)
ncores <- getOption("mc.cores", default=1)
registerDoParallel(cores=ncores)
library(BiocParallel)
register(DoparParam())

tsmsg("Using ", ncores, " cores.")

windowCountsParallel <- function(bam.files, ..., filter=10) {
    reslist <- bplapply(bam.files, windowCounts, ..., filter=0)
    res <- do.call(cbind, reslist)
    keep <- rowSums(assay(res)) >= filter
    res[keep,]
}

tsmsg("Loading sample data")

sample.table <- readRDS("saved_data/samplemeta-ChIPSeq.RDS") %>%
    ## Compute full path to BAM file
    mutate(bam_file=sprintf("aligned/chipseq_bowtie2_hg38.analysisSet/%s/Aligned.bam", SRA_run)) %>%
    ## Ensure that days_after_activation is a factor and can't be
    ## interpreted as a numeric
    mutate(days_after_activation=days_after_activation %>%
               factor %>% `levels<-`(str_c("Day", levels(.)))) %>%
    rename(time_point=days_after_activation)

assert_that(all(file.exists(sample.table$bam_file)))

tsmsg("Loading blacklist regions")
blacklist <- import("saved_data/ChIPSeq-merged-blacklist.bed", format="bed")

## Standard nuclear chromosomes only. (chrM is excluded because it is
## not located in the nucleus and is thus not subject to histone
## modification. The unplaced scaffolds are mostly not large enough to
## contain even a single typically-sized peak, so little is lost by
## excluding them for this analysis.)
std.chr <- extractSeqlevels("Homo sapiens", "UCSC") %>% setdiff("chrM")
param <- readParam(restrict=std.chr, discard=blacklist)

{
    tsmsg("Counting reads in 10kb bins in ", nrow(sample.table), " samples.")
    binned <- windowCountsParallel(sample.table$bam_file, bin=TRUE, width=10000, param=param)
    colData(binned) %<>% {cbind(sample.table, .[c("totals", "ext")])}
    saveRDS(binned, "saved_data/csaw-bigbin-counts-10kb.RDS")
    rm(binned)
}

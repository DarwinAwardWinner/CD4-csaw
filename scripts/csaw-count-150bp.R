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
    tsmsg("Counting reads in 150bp windows in ", nrow(sample.table), " samples.")
    ## Using extension=spacing+1 ensures that each read overlaps
    ## exactly 2 windows.
    wcounts <- windowCounts(sample.table$bam_file, ext=151, spacing=150, param=param)
    colData(wcounts) %<>% {cbind(sample.table, .[c("totals", "ext")])}
    saveRDS(wcounts, "saved_data/csaw-window-counts-150bp.RDS")
    rm(wcounts)
}

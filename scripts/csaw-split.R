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
library(SummarizedExperiment)
library(dplyr)

tsmsg("Splitting big bin counts")
bigbin.counts <- readRDS("saved_data/csaw-bigbin-counts-10kb.RDS")
chipgroups <- split(seq_len(ncol(bigbin.counts)), colData(bigbin.counts)$chip_antibody)
for (i in names(chipgroups)) {
    sel <- chipgroups[[i]]
    x <- bigbin.counts[,sel]
    x <- x[colSums(assay(x)) >= 10,]
    fname <- sprintf("saved_data/csaw-bigbin-counts-%s-10kb.RDS", i)
    tsmsg("Saving ", fname)
    saveRDS(x, fname)
}
rm(bigbin.counts)

tsmsg("Splitting window counts")
window.counts <- readRDS("saved_data/csaw-window-counts-150bp.RDS")
chipgroups <- split(seq_len(ncol(window.counts)), colData(window.counts)$chip_antibody)
for (i in names(chipgroups)) {
    sel <- chipgroups[[i]]
    x <- window.counts[,sel]
    x <- x[colSums(assay(x)) >= 10,]
    fname <- sprintf("saved_data/csaw-window-counts-%s-150bp.RDS", i)
    tsmsg("Saving ", fname)
    saveRDS(x, fname)
    gc()
}

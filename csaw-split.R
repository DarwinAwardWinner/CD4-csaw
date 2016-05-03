#!/usr/bin/env Rscript

library(stringr)
library(magrittr)
library(openxlsx)
library(doParallel)
options(mc.cores=parallel::detectCores())
registerDoParallel(cores=parallel::detectCores())
library(SummarizedExperiment)
library(dplyr)
library(edgeR)
library(limma)
library(csaw)
setwd("~/Projects/CD4-csaw")

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

tsmsg("Processing big bin counts")
bigbin.counts <- readRDS("saved_data/bigbin-counts.RDS")
chipgroups <- split(seq_len(ncol(bigbin.counts)), colData(bigbin.counts)$ChIP)
for (i in names(chipgroups)) {
    sel <- chipgroups[[i]]
    x <- bigbin.counts[,sel]
    x <- x[colSums(assay(x)) >= 10,]
    fname <- sprintf("saved_data/bigbin-counts-%s.RDS", i)
    tsmsg("Saving ", fname)
    saveRDS(x, fname)
}

tsmsg("Processing window counts")
window.counts <- readRDS("saved_data/window-counts-147bp.RDS")
chipgroups <- split(seq_len(ncol(window.counts)), colData(window.counts)$ChIP)
for (i in names(chipgroups)) {
    sel <- chipgroups[[i]]
    x <- window.counts[,sel]
    x <- x[colSums(assay(x)) >= 10,]
    fname <- sprintf("saved_data/window-counts-%s-147bp.RDS", i)
    tsmsg("Saving ", fname)
    saveRDS(x, fname)
    gc()
}

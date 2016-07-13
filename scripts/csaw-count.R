#!/usr/bin/env Rscript

library(stringr)
library(magrittr)
library(GenomicRanges)
library(Rsubread)
library(openxlsx)
library(annotate)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(doParallel)
options(mc.cores=parallel::detectCores())
registerDoParallel(cores=parallel::detectCores())
library(BiocParallel)
register(MulticoreParam(parallel::detectCores()))
library(GenomicRanges)
library(rtracklayer)
library(SummarizedExperiment)
library(dplyr)
library(purrr)
library(csaw)
library(Matrix)

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

windowCountsParallel <- function(bam.files, ..., filter=10) {
    reslist <- bplapply(bam.files, windowCounts, ..., filter=0)
    res <- do.call(cbind, reslist)
    keep <- rowSums(assay(res)) >= filter
    res[keep,]
}

tsmsg("Loading sample data")

sample.table <- read.xlsx("data_files/ChIP-Seq/sample-tables.xlsx", "Samples") %>%
    set_colnames(make.unique(colnames(.))) %>%
    ## Select/compute/rename desired columns
    transmute(
        Sample=title,
        ChIP=`characteristics:.sampletype`,
        Celltype=`characteristics:.celltype`,
        Activated=`characteristics:.activated`,
        Day=`characteristics:.days.after.activation`,
        Donor=`characteristics:.donor.ID`,
        file.name=raw.file) %>%
    ## Make sure no factorial variables can be accidentally
    ## numericized by prefixing them with letters
    mutate(
        Donor=sprintf("Dn%s", Donor),
        Day=sprintf("D%i", Day)
    ) %>%
    ## Compute full path to bam file
    mutate(
        bampath=file.path("data_files/ChIP-Seq", file.name)
    ) %>%
    ## Factor variables with proper levels (not ASCIIbetical order)
    mutate(
        ChIP=factor(ChIP),
        Celltype=factor(Celltype, levels=c("Naive", "Memory")),
        Day=factor(Day, levels=sprintf("D%s", c(0,1,5,14))),
        TreatmentGroup=interaction(Celltype, Day, sep=""),
        Group=interaction(ChIP, TreatmentGroup, sep="."),
        Donor=factor(Donor)
    ) %>% set_rownames(.$Sample)

stopifnot(all(file.exists(sample.table$bampath)))

tsmsg("Loading blacklist regions")
blacklist <- import("saved_data/wgEncodeDacMapabilityConsensusExcludable.bed.gz", format="bed")

## Standard nuclear chromosomes
std.chr <- extractSeqlevels("Homo sapiens", "UCSC") %>% setdiff("chrM")
param <- readParam(restrict=std.chr, discard=blacklist)

{
    tsmsg("Counting in 10 kb bins")
    binned <- windowCountsParallel(sample.table$bampath, bin=TRUE, width=10000, param=param)
    colData(binned) %<>% {cbind(sample.table, .[c("totals", "ext")])}
    saveRDS(binned, "saved_data/bigbin-counts-10kb.RDS")
    rm(binned)
    tsmsg("Counting in 147-bp windows")
    wcounts <- windowCounts(sample.table$bampath, ext=147, spacing=147, param=param)
    colData(wcounts) %<>% {cbind(sample.table, .[c("totals", "ext")])}
    saveRDS(wcounts, "saved_data/window-counts-147bp.RDS")
    rm(wcounts)
}

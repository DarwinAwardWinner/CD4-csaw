#!/usr/bin/env Rscript

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

library(here)
library(stringr)
library(GenomicRanges)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(purrr)
library(csaw)
library(rctutils)

use_futures("multicore")

tsmsg("Loading sample data")

sample.table <- readRDS(here("saved_data", "samplemeta-ChIPSeq.RDS")) %>%
    ## Compute full path to BAM file
    mutate(bam_file=here("aligned", "chipseq_bowtie2_hg38.analysisSet", SRA_run, "Aligned.bam")) %>%
    ## Ensure that days_after_activation is a factor and can't be
    ## interpreted as a numeric
    mutate(days_after_activation=days_after_activation %>%
               factor %>% `levels<-`(str_c("Day", levels(.)))) %>%
    rename(time_point=days_after_activation)

stopifnot(all(file.exists(sample.table$bam_file)))

tsmsg("Loading blacklist regions")
blacklist <- import(here("saved_data", "ChIPSeq-merged-blacklist.bed"), format="bed")

## Standard nuclear chromosomes only. (chrM is excluded because it is
## not located in the nucleus and is thus not subject to histone
## modification. The unplaced scaffolds are mostly not large enough to
## contain even a single typically-sized peak, so little is lost by
## excluding them for this analysis.)
std.chr <- extractSeqlevels("Homo sapiens", "UCSC") %>% setdiff("chrM")
param <- readParam(restrict=std.chr, discard=blacklist)
param.dedup.on <- reform(param, dedup=TRUE)
param.dedup.on.noBL <- reform(param.dedup.on, discard=GRanges())

## Determine fragment length using cross-correlation function, see
## csaw UG 2.4.1

tsmsg("Beginning CCF computation")

sample.ccf.noBL <-
    bplapply(sample.table$bam_file, function(bam) {
        tsmsg("Computing no-blacklist CCF for ", bam)
        correlateReads(bam, max.dist=1000, cross=TRUE,
                       param=param.dedup.on.noBL)

    })

sample.ccf <-
    bplapply(sample.table$bam_file, function(bam) {
        tsmsg("Computing no-blacklist CCF for ", bam)
        correlateReads(bam, max.dist=1000, cross=TRUE,
                       param=param.dedup.on)

    })

names(sample.ccf.noBL) <- names(sample.ccf) <- sample.table$SampleName
saveRDS(sample.ccf, here("saved_data", "chipseq-ccf.RDS"))
saveRDS(sample.ccf.noBL, here("saved_data", "chipseq-ccf-noBL.RDS"))
#!/usr/bin/env Rscript

library(here)
library(MASS)
library(magrittr)
library(dplyr)
library(reshape2)
library(assertthat)
library(csaw)
library(edgeR)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(rtracklayer)
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
    ## Use better names for things
    rename(Sample=SampleName,
           ChIP=chip_antibody,
           TimePoint=days_after_activation,
           CellType=cell_type,
           Donor=donor_id) %>%
    mutate(TreatmentGroup=interaction(CellType, TimePoint, sep="."))

assert_that(all(file.exists(sample.table$bam_file)))

input.sample.table <- sample.table %>% filter(ChIP == "input")

tsmsg("Preparing SeqInfo for HG38 standard chromosomes")
std.chr <- extractSeqlevels("Homo sapiens", "UCSC")
std.seqinfo <- BSgenome.Hsapiens.UCSC.hg38 %>%
    seqinfo %>% keepSeqlevels(std.chr, pruning.mode="coarse")

set.seed(1986)

tsmsg("Doing 1kb window counts")

param <- readParam(restrict=std.chr)
binned1kb <- windowCountsParallel(input.sample.table$bam_file, bin=TRUE, width=1000, param=param, filter=0)
binned1kb %<>% keepSeqlevels(std.chr)
seqinfo(binned1kb) <- std.seqinfo[seqlevels(binned1kb)]
colnames(binned1kb) <- input.sample.table$Sample

tsmsg("Fitting NB GLM to each sample")

nbfits <- bplapply(colnames(binned1kb), function(i) {
    tsmsg("Fitting NB model for ", i)
    x <- assay(binned1kb[,i])
    glm.nb(x ~ 1)
})
names(nbfits) <- colnames(binned1kb)

tsmsg("Saving NB GLM fits")
saveRDS(nbfits, here("saved_data", "ChIPSeq-input-depth-NBGLM-fits.RDS"))

mu <- sapply(nbfits, . %>% coef %>% exp %>% unname)
size <- sapply(nbfits, `[[`, "theta")

colData(binned1kb) %<>%
    cbind(data.frame(
        nbinom.mu=mu, nbinom.size=size))

tsmsg("Computing theoretical and empirical quantiles for each bin")

assay(binned1kb, "pnbinom") <-
    pnbinom(assay(binned1kb, "counts"),
            size=expandAsMatrix(size, dim(binned1kb)),
            mu=expandAsMatrix(mu, dim(binned1kb)))
assay(binned1kb, "pemp") <- assay(binned1kb, "counts") %>% apply(2, . %>% {ecdf(.)(.)})

tsmsg("Saving bin counts and quantiles")
saveRDS(binned1kb, here("saved_data", "window-counts-input-unfiltered-1kb.RDS"))

tsmsg("Generating greylist (bins with NB quantile > 0.99 in 2 or more samples)")

over99 <- assay(binned1kb, "pnbinom") > 0.99
greylisted <- rowSums(over99) >= 2
glranges <- rowRanges(binned1kb)[greylisted]

tsmsg("Saving greylist")

saveRDS(glranges, here("saved_data", "ChIPSeq-input-greylist.RDS"))
export(glranges, here("saved_data", "ChIPSeq-input-greylist.bed"))

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
library(reshape2)
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
param.dedup.on <- reform(param, dedup=TRUE)

## Determining fragment length using cross-correlation function, see
## csaw UG 2.4.1
sample.ccf <- lapply(sample.table$bampath,
                     . %T>% tsmsg("Computing CCF for ", .) %>%
                     correlateReads(max.dist=1000, cross=TRUE, param=param.dedup.on))
names(sample.ccf) <- sample.table$Sample
saveRDS(sample.ccf, "saved_data/csaw-ccf.RDS")
ccftable <- sample.ccf %>%
    lapply(function(x) data.frame(Delay=seq(from=0, by=1, length.out=length(x)), CCF=x)) %>%
    melt(measure.vars = "CCF") %>%
    rename(Sample=L1) %>%
    dcast(Sample + Delay ~ variable) %>%
    group_by(Sample) %>%
    mutate(RelCCF=CCF/max(CCF)) %>%
    ungroup %>%
    inner_join(sample.table, by="Sample")

refline.table <- data.frame(
    Reference=c("Read Length (100bp)", "Nucleosome Footprint (147bp)"),
    Xintercept=c(100, 147))

{
    baseplot <- ggplot(ccftable) +
        facet_wrap(~ChIP, scales="free") +
        aes(x=Delay, y=CCF, group=Sample, color=TreatmentGroup, linetype=NA) +
        ylim(0,NA) +
        geom_vline(data=refline.table,
                   aes(xintercept=Xintercept, linetype=Reference, color=NA),
                   color="black", alpha=0.5) +
        scale_color_hue(name="Group") +
        scale_linetype(name="Reference") +
        theme(legend.position="bottom")
    p <- list(
        Raw=baseplot +
            geom_line(size=0.25, linetype="solid") +
            ggtitle("Cross-Correlation Function, Raw"),
        loess_span0.05=baseplot +
            geom_smooth(fill=NA, method="loess", span=0.05, n=500, size=0.25, linetype="solid") +
            ggtitle("Cross-Correlation Function, Loess-Smoothed (span = 0.05)"),
        loess_span0.075=baseplot +
            geom_smooth(fill=NA, method="loess", span=0.075, n=500, size=0.25, linetype="solid") +
            ggtitle("Cross-Correlation Function, Loess-Smoothed (span = 0.075)"),
        loess_span0.1=baseplot +
            geom_smooth(fill=NA, method="loess", span=0.1, n=500, size=0.25, linetype="solid") +
            ggtitle("Cross-Correlation Function, Loess-Smoothed (span = 0.1)"))
    pdf("results/csaw/CCF-plots.pdf", width=16, height=16)
    print(p)
    dev.off()
    pdf("results/csaw/CCF-plots-relative.pdf", width=16, height=16)
    print(lapply(p, . %>% add(aes(y=RelCCF))))
    dev.off()
}

max.ccf <- ccftable %>% filter(Delay >= 100) %>% group_by(Sample) %>% summarize(CCF=max(CCF))
summary(max.ccf$Delay)

## Determining window size using findMaxima & profileSites, see csaw
## UG 2.5

sample.site.profiles <- lapply(sample.table$bampath, function(bam) {
    tsmsg("Profiling maxima for ", bam)
    windowed <- windowCounts(bam, spacing=50, width=50, ext=147, param=param.dedup.on, filter=20)
    rwsms <- rowSums(assay(windowed))
    maxed <- findMaxima(rowRanges(windowed), range=5000, metric=rwsms)
    profileSites(bam, rowRanges(windowed)[maxed], range=10000, param=param.dedup.on, weight=1/rwsms[maxed])
})
names(sample.site.profiles) <- sample.table$Sample
saveRDS(sample.site.profiles, "saved_data/csaw-siteprof.RDS")

profile.table <- sample.site.profiles %>%
    lapply(function(x) data.frame(Distance=as.numeric(names(x)), RelativeCoverage=x)) %>%
    melt(id.vars="Distance") %>%
    rename(Sample=L1) %>%
    dcast(Sample + Distance ~ variable) %>%
    inner_join(sample.table, by="Sample") %>%
    group_by(Sample)##  %>%
    ## mutate(Deriv=grad(approxfun(), x)) %$% unique(m)

{
    baseplot <- ggplot(profile.table) +
        facet_wrap(~ChIP, scales="free") +
        aes(x=Distance, y=RelativeCoverage, group=Sample, color=TreatmentGroup) +
        ylim(0,NA)
    p <- list(
        Raw10kb=baseplot +
            geom_line(size=0.25) +
            coord_cartesian(xlim=c(-10000, 10000)) +
            ggtitle("Relative Coverage Around Maxima, Raw, 10kb Radius"),
        Raw1kb=baseplot +
            geom_line(size=0.25) +
            coord_cartesian(xlim=c(-1000, 1000)) +
            ggtitle("Relative Coverage Around Maxima, Raw, 1kb Radius"),
        Raw300bp=baseplot +
            geom_line(size=0.25) +
            coord_cartesian(xlim=c(-300, 300)) +
            ggtitle("Relative Coverage Around Maxima, Raw, 300bp Radius"))
    pdf("results/csaw/site-profile-plots.pdf", width=16, height=16)
    print(p)
    dev.off()
}

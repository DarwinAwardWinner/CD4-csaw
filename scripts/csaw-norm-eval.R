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

positional_args <- function(argnames, args=commandArgs(trailingOnly = TRUE)) {
    argnames <- as.character(argnames)
    if (length(args) != length(argnames)) {
        stop(sprintf("Got %s arguments, but expected %s arguments: %s",
                     length(args), length(argnames), paste(argnames, collapse=", ")))
    }
    setNames(args, argnames)
}

library(stringr)
library(magrittr)
library(SummarizedExperiment)
library(readr)
library(edgeR)
library(limma)
library(csaw)
library(ggplot2)
library(scales)
library(rtracklayer)
library(parallel)
library(dplyr)

suppressPlot <- function(arg) {
    png("/dev/null")
    on.exit(dev.off())
    result <- arg
    result
}

## https://charlesjb.github.io/How_to_import_narrowPeak/
extraCols_narrowPeak <- c(signalValue = "numeric", log10pValue = "numeric",
                          log10qValue = "numeric", peak = "integer")
import.narrowPeak <- function(...) {
    import(..., format="BED", extraCols=extraCols_narrowPeak)
}

asDGEList <- function(...) {
    dge <- csaw::asDGEList(...)
    ## This seems to save a lot of memory
    rownames(dge$counts) <- NULL
    dge
}

## Version of cpm that uses an offset matrix instead of lib sizes
cpmWithOffset <- function(dge, offset=expandAsMatrix(getOffset(dge), dim(dge)),
                          log = FALSE, prior.count = 0.25, ...) {
    x <- dge$counts
    effective.lib.size <- exp(offset)
    if (log) {
        prior.count.scaled <- effective.lib.size/mean(effective.lib.size) * prior.count
        effective.lib.size <- effective.lib.size + 2 * prior.count.scaled
    }
    effective.lib.size <- 1e-06 * effective.lib.size
    if (log)
        log2((x + prior.count.scaled) / effective.lib.size)
    else x / effective.lib.size
}

withGC <- function(expr) {
    on.exit(gc())
    return(expr)
}

autoFactorize <- function(df) {
    for (i in colnames(df)) {
        if (is.character(df[[i]]) && anyDuplicated(df[[i]])) {
            df[[i]] %<>% factor
        }
    }
    df
}

chip <- positional_args("chip")

tsmsg("Reading data for ", chip)
bigbin.counts <- readRDS(sprintf("saved_data/csaw-bigbin-counts-%s-10kb.RDS", chip))
colnames(bigbin.counts) <- colData(bigbin.counts)$SampleName
window.counts <- readRDS(sprintf("saved_data/csaw-window-counts-%s-150bp.RDS", chip))
colnames(window.counts) <- colData(window.counts)$SampleName

chip.peaks <- import.narrowPeak(
    sprintf(fmt="peak_calls/epic_hg38.analysisSet/%s_condition.ALL_donor.ALL/peaks_noBL_IDR.narrowPeak", chip)) %>%
    .[.$log10qValue >= -log10(0.05)]

tsmsg("Filtering low-count bins")
## Filter to windows with at least 10 counts total
keep <- rowSums(assay(window.counts)) >= 10
window.counts %<>% .[keep,]
invisible(gc())

## Compute background normalization from big bins
tsmsg("Computing composition normalization from background")
colData(window.counts)$CompNormFactors <- normOffsets(bigbin.counts, type="scaling")

tsmsg("Computing average window abundances")
abundances <- withGC(aveLogCPM(asDGEList(window.counts)))

tsmsg("Computing window filter statistic")
## Compute enirchment of windows above average global background
filter.stat <- filterWindows(window.counts, bigbin.counts, type="global")

filter.stat$global.bg <- filter.stat %$% {abundances - filter} %>% mean
filter.stat.df <- filter.stat %$% {
    rbind(
        data.frame(BinType="ChIP", Abundance=abundances %>% quantile(seq(0, 1, length.out=1e6))),
        data.frame(BinType="BG", Abundance=back.abundances))
}
filter.stat.thresholds <- data.frame(Line=c("Background", "Threshold"),
                                     Abundance=filter.stat$global.bg + c(0, log2(3)))

## Pick a "high-abundance" threshold
filter.threshold <- filter.stat$global.bg + log2(3)

tsmsg("Computing efficiency normalization from high-abundance windows")
high.ab.window.counts <- window.counts[abundances >= filter.threshold,]
colData(window.counts)$HANormFactors <- normOffsets(high.ab.window.counts, type="scaling")

peak.overlap <- overlapsAny(window.counts, chip.peaks)
peak.window.counts <- window.counts[peak.overlap,]
tsmsg("Computing efficiency normalization from peak windows")
pnf <- normOffsets(peak.window.counts, type="scaling")
colData(window.counts)$PeakNormFactors <- pnf
colData(peak.window.counts)$PeakNormFactors <- pnf

withGC({
    rm(window.counts)
    rm(bigbin.counts)
})

## Create DGEList and populate $genes and $samples
dge <- asDGEList(peak.window.counts)
dge$genes <- rowRanges(peak.window.counts) %>% as.data.frame %>% lapply(Rle) %>% DataFrame
dge$samples %<>% cbind(colData(peak.window.counts)) %>% as.data.frame %>%
    autoFactorize %>%
    mutate(group = interaction(cell_type, str_replace(time_point, "Day", "D"), sep=""),
           nf.logratio = log2(PeakNormFactors / CompNormFactors))
colnames(dge) <- rownames(dge$samples) <- dge$samples$SampleName

tsmsg("Filtering")
ave.count.threshold <- 5
thresh <- aveLogCPM(ave.count.threshold, lib.size=mean(dge$samples$totals))
ab <- aveLogCPM(dge)
present <- ab >= thresh

dge <- dge[present,]
dge$offset <- normOffsets(dge$counts, lib.sizes=dge$samples$lib.size, type="loess")

## Try both normalizations
dge.comp <- dge
dge.eff <- dge
dge.loess <- dge

dge.comp$samples %<>% mutate(norm.factors=CompNormFactors)
dge.comp$offset <- NULL
dge.eff$samples %<>% mutate(norm.factors=PeakNormFactors)
dge.eff$offset <- NULL
dge.loess$samples %<>% mutate(norm.factors=1)

tsmsg("Creating MDS plots")
mdsdist.comp <- suppressPlot(plotMDS(dge.comp, top=15000)) %$% distance.matrix %>% as.dist
mdsdims.comp <- suppressWarnings(cmdscale(mdsdist.comp, k=ncol(dge.comp)-1)) %>%
    set_colnames(str_c("PC", seq_len(ncol(.)))) %>%
    data.frame %>%
    cbind(dge.comp$samples, .) %>%
    arrange(cell_type, time_point, donor_id)
pc.cols <- str_detect(colnames(mdsdims.comp),"PC[0-9]+")
flipsign <- mdsdims.comp %>% filter(time_point=="Day0") %>%
    .[pc.cols] %>% colMeans %>% sign %>% multiply_by(-1)
mdsdims.comp[names(flipsign)] %<>% scale(center=FALSE, scale=flipsign)
mdsplot.comp <- ggplot(mdsdims.comp) +
    aes(x=PC1, y=PC2, shape=cell_type, color=time_point, group=donor_id:cell_type, linetype=donor_id) +
    geom_point(size=5) +
    geom_path(aes(color=NA)) +
    coord_fixed() +
    ggtitle("MDS Plot, Composition Normalized")

mdsdist.eff <- suppressPlot(plotMDS(dge.eff, top=15000)) %$% distance.matrix %>% as.dist
mdsdims.eff <- suppressWarnings(cmdscale(mdsdist.eff, k=ncol(dge.eff)-1)) %>%
    set_colnames(str_c("PC", seq_len(ncol(.)))) %>%
    data.frame %>%
    cbind(., dge.eff$samples) %>%
    arrange(cell_type, time_point, donor_id)
pc.cols <- str_detect(colnames(mdsdims.eff),"PC[0-9]+")
flipsign <- mdsdims.eff %>% filter(time_point=="Day0") %>%
    .[pc.cols] %>% colMeans %>% sign %>% multiply_by(-1)
mdsdims.eff[names(flipsign)] %<>% scale(center=FALSE, scale=flipsign)
mdsplot.eff <- ggplot(mdsdims.eff) +
    aes(x=PC1, y=PC2, shape=cell_type, color=time_point, group=donor_id:cell_type, linetype=donor_id) +
    geom_point(size=5) +
    geom_path(aes(color=NA)) +
    coord_fixed() +
    ggtitle("MDS Plot, Efficiency Normalized")

yoffset <- cpmWithOffset(dge.loess, prior.count=2)
mdsdist.loess <- suppressPlot(plotMDS(yoffset, top=15000)) %$% distance.matrix %>% as.dist
mdsdims.loess <- suppressWarnings(cmdscale(mdsdist.loess, k=ncol(dge.loess)-1)) %>%
    set_colnames(str_c("PC", seq_len(ncol(.)))) %>%
    data.frame %>%
    cbind(., dge.loess$samples) %>%
    arrange(cell_type, time_point, donor_id)
pc.cols <- str_detect(colnames(mdsdims.loess),"PC[0-9]+")
flipsign <- mdsdims.loess %>% filter(time_point=="Day0") %>%
    .[pc.cols] %>% colMeans %>% sign %>% multiply_by(-1)
mdsdims.loess[names(flipsign)] %<>% scale(center=FALSE, scale=flipsign)
mdsplot.loess <- ggplot(mdsdims.loess) +
    aes(x=PC1, y=PC2, shape=cell_type, color=time_point, group=donor_id:cell_type, linetype=donor_id) +
    geom_point(size=5) +
    geom_path(aes(color=NA)) +
    coord_fixed() +
    ggtitle("MDS Plot, Loess Normalized")

tsmsg("Estimating dispersions")
design <- model.matrix(~0 + group + donor_id, dge$samples)
colnames(design) %<>%
    str_replace("^group", "") %>%
    str_replace("^donor_idD", "Dn")

dge.comp %<>% estimateDisp(design, robust=TRUE)
fit.comp <- glmQLFit(dge.comp, design, robust=TRUE)

dge.eff %<>% estimateDisp(design, robust=TRUE)
fit.eff <- glmQLFit(dge.eff, design, robust=TRUE)

dge.loess %<>% estimateDisp(design, robust=TRUE)
fit.loess <- glmQLFit(dge.loess, design, robust=TRUE)

{
    tsmsg("Printing plots")
    pdf(sprintf("plots/csaw/%s-norm-eval.pdf", chip), width=8, height=8)
    plotBCV(dge.comp, ylim=c(0, 1.5), main=sprintf("BCV Plot, Composition Normalized (prior df=%.2f)", max(dge.comp$prior.df)))
    plotBCV(dge.eff, ylim=c(0, 1.5), main=sprintf("BCV Plot, Efficiency Normalized (prior df=%.2f)", max(dge.eff$prior.df)))
    plotBCV(dge.loess, ylim=c(0, 1.5), main=sprintf("BCV Plot, Loess Normalized (prior df=%.2f)", max(dge.loess$prior.df)))
    plotQLDisp(fit.comp, ylim=c(0.5, 1.9), main=sprintf("QL Dispersion Plot, Composition Normalized (prior df=%.2f)", max(fit.comp$df.prior)))
    plotQLDisp(fit.eff, ylim=c(0.5, 1.9), main=sprintf("QL Dispersion Plot, Efficiency Normalized (prior df=%.2f)", max(fit.eff$df.prior)))
    plotQLDisp(fit.loess, ylim=c(0.5, 1.9), main=sprintf("QL Dispersion Plot, Loess Normalized (prior df=%.2f)", max(fit.loess$df.prior)))
    print(mdsplot.eff)
    print(mdsplot.comp)
    print(mdsplot.loess)
    dev.off()
}

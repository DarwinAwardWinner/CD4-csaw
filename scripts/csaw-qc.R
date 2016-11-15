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
library(openxlsx)
library(SummarizedExperiment)
library(readr)
library(edgeR)
library(limma)
library(csaw)
library(ggplot2)
library(scales)
library(rtracklayer)
library(RColorBrewer)
library(reshape2)
library(parallel)
library(ks)
library(dplyr)

power_trans <- function(pow) {
    name <- sprintf("^%s", pow)
    trans_new(name,
              transform=function(x) x ^ pow,
              inverse=function(x) x ^ (1/pow),
              domain =c(0,Inf))
}

clamp_trans <- function(lower_threshold=0, upper_threshold=1) {
    name <- sprintf("Clamp values outside of [%s, %s]", lower_threshold, upper_threshold)
    trans_new(name,
              transform=function(x) pmin(upper_threshold, pmax(lower_threshold, x)),
              ## transform is only invertible for part of the range
              inverse=identity)
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

chip <- positional_args("chip")

tsmsg("Reading data for ", chip)
bigbin.counts <- readRDS(sprintf("saved_data/csaw-bigbin-counts-%s-10kb.RDS", chip))
colnames(bigbin.counts) <- colData(bigbin.counts)$SampleName
window.counts <- readRDS(sprintf("saved_data/csaw-window-counts-%s-150bp.RDS", chip))
colnames(window.counts) <- colData(window.counts)$SampleName

chip.peaks <- import.narrowPeak(
    sprintf(fmt="peak_calls/epic_hg38.analysisSet/%s_condition.ALL_donor.ALL/peaks_noBL_IDR.narrowPeak", chip)) %>%
    .[.$log10qValue >= -log10(0.05)]

tsmsg("Computing fraction of reads in peaks (FRiP)")
## (Frag length - 1) / (bin width) + 1, assuming non-overlapping,
## non-gapped bins
count.duplication.factor <- (colData(window.counts)$ext - 1) / median(width(rowRanges(window.counts))) + 1
colData(window.counts)$RiP.approx <- {
    window.counts %>%
        subsetByOverlaps(chip.peaks) %>%
        assay("counts") %>% colSums %>%
        divide_by(count.duplication.factor)
}
colData(window.counts)$FRiP.approx <- {
    colData(window.counts)$RiP.approx / colData(window.counts)$totals
}

tsmsg("Counting number of non-zero bins")
colData(window.counts)$NonZeroBins <- withGC(colSums(assay(window.counts) > 0))

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

## Pick a threshold
filter.threshold <- filter.stat$global.bg + log2(3)

tsmsg("plotting abundance vs peak overlap")
## Plot abundance vs peak containment
keep.table <- data.frame(
    PeakKeep=overlapsAny(window.counts, GRangesList(list(chip.peaks))),
    AbundanceKeep=abundances >= filter.threshold,
    Abundance=abundances)

tsmsg("Computing efficiency normalization from high-abundance windows")
high.ab.window.counts <- window.counts[abundances >= filter.threshold,]
colData(window.counts)$HANormFactors <- normOffsets(high.ab.window.counts, type="scaling")

table(keep.table[c("PeakKeep", "AbundanceKeep")])

{
    p <- ggplot(keep.table) +
        aes(y=Abundance, x=PeakKeep) +
        geom_violin() +
        xlab("Peak Overlap") +
        ggtitle("Window Mean logCPM by Peak Status")
    pdf(sprintf("plots/csaw/%s-window-abundance-vs-peaks.pdf", chip))
    print(p)
    dev.off()
}

## Based on above plot, we decide to use the peaks as filter
## criterion
peak.overlap <- overlapsAny(window.counts, chip.peaks)
peak.window.counts <- window.counts[peak.overlap,]
tsmsg("Computing efficiency normalization from peak windows")
pnf <- normOffsets(peak.window.counts, type="scaling")
colData(window.counts)$PeakNormFactors <- pnf
colData(peak.window.counts)$PeakNormFactors <- pnf

offsets <- normOffsets(peak.window.counts, type="loess")
dimnames(offsets) <- dimnames(peak.window.counts)

## Create DGEList and populate $genes and $samples
dge <- asDGEList(window.counts)
dge$genes <- rowRanges(window.counts) %>% as.data.frame %>% lapply(Rle) %>% DataFrame
## aveLogCPM takes forever to run, so let's save the result here
dge$genes$Abundance <- abundances
dge$samples %<>% cbind(colData(window.counts)) %>% as.data.frame %>%
    mutate(group=interaction(cell_type, str_replace(time_point, "Day", "D"), sep=""))
dge$samples$nf.logratio <- dge$samples %$% log2(PeakNormFactors / CompNormFactors)
colnames(dge) <- rownames(dge$samples) <- dge$samples$SampleName

{
    tsmsg("Testing norm factors for association with experimental factors")
    xvars <- c("cell_type", "time_point", "donor_id", "group", "FRiP.approx")
    yvars <- c("lib.size", "CompNormFactors", "HANormFactors", "PeakNormFactors", "nf.logratio", "FRiP.approx")

    formulas <- list()
    for (xv in xvars) {
        for (yv in yvars) {
            if (xv != yv) {
                modname <- sprintf("%s.vs.%s", yv, xv)
                formulas[[modname]] <- as.formula(sprintf("%s ~ %s", yv, xv))
            }
        }
    }

    mods <- lapply(formulas, . %>% lm(data=dge$samples))
    tests <- mods %>% lapply(anova)
    results <- tests %>% sapply(. %$% `Pr(>F)` %>% .[1]) %>% data_frame(Test=names(.) , PValue=., FDR=p.adjust(., "BH")) %>% arrange(PValue, FDR)
    write.xlsx(results, sprintf("results/csaw/%s-normfactor-tests.xlsx", chip))
}

{
    tsmsg("Plotting norm factors vs experimental factors")
    ## Plot norm factors vs lib sizes
    nf.vs.ls <- dge$samples %>%
        select(LibSize=lib.size, HiAbNorm=HANormFactors,
               PeakNorm=PeakNormFactors, CompNorm=CompNormFactors) %>%
        melt(id.vars="LibSize", variable.name = "NormType", value.name="NormFactor")
    p <- list(ggplot(nf.vs.ls) +
              aes(x=LibSize, y=NormFactor, group=NormType, color=NormType, fill=NormType) +
              geom_point() + geom_smooth(alpha=0.15) +
              scale_x_continuous(trans=log_trans(2)) +
              scale_y_continuous(trans=log_trans(2)) +
              ggtitle("Norm Factors vs Library Sizes"))
    ## Plot norm factors & lib sizes vs FRiP
    nf.ls.vs.frip <- dge$samples %>%
        select(LibSize=lib.size, HiAbNorm=HANormFactors,
               PeakNorm=PeakNormFactors, CompNorm=CompNormFactors,
               FRiP=FRiP.approx, group=group, donor_id=donor_id) %>%
        melt(id.vars=c("FRiP", "group", "donor_id"), variable.name = "Var", value.name="Value")
    p <- c(p, list(ggplot(nf.ls.vs.frip) +
                   aes(x=FRiP, y=Value, shape=donor_id, color=group,
                       group=1) +
                   geom_smooth(method="lm", alpha=0.15, color="gray65", size=0.5) +
                   geom_point(size=2) +
                   scale_shape_manual(values=c(15:18)) +
                   xlab("Fraction of Reads In Peaks") +
                   theme(legend.position="bottom",
                         legend.direction="horizontal") +
                   facet_wrap(~Var, ncol=2, scales="free_y")))
    ## Plot lib sizes & norm factors vs experimental factors
    ls.nf.vs.exp <- dge$samples %>%
        select(
            SampleName, cell_type, time_point, donor_id,
            group=group, LibSize=lib.size,
            HighAbundanceNormFactor=HANormFactors,
            PeakNormFactor=PeakNormFactors,
            CompositionNormFactor=CompNormFactors,
            FRiP=FRiP.approx) %>%
        melt(id.vars=c("SampleName", "cell_type", "time_point", "donor_id",
                       "group"),
             variable.name="Factor", value.name="Value")
    p0 <- ggplot(ls.nf.vs.exp) +
        aes(y=Value) + scale_y_continuous(trans=log_trans(2)) +
        geom_boxplot(color="gray65", outlier.colour = NA, alpha=0.5) +
        geom_point(position=position_jitter(width=0.25)) +
        facet_wrap(~Factor, scales="free_y", nrow=3) +
        theme(axis.text.x=element_text(angle = 30, hjust = 1))
    vars.to.plot <- c("cell_type", "time_point", "donor_id", "group")
    for (i in vars.to.plot) {
        p <- c(p, list(p0 + aes_string(x=i) +
                       ggtitle(sprintf("Library Sizes and Norm Factors vs %s", i))))
    }
    pdf(sprintf("plots/csaw/%s-normfactors.pdf", chip), width=8, height=8)
    print(p)
    dev.off()
}

## Select the most and least extreme samples for plotting MA plots

middle.samples <- dge$samples$nf.logratio %>% abs %>% order
cn.higher.samples <- dge$samples$nf.logratio %>% order
pn.higher.samples <- rev(cn.higher.samples)

tsmsg("Computing logCPM")
logcpm <- cpm(dge, log=TRUE)
bigbin.logcpm <- cpm(asDGEList(bigbin.counts), log=TRUE)
peak.logcpm <- logcpm[peak.overlap,]

peak.logcpm.loess <- cpmWithOffset(dge[peak.overlap,], offset=offsets + getOffset(dge[peak.overlap,]), log=TRUE)

getLineData <- function(s1, s2) {
    c(Comp="CompNormFactors",
      Peaks="PeakNormFactors",
      HiAb="HANormFactors") %>%
        sapply(. %>% dge$samples[[.]] %>% log2 %>% {.[s2] - .[s1]}) %>%
        data.frame(NormFactor=., NormType=names(.))
}

getOffsetLineData <- function(s1, s2, n=1000) {
    x <- data.frame(A=dge$genes$Abundance[peak.overlap],
                    Offset=(offsets[,s2] - offsets[,s1]) / log(2))
    f <- splinefun(x$A, x$Offset)
    data.frame(A=seq(from=min(x$A), to=max(x$A), length.out = n)) %>%
        mutate(Offset=f(A))
}

doMAPlot <- function(logcpm.matrix, s1, s2, linedata=getLineData(s1, s2)) {
    pointdata <- data.frame(S1=logcpm.matrix[,s1], S2=logcpm.matrix[,s2]) %>%
        transmute(A=(S1+S2)/2, M=S2-S1) %>%
        filter(A >= -2)
    ## Compute bandwidth and kernel smooth surface
    H <- pointdata %>% Hbcv.diag(binned=TRUE) %>% divide_by(4)
    k <- pointdata %>%
        as.matrix %>%
        kde(gridsize=1024, bgridsize=rep(1024, 2), verbose=TRUE,
            H=H, binned=TRUE)
    ## Sometimes the estimate goes a bit negative, which is no good

    densdata <- melt(k$estimate) %>%
        transmute(
            A=k$eval.points[[1]][Var1],
            M=k$eval.points[[2]][Var2],
            Density=value %>% pmax(0),
            ## Part of a hack to make the alpha look less bad
            AlphaDens=value %>% pmax(1e-15))

    p <- ggplot(pointdata) +
        ## MA Plot density
        geom_raster(aes(x=A, y=M, fill=Density, alpha=AlphaDens),
                    data=densdata,
                    interpolate=TRUE) +
        scale_fill_gradientn(colors=suppressWarnings(brewer.pal(Inf, "Blues")),
                             trans=power_trans(1/8),
                             name="Density") +
        scale_alpha_continuous(trans=power_trans(1/40), guide=FALSE)
    if (!is.null(linedata) && nrow(linedata) > 0) {
        p <- p +
            ## Normalization lines
            geom_hline(data=linedata, aes(yintercept=NormFactor, color=NormType)) +
            scale_color_discrete(name="Norm Type")
    }
    p <- p +
        ## Loess curve
        geom_smooth(aes(x=A, y=M), span=0.3, fill=NA, color="black") +
        ## Scales
        scale_x_continuous(name="log2(CPM)", expand=c(0,0)) +
        scale_y_continuous(name="log2(FC)", expand=c(0,0)) +
        coord_fixed(0.5)
    p
}

p1 <- seq_len(floor(length(cn.higher.samples) / 2)) %>%
    lapply(function(i) {
        s1 <- cn.higher.samples[i]
        s2 <- cn.higher.samples[length(cn.higher.samples) - i + 1]
        title <- sprintf("MA Plot of 150bp Windows for \n %s vs %s",
                         colnames(dge)[s1], colnames(dge)[s2])
        tsmsg("Making ", title)
        withGC(doMAPlot(logcpm, s1, s2) +
               ggtitle(title) +
               theme(plot.title = element_text(hjust = 0.5)))
    })
p2 <- seq_len(floor(length(cn.higher.samples) / 2)) %>%
    lapply(function(i) {
        s1 <- cn.higher.samples[i]
        s2 <- cn.higher.samples[length(cn.higher.samples) - i + 1]
        title <- sprintf("10KB Bin MA Plot for \n %s vs %s",
                         colnames(dge)[s1], colnames(dge)[s2])
        tsmsg("Making ", title)
        withGC(doMAPlot(bigbin.logcpm, s1, s2) +
               ggtitle(title) +
               theme(plot.title = element_text(hjust = 0.5)))
    })
p3 <- seq_len(floor(length(cn.higher.samples) / 2)) %>%
    lapply(function(i) {
        s1 <- cn.higher.samples[i]
        s2 <- cn.higher.samples[length(cn.higher.samples) - i + 1]
        title <- sprintf("MA Plot of 150bp Windows Overlapping Peaks for \n %s vs %s",
                         colnames(dge)[s1], colnames(dge)[s2])
        tsmsg("Making ", title)
        withGC(doMAPlot(peak.logcpm, s1, s2) +
               ggtitle(title) +
               theme(plot.title = element_text(hjust = 0.5)))
    })
p4 <- seq_len(floor(length(cn.higher.samples) / 2)) %>%
    lapply(function(i) {
        s1 <- cn.higher.samples[i]
        s2 <- cn.higher.samples[length(cn.higher.samples) - i + 1]
        linedata <- getOffsetLineData(s1, s2)
        title <- sprintf("Loess-Normalized MA Plot of 150bp Windows Overlapping Peaks for \n %s vs %s",
                         colnames(dge)[s1], colnames(dge)[s2])
        tsmsg("Making ", title)
        withGC(doMAPlot(peak.logcpm.loess, s1, s2, linedata=NULL) +
               ggtitle(title) +
               theme(plot.title = element_text(hjust = 0.5)))
    })

{
    tsmsg("Printing MA plots")
    pdf(sprintf("plots/csaw/%s Selected Sample MA Plots.pdf", chip), width=10, height=10)
    withGC(print(p1))
    dev.off()
    pdf(sprintf("plots/csaw/%s Selected Sample 10KB Bin MA Plots.pdf", chip), width=10, height=10)
    withGC(print(p2))
    dev.off()
    pdf(sprintf("plots/csaw/%s Selected Sample Peak-Overlap MA Plots.pdf", chip), width=10, height=10)
    withGC(print(p3))
    dev.off()
    pdf(sprintf("plots/csaw/%s Selected Sample Peak-Overlap Normalized MA Plots.pdf", chip), width=10, height=10)
    withGC(print(p4))
    dev.off()
}

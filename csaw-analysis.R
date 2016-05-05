
library(stringr)
library(magrittr)
library(openxlsx)
library(doParallel)
options(mc.cores=parallel::detectCores())
registerDoParallel(cores=parallel::detectCores())
library(SummarizedExperiment)
library(readr)
library(dplyr)
library(edgeR)
library(limma)
library(csaw)
library(ggplot2)
library(scales)
library(rtracklayer)
library(RColorBrewer)
library(reshape2)
library(parallel)

## Ensure output directory exists
dir.create("results/csaw", FALSE, TRUE)

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

suppressPlot <- function(arg) {
    dev.new()
    result <- arg
    dev.off()
    result
}

add.numbered.colnames <- function(x, prefix="C") {
    x %>% set_colnames(sprintf("%s%i", prefix, seq(from=1, length.out=ncol(x))))
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

## Like limma's removeBatchEffect, but for DGEList. Modifies the
## offsets instead of the data.
offsetBatchEffect <- function (dge, batch = NULL, batch2 = NULL, covariates = NULL,
                               design = matrix(1, ncol(dge), 1), ...) {
    if (is.null(batch) && is.null(batch2) && is.null(covariates))
        return(dge)
    if (!is.null(batch)) {
        batch <- as.factor(batch)
        contrasts(batch) <- contr.sum(levels(batch))
        batch <- model.matrix(~batch)[, -1, drop = FALSE]
    }
    if (!is.null(batch2)) {
        batch2 <- as.factor(batch2)
        contrasts(batch2) <- contr.sum(levels(batch2))
        batch2 <- model.matrix(~batch2)[, -1, drop = FALSE]
    }
    if (!is.null(covariates))
        covariates <- as.matrix(covariates)
    X.batch <- cbind(batch, batch2, covariates)
    fit <- glmFit(dge, cbind(design, X.batch), ..., priot.count = 0)
    beta <- fit$coefficients[, -(1:ncol(design)), drop = FALSE]
    beta[is.na(beta)] <- 0
    dge$offset <- fit$offset + beta %*% t(X.batch)
    dge
}

## Version of voom that uses an offset matrix instead of lib sizes
voomWithOffset <- function (dge, design = NULL, offset=expandAsMatrix(getOffset(dge), dim(dge)),
                            normalize.method = "none", plot = FALSE, span = 0.5, ...)
{
    out <- list()
    out$genes <- dge$genes
    out$targets <- dge$samples
    if (is.null(design) && diff(range(as.numeric(counts$sample$group))) >
        0)
        design <- model.matrix(~group, data = counts$samples)
    counts <- dge$counts
    if (is.null(design)) {
        design <- matrix(1, ncol(counts), 1)
        rownames(design) <- colnames(counts)
        colnames(design) <- "GrandMean"
    }

    effective.lib.size <- exp(offset)

    y <- log2((counts + 0.5)/(effective.lib.size + 1) * 1e+06)
    y <- normalizeBetweenArrays(y, method = normalize.method)
    fit <- lmFit(y, design, ...)
    if (is.null(fit$Amean))
        fit$Amean <- rowMeans(y, na.rm = TRUE)
    sx <- fit$Amean + mean(log2(effective.lib.size + 1)) - log2(1e+06)
    sy <- sqrt(fit$sigma)
    allzero <- rowSums(counts) == 0
    if (any(allzero)) {
        sx <- sx[!allzero]
        sy <- sy[!allzero]
    }
    l <- lowess(sx, sy, f = span)
    if (plot) {
        plot(sx, sy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )",
            pch = 16, cex = 0.25)
        title("voom: Mean-variance trend")
        lines(l, col = "red")
    }
    f <- approxfun(l, rule = 2)
    if (fit$rank < ncol(design)) {
        j <- fit$pivot[1:fit$rank]
        fitted.values <- fit$coef[, j, drop = FALSE] %*% t(fit$design[,
            j, drop = FALSE])
    }
    else {
        fitted.values <- fit$coef %*% t(fit$design)
    }
    fitted.cpm <- 2^fitted.values
    ## fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
    fitted.count <- 1e-06 * fitted.cpm * (effective.lib.size + 1)
    fitted.logcount <- log2(fitted.count)
    w <- 1/f(fitted.logcount)^4
    dim(w) <- dim(fitted.logcount)
    out$E <- y
    out$weights <- w
    out$design <- design
    out$effective.lib.size <- effective.lib.size
    if (is.null(out$targets))
        out$targets <- data.frame(lib.size = exp(colMeans(offset)))
    else out$targets$lib.size <- exp(colMeans(offset))
    new("EList", out)
}

## Version of voom that uses an offset matrix instead of lib sizes
voomWithQualityWeightsAndOffset <-
    function (dge, design = NULL,
              offset=expandAsMatrix(getOffset(dge), dim(dge)),
              normalize.method = "none",
              plot = FALSE, span = 0.5, var.design = NULL, method = "genebygene",
              maxiter = 50, tol = 1e-10, trace = FALSE, replace.weights = TRUE,
              col = NULL, ...)
{
    counts <- dge$counts
    if (plot) {
        oldpar <- par(mfrow = c(1, 2))
        on.exit(par(oldpar))
    }
    v <- voomWithOffset(dge, design = design, offset = offset, normalize.method = normalize.method,
        plot = FALSE, span = span, ...)
    aw <- arrayWeights(v, design = design, method = method, maxiter = maxiter,
        tol = tol, var.design = var.design)
    v <- voomWithOffset(dge, design = design, weights = aw, offset = offset,
        normalize.method = normalize.method, plot = plot, span = span,
        ...)
    aw <- arrayWeights(v, design = design, method = method, maxiter = maxiter,
        tol = tol, trace = trace, var.design = var.design)
    wts <- asMatrixWeights(aw, dim(v)) * v$weights
    attr(wts, "arrayweights") <- NULL
    if (plot) {
        barplot(aw, names = 1:length(aw), main = "Sample-specific weights",
            ylab = "Weight", xlab = "Sample", col = col)
        abline(h = 1, col = 2, lty = 2)
    }
    if (replace.weights) {
        v$weights <- wts
        v$sample.weights <- aw
        return(v)
    }
    else {
        return(wts)
    }
}

estimateDispByGroup <- function(dge, group=as.factor(dge$samples$group), batch, ...) {
    stopifnot(nlevels(group) > 1)
    stopifnot(length(group) == ncol(dge))
    if (!is.list(batch)) {
        batch <- list(batch=batch)
    }
    batch <- as.data.frame(batch)
    stopifnot(nrow(batch) == ncol(dge))
    colnames(batch) %>% make.names(unique=TRUE)
    igroup <- seq_len(ncol(dge)) %>% split(group)
    lapply(igroup, function(i) {
        group.dge <- dge[,i]
        group.batch <- droplevels(batch[i,, drop=FALSE])
        group.batch <- group.batch[sapply(group.batch, . %>% unique %>% length %>% is_greater_than(1))]
        group.vars <- names(group.batch)
        if (length(group.vars) == 0)
            group.vars <- "1"
        group.batch.formula <- as.formula(str_c("~", str_c(group.vars, collapse="+")))
        des <- model.matrix(group.batch.formula, group.batch)
        estimateDisp(group.dge, des, ...)
    })
}

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

withGC <- function(expr) {
    on.exit(gc())
    return(expr)
}

## Eval expression in forked process in the foreground (not in
## parallel). This can be useful if the expression makes problematic
## irreversible changes to the environment, or causes memory usage to
## balloon irreversibly.
in.forked.process <- function(expr) {
  result <- mccollect(mcparallel(expr))
  if (length(result) == 0) {
    NULL
  } else {
    result <- result[[1]]
    if (is(result, "try-error")) {
      stop(result)
    } else {
      result
    }
  }
}

tsmsg("Reading ChIP input data")
input.bigbin.counts <- readRDS(sprintf("saved_data/bigbin-counts-%s.RDS", "input"))
input.window.counts <- readRDS(sprintf("saved_data/window-counts-%s.RDS", "input"))

comp.inputnorm <- normOffsets(input.bigbin.counts, type="scaling")

chips <- c("H3K4me3", "H3K4me2", "H3K27me3")

## Read peak files
peaks <- lapply(
    setNames(nm=chips),
    . %>% sprintf(fmt="data_files/ChIP-Seq/%s_peaks_IDR_filtered.bed") %>%
    import.narrowPeak)

peaks %>% sapply(. %>% width %>% sum %>% divide_by(73) %>% round)

for (chip in chips) { in.forked.process({
    tsmsg("Reading data for ", chip)
    bigbin.counts <- readRDS(sprintf("saved_data/bigbin-counts-%s.RDS", chip))
    window.counts <- readRDS(sprintf("saved_data/window-counts-%s.RDS", chip))
    chip.peaks <- peaks[[chip]]
    ## Filter to windows with at least 10 counts total
    keep <- rowSums(assay(window.counts)) >= 10
    window.counts %<>% .[keep,]
    invisible(gc())

    tsmsg("Counting number of non-zero bins")
    colData(window.counts)$NonZeroBins <- withGC(colSums(assay(window.counts) > 0))


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
        pdf(sprintf("results/csaw/%s-window-abundance-vs-peaks.pdf", chip))
        print(p)
        dev.off()
    }

    ## Based on above plot, we decide to use the peaks as filter
    ## criterion
    peak.overlap <- overlapsAny(window.counts, chip.peaks)
    peak.window.counts <- window.counts[peak.overlap,]
    tsmsg("Computing efficiency normalization from peak windows")
    pnf <- normOffsets(peak.window.counts, type="scaling")
    colData(window.counts)$
    PeakNormFactors <- pnf
    colData(peak.window.counts)$PeakNormFactors <- pnf

    ## Create DGEList and populate $genes and $samples
    dge <- asDGEList(window.counts)
    ## Saves memory
    rownames(dge$counts) <- NULL
    gc()
    dge$genes <- rowRanges(window.counts) %>% as.data.frame %>% lapply(Rle) %>% DataFrame
    dge$samples %<>% cbind(colData(window.counts)) %>% as.data.frame
    dge$samples$nf.logratio <- dge$samples %$% log2(PeakNormFactors / CompNormFactors)

    {
        tsmsg("Testing norm factors for association with experimental factors")
        xvars <- c("Celltype", "Day", "Donor", "TreatmentGroup")
        yvars <- c("lib.size", "CompNormFactors", "HANormFactors", "PeakNormFactors", "nf.logratio")

        formulas <- list()
        for (xv in xvars) {
            for (yv in yvars) {
                modname <- sprintf("%s.vs.%s", yv, xv)
                formulas[[modname]] <- as.formula(sprintf("%s ~ %s", yv, xv))
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
        nf.vs.ls <- dge$samples %>% select(LibSize=lib.size, HiAb=HANormFactors, Peak=PeakNormFactors, Comp=CompNormFactors) %>%
            melt(id.vars="LibSize", variable.name = "NormType", value.name="NormFactor")
        p <- list(ggplot(nf.vs.ls) +
                  aes(x=LibSize, y=NormFactor, group=NormType, color=NormType, fill=NormType) +
                  geom_point() + geom_smooth(alpha=0.15) +
                  scale_x_continuous(trans=log_trans(2)) +
                  scale_y_continuous(trans=log_trans(2)) +
                  ggtitle("Norm Factors vs Library Sizes"))
        ## Plot lib sizes & norm factors vs experimental factors
        ls.nf.vs.exp <- dge$samples %>%
            select(
                Sample, Celltype, Day, Donor,
                Group=TreatmentGroup, LibSize=lib.size,
                HighAbundanceNormFactor=HANormFactors,
                PeakNormFactor=PeakNormFactors,
                CompositionNormFactor=CompNormFactors) %>%
            melt(id.vars=c("Sample", "Celltype", "Day", "Donor",
                           "Group"),
                 variable.name="Factor", value.name="Value")
        p0 <- ggplot(ls.nf.vs.exp) +
            aes(x=Value) + scale_x_continuous(trans=log_trans(2)) +
            geom_point(position=position_jitter(height=0.25)) +
            facet_wrap(~Factor, scales="free_x")
        vars.to.plot <- c("Celltype", "Day", "Donor", "Group")
        for (i in vars.to.plot) {
            p <- c(p, list(p0 + aes_string(y=i) +
                           ggtitle(sprintf("Library Sizes and Norm Factors vs %s", i))))
        }
        pdf(sprintf("results/csaw/%s-normfactors.pdf", chip))
        print(p)
        dev.off()
    }

    ## Select the most and least extreme samples for plotting MA plots

    middle.samples <- dge$samples$nf.logratio %>% abs %>% order
    cn.higher.samples <- dge$samples$nf.logratio %>% order
    pn.higher.samples <- rev(cn.higher.samples)

    tsmsg("Computing logCPM")
    logcpm <- cpm(dge, log=TRUE)

    doMAPlot <- function(s1, s2) {
        pointdata <- data.frame(S1=logcpm[,s1], S2=logcpm[,s2]) %>%
            transmute(A=(S1+S2)/2, M=S2-S1) %>%
            filter(A >= -2)
        bw <- sapply(pointdata[c("A", "M")], . %>% (MASS::bandwidth.nrd))
        dataranges <- sapply(pointdata[c("A", "M")], . %>% range %>% diff)

        ## Use the same perceptual bandwidth in both dimensions
        ## (yielding circular blobs)
        bw <- min(bw / c(1,2)) * c(1,2)

        linedata <- c(Comp="CompNormFactors",
                      Peaks="PeakNormFactors",
                      HiAb="HANormFactors") %>%
            sapply(. %>% dge$samples[[.]] %>% log2 %>% {.[s2] - .[s1]}) %>%
            data.frame(NormFactor=., NormType=names(.))

        ggplot(pointdata) +
            aes(x=A, y=M) +
            stat_density2d(contour=FALSE, geom="raster", n=512, interpolate=TRUE,
                           aes(fill=..density.., alpha=..density..)) +
            scale_fill_gradientn(colors=suppressWarnings(brewer.pal(Inf, "Blues")),
                                 trans=power_trans(1/8),
                                 name="Density") +
            scale_alpha_continuous(trans=power_trans(1/20), guide=FALSE) +
            geom_hline(data=linedata, aes(yintercept=NormFactor, color=NormType)) +
            scale_color_discrete(name="Norm Type") +
            scale_x_continuous(limits=pointdata$A %>% range %>% expand_range(add=bw[1] * 3),
                               expand=c(0,0)) +
            scale_y_continuous(limits=pointdata$M %>% range %>% expand_range(add=bw[2] * 3),
                               expand=c(0,0)) +
            coord_fixed(0.5)
    }
    p <- seq_len(floor(length(cn.higher.samples) / 2)) %>%
        lapply(function(i) {
            s1 <- cn.higher.samples[i]
            s2 <- cn.higher.samples[length(cn.higher.samples) - i + 1]
            title <- sprintf("MA Plot for %s vs %s",
                             colnames(dge)[s1], colnames(dge)[s2])
            tsmsg("Making ", title)
            withGC(doMAPlot(s1, s2) +
                   ggtitle(title))
        })
    {
        tsmsg("Printing MA plots")
        pdf(sprintf("results/csaw/%s Selected Sample MA Plots.pdf", chip), width=10, height=10)
        withGC(print(p))
        dev.off()
    }
    tsmsg("Saving image")
    save.image(sprintf("csaw-%s.rda", chip))
    NULL
})}

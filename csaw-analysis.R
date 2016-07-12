#!/usr/bin/env Rscript

library(stringr)
library(magrittr)
library(openxlsx)
library(doParallel)
options(mc.cores=parallel::detectCores())
registerDoParallel(cores=parallel::detectCores())
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

## Ensure output directory exists
dir.create("results/csaw", FALSE, TRUE)

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

suppressPlot <- function(arg) {
    png("/dev/null")
    on.exit(dev.off())
    result <- arg
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

chips <- c("H3K4me3", "H3K4me2", "H3K27me3")

ave.count.threshold <- 5

for (chip in chips) { in.forked.process({
    tsmsg("Reading data for ", chip)
    ## Let's avoid re-doing all the normalization and filtering by
    ## loading the saved data from csaw-qc.R
    dge <- readRDS(sprintf("saved_data/csaw-DGEList-%s.RDS", chip))
    dge <- dge[,dge$samples %$% order(ChIP, Celltype, Day, Donor, Sample)]
    chip.peaks <- import.narrowPeak(
        sprintf(fmt="data_files/ChIP-Seq/%s_peaks_IDR_filtered.bed", chip))

    tsmsg("Filtering")
    ab <- dge$genes$Abundance
    thresh <- aveLogCPM(5, lib.size=mean(dge$samples$totals))
    present <- ab >= thresh
    dge <- dge[present,]

    ## Use efficiency normalization
    dge$samples %<>% mutate(norm.factors=PeakNormFactors)

    design <- model.matrix(~0 + TreatmentGroup + Donor, dge$samples)
    colnames(design) %<>%
        str_replace("^TreatmentGroup", "") %>%
        str_replace("^DonorDn", "Donor")

    tsmsg("Estimating dispersions")
    dge %<>% estimateDisp(design, robust=TRUE)
    fit <- glmQLFit(dge.comp, design, robust=TRUE)

    tsmsg("Saving image")
    save.image(sprintf("saved_data/csaw-norm-eval-%s.rda", chip))
    NULL
})}

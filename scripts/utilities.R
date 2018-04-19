## TODO: Convert to custom utils package, and make all my scripts use it

suppressMessages({
    library(stringr)
    library(rex)
    library(glue)
    library(sitools)
    library(glue)
    library(sitools)
    library(rex)
    library(magrittr)
    library(dplyr)
    library(readr)
    library(tidyr)
    library(forcats)
    library(rlang)
    library(assertthat)
    library(BiocParallel)
    library(lazyeval)
    library(future)
    library(Rtsne)
    library(qvalue)
    library(fdrtool)
    library(edgeR)
    library(csaw)
    library(rtracklayer)
})

# Inverse of sitools::f2si
si2f <- function(string, unit="") {
    if (length(string) == 0) {
        return(numeric(0))
    }
    sifactor <- c(1e-24, 1e-21, 1e-18, 1e-15, 1e-12, 1e-09, 1e-06,
        0.001, 1, 1000, 1e+06, 1e+09, 1e+12, 1e+15, 1e+18, 1e+21,
        1e+24)
    pre <- c("y", "z", "a", "f", "p", "n", "u", "m",
        "", "k", "M", "G", "T", "P", "E", "Z", "Y")

    rx <- rex(
        # Leading whitespace
        start,
        zero_or_more(space),

        # Capture a floating point number
        capture(
            # Sign
            maybe(one_of("+", "-")),
            # Integer part
            zero_or_more(digit),
            # Decimal point
            maybe("."),
            # Fractional part (or integer part when decimal is not
            # present)
            one_or_more(digit),
            # Exponential notation
            maybe(
                one_of("e", "E"),
                maybe(one_of("+", "-")),
                one_or_more(digit)
            )
        ),

        # Space between number and unit
        zero_or_more(space),

        # Capture SI prefix
        capture(maybe(one_of(pre))),
        # User-specified unit suffix
        unit,

        # Trailing whitespace
        zero_or_more(space),
        end
    )

    m <- str_match(string, rx)
    base <- as.numeric(m[,2])
    p <- m[,3]
    fac <- sifactor[match(p, pre)]
    base * fac
}

# Convert e.g. "2.5kbp" to 2500
parse.bp <- function(size) {
    suppressWarnings({
        result <- si2f(size, unit="bp")
        # Fall back to just parsing a number without the "bp" suffix
        result[is.na(result)] <- si2f(size[is.na(result)])
    })
    assert_that(!any(is.na(result)))
    result
}

# Convert e.g. 2500 to "2.5kbp"
format.bp <- function(x) {
    x %>% round %>% f2si(unit="bp") %>% str_replace_all(rex(one_or_more(space)), "")
}

withGC <- function(expr) {
    on.exit(gc())
    return(expr)
}

# Use to assign to complex sub-expressions in the middle of a dplyr pipeline.
# For example, you can't easily do the following in the middle of a pipeline:
# "assays(x[[1]])$counts[3,5] <- 45". But now you can do it like: "x %>%
# assign_into(assays(.[[1]])$counts[3,5], 45) %>% another_fun() %>% ..."
assign_into <- function(x, expr, value) {
    expr <- lazyeval::lazy(expr)$expr
    f_eval(f_interp(~ x %>% { uq(expr) <- uq(value); . }))
}

# Useful to wrap functions that both produce a plot and return a useful value,
# when you only want the return value and not the plot.
suppressPlot <- function(arg) {
    png("/dev/null")
    result <- arg
    dev.off()
    result
}

# Always returns a list of ggplot objects. Flattens nested lists,
# encapsulates single plots into a 1-element list, ensures that all
# elements are ggplots.
get.ggplots <- function(plots) {
    UseMethod("get.ggplots")
}

get.ggplots.default <- function(plots) {
    stop(glue("Don't know how to get ggplots from an object of class {deparse(class(plots)[1])}"))
}

get.ggplots.gg <- function(plots) {
    list(plots)
}

get.ggplots.list <- function(plots) {
    plotlists <- lapply(plots, get.ggplots)
    do.call(c, plotlists)
}

## Returns TRUE if x refers to the device number of a currently active
## graphics device.
is_dev <- function(x) {
    is_scalar_integer(x) && x %in% dev.list()
}

with_dev <- function(dev, code, closedev) {
    orig.device <- dev.cur()
    new.device <- force(dev)
    # Functions that create devices don't generally return them, they
    # just set them as the new current device, so get the actual
    # device from dev.cur() instead.
    if (is.null(new.device)) {
        new.device <- dev.cur()
    }
    assert_that(is_dev(new.device) || new.device == 1)
    if (missing(closedev)) {
         closedev <- new.device != orig.device
    }
    on.exit({
        if (closedev) {
            dev.off(new.device)
        }
        if (is_dev(orig.device)) {
            dev.set(orig.device)
        }
    })
    force(code)
}

ggprint <- function(plots, device=dev.cur(), closedev, printfun=print) {
    p <- get.ggplots(plots)
    with_dev(device, lapply(p, printfun), closedev)
    invisible(p)
}

# Printer function for ggplotly, to be passed as the prinfun for ggprint
ggplotly.printer <- function(...) {
    dots <- list(...)
    function(p) {
        args <- c(list(p=p), dots)
        print(do.call(plotly::ggplotly, args))
    }
}

rasterpdf <- function(pdffile, outfile=pdffile, resolution=600) {
    tempf <- tempfile(pattern="raster", fileext=".pdf")
    on.exit(unlink(tempf))
    exitcode <- system2("convert", args=c("-density", resolution, pdffile, tempf),
        stdout=FALSE, stderr=FALSE)
    assert_that(exitcode == 0)
    assert_that(file.exists(tempf))
    suppressWarnings(file.rename(tempf, outfile))
    # If file still exists, then the rename failed because it's a
    # cross-device move, so copy and delete instead.
    if (file.exists(tempf)) {
        file.copy(tempf, outfile)
    }
}

add.numbered.colnames <- function(x, prefix="C") {
    x %>% set_colnames(glue("{prefix}{num}", num=seq(from=1, length.out=ncol(x))))
}

# For each column of a data frame, if it is a character vector with at least one
# repeated value, convert that column to a factor
autoFactorize <- function(df) {
    for (i in colnames(df)) {
        if (is.character(df[[i]]) && anyDuplicated(df[[i]])) {
            df[[i]] %<>% factor
        }
    }
    df
}

library(GGally)
# Variant of ggduo that accepts dataX and dataY as separate data frames
ggduo.dataXY <- function(dataX, dataY, extraData=NULL, ...) {
    assert_that(ncol(dataX) > 0)
    assert_that(ncol(dataY) > 0)
    alldata <- cbind(dataX, dataY)
    if (!is.null(extraData)) {
        alldata <- cbind(alldata, extraData)
    }
    ggduo(alldata, columnsX=colnames(dataX), columnsY=colnames(dataY), ...)
}

# Make cairo_pdf use onefile=TRUE by default
cairo_pdf <- function(..., onefile=TRUE) {
    grDevices::cairo_pdf(..., onefile = onefile)
}

library(edgeR)

# Split dge into groups and estimate dispersion for each group. Returns a list
# of DGELists.
estimateDispByGroup <- function(dge, group=as.factor(dge$samples$group), batch, ...) {
    assert_that(nlevels(group) > 1)
    assert_that(length(group) == ncol(dge))
    if (!is.list(batch)) {
        batch <- list(batch=batch)
    }
    batch <- as.data.frame(batch)
    assert_that(nrow(batch) == ncol(dge))
    colnames(batch) %>% make.names(unique=TRUE)
    igroup <- seq_len(ncol(dge)) %>% split(group)
    bplapply(igroup, function(i) {
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

# Versions of cpm and aveLogCPM that use an offset matrix instead of
# lib sizes
cpmWithOffset <- function(dge, offset=getOffset(dge),
                         log = FALSE, prior.count = 0.25, preserve.mean=TRUE, ...) {
    if (preserve.mean) {
        dge <- scaleOffset(dge, offset)
    } else {
        dge$offset <- offset
    }
    cpm(dge$counts, lib.size=exp(getOffset(dge)), log=log, prior.count=prior.count, ...)
}

aveLogCPMWithOffset <- function(y, ...) {
    UseMethod("aveLogCPMWithOffset")
}

aveLogCPMWithOffset.default <- function (y, offset = NULL, prior.count = 2,
                                         dispersion = NULL, weights = NULL, ...) {
    aveLogCPM(y, lib.size = NULL, offset = offset, prior.count = prior.count,
              dispersion = dispersion, weights = weights, ...)
}

aveLogCPMWithOffset.DGEList <- function (y, offset = expandAsMatrix(getOffset(y), dim(y)),
                                         prior.count = 2, dispersion = NULL, ...) {
    if (is.null(dispersion)) {
        dispersion <- y$common.dispersion
    }
    offsetMat <- offset
    aveLogCPMWithOffset(
        y$counts, offset = offsetMat, prior.count = prior.count,
        dispersion = dispersion, weights = y$weights)
}

library(limma)

# Version of voom that uses an offset matrix instead of lib sizes
voomWithOffset <-
    function (dge, design = NULL, offset=expandAsMatrix(getOffset(dge), dim(dge)),
              normalize.method = "none", plot = FALSE, span = 0.5, ...)
{
    out <- list()
    out$genes <- dge$genes
    out$targets <- dge$samples
    if (is.null(design) && diff(range(as.numeric(counts$sample$group))) > 0)
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
        fitted.values <- fit$coef[, j, drop = FALSE] %*%
            t(fit$design[, j, drop = FALSE])
    } else {
        fitted.values <- fit$coef %*% t(fit$design)
    }
    fitted.cpm <- 2^fitted.values
    # fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
    fitted.count <- 1e-06 * fitted.cpm * (effective.lib.size + 1)
    fitted.logcount <- log2(fitted.count)
    w <- 1/f(fitted.logcount)^4
    dim(w) <- dim(fitted.logcount)
    out$E <- y
    out$weights <- w
    out$design <- design
    out$effective.lib.size <- effective.lib.size
    if (is.null(out$targets)) {
        out$targets <- data.frame(lib.size = exp(colMeans(offset)))
    } else {
        out$targets$lib.size <- exp(colMeans(offset))
    }
    new("EList", out)
}

# Version of voom that uses an offset matrix instead of lib sizes
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
    } else {
        return(wts)
    }
}

# Convenience function for alternating dupCor and voom until convergence
voomWithDuplicateCorrelation <-
    function(counts, design = NULL, plot = FALSE, block = NULL, trim = 0.15,
             voom.fun=voom, dupCor.fun=duplicateCorrelation, initial.correlation=0,
             niter=5, tol=1e-6, verbose=TRUE, ...)
{
    assert_that(niter >= 1)
    assert_that(is.finite(niter))

    if (niter < 2) {
        warning("Using less than 2 iterations of voom and duplicateCorrelation is not recommended.")
    }
    # Do first iteration
    prev.cor <- 0
    iter.num <- 0
    if (initial.correlation != 0 && !is.null(block)) {
        elist <- voom.fun(counts, design=design, plot=FALSE, block=block, correlation=initial.correlation, ...)
    } else {
        elist <- voom.fun(counts, design=design, plot=FALSE, ...)
    }
    if (is.null(block)) {
        warning("Not running duplicateCorrelation because block is NULL.")
    }
    if (verbose) {
        message(glue("Initial guess for duplicate correlation before 1st iteration: {initial.correlation}"))
    }
    dupcor <- dupCor.fun(elist, design, block=block, trim=trim)
    iter.num <- iter.num + 1
    if (verbose) {
        message(glue("Duplicate correlation after {toOrdinal(iter.num)} iteration: {dupcor$consensus.correlation}"))
    }
    while (iter.num < niter) {
        if (!is.null(tol) && is.finite(tol)) {
            if (abs(dupcor$consensus.correlation - prev.cor) <= tol) {
                if (verbose) {
                    message(glue("Stopping after {toOrdinal(iter.num)} iteration because tolerance threshold was reached."))
                }
                break
            }
        }
        prev.cor <- dupcor$consensus.correlation
        elist <- voom.fun(counts, design=design, plot=FALSE, block=block, correlation=prev.cor, ...)
        dupcor <- dupCor.fun(elist, design, block=block, trim=trim)
        iter.num <- iter.num + 1
        if (verbose) {
            message(glue("Duplicate correlation after toOrdinal(iter.num) iteration: dupcor$consensus.correlation"))
        }
    }
    elist <- voom.fun(counts, design=design, plot=plot, block=block, correlation=dupcor$consensus.correlation, ...)
    for (i in names(dupcor)) {
        elist[[i]] <- dupcor[[i]]
    }
    return(elist)
}

library(codingMatrices)
# Variant of code_control that generates more verbose column names
code_control_named <- function (n, contrasts = TRUE, sparse = FALSE) {
    if (is.numeric(n) && length(n) == 1L) {
        if (n > 1L)
            levels <- .zf(seq_len(n))
        else stop("not enough degrees of freedom to define contrasts")
    }
    else {
        levels <- as.character(n)
        n <- length(n)
    }
    B <- diag(n)
    dimnames(B) <- list(1:n, levels)
    if (!contrasts) {
        if (sparse)
            B <- .asSparse(B)
        return(B)
    }
    B <- B - 1/n
    B <- B[, -1, drop=FALSE]
    colnames(B) <- paste(levels[1], levels[-1], sep = ".vs.")
    if (sparse) {
        .asSparse(B)
    }
    else {
        B
    }
}
environment(code_control_named) <- new.env(parent = environment(code_control))

# If your factor levels are globally unique, then it's safe to strip the factor
# name prefixes from column names of the design matrix. This function does that.
strip.design.names <- function(design, prefixes=names(attr(design, "contrasts"))) {
    if (!is.null(prefixes)) {
        regex.to.delete <- rex(or(prefixes) %if_prev_is% or(start, ":") %if_next_isnt% end)
        colnames(design) %<>% str_replace_all(regex.to.delete, "")
        if (anyDuplicated(colnames(design))) {
            warning("Stripping design names resulted in non-unique names")
        }
    }
    design
}

# Variant of model.matrix with an additional option to strip the prefixes as
# described above.
model.matrix <- function(..., strip.prefixes=FALSE) {
    design <- stats::model.matrix(...)
    if (strip.prefixes) {
        design <- strip.design.names(design)
    }
    return(design)
}

# Helper function for ensureAtomicColumns
collapseToAtomic <- function(x, sep=",") {
    if (is.atomic(x)) {
        return(x)
    } else {
        y <- lapply(x, str_c, collapse=sep)
        y[lengths(y) == 0] <- NA
        assert_that(all(lengths(y) == 1))
        y <- unlist(y)
        assert_that(length(y) == length(x))
        return(y)
    }
}

# Function to ensure that all columns of a data frame are atomic vectors.
# Columns that are lists have each of their elements collapsed into strings
# using the sepcified separator.
ensureAtomicColumns <- function(df, sep=",") {
    df[] %<>% lapply(collapseToAtomic, sep=sep)
    df
}

library(ggplot2)
# Function to create an annotatied p-value histogram
plotpvals <- function(pvals, ptn=propTrueNull(pvals)) {
    df <- data.frame(p=pvals)
    linedf <- data.frame(y=c(1, ptn), Line=c("Uniform", "Est. Null") %>% factor(levels=unique(.)))
    ggplot(df) + aes(x=p) +
        geom_histogram(aes(y = ..density..), binwidth=0.01, boundary=0) +
        geom_hline(aes(yintercept=y, color=Line),
            data=linedf, alpha=0.5, show.legend=TRUE) +
        scale_color_manual(name="Ref. Line", values=c("blue", "red")) +
        xlim(0,1) + ggtitle(glue("P-value distribution (Est. {format(100 * (1-ptn), digits=3)}% signif.)")) +
        expand_limits(y=c(0, 1.25)) +
        xlab("p-value") + ylab("Relative frequency") +
        theme(legend.position=c(0.95, 0.95),
            legend.justification=c(1,1))
}

# Misc functions to facilitate alternative FDR calculations for limma/edgeR
# results

# Variant of eBayes that uses propTrueNull (or a custom function) to set the
# proportion argument
eBayes_autoprop <- function(..., prop.method="lfdr") {
    eb <- eBayes(...)
    if (is.function(prop.method)) {
        ptn <- prop.method(eb$p.value)
    } else {
        ptn <- propTrueNull(eb$p.value, method=prop.method)
    }
    eBayes(..., proportion=1-ptn)
}

# Compute posterior probabilities and Bayesian FDR values from limma's B
# statistics
bfdr <- function(B) {
    o <- order(B, decreasing = TRUE)
    ro <- order(o)
    B <- B[o]
    positive <- which(B > 0)
    PP <- exp(B)/(1+exp(B))
    # Computing from 1-PP gives better numerical precision for the
    # most significant genes (large B-values)
    oneMinusPP <- 1/(1+exp(B))
    BayesFDR <- cummean(oneMinusPP)
    data.frame(B, PP, BayesFDR)[ro,]
}

# Add Bayesian FDR to a limma top table
add.bfdr <- function(ttab) {
    B <- ttab[["B"]]
    if (is.null(B)) {
        warning("Cannot add BFDR to table with no B statistics")
        return(ttab)
    }
    btab <- bfdr(B)[c("PP", "BayesFDR")]
    for (i in names(btab)) {
        ttab[[i]] <- btab[[i]]
    }
    ttab
}

# limma uses "P.Value", edgeR uses "PValue", so we need an abstraction
get.pval.colname <- function(ttab) {
    if (is.character(ttab)) {
        cnames <- ttab
    } else {
        cnames <- colnames(ttab)
    }
    pcol <- match(c("p.value", "pvalue", "pval", "p"),
        tolower(cnames)) %>%
        na.omit %>% .[1]
    pcolname <- cnames[pcol]
    if (length(pcolname) != 1)
        stop("Could not determine p-value column name")
    return(pcolname)
}

# Add a q-value column to any table with a p-value column
add.qvalue <- function(ttab, ...) {
    tryCatch({
        P <- ttab[[get.pval.colname(ttab)]]
        qobj <- qvalue(P, ...)
        ttab$QValue <- qobj$qvalues
        ttab$LocFDR <- qobj$lfdr
        attr(ttab, "qvalue") <- qobj
    }, error=function(e) {
        warning(str_c("Failed to compute q-values: ", e$message))
    })
    ttab
}

# Variant that does not restrict to the range [0,1]. Useful for identifying
# potential atypical p-value distributions (too many large p-values).
propTrueNullByLocalFDR.unrestricted <- function (p) {
    n <- length(p)
    i <- n:1L
    p <- sort(p, decreasing = TRUE)
    q <- n/i * p
    n1 <- n + 1L
    sum(i * q)/n/n1 * 2
}

# Add corrected p-value, q-value, and locfdr columns to any table with a p-value
# column. Works by the possibly questionable method of converting p-values to
# equivalent z-scores and running fdrtool on the z-scores.
add.fdrtool <- function(ttab, verbose=FALSE, plot=TRUE, convert.to.zscores, cutoff.method, ...) {
    P <- ttab[[get.pval.colname(ttab)]]
    assert_that(!is.null(P), length(P) > 0, !any(is.na(P)))
    if (missing(convert.to.zscores)) {
        # If pval dist is high-biased, use normal modelling instead.
        ptn <- propTrueNullByLocalFDR.unrestricted(P)
        convert.to.zscores <- ptn > 1
    }
    if (convert.to.zscores) {
        # Try normal modelling instead
        Zscore <- qnorm(1-(P/2))
        if (missing(cutoff.method)) {
            co <- fndr.cutoff(Zscore, statistic="normal")
            if (co >= 0.3) {
                cutoff.method <- "fndr"
            } else {
                cutoff.method <- "locfdr"
            }
        }
        fdrmod <- fdrtool(Zscore, statistic="normal", verbose=verbose,
            plot=plot, cutoff.method=cutoff.method, ...)

    } else {
        if (missing(cutoff.method)) {
            cutoff.method <- "fndr"
        }
        fdrmod <- fdrtool(P, statistic="pvalue", verbose=verbose,
            plot=plot, cutoff.method=cutoff.method,...)
    }
    fdrdf <- do.call(data.frame, fdrmod[c("pval", "qval", "lfdr")])
    for (i in names(fdrdf)) {
        ttab[[str_c("fdrtool.", i)]] <- fdrdf[[i]]
    }
    attr(ttab, "fdrtool") <- fdrmod
    ttab
}

# Functions for reading and writing narrowPeak files
read.narrowPeak <- function(file, ...) {
    peaks.df <- read.table(file, sep="\t", row.names=NULL, ...)
    names(peaks.df) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "summit")
    peaks.df$name <- as.character(peaks.df$name)
    peaks.df
}

write.narrowPeak <- function(x, file, ...) {
    x <- as(x, "data.frame")
    if("seqnames" %in% names(x))
        names(x)[names(x) == "seqnames"] <- "chr"
    x <- x[c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "summit")]
    write.table(x, file, sep="\t", row.names=FALSE, col.names=FALSE, ...)
}

# Like lapply, but returns future objects. The results can be fetched all at
# once with values().
future.lapply <- function(X, FUN, ...) {
    FUTUREFUN <- function(x) future(FUN(x, ...))
    lapply(X, FUTUREFUN, ...)
}

# Additional args passed to estimateDisp(), except for rawdisp, which can either be a vector of raw dispersions or
getBCVTable <- function(y, design, ..., rawdisp) {
    assert_that(is(y, "DGEList"))
    design.passed <- !missing(design)
    if (!design.passed) {
        # Can be NULL
        design <- y$design
    }
    # Estimate dispersions now if they are not already present
    if (! all(c("common.dispersion", "trended.dispersion", "tagwise.dispersion") %in%
                  names(y))) {
        if (is.null(design) && !design.passed) {
            warning("Estimating dispersions with no design matrix")
        }
        y %<-% estimateDisp(y, design=design, ...)
    }
    assert_that(!is.null(y$prior.df))
    if (all(y$prior.df == 0)) {
        if (missing(rawdisp) || is.null(rawdisp)) {
            rawdisp <- y
        }
        y %<-% estimateDisp(y, design=design, ...)
    }
    # Get raw (unsqueezed dispersions)
    if (missing(rawdisp) || is.na(rawdisp)) {
        # Estimate raw disperions using given design
        y.raw %<-% estimateDisp(y, design=design, prior.df=0)
    } else if (is.null(rawdisp)) {
        # Explicitly passing NULL means no raw disperions are desired.
        y.raw <- NULL
    } else if (is(rawdisp, "DGEList")) {
        # Assume DGEList already contains raw dispersions
        assert_that(all(dim(rawdisp) == dim(y)))
        y.raw <- rawdisp
    } else {
        # Assume anything else is a numeric vector of raw dispersions
        assert_that(is.numeric(rawdisp),
            length(rawdisp) == nrow(y))
        y.raw <- y
        y.raw$tagwise.dispersion <- rawdisp
        y.raw$prior.df <- y.raw$prior.n <- rep(0, length(rawdisp))
    }
    assert_that(all(y.raw$prior.df == 0))
    disptable <- y %>% as.list %$% data.frame(
        logCPM=AveLogCPM,
        CommonBCV=common.dispersion %>% sqrt,
        TrendBCV=trended.dispersion %>% sqrt,
        PriorDF=prior.df,
        eBayesBCV=tagwise.dispersion %>% sqrt)
    if (!is.null(y.raw)) {
        disptable$RawBCV <- y.raw$tagwise.dispersion %>% sqrt
    }
    return(disptable)
}

# ggplot version of edgeR::plotBCV. Additional arguments passed to getBCVTable
ggplotBCV <- function(y, xlab="Average log CPM", ylab="Biological coefficient of variation", rawdisp=NULL, ...) {
    if (is(y, "DGEList")) {
        disptable <- getBCVTable(y, ..., rawdisp=rawdisp)
    } else {
        disptable <- as.data.frame(y)
        assert_that(all(c("logCPM", "CommonBCV", "TrendBCV", "eBayesBCV") %in% names(disptable)))
    }

    # Reduce the number of points to plot for each line for performance
    # reasons
    npoints <- c(Common=2, Trend=500)
    disp.line.table <-
        disptable %>%
        select(logCPM, TrendBCV, CommonBCV) %>%
        melt(id.vars="logCPM", variable.name="DispType", value.name = "BCV") %>%
        mutate(DispType=str_replace(DispType, "BCV$", "")) %>%
        group_by(DispType) %>%
        do({
            approx(x=.$logCPM, y=.$BCV, n=npoints[.$DispType[1]]) %>%
                data.frame(logCPM=.$x, BCV=.$y)
        })

    p <- ggplot(disptable) +
        aes(x=logCPM)
    if ("RawBCV" %in% names(disptable)) {
        p <- p +
            geom_point(aes(y=RawBCV), size=0.4, color="black") +
            geom_density2d(aes(y=RawBCV), color="gray30", n=512) +
            labs(subtitle="Raw BCV (black) and eBayes-squeezed (blue)")
    }
    p <- p +
        geom_point(aes(y=eBayesBCV), size=0.1, color="darkblue") +
        geom_density2d(aes(y=eBayesBCV), color="blue", n=512) +
        geom_line(data=disp.line.table, aes(x=logCPM, y=BCV, group=DispType), color="white", size=1.5, alpha=0.5) +
        geom_line(data=disp.line.table, aes(x=logCPM, y=BCV, linetype=DispType), color="darkred", size=0.5) +
        scale_linetype_manual(name="Dispersion Type", values=c(Trend="solid", Common="dashed")) +
        labs(title="BCV plot", x=xlab, y=ylab)
    p
}

# Utilities for ggplot2 corrdinate transformation
power_trans <- function(pow) {
    name <- glue("^{pow}")
    trans_new(name,
        transform=function(x) x ^ pow,
        inverse=function(x) x ^ (1/pow),
        domain =c(0,Inf))
}

clamp_trans <- function(lower_threshold=0, upper_threshold=1) {
    name <- glue("Clamp values outside of [{lower_threshold}, {upper_threshold}]")
    trans_new(name,
        transform=function(x) pmin(upper_threshold, pmax(lower_threshold, x)),
        # transform is only invertible for part of the range
        inverse=identity)
}

# Do multiple runs of t-SNE and pick the one with the lowest final
# cost. (TODO: Make reproducible)
Rtsne.multi <- function(..., num.repeats=10, BPPARAM=bpparam()) {
    args <- list(...)
    results <- bplapply(seq_len(num.repeats), function(i) {
        do.call(Rtsne, args)
    })
    final.costs <- sapply(results, . %$% costs %>% tail(1))
    results[[which.min(final.costs)]]
}

# Alternate interface to removeBatchEffect: provide a single design
# matrix and specify the coefficient effects to remove. Note: use
# sum-to-zero contrasts for factors and mean-centered transformations
# for numerics if you want the grand means to remain unaffected by
# coefficient subtraction.
subtractCoefs <- function(x, design, coefsToSubtract, ...) {
    assert_that(!anyDuplicated(colnames(design)))
    subtract.design <- design[,coefsToSubtract]
    keep.design <- design[,setdiff(colnames(design), colnames(subtract.design))]
    removeBatchEffect(x, design=keep.design, covariates=subtract.design, ...)
}

BPselectModel <- function (y, designlist, criterion = "aic", df.prior = 0, s2.prior = NULL,
                           s2.true = NULL, ..., BPPARAM=bpparam())
{
    ym <- as.matrix(dge)
    if (any(is.na(ym)))
        stop("NAs not allowed")
    narrays <- ncol(ym)
    rm(ym)
    nmodels <- length(designlist)
    models <- names(designlist)
    if (is.null(models))
        models <- as.character(1:nmodels)
    if (df.prior > 0 & is.null(s2.prior))
        stop("s2.prior must be set")
    if (df.prior == 0)
        s2.prior <- 0
    criterion <- match.arg(criterion, c("aic", "bic", "mallowscp"))
    if (criterion == "mallowscp") {
        if (is.null(s2.true))
            stop("Need s2.true values")
        fits <- bplapply(designlist, lmFit, object=y, BPPARAM=BPPARAM)
        for (i in 1:nmodels) {
            fit <- fits[[i]]
            npar <- narrays - fit$df.residual[1]
            if (i == 1) {
                IC <- matrix(nrow = nrow(fit), ncol = nmodels,
                    dimnames = list(Probes = rownames(fit), Models = models))
                if (length(s2.true) != nrow(fit) && length(s2.true) != 1)
                    stop("s2.true wrong length")
            }
            IC[, i] <- fit$df.residual * fit$sigma^2/s2.true +
                npar * 2 - narrays
        }
    }
    else {
        ntotal <- df.prior + narrays
        penalty <- switch(criterion, bic = log(narrays), aic = 2)
        fits <- bplapply(designlist, lmFit, object=y, BPPARAM=BPPARAM)
        for (i in 1:nmodels) {
            fit <- fits[[i]]
            npar <- narrays - fit$df.residual[1] + 1
            s2.post <- (df.prior * s2.prior + fit$df.residual *
                            fit$sigma^2)/ntotal
            if (i == 1)
                IC <- matrix(nrow = nrow(fit), ncol = nmodels,
                    dimnames = list(Probes = rownames(fit), Models = models))
            IC[, i] <- ntotal * log(s2.post) + npar * penalty
        }
    }
    pref <- factor(apply(IC, 1, which.min), levels = 1:nmodels,
        labels = models)
    list(IC = IC, pref = pref, criterion = criterion)
}

# Return TRUE if namefun(obj) is a non-NULL character vector with no missing
# elements.
is.fully.named <- function(obj, namefun=names) {
    the.names <- namefun(obj)
    if (is.null(the.names)) {
        return(FALSE)
    } else if (!is.character(the.names)) {
        return(FALSE)
    } else if (any(is.na(the.names))) {
        return(FALSE)
    } else {
        return(TRUE)
    }
}

# Get log2-fold-change y-axis position of normalization lines between two
# samples for a list of DGEList objects.
getNormLineData <- function(dgelists, s1, s2) {
    assert_that(is.fully.named(dgelists))
    dgelists %>%
        sapply(. %>% {.$samples$norm.factors} %>% log2 %>% {.[s2] - .[s1]}) %>%
        data.frame(NormFactor=., NormType=names(dgelists))
}

# Get a curve representing the loess-based normaliation inplied by the offsets
# of two samples in a DGEList object. N is the number of points to interpolate the curve at.
getOffsetNormCurveData <- function(dge, s1, s2, n=1000) {
    assert_that(is.numeric(dge$offset))
    a <- aveLogCPM(dge, dispersion=0.05, prior.count=0.5)
    # Need to subtract the library size difference out of the offset
    raw.offset <- dge$offset %>% {.[,s2] - .[,s1]} %>% divide_by(log(2))
    lib.size.offset <- dge$samples$lib.size %>% {.[s2] / .[s1]} %>% log2
    x <- data.frame(A=a, Offset=raw.offset - lib.size.offset)
    f <- approxfun(x$A, x$Offset)
    data.frame(A=seq(from=min(x$A), to=max(x$A), length.out = n)) %>%
        mutate(M=f(A))
}

## Read MotifMap-provided BED file into a GRanges object. We can't use
## rtracklayer::import.bed because it chokes on spaces in fields,
## which MotifMap contains.
read.motifmap <- function(x, parse_name=TRUE) {
    tab <- read_tsv(x, col_names=c("chr", "start", "end", "name", "score", "strand"),
                    col_types="ciicdc", progress=FALSE)
    if (parse_name) {
        tab %<>% separate(name, into=c("motif_ID", "TF_name"), sep="=")
    }
    gr <- makeGRangesFromDataFrame(tab, starts.in.df.are.0based=TRUE)
    gr
}

write.motifmap <- function(x, file) {
    assert_that(is(x, "GRanges"))
    if (! "name" %in% names(mcols(x))) {
        assert_that(all(c("motif_ID", "TF_name") %in% names(mcols(x))))
        mcols(x) %<>% as.data.frame %>% unite(name, c(motif_ID, TF_name), sep="=") %>% as("DataFrame")
    }
    export(x, file, format="BED")
}

## Like rtracklayer::liftOver but "fills in" small gaps induced by the
## liftOver process (i.e. no larger than allow.gap). If allow.gap is
## zero, this is equivalent to liftOver.
liftOverLax <- function(x, chain, ..., allow.gap=0) {
    newx <- liftOver(x, chain)
    if (allow.gap > 0) {
        gapped <- which(lengths(newx) > 1)
        newx.gapped.reduced <- reduce(newx[gapped], min.gapwidth = allow.gap + 1, with.revmap=TRUE)
        mcols(newx.gapped.reduced@unlistData) <- rep(mcols(x[gapped]), lengths(newx.gapped.reduced))
        newx[gapped] <- newx.gapped.reduced
    }
    return(newx)
}

# Relevel multple columns using fct_relevel
relevel_columns <- function(df, ...) {
    relevel_specs <- list(...)
    assert_that(is_named(relevel_specs))
    for (i in names(relevel_specs)) {
        df[[i]] %<>% fct_relevel(relevel_specs[[i]])
    }
    df
}

# Variant of save.image that allows excluding specific names
save.image.filtered <- function (file = ".RData", version = NULL, ascii = FALSE, compress = !ascii,
                                 safe = TRUE, exclude = NULL)
{
    if (!is.character(file) || file == "")
        stop("'file' must be non-empty string")
    opts <- getOption("save.image.defaults")
    if (is.null(opts))
        opts <- getOption("save.defaults")
    if (missing(safe) && !is.null(opts$safe))
        safe <- opts$safe
    if (missing(ascii) && !is.null(opts$ascii))
        ascii <- opts$ascii
    if (missing(compress) && !is.null(opts$compress))
        compress <- opts$compress
    if (missing(version))
        version <- opts$version
    if (safe) {
        outfile <- paste0(file, "Tmp")
        i <- 0
        while (file.exists(outfile)) {
            i <- i + 1
            outfile <- paste0(file, "Tmp", i)
        }
    }
    else outfile <- file
    on.exit(file.remove(outfile))
    vars.to.save <- setdiff(names(.GlobalEnv), exclude)
    save(list = vars.to.save, file = outfile, version = version,
         ascii = ascii, compress = compress, envir = .GlobalEnv,
         precheck = FALSE)
    if (safe)
        if (!file.rename(outfile, file)) {
            on.exit()
            stop(gettextf("image could not be renamed and is left in %s",
                          outfile), domain = NA)
        }
    on.exit()
}

# Unlike load, returns the environment itself
load.in.new.env <- function(file, envir=new.env(), ...) {
    load(file, envir, ...)
    return(envir)
}

load.filtered <- function(file, envir = parent.frame(), ..., exclude=NULL) {
    if (!length(exclude)) {
        return(load(file, envir, ...))
    }
    tempenv <- load.in.new.env(file=file, ...)
    for (i in setdiff(names(tempenv), exclude)) {
        envir[[i]] <- tempenv[[i]]
    }
}


## Custom runMOFA function (see https://github.com/PMBio/MOFA/pull/6 )

##############################################
## Functions to run MOFA from the R package ##
##############################################

#' @title runMOFA:
#' @name runMOFA
#' @description train a \code{\link{MOFAmodel}}
#' @param object an untrained \code{\link{MOFAmodel}}
#' @param DirOptions list with I/O options, should contain at least 'dataDir' where the input matrices as stored as .txt files and 'outFile' where the model is going to be stored as a .hdf5 file
#' @param ... Extra options to add to the mofa command
#' @return a trained \code{\link{MOFAmodel}}
#' @export
runMOFA <- function(object, DirOptions, ..., mofaPath="mofa") {

    # Sanity checks
    if (! is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
    stopifnot(all(c("dataDir","outFile") %in% names(DirOptions)))

    arglist <- list(
        inFiles = paste0(DirOptions$dataDir, "/", viewNames(object), ".txt"),
        header_cols = TRUE,
        header_rows = TRUE,
        delimiter = object@DataOpts$delimiter,
        outFile = DirOptions$outFile,
        views = viewNames(object),
        likelihoods = object@ModelOpts$likelihood,
        factors = object@ModelOpts$numFactors,
        iter = object@TrainOpts$maxiter,
        dropR2 =  object@TrainOpts$DropFactorThreshold,
        tolerance = object@TrainOpts$tolerance
    )

    # Setting the below arguments to NULL doesn't actually add them to
    # the argument list, but reserves that argument name to prevent
    # extra.arglist from using it.
    if (!is.null(object@ModelOpts$covariates)) {
        arglist$covariatesFile <- file.path(DirOptions$dataDir, "covariates.txt")
        arglist$scale_covariates <- rep(1,ncol(object@ModelOpts$covariates))
    } else {
        arglist$covariatesFile <- NULL
        arglist$scale_covariates <- NULL
    }
    arglist$learnIntercept <- as.logical(object@ModelOpts$learnIntercept)
    if (! object@ModelOpts$sparsity) {
        arglist$learnTheta <- 0
    } else {
        arglist$learnTheta <- NULL
    }

    arglist$center_features <- as.logical(object@DataOpts$centerFeatures)
    arglist$scale_views <- as.logical(object@DataOpts$scaleViews)
    if (object@DataOpts$removeIncompleteSamples == T) { command <- paste(command, "--RemoveIncompleteSamples", sep=" ") }
    arglist$verbose <- as.logical(object@TrainOpts$verbose)

    extra.arglist <- list(...)
    if (any(is.na(names(extra.arglist)) | names(extra.arglist) == "")) {
        stop("All extra options must be named")
    }

    # Remove leading "--" from extra arg names if present (it will be
    # added back later)
    names(arglist) <- sub("^--", "", names(arglist))
    conflicting.argnames <- intersect(names(extra.arglist), names(arglist))
    if (length(conflicting.argnames) > 0)
        stop(paste0("You cannot pass the following arguments as extra options to runMOFA: ",
                    deparse(conflicting.argnames)))

    # No conflicting argument names,
    arglist <- c(arglist, extra.arglist)

    argv <- character(0)
    for (argname in names(arglist)) {
        argval <- arglist[[argname]]
        argname <- paste0("--", argname)

        if (is.null(argval)) {
            # Placeholder option; don't add it
        }
        if (is.logical(argval)) {
            # Flag option
            if (length(argval) != 1) {
                stop(paste("Invalid argument value:", deprase(argval)))
            } else if (argval == FALSE || is.na(argval)) {
                # Unset flag: don't add it
            } else if (argval == TRUE) {
                # Set flag: add it
                argv <- c(argv, argname)
            }
        } else {
            # Option with arguments: add the option followed by it args
            argv <- c(argv, argname, argval)
        }
    }
    argv <- unlist(argv)

    if (length(mofaPath) != 1) stop("Invalid mofaPath")

    # If output already exists, remove it
    if (file.exists(DirOptions$outFile)) {
        if (arglist$verbose) {
            message("Deleting old output file")
        }
        file.remove(DirOptions$outFile)
    }

    if (arglist$verbose) {
        message("Running MOFA command: ", paste(collapse=" ", shQuote(c(mofaPath, argv))))
    }
    # Run!
    exitcode <- system2(command=mofaPath, args=shQuote(argv), wait=T)
    if (exitcode != 0) {
        stop(paste("mofa command failed with exit code", exitcode))
    }

    # Load trained model
    object <- loadModel(DirOptions$outFile, object)

    return(object)
}

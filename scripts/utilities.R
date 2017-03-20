library(magrittr)
library(dplyr)
library(assertthat)
library(BiocParallel)
library(rex)
library(lazyeval)
library(future)
library(Rtsne)

withGC <- function(expr) {
    on.exit(gc())
    return(expr)
}

# Use to assign to complex sub-expressions in the middle of a dplyr pipeline.
# For example, you can't easily do the following in the middle of a pipeline:
# "assays(x[[1]])$counts[3,5] <- 45". But now you can do it like: "x %>%
# assign_into(assays(.[[1]])$counts[3,5], 45) %>% another_fun() %>% ..."
assign_into <- function(x, expr, value) {
    expr <- lazy(expr)$expr
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

# Somewhat surprisingly, this function doesn't actually rely on ggplot2 being
# loaded.
ggprint <- function(plots, device=dev.cur(), closedev, printfun=print) {
    orig.device <- dev.cur()
    new.device <- device
    # Functions that create devices don't generally return them, they just set
    # them as the new current device, so get the actual device from dev.cur()
    # instead.
    if (is.null(new.device)) {
        new.device <- dev.cur()
    }
    if (missing(closedev)) {
        closedev = orig.device != new.device && new.device != 1
    }
    if (!is.null(device) && !is.na(device)) {
        dev.set(new.device)
        on.exit({
            if (closedev) {
                dev.off(new.device)
            }
            if (new.device != orig.device) {
                dev.set(orig.device)
            }
        })
    }
    assertthat::assert_that(is(plots, "gg") ||
                                is.list(plots))
    if (is(plots, "gg")) {
        printfun(plots)
    } else if (is.list(plots)) {
        lapply(plots, ggprint, device=NA, closedev=FALSE, printfun=printfun)
    } else {
        stop("Argument is not a ggplot or list of ggplots")
    }
    invisible(NULL)
}

# Printer function for ggplotly, to be passed as the prinfun for ggprint
ggplotly.printer <- function(...) {
    dots <- list(...)
    function(p) {
        args <- c(list(p=p), dots)
        print(do.call(plotly::ggplotly, args))
    }
}

## Based on this: https://www.r-bloggers.com/shrinking-rs-pdf-output/
## Requires some command-line image and PDF manip tools.
rasterpdf <- function(pdffile, outfile=pdffile, resolution=600, verbose=FALSE) {
    vmessage <- function(...) {
        if (verbose) message(...)
    }
    require(parallel)
    wd=getwd()
    td=tempfile(pattern="rasterpdf")
    dir.create(td, recursive = TRUE)
    on.exit(unlink(td, recursive=TRUE))
    file.copy(pdffile, file.path(td,"toraster.pdf"))
    setwd(td)
    on.exit(setwd(wd), add=TRUE)
    assert_that(file.exists("toraster.pdf"))
    system2("pdftk", args=c("toraster.pdf", "burst"))
    files=list.files(pattern="pg_")

    vmessage(paste0("Rasterizing ",length(files)," pages:  (",paste(files,collapse=","),")"))
    bplapply(files,function(f){
        system2("gs", args=c("-dBATCH", "-dTextAlphaBits=4", "-dNOPAUSE",
                             paste0("-r", resolution), "-q", "-sDEVICE=png16m",
                             paste0("-sOutputFile=",f,".png"),f))
        system2("convert", args=c("-quality", "100", "-density", resolution,
                                  paste0(f,".png"),
                                  paste0(strsplit(f,".",fixed=T)[[1]][1],".pdf")))
        vmessage(paste0("Finished page ",f))
        return()
    })
    vmessage("Compiling the final pdf")
    file.remove("toraster.pdf")
    file.remove(list.files(pattern="png"))
    setwd(wd)
    system2("gs", args=c("-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=pdfwrite",
                         paste0("-sOutputFile=",outfile),
                         list.files(path=td, pattern=glob2rx("*.pdf"), full.names=TRUE)))
    vmessage("Finished!!")
}

add.numbered.colnames <- function(x, prefix="C") {
    x %>% set_colnames(sprintf("%s%i", prefix, seq(from=1, length.out=ncol(x))))
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

## Versions of cpm and aveLogCPM that use an offset matrix instead of lib sizes
cpmWithOffset <- function(dge, offset=expandAsMatrix(getOffset(dge), dim(dge)),
                          log = FALSE, prior.count = 0.25, preserve.mean=TRUE, ...) {
    x <- dge$counts
    if (preserve.mean) {
        # Ensure that the mean logcpm is not changed by the offsets, by setting
        # each row mean to the mean offset derived from the library sizes
        mean.lib.size.offset <- dge %>% assign_into(.$offset, NULL) %>% getOffset %>% mean
        offset %<>% subtract(rowMeans(.)) %>% add(mean.lib.size.offset)
    }
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

aveLogCPMWithOffset <- function(y, ...) {
    UseMethod("aveLogCPM")
}

aveLogCPMWithOffset.default <- function (y, offset = NULL, prior.count = 2,
                                         dispersion = NULL, weights = NULL, ...)
{
    aveLogCPM(y, lib.size = NULL, offset = offset, prior.count = prior.count, dispersion = dispersion, weights = weights, ...)
}

aveLogCPMWithOffset.DGEList <- function (
    y, offset = expandAsMatrix(getOffset(y), dim(y)),
    prior.count = 2, dispersion = NULL, ...) {
    if (is.null(dispersion)) {
        dispersion <- y$common.dispersion
    }
    aveLogCPMWithOffset(
        y$counts, offset = y$offset, prior.count = prior.count,
        dispersion = dispersion, weights = y$weights)
}

library(limma)

## Version of voom that uses an offset matrix instead of lib sizes
voomWithOffset <- function (
    dge, design = NULL, offset=expandAsMatrix(getOffset(dge), dim(dge)),
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
voomWithQualityWeightsAndOffset <-function (
    dge, design = NULL,
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

# Convenience function for alternating dupCor and voom until convergence
voomWithDuplicateCorrelation <- function(
    counts, design = NULL, plot = FALSE, block = NULL, trim = 0.15,
    voom.fun=voom, dupCor.fun=duplicateCorrelation, initial.correlation=0,
    niter=5, tol=1e-6, verbose=TRUE, ...) {
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
        message(sprintf("Initial guess for duplicate correlation before 1st iteration: %s", initial.correlation))
    }
    dupcor <- dupCor.fun(elist, design, block=block, trim=trim)
    iter.num <- iter.num + 1
    if (verbose) {
        message(sprintf("Duplicate correlation after %s iteration: %s", toOrdinal(iter.num), dupcor$consensus.correlation))
    }
    while (iter.num < niter) {
        if (!is.null(tol) && is.finite(tol)) {
            if (abs(dupcor$consensus.correlation - prev.cor) <= tol) {
                if (verbose) {
                    message(sprintf("Stopping after %s iteration because tolerance threshold was reached.", toOrdinal(iter.num)))
                }
                break
            }
        }
        prev.cor <- dupcor$consensus.correlation
        elist <- voom.fun(counts, design=design, plot=FALSE, block=block, correlation=prev.cor, ...)
        dupcor <- dupCor.fun(elist, design, block=block, trim=trim)
        iter.num <- iter.num + 1
        if (verbose) {
            message(sprintf("Duplicate correlation after %s iteration: %s", toOrdinal(iter.num), dupcor$consensus.correlation))
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
code_control_named <- function (n, contrasts = TRUE, sparse = FALSE)
{
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
        xlim(0,1) + ggtitle(sprintf("P-value distribution (Est. %0.2f%% signif.)",
                                    100 * (1-ptn))) +
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
    ## Computing from 1-PP gives better numerical precision for the
    ## most significant genes (large B-values)
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
    cbind(ttab, bfdr(B)[c("PP", "BayesFDR")])
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
        qobj %$%
            cbind(ttab,
                  QValue=qvalues,
                  LocFDR=lfdr)
    }, error=function(e) {
        warning(str_c("Failed to compute q-values: ", e$message))
        ttab
    })
}

# Functions for reading and writing narrowPeak files
read.narrowPeak <- function(file, ...) {
    peaks.df <- read.table(file, sep="\t", row.names=NULL, ...)
    names(peaks.df) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "summit")
    peaks.df$name <- as.character(peaks.df$name)
    ## havenames <- !any(peaks.df$name == ".")
    ## res <- data.frame2GRanges(peaks.df, keepColumns=TRUE, startOffset=1, endOffset=0)
    ## ## Eliminate the dummy row names from the data
    ## if (havenames)
    ##     names(res) <- res$name
    ## else
    ##     names(res) <- NULL
    ## res
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

    ## Reduce the number of points to plot for each line for performance
    ## reasons
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

## Do multiple runs of t-SNE and pick the one with the lowest final
## cost. (TODO: Make reproducible)
Rtsne.multi <- function(..., num.repeats=10, BPPARAM=BPPARAM()) {
    args <- list(...)
    results <- bplapply(seq_len(num.repeats), function(i) {
        do.call(Rtsne, args)
    })
    final.costs <- sapply(results, . %$% costs %>% tail(1))
    results[[which.min(final.costs)]]
}

## Alternate interface to removeBatchEffect: provide a single design
## matrix and specify the coefficient effects to remove. Note: use
## sum-to-zero contrasts for factors and mean-centered transformations
## for numerics if you want the grand means to remain unaffected by
## coefficient subtraction.
subtractCoefs <- function(x, design, coefsToSubtract, ...) {
    assert_that(!anyDuplicated(colnames(design)))
    subtract.design <- design[,coefsToSubtract]
    keep.design <- design[,setdiff(colnames(design), colnames(subtract.design))]
    removeBatchEffect(x, design=keep.design, covariates=subtract.design, ...)
}

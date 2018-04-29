#!/usr/bin/env Rscript

library(getopt)
library(optparse)
library(rctutils)

## Don't default to more than 4 cores
num.cores <- min(availableCores(), 4)

get_options <- function(opts) {

    ## Do argument parsing early so the script exits quickly if arguments are invalid
    optlist <- list(
        make_option(c("-s", "--samplemeta-file"), metavar="FILENAME.RDS", type="character",
                    help="(REQUIRED) RDS/RData/xlsx/csv file containing a table of sample metadata. Any existing rownames will be replaced with the values in the sample ID  column (see below)."),
        make_option(c("-c", "--sample-id-column"), type="character", default="Sample",
                    help="Sample metadata column name that holds the sample IDs. These will be substituted into '--bam-file-pattern' to determine the BAM file names."),
        make_option(c("-f", "--filter-sample-ids"), type="character",
                    help="Comma-separated list of sample IDs. If this options is provided, only the specified sample IDs will be used."),
        make_option(c("-p", "--bam-file-pattern"), metavar="PATTERN", type="character",
                    help="(REQUIRED) Format string to convert sample IDs into BAM file paths. This should contain the string '{SAMPLE}' wherever the sample ID should be substituted (this can occur multiple times),. Example: 'bam_files/Sample_{SAMPLE}/Aligned.bam"),
        make_option(c("-o", "--output-file"), metavar="FILENAME.RDS", type="character",
                    help="(REQUIRED) Output file name. The SummarizedExperiment object containing the counts will be saved here using saveRDS, so it should end in '.RDS'."),
        make_option(c("-b", "--expected-bam-files"), metavar="BAMFILE1,BAMFILE2,...", type="character",
                    help="Comma-separated list of bam file names expected to be used as input. This argument is optional, but if it is provided, it will be checked against the list of files determined from '--samplemeta-file' and '--bam-file-pattern', and an error will be raised if they don't match exactly."),
        make_option(c("-j", "--threads"), metavar="N", type="integer", default=num.cores,
                    help="Number of threads to use while counting reads"),
        make_option(c("-t", "--annotation-txdb"), metavar="PACKAGENAME", type="character",
                    help="Name of TxDb package, or the name of a database file, to use for gene annotation"),
        make_option(c("-g", "--annotation-gff"), metavar="FILENAME", type="character",
                    help="File Name of GFF3 file to use for gene annotation."),
        make_option(c("-f", "--gff-exon-featuretype"), metavar="FEATURETYPE", type="character", default="exon",
                    help="GFF feature type to use"),
        make_option(c("-i", "--gff-geneid-attr"), metavar="ATTRNAME", type="character", default="gene_id",
                    help="GFF feature attribute to use as a feature's Gene ID."),
        make_option(c("-e", "--gff-gene-featuretype"), metavar="FEATURETYPE", type="character", default="gene",
                    help="GFF feature type from which gene metadata should be extracted."),
        make_option(c("-r", "--annotation-rds"), metavar="FILENAME", type="character",
                    help="File Name of RDS or RData file to use for gene annotation. It should contain a single GRanges or GRangesList object (or something that can be coereced into one), with each element representing a gene/feature to be counted. Any metadata columns on the object will be carried through to the output SummarizedExperiment."),
        make_option(c("--annotation-saf"), metavar="FILENAME", type="character",
                    help="File Name of RDS/RData/xlsx/csv file containing gene annotations in SAF format (i.e. 5 columns named GeneID, Chr, Start, End, Strand). Additional columns beyond the first 5 will be retained as metadata columns on the genes/exons."),
        make_option(c("-a", "--additional-gene-info"), metavar="FILENAME", type="character",
                    help="RDS/RData/xlsx/csv file containing a table of gene metadata. Row names (or the first column of the file if there are no row names) should be gene/feature IDs that match the ones used in the main annotation, and these should be unique. This can also be a GFF3 file where the metadata is in the attributes of elements of type specified by '--gff-gene-featuretype' ('gene' by default), where the 'ID' attribute specifies the gene ID."))
    progname <- na.omit(c(get_Rscript_filename(), "rnaseq-count.R"))[1]
    parser <- OptionParser(
        usage="Usage: %prog [options] -s SAMPLEMETA.RDS -p PATTERN -o SUMEXP.RDS [ -t TXDB.PACKAGE.NAME | -g ANNOTATION.GFF3 | -r ANNOTATION.RDS ]",
        description="Count reads in genes using Rsubread::featureCounts.

Counts are performed for stranded, non-stranded, and reverse-stranded modes, and are stored along with the sample and gene metadata in a SummarizedExperiment object. Note that the '-s', '-p', and '-o' arguments are all required, since they specify the input and output files. Also, exactly one of '-t', '-g', or '-r' is required to specify the annotation.",
option_list = optlist,
add_help_option = TRUE,
prog = progname,
epilogue = "")

    cmdopts <- parse_args(parser, opts)
    ## Ensure that all required arguments were provided
    required.opts <- c("samplemeta-file", "bam-file-pattern", "output-file")
    missing.opts <- setdiff(required.opts, names(cmdopts))
    if (length(missing.opts) > 0) {
        stop(str_c("Missing required arguments: ", deparse(missing.opts)))
    }

    ## Split list arguments
    for (i in c("filter-sample-ids", "expected-bam-files")) {
        if (i %in% names(cmdopts)) {
            cmdopts[[i]] %<>% str_split(",") %>% unlist
        }
    }

    ## Ensure that exactly one annotation was provided
    annot.opts <- c("annotation-txdb", "annotation-gff", "annotation-rds", "annotation-saf")
    provided.annot.opts <- intersect(annot.opts, names(cmdopts))
    if (length(provided.annot.opts) < 1) {
        stop("No annotation provided")
    } else if (length(provided.annot.opts) > 1) {
        stop("Multiple annotations provided. Please provide only one annotation.")
    }

    ## Replace dashes with underscores so that all options can easily
    ## be accessed by "$"
    cmdopts %>% setNames(chartr("-", "_", names(.)))
}

## Do this early to handle "--help" before wasting time loading
## pacakges & stuff
invisible(get_options(commandArgs(TRUE)))

library(assertthat)
library(dplyr)
library(magrittr)
library(openxlsx)
library(stringr)

library(annotate)
library(GenomicRanges)
library(Rsubread)
library(SummarizedExperiment)
library(withr)
library(org.Hs.eg.db)

## Guess type of ID (currenly unused)
identify.ids <- function(ids, db="org.Hs.eg.db", idtypes=c("ENTREZID", "ENSEMBL", "UNIGENE"), threshold=0.5) {
    if (is.character(db)) {
        library(db, character.only=TRUE)
        pos <- str_c("package:", db)
        db <- get(db, pos)
    }
    assert_that(is(db, "AnnotationDb"))
    idtypes %<>% intersect(keytypes(db))
    assert_that(length(idtypes) > 0)
    idcounts <- sapply(idtypes, function(idtype) {
        sum(ids %in% keys(db, idtype))
    }) %>% setNames(idtypes)
    idcounts %<>% sort(decreasing = TRUE)
    result <- names(idcounts)[1]
    if (idcounts[result] / length(ids) < threshold) {
        stop(glue("Could not identify more than {format(threshold * 100, digits=2)}%% of given IDs as any of {deparse(idtypes)}"))
    }
    tsmsg("Detected gene IDs as ", result)
    attr(result, "counts") <- idcounts
    return(result)
}

{
    cmdopts <- get_options(commandArgs(TRUE))
    ## myargs <- c("-s", "./saved_data/samplemeta-RNASeq.RDS", "-c", "SRA_run", "-p", "aligned/rnaseq_star_hg38.analysisSet_gencode.v25/%s/Aligned.sortedByCoord.out.bam", "-o", "sexp.rds", "-j", "2", "-g", "~/references/hg38/gencode.v25.gff3")
    ## cmdopts <- get_options(myargs)
    cmdopts$help <- NULL

    cmdopts$threads %<>% round %>% max(1)
    tsmsg("Running with ", cmdopts$threads, " threads")
    use_multicore_futures(cmdopts$threads)

    tsmsg("Args:")
    print_var_vector(cmdopts)

    ## Delete the output file if it exists
    suppressWarnings(file.remove(cmdopts$output_file))
    assert_that(!file.exists(cmdopts$output_file))

    tsmsg("Loading sample metadata")
    samplemeta <- read_table_general(cmdopts$samplemeta_file)

    tsmsg("Got metadata for ", nrow(samplemeta), " samples")

    assert_that(cmdopts$sample_id_column %in% colnames(samplemeta))
    assert_that(!anyDuplicated(samplemeta[[cmdopts$sample_id_column]]))

    rownames(samplemeta) <- samplemeta[[cmdopts$sample_id_column]]

    samplemeta$bam_file <- glue(cmdopts$bam_file_pattern, SAMPLE=samplemeta[[cmdopts$sample_id_column]], .envir=emptyenv())

    if (!is.null(cmdopts$filter_sample_ids)) {
        tsmsg("Selecting only ", length(cmdopts$filter_sample_ids), " specified samples.")
        assert_that(all(cmdopts$filter_sample_ids %in% samplemeta[[cmdopts$sample_id_column]]))
        samplemeta %<>% .[.[[cmdopts$sample_id_column]] %in% cmdopts$filter_sample_ids,]
    }

    if ("expected_bam_files" %in% names(cmdopts)) {
        tryCatch({
            assert_that(setequal(samplemeta$bam_file, cmdopts$expected_bam_files))
            tsmsg("Sample metadata contains all expected bam files")
        }, error=function(...) {
            unexpected_existing <- setdiff(samplemeta$bam_file, cmdopts$expected_bam_files)
            expected_but_missing <- setdiff(cmdopts$expected_bam_files, samplemeta$bam_file)
            if (length(unexpected_existing) > 0) {
                tsmsg(glue("Got unexpected bam files: {deparse(unexpected_existing)}"))
            }
            if (length(expected_but_missing) > 0) {
                tsmsg(glue("Didn't find expected bam files: {deparse(expected_but_missing)}"))
            }
            stop("Bam file list was not as expected")
        })
    }

    assert_that(all(file.exists(samplemeta$bam_file)))

    tsmsg("Reading annotation data")

    if ("annotation_txdb" %in% names(cmdopts)) {
        txdb <- get_txdb(cmdopts$annotation_txdb)
        annot <- exonsBy(txdb, "gene")
    } else if ("annotation_gff" %in% names(cmdopts)) {
        annot <- cmdopts %$%
            read_annotation_from_gff(
                annotation_gff,
                exonFeatureType=gff_exon_featuretype,
                geneIdAttr=gff_geneid_attr,
                geneFeatureType=gff_gene_featuretype)
    } else if ("annotation_rds" %in% names(cmdopts)) {
        annot <- read_annotation_from_rdata(cmdopts$annotation_rds)
    } else if ("annotation_saf" %in% names(cmdopts)) {
        annot <- read_annotation_from_saf(cmdopts$annotation_saf)
    }

    assert_that(is(annot, "GRangesList"))
    tsmsg("Annotation has ", length(annot), " features")

    if ("additional_gene_info" %in% names(cmdopts)) {
        tsmsg("Reading additional gene annotation metadata")
        additional_gene_info <- read_additional_gene_info(cmdopts$additional_gene_info)
        genes_without_info <- setdiff(names(annot), rownames(additional_gene_info))
        if (length(genes_without_info) > 0) {
            empty_row <- list(character(0)) %>% rep(ncol(additional_gene_info)) %>% setNames(colnames(additional_gene_info))
            single.val.cols <- sapply(additional_gene_info, function(x) all(lengths(x) == 1))
            for (i in seq_along(empty_row)) {
                if (single.val.cols[i]) {
                    empty_row[[i]] <- NA
                } else {
                    empty_row[[i]] <- list(logical(0)) %>% as(class(additional_gene_info[[i]]))
                }
            }
            empty_row %<>% DataFrame
            empty_gene_table <- empty_row[rep(1, length(genes_without_info)),] %>%
                set_rownames(genes_without_info)
            additional_gene_info %<>% rbind(empty_gene_table)
        }
        assert_that(all(names(annot) %in% rownames(additional_gene_info)))
        mcols(annot)[colnames(additional_gene_info)] <- additional_gene_info[names(annot),]
        metadata(annot) %<>% c(metadata(additional_gene_info))
    }

    saf <- grl_to_saf(annot)

    if (all(lengths(annot) == 1)) {
        annot.mcols <- mcols(annot)
        annot <- unlist(annot)
        mcols(annot) <- annot.mcols
        rm(annot.mcols)
    }

    empty.counts <- matrix(NA, nrow=length(annot), ncol=nrow(samplemeta))
    sexp <- SummarizedExperiment(
        assays=List(
            counts=empty.counts,
            sense.counts=empty.counts,
            antisense.counts=empty.counts
        ),
        colData=as(samplemeta, "DataFrame"),
        rowRanges=annot,
        metadata=list()
    )
    colnames(sexp) <- colData(sexp)[[cmdopts$sample_id_column]]
    rownames(sexp) <- names(annot)

    tsmsg("Computing sense counts")
    sense.fc <- featureCountsParallel(
        samplemeta$bam_file, annot.ext=saf,
        strandSpecific=1)
    tsmsg("Computing antisense counts")
    antisense.fc <- featureCountsParallel(
        samplemeta$bam_file, annot.ext=saf,
        strandSpecific=2)
    tsmsg("Computing unstranded counts")
    unstranded.fc <- featureCountsParallel(
        samplemeta$bam_file, annot.ext=saf,
        strandSpecific=0)
    assay(sexp, "counts")[,] <- unstranded.fc$counts
    assay(sexp, "sense.counts")[,] <- sense.fc$counts
    assay(sexp, "antisense.counts")[,] <- antisense.fc$counts

    count.stats <- List(counts=unstranded.fc$stat,
                        sense.counts=sense.fc$stat,
                        antisense.counts=antisense.fc$stat) %>%
        endoapply(function(x) {
            x <- data.frame(set_colnames(t(x[,-1]), x[[1]]))
            rownames(x) <- colnames(sexp)
            x
        })
    metadata(sexp)$stat <- count.stats

    tsmsg("Saving SummarizedExperiment")
    save_RDS_or_RDA(sexp, cmdopts$output_file)
    invisible(NULL)
}

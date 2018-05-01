#!/usr/bin/env Rscript

library(getopt)
library(optparse)
library(magrittr)
library(assertthat)
library(glue)
library(rctutils)

get_options <- function(opts) {
    optlist <- list(
        make_option(c("-s", "--samplemeta-file"), metavar = "FILENAME.RDS", type = "character",
                    help = "(REQUIRED) RDS/RData/xlsx/csv file containing a table of sample metadata. Any existing rownames will be replaced with the values in the sample ID  column (see below)."),
        make_option(c("-c", "--sample-id-column"), type = "character", default = "Sample",
                    help = "Sample metadata column name that holds the sample IDs. These will be substituted into '--abundance-file-pattern' to determine the abundance file names."),
        make_option(c("-a", "--abundance-file-pattern"), metavar = "PATTERN", type = "character",
                    help = "(REQUIRED) Format string to convert sample IDs into file paths to the abundance.h5 file for each sample. This should contain the string '{SAMPLE}' wherever the sample ID should be substituted (this can occur multiple times),. Example: 'kallisto_quant/Sample_{SAMPLE}/abundance.h5'"),
        make_option(c("-l", "--aggregate-level"), metavar = "LEVEL", type = "character", default = "auto",
                    help = "Whether to save aggregated gene counts or transcript counts in the output file. By default, aggregated gene counts are saved if a gene annotation is provided, and transcript counts are saved otherwise. You can force one or the other by specifying 'gene' or 'transcript' for this option."),
        make_option(c("-o", "--output-file"), metavar = "FILENAME.RDS", type = "character",
                    help = "(REQUIRED) Output file name. The SummarizedExperiment object containing the counts will be saved here using saveRDS, so it should end in '.RDS'."),
        make_option(c("-b", "--expected-abundance-files"), metavar = "FILE1.h5,FILE2.h5,...", type = "character",
                    help = "Comma-separated list of file names expected to be used as input. This argument is optional, but if it is provided, it will be checked against the list of files determined from '--samplemeta-file' and '--abundance-file-pattern', and an error will be raised if they don't match exactly."),
        make_option(c("-m", "--genemap-file"), metavar = "FILENAME", type = "character",
                    help = "Genemap file in the Salmon simple gene map format (see 'salmon quant --help-reads')"),
        make_option(c("-d", "--annotation-txdb"), metavar = "PACKAGE_OR_FILE_NAME", type = "character",
                    help = "Name of TxDb package, or the name of a database file, to use for gene annotation"),
        make_option(c("-g", "--gene-info"), metavar = "FILENAME", type = "character",
                    help = "RDS/RData/xlsx/csv file containing a table of gene metadata. Row names (or the first column of the file if there are no row names) should be gene/feature IDs that match the ones used in the main annotation, and these should be unique. This option is ignored when not aggregating counts to the gene level."),
        make_option(c("--transcript-info"), metavar = "FILENAME", type = "character",
                    help = "RDS/RData/xlsx/csv file containing a table of transcript metadata. Row names (or the first column of the file if there are no row names) should be transcript IDs that match the ones used in the quantification files, and these should be unique. This option is ignored when aggregating counts to the gene level."))
    progname <- na.omit(c(get_Rscript_filename(), "convert-quant-to-sexp.R"))[1]
    parser <- OptionParser(
        usage = "Usage: %prog [ -d TXDB | -m GENEMAP ] [ -g GENEINFO | -t TXINFO ] -s SAMPLEMETA.RDS -a PATTERN -l (gene|transcript) -o SUMEXP.RDS",
        description = "Collect RNA-seq quantification results into a SummarizedExperiment object.

TODO UPDATE Counts are stored along with the sample and gene metadata in a SummarizedExperiment object. Note that the '-s', '-a', '-t', and '-o' arguments are all required, since they specify the essential input and output files and formats.",
option_list = optlist,
add_help_option = TRUE,
prog = progname,
epilogue = "")

    cmdopts <- parse_args(parser, opts)
    ## Ensure that all required arguments were provided
    required.opts <- c("samplemeta-file", "abundance-file-pattern", "output-file")
    missing.opts <- setdiff(required.opts, names(cmdopts))
    if (length(missing.opts) > 0) {
        stop(str_c("Missing required arguments: ", deparse_onestring(missing.opts)))
    }

    ## Ensure that no more than one annotation was provided, and that
    ## exactly one annotation was provided if gene-level aggregation
    ## was requested.
    annot.opts <- c("annotation-txdb", "genemap-file")
    provided.annot.opts <- intersect(annot.opts, names(cmdopts))
    if (length(provided.annot.opts) > 1) {
        stop("Multiple gene annotations were provided. Please provide only one.")
    }
    quant.level.options <- c("auto", "gene", "transcript", "tx")
    cmdopts[['aggregate-level']] %<>% tolower %>% match_arg(choices = quant.level.options, argname = "--aggregate-level", ignore.case = TRUE)
    if (cmdopts[['aggregate-level']] == "auto") {
        cmdopts[['aggregate-level']] = ifelse(length(provided.annot.opts) == 1, "gene", "transcript")
    }
    ## "tx" is an undocumented shortcut for "transcript", for
    ## consistency with tximport, TxDb, etc.
    if (cmdopts[['aggregate-level']] == "tx") {
        cmdopts[['aggregate-level']] <- "transcript"
    }
    assert_that(cmdopts[['aggregate-level']] %in% c("gene", "transcript"))
    if (cmdopts[['aggregate-level']] == "gene" && length(provided.annot.opts) < 1) {
        stop("Gene-level quantification was requested but no gene annotations were provided")
    }

    ## Replace dashes with underscores so that all options can easily
    ## be accessed by "$"
    cmdopts %>% setNames(chartr("-", "_", names(.)))
}

## Do argument parsing early so the script exits quickly if arguments are invalid
get_options(commandArgs(TRUE))

library(assertthat)
library(future)
library(magrittr)
library(stringr)

library(S4Vectors)
library(SummarizedExperiment)
library(tximport)

{

    cmdopts <- get_options(commandArgs(TRUE))

    ## cmdopts <- list(
    ##     samplemeta_file = "saved_data/samplemeta-RNASeq.RDS",
    ##     sample_id_column = "SRA_run",
    ##     abundance_file_pattern = "salmon_quant/hg38.analysisSet_knownGene/{SAMPLE}/abundance.h5",
    ##     aggregate_level = "transcript",
    ##     output_file = "saved_data/SummarizedExperiment_rnaseq_transcript_salmon_hg38.analysisSet_knownGene.RDS",
    ##     annotation_txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene")

    cmdopts$help <- NULL

    ## Expand expected_abundance_files into vector
    if ("expected_abundance_files" %in% names(cmdopts)) {
        cmdopts$expected_abundance_files %<>% str_split(",") %>% unlist
    }

    tsmsg("Args:")
    print_var_vector(cmdopts)

    ## Delete the output file if it exists
    suppressWarnings(file.remove(cmdopts$output_file))
    assert_that(!file.exists(cmdopts$output_file))

    tsmsg("Loading sample metadata")
    samplemeta <- read.table.general(cmdopts$samplemeta_file)

    tsmsg("Got metadata for ", nrow(samplemeta), " samples")

    assert_that(cmdopts$sample_id_column %in% colnames(samplemeta))
    assert_that(!anyDuplicated(samplemeta[[cmdopts$sample_id_column]]))

    rownames(samplemeta) <- samplemeta$sample <- samplemeta[[cmdopts$sample_id_column]]

    ## Use glue(), but only allow interpolation of the SAMPLE variable
    ## and no other code
    safe_env <- new.env(parent = emptyenv())
    assign("SAMPLE", samplemeta[[cmdopts$sample_id_column]], envir = safe_env)
    samplemeta$path <- glue(cmdopts$abundance_file_pattern, .envir = safe_env)

    if ("expected_abundance_files" %in% names(cmdopts)) {
        tryCatch({
            assert_that(setequal(samplemeta$path, cmdopts$expected_abundance_files))
            tsmsg("Sample metadata contains all expected abundance files")
        }, error = function(...) {
            unexpected_existing <- setdiff(samplemeta$path, cmdopts$expected_abundance_files)
            expected_but_missing <- setdiff(cmdopts$expected_abundance_files, samplemeta$path)
            if (length(unexpected_existing) > 0) {
                tsmsg(glue("Got unexpected abundance files: {deparse_onestring(unexpected_existing)}"))
            }
            if (length(expected_but_missing) > 0) {
                tsmsg(glue("Didn't find expected abundance files: {deparse_onestring(expected_but_missing)}"))
            }
            stop("Abundance file list was not as expected")
        })
    }

    assert_that(all(file.exists(samplemeta$path)))

    annot <- NULL
    annot_ranges <- NULL
    tx2gene <- NULL
    if (cmdopts$aggregate_level == "gene") {
        if ("annotation_txdb" %in% names(cmdopts)) {
            tsmsg("Reading transcript ranges from TxDb")
            txdb <- get_txdb(cmdopts$annotation_txdb)
            annot_ranges <- transcriptsBy(txdb, "gene")
            tsmsg("Reading transcript gene mappings from TxDb")
            tx2gene <- get_tx2gene_from_txdb(txdb)
        } else if ("genemap_file" %in% names(cmdopts)) {
            tsmsg("Reading transcript gene mappings from genemap file")
            tx2gene <- read_tx2gene_from_genemap(cmdopts$genemap_file)
        } else {
            stop("Need a gene annotation to aggregate at the gene level.")
        }
        if ("gene_info" %in% names(cmdopts)) {
            tsmsg("Reading gene annotations")
            annot <- read_table_general(cmdopts$gene_info, dataframe.class = "DataFrame")
            ## Nonexistent or automatic row names
            if (.row_names_info(annot) <= 0) {
                row.names(annot) <- annot[[1]]
            }
        }
    } else {
        if ("annotation_txdb" %in% names(cmdopts)) {
            tsmsg("Reading transcript ranges from TxDb")
            txdb <- get_txdb(cmdopts$annotation_txdb)
            annot_ranges <- transcripts(txdb)
            names(annot_ranges) <- annot_ranges$tx_name
        }
        if ("transcript_info" %in% names(cmdopts)) {
            tsmsg("Reading transcript annotations")
            annot <- read_table_general(cmdopts$transcript_info, dataframe.class = "DataFrame")
            ## Nonexistent or automatic row names
            if (.row_names_info(annot) <= 0) {
                row.names(annot) <- annot[[1]]
            }
        }
    }

    tsmsg("Reading quantification files")
    txi <- tximport(samplemeta$path, type = "kallisto", txOut = TRUE)
    if (cmdopts$aggregate_level == "gene") {
        txi %<>% summarizeToGene(tx2gene)
    }

    tsmsg("Matching annotation to quantification tables")
    txi_assayNames <- c("counts", "abundance", "length")
    txi_featureNames <- rownames(txi[[txi_assayNames[1]]])
    if (is.null(annot)) {
        annot <- DataFrame(row.names = txi_featureNames)
        annot[[cmdopts$aggregate_level]] <- txi_featureNames
    } else {
        annot %<>% .[txi_featureNames,] %>% set_rownames(txi_featureNames)
    }

    ## If we have ranges (from a TxDb), we make a
    ## RangedSummarizedExperiment, otherwise we make a regular one
    if (!is.null(annot_ranges)) {
        tsmsg("Constructing the RangedSummarizedExperiment object")
        annot_ranges %<>% .[txi_featureNames]
        mcols(annot_ranges) <- as(annot, "DataFrame")
        sexp <- SummarizedExperiment(
            assays = List(txi[txi_assayNames]),
            colData = as(samplemeta, "DataFrame"),
            rowRanges = annot_ranges,
            ## Put non-assay elements of txi into the metadata
            metadata = SimpleList(txi[!names(txi) %in% txi_assayNames]))
    } else {
        tsmsg("Constructing the SummarizedExperiment object")
        sexp <- SummarizedExperiment(
            assays = List(txi[txi_assayNames]),
            colData = as(samplemeta, "DataFrame"),
            rowData = as(annot, "DataFrame"),
            ## Put non-assay elements of txi into the metadata
            metadata = SimpleList(txi[!names(txi) %in% txi_assayNames]))
    }

    tsmsg("Saving SummarizedExperiment")
    save_RDS_or_RDA(sexp, cmdopts$output_file)
    invisible(NULL)
}

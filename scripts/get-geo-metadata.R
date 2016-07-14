#!/usr/bin/env Rscript

getScriptPath <- function() {
    argv <-commandArgs()
    dir <- na.omit(stringr::str_match(argv, "^--file=(.*)$")[,2])[1]
    if (!is.na(dir) && !is.null(dir))
        return(dir)
}
setwd(file.path(dirname(getScriptPath()), ".."))

library(GEOquery)
library(SRAdb)
library(stringr)
library(magrittr)
library(dplyr)
library(assertthat)
library(lubridate)

tsmsg <- function(...) {
    message(base::date(), ": ", ...)
}

fac2char <- function(df) {
    df[sapply(df, class) == "factor"] %<>% lapply(as.character)
    df
}

getGEO <- function(...) {
    GEOquery::getGEO(..., destdir="saved_data")
}

sra_con <- {
    sqlfile <- file.path("saved_data", "SRAmetadb.sqlite")
    if(!file.exists(sqlfile)) {
        getSRAdbFile(destdir=dirname(sqlfile), destfile=str_c(basename(sqlfile), ".gz"))
    }
    stopifnot(file.exists(sqlfile))
    dbConnect(SQLite(),sqlfile)
}

get_channel_meta <- function (pdata) {
    num.channels <- max(as.numeric(pdata$channel_count))
    assert_that(num.channels >= 1)
    lapply(seq_len(num.channels), function(chan) {
        chan_regexp <- sprintf("ch%s(\\.\\d+)?$", chan)
        pdata[str_detect(colnames(pdata), chan_regexp)]
    })
}

get_characteristics <- function(pdata) {
    pdata %<>% .[str_detect(colnames(.), "^characteristics_ch")]
    cnames <- lapply(pdata, . %>% str_extract("^[^:]+") %>% unique)
    assert_that(all(lengths(cnames) == 1))
    cnames %<>% unlist %>% str_replace_all("\\s+", "_")
    colnames(pdata) <- cnames
    pdata[] <- lapply(pdata, . %>% str_replace("^[^:]+: ", ""))
    pdata
}

get_relations <- function(pdata) {
    pdata %<>% .[str_detect(colnames(.), "^relation")]
    cnames <- lapply(pdata, . %>% str_extract("^[^:]+") %>% unique)
    assert_that(all(lengths(cnames) == 1))
    cnames %<>% unlist %>% str_replace_all("\\s+", "_")
    colnames(pdata) <- cnames
    pdata[] <- lapply(pdata, . %>% str_replace("^[^:]+: ", ""))
    pdata
}

get_samplemeta_from_geo_pdata <- function(eset) {
    pdata <- pData(eset)
    front_data <- pdata %>% select(SampleName=title, geo_accession)
    back_data  <- pdata %>%
        select(submission_date, last_update_date)
    samplemeta <- get_characteristics(get_channel_meta(pdata)[[1]])
    relations <- get_relations(pdata)
    cbind(front_data, samplemeta, relations, back_data) %>% fac2char
}

## Only performs mutations if all of given names are present in .data
mutate_if_present <- function(.data, names, ...) {
    if (all(names %in% base::names(.data))) {
        mutate_(.data, .dots = lazyeval::lazy_dots(...))
    } else {
        .data
    }
}

geo_ids <- c(RNASeq="GSE73213",
             ChIPSeq="GSE73212")

tsmsg("Loading GEO data")
esets <- lapply(geo_ids, . %>% getGEO %>% .[[1]])

tsmsg("Extracting sample metadata")
samplemeta <- lapply(esets, function(eset) {
    eset %>% get_samplemeta_from_geo_pdata %>%
        ## Extract just the IDs from the URLs
        mutate(BioSample=str_replace(BioSample, "^\\Qhttp://www.ncbi.nlm.nih.gov/biosample/", ""),
               SRA=str_replace(SRA, "^\\Qhttp://www.ncbi.nlm.nih.gov/sra?term=", "")) %>%
        ## Get all relevant SRA IDs
        cbind(sraConvert(.$SRA, out_type=as.list(args(sraConvert))$out_type, sra_con) %>% setNames(str_c("SRA_", names(.)))) %>%
        ## Get the basename of the download URL
        mutate(SRA=NULL,
               SRA_file=listSRAfile (in_acc = SRA_run, sra_con = sra_con, fileType ="sra", srcType="fasp") %$%
                   fasp[match(SRA_run, run)] %>% basename) %>%
        ## Fix column types and ensure that non-numeric variables
        ## never start with numbers
        mutate(cell_type=factor(cell_type, levels=c("Naive", "Memory")),
               activated=as.logical(activated),
               days_after_activation=as.numeric(days_after_activation),
               donor_id=sprintf("D%s", donor_id),
               submission_date=mdy(submission_date),
               last_update_date=mdy(last_update_date)) %>%
        ## Only in RNA-seq
        mutate_if_present(
            "technical_batch",
            technical_batch=sprintf("B%s", technical_batch)) %>%
        ## Only in ChIP-seq
        mutate_if_present(
            "chip_antibody",
            chip_antibody=factor(chip_antibody,
                                 levels=c("input", "H3K4me2", "H3K4me3", "H3K27me3")))
})
assert_that(all(!is.na(unlist(samplemeta))))

tsmsg("Saving sample metadata")
for (i in names(samplemeta)) {
    saveRDS(samplemeta[[i]], file.path("saved_data", sprintf("samplemeta-%s.RDS", i)))
}
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

library(GEOquery)
library(SRAdb)
library(stringr)
library(glue)
library(magrittr)
library(dplyr)
library(assertthat)
library(lubridate)
library(rex)

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
    assert_that(file.exists(sqlfile))
    dbConnect(SQLite(),sqlfile)
}

get_channel_meta <- function (pdata) {
    num.channels <- max(as.numeric(pdata$channel_count))
    assert_that(num.channels >= 1)
    lapply(seq_len(num.channels), function(chan) {
        chan_regexp <- rex(
            "ch", chan,
            maybe(
                ".",
                one_or_more(digit)
            ),
            end
        )
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
        mutate(.data, !!!quos(...))
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
        rename(SRA_experiment=SRA) %>%
        ## Get all relevant SRA IDs
        inner_join(sraConvert(.$SRA, out_type=as.list(args(sraConvert))$out_type, sra_con) %>% setNames(str_c("SRA_", names(.)))) %>%
        ## Fix column types and ensure that non-numeric variables
        ## never start with numbers
        mutate(cell_type=factor(cell_type, levels=c("Naive", "Memory")),
               activated=as.logical(activated),
               days_after_activation=as.numeric(days_after_activation),
               donor_id=glue("D{donor_id}"),
               submission_date=mdy(submission_date),
               last_update_date=mdy(last_update_date)) %>%
        ## Only in RNA-seq
        mutate_if_present(
            "technical_batch",
            technical_batch=glue("B{technical_batch}"),
            libType=ifelse(technical_batch == "B1", "SF", "SR")) %>%
        ## Only in ChIP-seq
        mutate_if_present(
            "chip_antibody",
            chip_antibody=factor(chip_antibody,
                                 levels=c("input", "H3K4me2", "H3K4me3", "H3K27me3")))
})
assert_that(all(!is.na(unlist(samplemeta))))

tsmsg("Saving sample metadata")
for (i in names(samplemeta)) {
    saveRDS(samplemeta[[i]], file.path("saved_data", glue("samplemeta-{i}.RDS")))
}

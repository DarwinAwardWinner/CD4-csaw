#!/usr/bin/env Rscript

library(GEOquery)
library(SRAdb)
library(stringr)
library(magrittr)
library(dplyr)
library(assertthat)

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

auto_fac2char <- function(df) {
    for (i in names(df)) {
        if (is.factor(df[[i]]) && length(unique(df[[i]])) == nrow(df)) {
            df[[i]] %<>% as.character
        }
    }
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
        chan_regexp <- sprintf("ch%s\\.\\d+$", chan)
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
        select(status, submission_date, last_update_date)
    samplemeta <- get_characteristics(get_channel_meta(pdata)[[1]])
    relations <- get_relations(pdata)
    cbind(front_data, samplemeta, relations, back_data) %>% auto_fac2char
}

first.accessible <- function(paths, mode=0) {
    for (path in paths) {
        if (file.access(path, mode) == 0) {
            return(path)
        }
    }
    return(NA_character_)
}

## Check if we can use ascp download method, or fall back to FTP if
## not
getSRAfile <- function(...) {
    SRAdb::getSRAfile(..., srcType="ftp")
}
ascp.path <- first.accessible(c(Sys.which("ascp"), path.expand("~/.aspera/connect/bin/ascp")), mode=1)
if (!is.na(ascp.path)) {
    ascp.key.file <- normalizePath(first.accessible(
        file.path(dirname(ascp.path), c("asperaweb_id_dsa.openssh", "../etc/asperaweb_id_dsa.openssh"))))
    if (!is.na(ascp.key.file)) {
        cmd <- sprintf(
            "%s -T -k1 -i %s",
            shQuote(ascp.path),
            shQuote(ascp.key.file))
        ## cmd <- sprintf(
        ##     "%s -q -k2",
        ##     shQuote(ascp.path),
        ##     shQuote(ascp.key.file))
        getSRAfile <- function(...) {
            SRAdb::getSRAfile(..., srcType="fasp",
                              ascpCMD = cmd)
        }
    }
}
sra_dir <- "sra_files"

rnaseq_geo <- "GSE73213"
chipseq_geo <- "GSE73212"

rnaseq_eset <- getGEO(rnaseq_geo)[[1]]
chipseq_eset <- getGEO(chipseq_geo)[[1]]

rnaseq_samplemeta <- get_samplemeta_from_geo_pdata(rnaseq_eset) %>%
    ## Extract just the IDs from the URLs
    mutate(BioSample=str_replace(BioSample, "^\\Qhttp://www.ncbi.nlm.nih.gov/biosample/", ""),
           SRA=str_replace(SRA, "^\\Qhttp://www.ncbi.nlm.nih.gov/sra?term=", "")) %>%
    ## Get all relevant SRA IDs
    cbind(sraConvert(.$SRA, out_type=as.list(args(sraConvert))$out_type, sra_con) %>% setNames(str_c("SRA_", names(.)))) %>%
    ## Get the basename of the download URL
    mutate(SRA=NULL,
           download_path=listSRAfile (in_acc = SRA_run, sra_con = sra_con, fileType ="sra", srcType="fasp") %$%
               fasp[match(SRA_run, run)] %>% basename %>% file.path(sra_dir, .))
chipseq_samplemeta <- get_samplemeta_from_geo_pdata(chipseq_eset) %>%
    ## Extract just the IDs from the URLs
    mutate(BioSample=str_replace(BioSample, "^\\Qhttp://www.ncbi.nlm.nih.gov/biosample/", ""),
           SRA=str_replace(SRA, "^\\Qhttp://www.ncbi.nlm.nih.gov/sra?term=", "")) %>%
    ## Get all relevant SRA IDs
    cbind(sraConvert(.$SRA, out_type=as.list(args(sraConvert))$out_type, sra_con) %>% setNames(str_c("SRA_", names(.)))) %>%
    ## Get the basename of the download URL
    mutate(SRA=NULL,
           download_path=listSRAfile (in_acc = SRA_run, sra_con = sra_con, fileType ="sra", srcType="fasp") %$%
               fasp[match(SRA_run, run)] %>% basename %>% file.path(sra_dir, .))

getSRAfile(c(rnaseq_samplemeta$SRA_run, chipseq_samplemeta$SRA_run),
           sra_con, destDir=sra_dir, makeDirectory = TRUE)
assert_that(all(file.exists(rnaseq_samplemeta$download_path)))
assert_that(all(file.exists(chipseq_samplemeta$download_path)))

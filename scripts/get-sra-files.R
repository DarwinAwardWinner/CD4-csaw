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
library(glue)
library(magrittr)
library(dplyr)
library(assertthat)

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

first.accessible <- function(paths, mode=0) {
    for (path in paths) {
        if (file.access(path, mode) == 0) {
            return(path)
        }
    }
    return(NA_character_)
}

sra_con <- {
    sqlfile <- file.path("saved_data", "SRAmetadb.sqlite")
    if(!file.exists(sqlfile)) {
        getSRAdbFile(destdir=dirname(sqlfile), destfile=str_c(basename(sqlfile), ".gz"))
    }
    stopifnot(file.exists(sqlfile))
    dbConnect(SQLite(),sqlfile)
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
        cmd <- glue(
            "{ascp} -T -k1 -i {keyfile}",
            ascp=shQuote(ascp.path),
            keyfile=shQuote(ascp.key.file))
        getSRAfile <- function(...) {
            SRAdb::getSRAfile(..., srcType="fasp",
                              ascpCMD = cmd)
        }
    }
}
sra_dir <- "sra_files"

samplemeta.files <- list.files(path="saved_data", pattern="^samplemeta-(.*)\\.RDS$", full.names = TRUE)
samplemeta <- lapply(samplemeta.files, readRDS) %>% do.call(what=rbind.fill)

getSRAfile(samplemeta$SRA_run,
           sra_con, destDir=sra_dir, makeDirectory = TRUE)
assert_that(all(file.exists(file.path(sra_dir, samplemeta$SRA_file))))

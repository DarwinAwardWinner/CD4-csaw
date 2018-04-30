#!/usr/bin/env Rscript

suppressPackageStartupMessages(suppressMessages({
    library(stringr)
    library(glue)
    library(SRAdb)
    library(assertthat)
    library(rctutils)

    sra_con <- {
        sqlfile <- here("saved_data", "SRAmetadb.sqlite")
        if(!file.exists(sqlfile)) {
            getSRAdbFile(destdir=dirname(sqlfile),
                         destfile=str_c(basename(sqlfile), ".gz"))
        }
        stopifnot(file.exists(sqlfile))
        dbConnect(SQLite(),sqlfile)
    }

    ## Check if we can use ascp download method, otherwise fall back
    ## to FTP
    getSRAfile <- function(...) {
        SRAdb::getSRAfile(..., srcType="ftp")
    }

    ascp.path <- first_accessible_path(
        c(Sys.which("ascp"),
          path.expand("~/.aspera/connect/bin/ascp")), mode=1)
    if (!is.na(ascp.path)) {
        ascp.key.file <- normalizePath(first.accessible(
            file.path(dirname(ascp.path),
                      c("asperaweb_id_dsa.openssh",
                        "../etc/asperaweb_id_dsa.openssh"))))
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

}))

{
    sra_dir <- "sra_files"

    sra_runs <- commandArgs(TRUE)
    assert_that(all(str_detect(sra_runs, "^SRR")))

    getSRAfile(sra_runs, sra_con, destDir=sra_dir, makeDirectory = TRUE)

    expected_files <- here(sra_dir, str_c(sra_runs, ".sra"))
    assert_that(all(file.exists(expected_files)))
    invisible(NULL)                     # Avoid output on console
}


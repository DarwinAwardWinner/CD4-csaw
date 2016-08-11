#!/usr/bin/env Rscript

library(assertthat)
library(SRAdb)

sra_con <- {
    sqlfile <- file.path("saved_data", "SRAmetadb.sqlite")
    if(!file.exists(sqlfile)) {
        getSRAdbFile(destdir=dirname(sqlfile), destfile=str_c(basename(sqlfile), ".gz"))
    }
    assert_that(file.exists(sqlfile))
}

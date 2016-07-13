#!/usr/bin/env Rscript

options(download.file.method="wget", download.file.extra="--continue --tries=2")

## Remember to cite: https://sites.google.com/site/anshulkundaje/projects/blacklists
urls <- c(
    "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz",
    "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz")

for (url in urls) {
    download.file(url, destfile = file.path("saved_data", basename(url)))
}

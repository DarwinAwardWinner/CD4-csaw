#!/bin/bash

perl -lane 'print $1 if m/library\((\w+)\)/' *.R | sort -u | R --no-save <<EOF
source("https://bioconductor.org/biocLite.R")
needed_pkgs <- setdiff(readLines(stdin()), c("", rownames(installed.packages())))
if (length(needed_pkgs)) {
    message("Installing the following packages:")
    dput(pkgs)
    biocLite(needed_pkgs, suppressUpdates=TRUE)
} else {
    message("All needed packages are already installed.")
}
EOF

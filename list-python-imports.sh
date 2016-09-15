#!/bin/bash
find . -maxdepth 2 -name '*.py' -o -name '*Snakefile' | parallel cat | \
    perl -lane '
s/#.*//;
if (m/from \s+ (\S+) \s+ import \s+ (\S+)/x) {
    print "$1";
} elsif (m/^import \s+ (\S+)/x) {
    print "$1";
}' | \
    sort -u | \
    grep -v tool_versions

#!/bin/bash
find . -maxdepth 2 -name '*.R' | parallel cat | \
    perl -lane 'print $1 if m/library\(("?[^\s,"]+)"?\)/' | \
    sort -u

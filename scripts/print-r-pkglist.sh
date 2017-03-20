#!/bin/bash

scripts_dir="$(dirname "$0")"
cat "$scripts_dir"/*.R{,md} | \
    grep -v 'character\.only' | \
    perl -lane 'print $1 if m{library\(([A-Za-z0-9._]+)}' | \
    sort -u

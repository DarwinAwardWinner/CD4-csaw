#!/usr/bin/env python

import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def read_fastq(file):
    while True:
        record = (file.readline(), file.readline(),
                  file.readline(), file.readline(),)
        if len(record[0]) == 0:
            break
        if len(record[3]) == 0:
            raise IOError("Invalid FASTQ: Input did not contain a multiple of 4 lines")
        if not record[0].startswith("@"):
            raise IOError("Invalid FASTQ: Read's first line did not start with '@'")
        if not record[2].startswith("+"):
            raise IOError("Invalid FASTQ: Read's third line did not start with '+'")
        title = record[0][1:-1]
        seq = record[1][:-1]
        qual = record[3][:-1]
        yield (title, seq, qual)

for title, seq, qual in read_fastq(sys.stdin) :
    # Replace missing quality with all minimum quality
    if qual == "*":
        qual = "!" * len(seq)
    sys.stdout.write("@{title}\n{seq}\n+{title}\n{qual}\n".format(
        title=title, seq=seq, qual=qual))

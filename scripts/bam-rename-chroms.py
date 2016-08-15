#!/usr/bin/env python

import logging
import plac
import pysam
import sys

from copy import deepcopy

def tqdm_fake(object, *args, **kwargs):
    return object

def read_chrom_mapping(filename):
    mapping = {}
    with open(filename, 'r') as infile:
        for line in infile:
            (from_id, to_id) = line.split('\t')
            to_id = to_id.strip()
            if to_id:
                mapping[from_id] = to_id
    return mapping

def fixup_header(header, chrom_mapping, _do_copy=True):
    if _do_copy:
        header = deepcopy(header)
    sq = header['SQ']
    for sqline in sq:
        try:
            sqline['SN'] = chrom_mapping[sqline['SN']]
        except KeyError:
            # Need to leave unmatched SN values unchanged so that
            # reads mapped to those references will not become invalid
            pass
    return header

@plac.annotations(
    # arg=(helptext, kind, abbrev, type, choices, metavar)
    chrom_mapping_file=('File containing chromsome mapping', 'positional'),
    input_file=('Input bam file (stdin by default)', 'positional'),
    output_file=('Output bam file (stdout by default)', 'positional'),
    quiet=('Do not print informational messages.', 'flag', 'q'),
    verbose=('Print debug messages that are probably only useful if something is going wrong.', 'flag', 'v'),
    )
def main(chrom_mapping_file, input_file="-", output_file="-",
         quiet=False, verbose=False,
         ):
    '''Rename reference sequences in BAM file.'''
    if quiet:
        logging.basicConfig(level=logging.WARN)
    elif verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    if quiet:
        tqdm = tqdm_fake
    else:
        try:
            from tqdm import tqdm
        except ImportError:
            tqdm = tqdm_fake

    logging.info("Reading chromosome mapping file")
    cmap = read_chrom_mapping(chrom_mapping_file)
    with pysam.Samfile(input_file) as bam_input:
        logging.info("Modifying header for output file")
        output_header = fixup_header(bam_input.header, cmap)
        with pysam.AlignmentFile(output_file, 'wb', header=output_header) as bam_output:
            logging.info("Copying reads into output file")
            for ar in tqdm(bam_input.fetch(until_eof=True), desc="Copying"):
                bam_output.write(ar)

if __name__ == '__main__':
    plac.call(main)

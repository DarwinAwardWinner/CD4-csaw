import csv
import os
import os.path
import rpy2.rinterface
import re
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

from atomicwrites import atomic_write, AtomicWriter
from math import log10, ceil
from itertools import product, count, chain
from subprocess import check_call, Popen, PIPE, CalledProcessError, list2cmdline
from rpy2 import robjects
from rpy2.robjects import r, pandas2ri
from rpy2.robjects import globalenv as r_env
from warnings import warn
from tempfile import TemporaryDirectory

from snakemake.io import expand
from snakemake.utils import min_version
min_version('3.7.1')

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

FTP = FTPRemoteProvider()
HTTP = HTTPRemoteProvider()

pandas2ri.activate()
rpy2.rinterface.set_writeconsole_warnerror(lambda x: sys.stderr.write(x))

# Commands to compress and decompress a variety of fastq compression
# methods
fastq_compression_cmds = {
    'fq.gz': {
        'compress': ['gzip', '-c'],
        'decompress': ['gzip', '-d', '-c'],
    },
    'fq.bz2': {
        'compress': ['bzip2', '-z', '-c'],
        'decompress': ['bzip2', '-d', '-c'],
    },
    'fq.qp': {
        'compress': ['quip', '-i', 'fastq', '-o', 'quip', '-c' ],
        'decompress': ['quip', '-i', 'quip', '-o', 'fastq', '-c' ],
    },
    # i.e. no compression
    'fq': {
        'compress': ['cat'],
        'decompress': ['cat'],
    },
}

def Popen_pipeline(cmds, stdin=None, stdout=None, *args, **kwargs):
    '''Popen a pipeline of several commands.

    Returns a list with all the process objects returned by Popen.

    Each command's stdout becomes the next command's stdin. The stdin
    argument becomes the stdin of the first command, while the stdout
    argument becomes the stdout of the last command. All other
    arguments are passed to every invocation of Popen(), so ensure
    that they make sense in that context.

    '''
    cmds = list(cmds)
    if len(cmds) == 0:
        raise ValueError('Cannot run a pipeline with zero commands')
    if len(cmds) == 1:
        return [Popen(cmds[0], stdin=stdin, stdout=stdout, *args, **kwargs)]
    else:
        first_cmd = cmds[0]
        last_cmd = cmds[-1]
        middle_cmds = cmds[1:-1]
        proclist = [Popen(first_cmd, stdin=stdin, stdout=PIPE, *args, **kwargs)]
        for cmd in middle_cmds:
            prev_cmd_stdout = proclist[-1].stdout
            proclist.append(Popen(cmd, stdin=prev_cmd_stdout, stdout=PIPE, *args, **kwargs))
        prev_cmd_stdout = proclist[-1].stdout
        proclist.append(Popen(last_cmd, stdin=prev_cmd_stdout, stdout=stdout, *args, **kwargs))
        return proclist

def wait_for_subprocs(proclist, expected_exitcodes=0, wait_for_all=True):
    '''Wait for a list of subprocesses to exit.

    If any subprocess returns an exit code that is not in
    expected_exitcodes (default: only 0), a CalledProcessError is
    raised (as if subprocess.check_call was run). Even if an exception
    is raised, this function will still wait for all remaining
    processes to finish unless wait_for_all is False.

    '''
    # Allow supplying a single number or a list of numbers
    try:
        0 in expected_exitcodes
    except TypeError:
        expected_exitcodes = [ expected_exitcodes ]
    if wait_for_all:
        for proc in proclist:
            proc.wait()
    for proc in proclist:
        exitcode = proc.wait()
        if exitcode not in expected_exitcodes:
            raise CalledProcessError(exitcode, proc.args)

def read_R_dataframe(rdsfile):
    '''Read an R data frame stored in an RDS file.

    Returns the result as a Pandas DataFrame.

    '''
    readRDS = robjects.r('readRDS')
    df = readRDS((robjects.StrVector([rdsfile])))
    return(pandas2ri.ri2py(df))

# TODO: Does a function like this already exist in pandas?
def dfselect(dframe, what=None, where=None, **where_kwargs):
    '''Filter DataFrame and select columns.

    First, the DataFrame is filtered according to the 'where'
    arguments. Each key (or keyword argument) in 'where' should be a
    column name, and its value should be the allowed value or values
    for that column. Rows will be selected if they match all the
    requirements in 'where'. Each filter can also be a callable, which
    should take a single value and return True for including that
    value and False for excluding it.

    Second, if 'what' is not None, it will be used to select one or
    more columns. As normal, using a single string returns a Series
    object, while using a list of strings returns a DataFrame with a
    subset of columns selected.

    '''
    if where is None:
        where = dict()
    where.update(where_kwargs)
    if where:
        selected = pd.Series(True, index=dframe.index)
        for (colname, allowed_vals) in where.items():
            if callable(allowed_vals):
                allow_func = allowed_vals
                allowed_rows = [allow_func(x) for x in dframe[colname]]
            else:
                try:
                    allowed_vals = pd.Series(allowed_vals)
                except TypeError:
                    allowed_vals = pd.Series(list(allowed_vals))
                allowed_rows = dframe[colname].isin(allowed_vals)
            selected &= allowed_rows
        dframe = dframe[selected]
    if what is None:
        return dframe
    else:
        return dframe[what]

def df_cartesian_product(*dfs):
    '''Return the cartesian product of 2 or more DataFrames.

    None of dfs should share a column name with any other, or else the
    column names will have arbitrary suffixes.

    '''
    if len(dfs) == 0:
        raise ValueError("Cannot generate empty Cartesian product")
    elif len(dfs) == 1:
        return dfs[0]
    else:
        all_colnames = list(chain.from_iterable(df.columns for df in dfs))
        # Get an unused column name
        merge_key = 'key_to_merge_on_'
        while merge_key in all_colnames:
            merge_key += 'xxxxx'
        merged_df = dfs[0].copy()
        merged_df[merge_key] = 1
        for next_df in (df.copy() for df in dfs[1:]):
            next_df[merge_key] = 1
            merged_df = merged_df.merge(next_df, on=merge_key)
        for i in merged_df:
            if i.startswith(merge_key):
                merged_df.drop(i, 1, inplace=True)
        return merged_df

def recycled(it, length=None):
    '''Recycle iterable to specified length.

    Items are returned in round-robin order until the specified length
    is reached. If length is None, the iterable is recycled
    indefinitely.

    Cannot handle infinite-length iterators as inputs, since it needs
    to enumerate all the elements of an iterable in order to recycle
    them.

    '''
    if length is None:
        x = list(it)
        while True:
            yield from iter(x)
    else:
        need_recycling = True
        try:
            if len(it) >= length:
                need_recycling = False
        except TypeError:
            pass
        if need_recycling:
            it = recycled(it)
            i = 0
            while i < length:
                i += 1
                yield next(it)
        else:
            yield from iter(it[:length])
    raise StopIteration

def zip_recycled(*args, length=None):
    '''Like zip(), but recycles all iterables to the specified length.

    If length is None,

    Cannot handle infinite-length iterators as inputs, since it needs
    to enumerate all the elements of an iterable in order to recycle
    them.

    '''
    return zip(*(recycled(arg, length) for arg in args))

def zip_longest_recycled(*args, warn_on_uneven=True):
    '''Like itertools.zip_longest(), but recycles shorter iterables.

    Unless kwarg warn_on_mismatch is set to False, a warning will be
    raised if all the iterable lengths do not divide evenly into the
    length of the longest iterable.

    '''
    args = [list(arg) for arg in args]
    maxlen = max(map(len, args))
    if warn_on_uneven:
        if max(maxlen % len(arg) for arg in args) > 0:
            warn("Longest iterable's length is not a multiple of shorter.")
    return zip_recycled(*args, length=maxlen)

def list_salmon_output_files(outdirs, alignment=False):
    file_list = [
        'aux_info/bootstrap/bootstraps.gz',
        'aux_info/bootstrap/names.tsv.gz',
        'aux_info/eq_classes.txt',
        'aux_info/exp3_seq.gz',
        'aux_info/exp5_seq.gz',
        'aux_info/expected_bias.gz',
        'aux_info/fld.gz',
        'aux_info/meta_info.json',
        'aux_info/obs3_seq.gz',
        'aux_info/obs5_seq.gz',
        'aux_info/observed_bias.gz',
        'aux_info/observed_bias_3p.gz',
        'cmd_info.json',
        'quant.genes.sf',
        'quant.sf',
    ]
    if alignment:
        file_list += ['logs/salmon.log',]
    else:
        file_list += ['libParams/flenDist.txt', 'logs/salmon_quant.log',]
    if isinstance(outdirs, str):
        outdirs = [outdirs]
    return [ os.path.join(od, f) for od in outdirs for f in file_list ]

def list_kallisto_output_files(outdir):
    file_list = [
        'abundance.h5', 'abundance.tsv', 'run_info.json',
    ]
    return [ os.path.join(outdir, f) for f in file_list ]

def list_macs_callpeak_output_files(dirname):
    file_list = [
        'peaks.narrowPeak',
        'peaks.xls',
        'summits.bed',
        # 'control_lambda.bdg',
        # 'treat_pileup.bdg',
    ]
    return [ os.path.join(dirname, fname) for fname in file_list ]

def read_narrowpeak(infile):
    peaks = pd.DataFrame.from_csv(infile, header=None, sep='\t', index_col=None)
    peaks.columns = ('chr', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'summit')
    return peaks

def write_narrowpeak(peaks, outfile):
    peaks.to_csv(outfile, sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)

def pick_top_peaks(infile, outfile, by='score', ascending=False, number=150000, *args, **kwargs):
    '''Copy the top N peaks from infile to outfile.

    Peaks are read from 'infile', sorted, and then the top 'number'
    are written to 'outfile'. Peaks are read and written in narrowPeak
    format. Arguments 'by', 'ascending', and any other arguments are
    passed to pandas.DataFrame.sort_values to determine how to sort.

    Reasonable values for 'by', include 'score', 'signalValue', and
    'pValue'. Typically you want 'ascending=False' for all of these,
    including 'pValue', which is typically on a negative log10 scale,
    so higher values are more significant.

    '''
    peaks = read_narrowpeak(infile)
    peaks.sort_values(by=by, axis=0, ascending=ascending, inplace=True, *args, **kwargs)
    write_narrowpeak(peaks.head(number), outfile)

def call_R_external(f, *args, **kwargs):
    arglist_string = r['paste'](r['deparse'](r['list'](*args, **kwargs), backtick=True, nlines=-1), collapse=" ")[0]
    rcode = "do.call(%s, %s)" % (f, arglist_string)
    check_call(['Rscript', '-e', rcode])

def dict_to_R_named_list(d):
    return r['list'](**d)

rmd_default_formats = {
    # The notebook format for html has additional bells & whistles
    # that are useful even outside the context of interactive
    # operation, so we use that format for html output.
    'html': 'html_notebook',
    'pdf': 'pdf_document',
}

def rmd_render(input, output_file, output_format=None, **kwargs):
    if output_format is None:
        if output_file is not None:
            ext = os.path.splitext(output_file)[1][1:]
            try:
                output_format = rmd_default_formats[ext]
            # If no specific output format is specified, just append
            # "_document" and hope that works
            except KeyError:
                if ext == '':
                    raise ValueError("Cannot determine output format from file name.")
                else:
                    output_format = ext + '_document'
        # Output file will not be saved, so just pick something
        # arbitrarily.
        else:
            output_format = 'html_document'
    arg_converters = {
        'params': dict_to_R_named_list,
        'output_options': dict_to_R_named_list,
    }
    for (k, convfun) in arg_converters.items():
        if k in kwargs:
            kwargs[k] = convfun(kwargs[k])
    with TemporaryDirectory() as tmpdir:
        # The output name must not have a file extension because of
        # https://github.com/rstudio/rmarkdown/issues/1180
        tmp_output_file = os.path.join(tmpdir, "output_file")
        call_R_external('rmarkdown::render', input=input, output_file=tmp_output_file, output_format=output_format, **kwargs)
        if output_file is not None:
            shutil.move(tmp_output_file, output_file)

def rmd_run_without_rendering(input, **kwargs):
    '''Run the code in an Rmd file but don't produce a report.'''
    rmd_render(input, output_file=None, output_format=None, **kwargs)

# Run a separate Snakemake workflow (if needed) to fetch the sample
# metadata, which must be avilable before evaluating the rules below.
# Without this two-step workflow, the below rules would involve quite
# complex use of multiple dynamic() inputs and outputs.
try:
    rnaseq_samplemeta = read_R_dataframe('saved_data/samplemeta-RNASeq.RDS')
    chipseq_samplemeta = read_R_dataframe('saved_data/samplemeta-ChIPSeq.RDS')
except Exception:
    from processify import processify
    from snakemake import snakemake
    snakemake = processify(snakemake)
    result = snakemake(
        'pre.Snakefile',
        targets=['saved_data/samplemeta-RNASeq.RDS', 'saved_data/samplemeta-ChIPSeq.RDS'],
        lock=False,
        quiet=True)
    if not result:
        raise Exception('Could not retrieve experiment metadata from GEO')
    rnaseq_samplemeta = read_R_dataframe('saved_data/samplemeta-RNASeq.RDS')
    chipseq_samplemeta = read_R_dataframe('saved_data/samplemeta-ChIPSeq.RDS')

rnaseq_samplemeta['time_point'] = rnaseq_samplemeta['days_after_activation'].apply(lambda x: 'Day{:.0f}'.format(x))
chipseq_samplemeta['time_point'] = chipseq_samplemeta['days_after_activation'].apply(lambda x: 'Day{:.0f}'.format(x))

promoter_radii = {
    'H3K4me3': '1kbp',
    'H3K4me2': '1kbp',
    'H3K27me3': '2.5kbp',
    'input': None,
}

chipseq_samplemeta['promoter_radius'] = chipseq_samplemeta['chip_antibody'].apply(lambda x: promoter_radii[x])

rnaseq_sample_libtypes = dict(zip(rnaseq_samplemeta['SRA_run'], rnaseq_samplemeta['libType']))

rnaseq_star_outdirs = [
    'rnaseq_star_hg38.analysisSet_knownGene',
    'rnaseq_star_hg38.analysisSet_ensembl.85',
]
rnaseq_hisat_outdir = 'rnaseq_hisat2_grch38_snp_tran'

aligned_rnaseq_star_bam_files = expand(
    'aligned/{dirname}/{samp}/Aligned.sortedByCoord.out.bam',
    dirname=rnaseq_star_outdirs, samp=rnaseq_samplemeta['SRA_run'])

aligned_rnaseq_hisat_bam_files = expand(
    'aligned/{dirname}/{samp}/Aligned.bam',
    dirname=rnaseq_hisat_outdir, samp=rnaseq_samplemeta['SRA_run'])

aligned_rnaseq_bam_files = aligned_rnaseq_star_bam_files + aligned_rnaseq_hisat_bam_files
aligned_rnaseq_bai_files = [ bam + '.bai' for bam in aligned_rnaseq_bam_files ]

chipseq_samplemeta_noinput = dfselect(chipseq_samplemeta, chip_antibody=lambda x: x != 'input')

promoter_radius_table = dfselect(chipseq_samplemeta_noinput, what=['chip_antibody', 'promoter_radius']).drop_duplicates()

def all_pairs(v, *, include_equal=False, include_reverse=False):
    '''Return iterator over 2-tuples of input elements.

    Tuples are yielded in arbitrary order. If 'include_equal' is True,
    tuples with both elements the same will be allowed. If
    'include_reverse' is True, tuples with the same elements in
    opposite order (e.g. (1,2) and (2,1)) will both be yielded.

    '''
    try:
        len(v)
    except Exception:
        v = list(v)
    for i,j in product(range(len(v)), repeat=2):
        if not include_equal and i == j:
            continue
        if not include_reverse and i > j:
            continue
        yield (v[i], v[j])

idr_sample_pairs = chipseq_samplemeta_noinput.\
                   groupby(['chip_antibody', 'cell_type', 'time_point']).\
                   apply(lambda x: pd.DataFrame.from_items(zip(count(), list(all_pairs(list(x['donor_id'])))),
                                                      orient='index', columns=['donorA', 'donorB'])).\
                   reset_index(level=3, drop=True).reset_index()

aligned_chipseq_input_bam_files = expand(
    'aligned/chipseq_bowtie2_hg38.analysisSet/{SRA_run}/Aligned.bam',
    SRA_run=list(chipseq_samplemeta['SRA_run'][chipseq_samplemeta['chip_antibody'] == 'input']))
aligned_chipseq_bam_files = expand(
    'aligned/chipseq_bowtie2_hg38.analysisSet/{SRA_run}/Aligned.bam',
    SRA_run=list(chipseq_samplemeta['SRA_run'][chipseq_samplemeta['chip_antibody'] != 'input']))

include: 'rulegraph.Snakefile'
include: 'tool_versions.py'
include: 'mem_requirements.py'

shell.executable('bash')

# TODO: Put this in a config file
HG38_REF_PATH='~/references/hg38'

subworkflow hg38_ref:
    workdir: os.path.expanduser(HG38_REF_PATH)

targets = {
    'rnaseq_counts': expand(
        'saved_data/SummarizedExperiment_rnaseq_{aligner_and_genome}_{transcriptome}.RDS',
        aligner_and_genome=['star_hg38.analysisSet', 'hisat2_grch38_snp_tran'],
        transcriptome=['ensembl.85','knownGene']),
    'rnaseq_quant': expand(
        'saved_data/SummarizedExperiment_rnaseq_{quantifier}_{genome}_{transcriptome}.RDS',
        quantifier=['kallisto','salmon','shoal'],
        genome="hg38.analysisSet",
        transcriptome=['ensembl.85','knownGene']),
    'rnaseq_eda' : expand(
        'reports/RNA-seq/{tool_and_genome}_{transcriptome}-exploration.html',
        tool_and_genome=[
            'star_hg38.analysisSet', 'hisat2_grch38_snp_tran',
            'salmon_hg38.analysisSet', 'kallisto_hg38.analysisSet',
            'shoal_hg38.analysisSet',
        ],
        transcriptome=['ensembl.85', 'knownGene']),
    'rnaseq_compare' : 'reports/RNA-seq/rnaseq-compare.html',
    'rnaseq_diffexp' : expand(
        'reports/RNA-seq/{tool_and_genome}_{transcriptome}-diffexp.html',
        tool_and_genome=[
            'star_hg38.analysisSet', 'hisat2_grch38_snp_tran',
            'salmon_hg38.analysisSet', 'kallisto_hg38.analysisSet',
            'shoal_hg38.analysisSet',
        ],
        transcriptome=['ensembl.85', 'knownGene']),
    'rnaseq_gst' : 'saved_data/CAMERA-results-RNA.RDS',
    'chipseq_eda' : expand(
        'reports/ChIP-seq/{chip_antibody}-exploration.html',
        chip_antibody=chipseq_samplemeta_noinput['chip_antibody'].unique(),
    ),
    # TODO: Output tables
    'chipseq_diffmod' : expand(
        'reports/ChIP-seq/{chip_antibody}-diffmod.html',
        chip_antibody=chipseq_samplemeta_noinput['chip_antibody'].unique(),
    ),
    'chipseq_nvm_diminish': expand(
        'reports/ChIP-seq/{chip_antibody}-NvM-diminish-analysis.html',
        chip_antibody=chipseq_samplemeta_noinput['chip_antibody'].unique(),
    ),
    'promoter_eda': expand(
        'reports/ChIP-seq/{genome}_{transcriptome}_{chip_antibody}-{promoter_radius}-promoter-exploration.html',
        zip_longest_recycled, genome=["hg38.analysisSet"],
        **df_cartesian_product(pd.DataFrame({'transcriptome': ["knownGene", "ensembl.85"]}),
                               promoter_radius_table),
    ),
    'promoter_diffmod': expand(
        'reports/ChIP-seq/{genome}_{transcriptome}_{chip_antibody}-{promoter_radius}-promoter-diffmod.html',
        zip_longest_recycled, genome=["hg38.analysisSet"],
        **df_cartesian_product(pd.DataFrame({'transcriptome': ["knownGene", "ensembl.85"]}),
                               promoter_radius_table),
    ),
    'promoter_gst': expand(
        'saved_data/CAMERA-results-{chip_antibody}-{promoter_radius}-promoter.RDS',
        zip_longest_recycled, **promoter_radius_table
    ),
    'macs_predictd' : 'results/macs_predictd/output.log',
    'idr_peaks_epic' :expand(
        'peak_calls/epic_{genome_build}/{chip_antibody}_condition.{condition}_donor.ALL/peaks_noBL_IDR.narrowPeak',
        genome_build='hg38.analysisSet',
        chip_antibody=chipseq_samplemeta_noinput['chip_antibody'].unique(),
        condition = list(chipseq_samplemeta_noinput\
                         .apply(lambda x: '%s.%s' % (x['cell_type'], x['time_point']), axis=1).unique()) + ['ALL']),
    'idr_peaks_macs': expand(
        'peak_calls/macs_{genome_build}/{chip_antibody}_condition.{condition}_donor.ALL/peaks_noBL_IDR.narrowPeak',
        genome_build='hg38.analysisSet',
        chip_antibody=chipseq_samplemeta_noinput['chip_antibody'].unique(),
        condition = list(chipseq_samplemeta_noinput\
                         .apply(lambda x: '%s.%s' % (x['cell_type'], x['time_point']), axis=1).unique()) + ['ALL']),
    'idr_plots_one_cond': set(expand(
        expand('plots/IDR/{{peak_caller}}_{{genome_build}}/{chip_antibody}/condition.{cell_type}.{time_point}/{donorA}vs{donorB}_idrplots.pdf',
               zip_longest_recycled,
               **dict(idr_sample_pairs.iteritems())),
        peak_caller=['macs', 'epic'], genome_build='hg38.analysisSet')),
    'idr_plots_all_cond': set(expand(
        expand('plots/IDR/{{peak_caller}}_{{genome_build}}/{chip_antibody}/condition.ALL/{donorA}vs{donorB}_idrplots.pdf',
               zip_longest_recycled,
               **dict(idr_sample_pairs.iteritems())),
        peak_caller=['macs', 'epic'], genome_build='hg38.analysisSet')),
    'ccf_plots': expand('plots/csaw/CCF-plots{suffix}.pdf',
                        suffix=('', '-relative', '-noBL', '-relative-noBL')),
    'site_profile_plot': 'plots/csaw/site-profile-plots.pdf',
    'macs_peaks_allcond_alldonor': set(expand(
        'peak_calls/macs_{genome_build}/{chip_antibody}_condition.ALL_donor.ALL/peaks.narrowPeak', zip_longest_recycled,
        genome_build='hg38.analysisSet',
        chip_antibody=chipseq_samplemeta_noinput['chip_antibody'])),
    'macs_peaks_allcond_onedonor': set(expand(
        'peak_calls/macs_{genome_build}/{chip_antibody}_condition.ALL_donor.{donor}/peaks.narrowPeak', zip_longest_recycled,
        genome_build='hg38.analysisSet',
        chip_antibody=chipseq_samplemeta_noinput['chip_antibody'],
        donor=chipseq_samplemeta_noinput['donor_id'])),
    'macs_peaks_onecond_alldonor': set(expand(
        'peak_calls/macs_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point}_donor.ALL/peaks.narrowPeak', zip_longest_recycled,
        genome_build='hg38.analysisSet',
        chip_antibody=chipseq_samplemeta_noinput['chip_antibody'],
        cell_type=chipseq_samplemeta_noinput['cell_type'],
        time_point=chipseq_samplemeta_noinput['time_point'])),
    'macs_peaks_onecond_onedonor': set(expand(
        'peak_calls/macs_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point}_donor.{donor}/peaks.narrowPeak', zip_longest_recycled,
        genome_build='hg38.analysisSet',
        chip_antibody=chipseq_samplemeta_noinput['chip_antibody'],
        cell_type=chipseq_samplemeta_noinput['cell_type'],
        time_point=chipseq_samplemeta_noinput['time_point'],
        donor=chipseq_samplemeta_noinput['donor_id'])),
    'epic_peaks_allcond_alldonor': set(expand(
        'peak_calls/epic_{genome_build}/{chip_antibody}_condition.ALL_donor.ALL/peaks.tsv', zip_longest_recycled,
        genome_build='hg38.analysisSet',
        chip_antibody=chipseq_samplemeta_noinput['chip_antibody'])),
    'epic_peaks_allcond_onedonor': set(expand(

        'peak_calls/epic_{genome_build}/{chip_antibody}_condition.ALL_donor.{donor}/peaks.tsv', zip_longest_recycled,
        genome_build='hg38.analysisSet',
        chip_antibody=chipseq_samplemeta_noinput['chip_antibody'],
        donor=chipseq_samplemeta_noinput['donor_id'])),
    'epic_peaks_onecond_alldonor': set(expand(
        'peak_calls/epic_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point}_donor.ALL/peaks.tsv', zip_longest_recycled,
        genome_build='hg38.analysisSet',
        chip_antibody=chipseq_samplemeta_noinput['chip_antibody'],
        cell_type=chipseq_samplemeta_noinput['cell_type'],
        time_point=chipseq_samplemeta_noinput['time_point'])),
    'epic_peaks_onecond_onedonor': set(expand(
        'peak_calls/epic_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point}_donor.{donor}/peaks.tsv', zip_longest_recycled,
        genome_build='hg38.analysisSet',
        chip_antibody=chipseq_samplemeta_noinput['chip_antibody'],
        cell_type=chipseq_samplemeta_noinput['cell_type'],
        time_point=chipseq_samplemeta_noinput['time_point'],
        donor=chipseq_samplemeta_noinput['donor_id'])),
    'all_idr_one_cond': set(expand(
        expand('plots/IDR/{{peak_caller}}_{{genome_build}}/{chip_antibody}/condition.{cell_type}.{time_point}/{donorA}vs{donorB}_idrplots.pdf',
               zip_longest_recycled,
               **dict(idr_sample_pairs.iteritems())),
        peak_caller=['macs', 'epic'], genome_build='hg38.analysisSet')),
    'all_idr_all_cond': set(expand(
        expand('plots/IDR/{{peak_caller}}_{{genome_build}}/{chip_antibody}/condition.ALL/{donorA}vs{donorB}_idrplots.pdf',
               zip_longest_recycled,
               **dict(idr_sample_pairs.iteritems())),
        peak_caller=['macs', 'epic'], genome_build='hg38.analysisSet')),
    'all_idr_filtered_peaks_epic': expand(
        'peak_calls/epic_{genome_build}/{chip_antibody}_condition.{condition}_donor.ALL/peaks_noBL_IDR.narrowPeak',
        genome_build='hg38.analysisSet',
        chip_antibody=chipseq_samplemeta_noinput['chip_antibody'].unique(),
        condition = list(chipseq_samplemeta_noinput.apply(lambda x: '%s.%s' % (x['cell_type'], x['time_point']), axis=1).unique()) + ['ALL']),
    'all_idr_filtered_peaks_macs': expand(
        'peak_calls/macs_{genome_build}/{chip_antibody}_condition.{condition}_donor.ALL/peaks_noBL_IDR.narrowPeak',
        genome_build='hg38.analysisSet',
        chip_antibody=chipseq_samplemeta_noinput['chip_antibody'].unique(),
        condition = list(chipseq_samplemeta_noinput.apply(lambda x: '%s.%s' % (x['cell_type'], x['time_point']), axis=1).unique()) + ['ALL']),
    'mofa': [
        'reports/promoter-mofa-analyze.html',
        'reports/peak-mofa-analyze.html',
    ],
}


rule all:
    '''This rule aggregates all the final outputs of the pipeline.'''
    input:
        targets['rnaseq_eda'],
        targets['rnaseq_compare'],
        targets['rnaseq_diffexp'],
        targets['rnaseq_gst'],
        targets['chipseq_eda'],
        targets['chipseq_diffmod'],
        targets['chipseq_nvm_diminish'],
        targets['promoter_eda'],
        targets['promoter_diffmod'],
        targets['promoter_gst'],
        targets['macs_predictd'],
        targets['idr_peaks_epic'],
        targets['idr_peaks_macs'],
        targets['idr_plots_one_cond'],
        targets['idr_plots_all_cond'],
        targets['ccf_plots'],
        targets['site_profile_plot'],
        targets['mofa'],
        'reports/lamere_2016_fig7.html',

rule all_rnaseq:
    '''This rule aggregates all the final outputs of the pipeline.'''
    input:
        targets['rnaseq_eda'],
        targets['rnaseq_compare'],
        targets['rnaseq_diffexp'],
        targets['rnaseq_gst'],

rule all_chipseq:
    '''This rule aggregates all the final outputs of the pipeline.'''
    input:
        targets['chipseq_eda'],
        targets['chipseq_diffmod'],
        targets['chipseq_nvm_diminish'],
        targets['promoter_eda'],
        targets['promoter_diffmod'],
        targets['promoter_gst'],
        targets['macs_predictd'],
        targets['idr_peaks_epic'],
        targets['idr_peaks_macs'],
        targets['idr_plots_one_cond'],
        targets['idr_plots_all_cond'],
        targets['ccf_plots'],
        targets['site_profile_plot'],
        'reports/lamere_2016_fig7.html',

rule all_rnaseq_counts:
    input: targets['rnaseq_counts']

rule all_rnaseq_eda:
    input: targets['rnaseq_eda']

rule all_rnaseq_quant:
    input: targets['rnaseq_quant']

rule all_rnaseq_diffexp:
    input: targets['rnaseq_diffexp']

rule all_macs_callpeak:
    input:
        targets['macs_peaks_allcond_alldonor'],
        targets['macs_peaks_allcond_onedonor'],
        targets['macs_peaks_onecond_alldonor'],
        targets['macs_peaks_onecond_onedonor'],

rule all_epic_callpeak:
    input:
        targets['epic_peaks_allcond_alldonor'],
        targets['epic_peaks_allcond_onedonor'],
        targets['epic_peaks_onecond_alldonor'],
        targets['epic_peaks_onecond_onedonor'],

rule all_idr:
    input:
        targets['all_idr_one_cond'],
        targets['all_idr_all_cond'],

rule all_idr_filtered_peaks:
    input:
        targets['all_idr_filtered_peaks_epic'],
        targets['all_idr_filtered_peaks_macs'],

rule all_mofa:
    input:
        'reports/promoter-mofa-analyze.html',
        'reports/peak-mofa-analyze.html',

rule fetch_sra_run:
    '''Script to fetch the .sra file for an SRA run.

    (An SRA run identifier starts with SRR.)

    https://www.ncbi.nlm.nih.gov/sra

    '''
    output: 'sra_files/{sra_run,SRR.*}.sra'
    version: SOFTWARE_VERSIONS['ASCP']
    resources: concurrent_downloads=1
    shell: 'scripts/get-sra-run-files.R {wildcards.sra_run:q}'

rule extract_fastq:
    '''Extract FASTQ from SRA files.

    Because the SRA files were originally generated from
    coordinate-sorted BAM files, the reads in the SRA files are likely
    also sorted. Hence, during extraction, the reads are also shuffled
    deterministically (i.e. using a fixed seed). This ensures that
    downstream tools expecting the reads in random order with respect
    to their mapping position (e.g. Salmon) will be satisfied.

    https://ncbi.github.io/sra-tools/

    http://homes.cs.washington.edu/~dcjones/fastq-tools/

    '''
    input: 'sra_files/{sra_run}.sra'
    output:
        fqfile='fastq_files/{sra_run}.{fqext,fq(|\\.gz|\\.bz2|\\.qp)}',
        temp_unshuffled=temp('fastq_files/{sra_run}_unshuffled.{fqext,fq(|\\.gz|\\.bz2|\\.qp)}_temp'),
        temp_shuffled=temp('fastq_files/{sra_run}_shuffled.{fqext,fq(|\\.gz|\\.bz2|\\.qp)}_temp'),
    version: (SOFTWARE_VERSIONS['SRATOOLKIT'], SOFTWARE_VERSIONS['FASTQ_TOOLS'])
    resources: diskio=1
    params:
        compress_cmd = lambda wildcards: fastq_compression_cmds[wildcards.fqext]['compress']
    shell:'''
    echo "Dumping fastq for {wildcards.sra_run:q}..."
    fastq-dump --stdout {input:q} | \
      scripts/fill-in-empty-fastq-qual.py \
      > {output.temp_unshuffled:q}
    echo "Shuffling fastq for {wildcards.sra_run:q}..."
    fastq-sort --random --seed=1986 {output.temp_unshuffled:q} > {output.temp_shuffled:q}
    echo "Compressing fastq for {wildcards.sra_run:q}..."
    {params.compress_cmd} < {output.temp_shuffled:q} > {output:q}
    rm -f {output.temp_unshuffled:q} {output.temp_shuffled:q}
    '''

rule align_rnaseq_with_star_single_end:
    '''Align fastq file with STAR.

    https://github.com/alexdobin/STAR

    '''
    input:
        fastq='fastq_files/{samplename}.fq.gz',
        index_sa=hg38_ref('STAR_index_{genome_build}_{transcriptome}/SA'),
        transcriptome_gff=hg38_ref('{transcriptome}.gff3'),
    output:
        sam=temp('aligned/rnaseq_star_{genome_build}_{transcriptome}/{samplename}/Aligned.out.sam'),
        sj='aligned/rnaseq_star_{genome_build}_{transcriptome}/{samplename}/SJ.out.tab',
        logs=[ os.path.join('aligned/rnaseq_star_{genome_build}_{transcriptome}/{samplename}', fname)
               for fname in ['Log.final.out', 'Log.out', 'Log.progress.out'] ],
    params:
        # Note: trailing slash is significant here
        outdir='aligned/rnaseq_star_{genome_build}_{transcriptome}/{samplename}/',
        index_genomedir=hg38_ref('STAR_index_{genome_build}_{transcriptome}'),
        read_cmd=list2cmdline(fastq_compression_cmds['fq.gz']['decompress']),
    version: SOFTWARE_VERSIONS['STAR']
    threads: 8
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['star']
    shell: '''
    STAR \
      --runThreadN {threads:q} \
      --runMode alignReads \
      --genomeDir {params.index_genomedir:q} \
      --sjdbGTFfile {input.transcriptome_gff:q} \
      --sjdbGTFfeatureExon CDS \
      --sjdbGTFtagExonParentTranscript Parent \
      --sjdbGTFtagExonParentGene gene_id \
      --sjdbOverhang 100 \
      --readFilesIn {input.fastq:q} \
      --readFilesCommand {params.read_cmd:q} \
      --outSAMattributes Standard \
      --outSAMunmapped Within \
      --outFileNamePrefix {params.outdir:q} \
      --outSAMtype SAM
    '''

rule convert_star_sam_to_bam:
    input:
        sam='aligned/rnaseq_star_{genome_build}_{transcriptome}/{samplename}/Aligned.out.sam',
    output:
        bam='aligned/rnaseq_star_{genome_build}_{transcriptome}/{samplename}/Aligned.sortedByCoord.out.bam',
    shell: '''
    picard-tools SortSam \
      I={input.sam:q} O={output.bam:q} \
      SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT
    '''

rule align_rnaseq_with_hisat2_single_end:
    '''Align fastq file with HISAT2.

    https://ccb.jhu.edu/software/hisat2/index.shtml

    '''
    input: fastq='fastq_files/{samplename}.fq.gz',
           index_f1=hg38_ref('HISAT2_index_grch38_snp_tran/index.1.ht2'),
           transcriptome_gff=hg38_ref('knownGene.gff3'),
           chrom_mapping=hg38_ref('chrom_mapping_GRCh38_ensembl2UCSC.txt'),
    output: bam='aligned/rnaseq_hisat2_grch38_snp_tran/{samplename}/Aligned.bam',
    log: 'aligned/rnaseq_hisat2_grch38_snp_tran/{samplename}/hisat2.log'
    version: SOFTWARE_VERSIONS['HISAT2']
    threads: 8
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['hisat2']
    run:
        index_basename = re.sub('\\.1\\.ht2', '', input.index_f1)
        outdir = os.path.dirname(output.bam)
        cmds = [
            [
                'hisat2',
                '--threads', str(threads),
                '-q', '--phred33',
                '--very-sensitive',
                '--dta-cufflinks',
                '-x', index_basename,
                '-U', input.fastq,
                '-k', '20',
                '--time',
            ],
            [
                # Convert to UCSC chromosome names
                'scripts/bam-rename-chroms.py', input.chrom_mapping,
            ],
            [
                'picard-tools', 'SortSam', 'I=/dev/stdin', 'O=/dev/stdout',
                'SORT_ORDER=coordinate', 'VALIDATION_STRINGENCY=LENIENT',
            ]
        ]
        with atomic_write(output.bam, mode='wb', overwrite=True) as outfile, \
             open(log[0], mode='wb') as logfile:
            pipeline = Popen_pipeline(cmds, stdout=outfile, stderr=logfile)
            wait_for_subprocs(pipeline)

# There are multiple index_bam rules each restricted to a subset of
# bam files in order to improve the rulegraph appearance.
rule index_bam_rnaseq:
    '''Create .bai file for a bam file.

    This rule is identical to index_bam_chipseq. They are only
    separated in order to yield a less-confusing rule graph
    visualization.

    https://broadinstitute.github.io/picard/

    '''
    input: '{basename}.bam'
    output: '{basename,aligned/rnaseq_.*}.bam.bai'
    shell: '''
    picard-tools BuildBamIndex I={input:q} O={output:q} \
        VALIDATION_STRINGENCY=LENIENT
    '''

rule index_bam_chipseq:
    '''Create .bai file for a bam file.

    This rule is identical to index_bam_rnaseq. They are only
    separated in order to yield a less-confusing rule graph
    visualization.

    https://broadinstitute.github.io/picard/

    '''
    input: '{basename}.bam'
    output: '{basename,aligned/chipseq_.*}.bam.bai'
    shell: '''
    picard-tools BuildBamIndex I={input:q} O={output:q} \
        VALIDATION_STRINGENCY=LENIENT
    '''

rule bam2bed:
    '''Convert a bam file to bed using bedtools.

    http://bedtools.readthedocs.io/en/latest/

    '''
    input: '{basename}.bam'
    output: '{basename}_reads.bed'
    version: SOFTWARE_VERSIONS['BEDTOOLS']
    shell: '''
    bedtools bamtobed -i {input:q} > {output:q}
    '''

rule bam2bed_macs_filterdup:
    '''Convert a bam file to bed, filtering duplicates.

    This rule uses the 'macs2 filterdup' command to remove excess
    duplicate reads while converting to bed format. Note that this
    command does not necessarily remove all duplicates, only those in
    excess of what would be expected by chance. The log file reports
    the maximum number of duplicates allowed at each locus for a given
    sample.

    https://github.com/taoliu/MACS

    '''
    input: '{basename}.bam'
    output: bed='{basename}_reads_macs_filterdup.bed',
    log: '{basename}_macs_filterdup.log'
    version: SOFTWARE_VERSIONS['MACS']
    shell: '''
    macs2 filterdup --ifile {input:q} --format BAM \
      --gsize hs --keep-dup auto \
      --ofile {output.bed:q} \
      2>&1 | tee {log:q} 1>&2
    '''

# The hisat2 documentation doesn't specify which version of Ensembl
# they used to build the prebuilt index. Hopefully it doesn't matter
# too much.
rule count_rnaseq_hisat2_ensembl:
    '''Assign & count reads reads aligned to Ensembl genes by HISAT2.

    https://bioconductor.org/packages/release/bioc/html/Rsubread.html

    '''
    input:
        samplemeta='saved_data/samplemeta-RNASeq.RDS',
        bam_files=expand(
            'aligned/rnaseq_hisat2_grch38_snp_tran/{SRA_run}/Aligned.bam',
            SRA_run=rnaseq_samplemeta['SRA_run']),
        bai_files=expand(
            'aligned/rnaseq_hisat2_grch38_snp_tran/{SRA_run}/Aligned.bam.bai',
            SRA_run=rnaseq_samplemeta['SRA_run']),
        txdb=hg38_ref('TxDb.Hsapiens.ensembl.hg38.v85.sqlite3'),
        genemeta=hg38_ref('genemeta.ensembl.85.RDS')
    output: sexp='saved_data/SummarizedExperiment_rnaseq_hisat2_grch38_snp_tran_ensembl.{release}.RDS'
    params:
        expected_bam_files=','.join(expand(
            'aligned/rnaseq_hisat2_grch38_snp_tran/{SRA_run}/Aligned.bam',
            SRA_run=rnaseq_samplemeta['SRA_run'])),
        bam_file_pattern='aligned/rnaseq_hisat2_grch38_snp_tran/{{SAMPLE}}/Aligned.bam',
    version: R_package_version('RSubread')
    threads: len(rnaseq_samplemeta)
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['rnaseq_count']
    shell: '''
    scripts/rnaseq-count.R \
      --samplemeta-file {input.samplemeta:q} \
      --sample-id-column SRA_run \
      --bam-file-pattern {params.bam_file_pattern:q} \
      --output-file {output.sexp:q} \
      --expected-bam-files {params.expected_bam_files:q} \
      --threads {threads:q} \
      --annotation-txdb {input.txdb:q} \
      --additional-gene-info {input.genemeta:q}
    '''

rule count_rnaseq_hisat2_knownGene:
    '''Assign & count reads reads aligned to UCSC known genes by HISAT2.

    https://bioconductor.org/packages/release/bioc/html/Rsubread.html

    '''
    input:
        samplemeta='saved_data/samplemeta-RNASeq.RDS',
        bam_files=expand(
            'aligned/rnaseq_hisat2_grch38_snp_tran/{SRA_run}/Aligned.bam',
            SRA_run=rnaseq_samplemeta['SRA_run']),
        bai_files=expand(
            'aligned/rnaseq_hisat2_grch38_snp_tran/{SRA_run}/Aligned.bam.bai',
            SRA_run=rnaseq_samplemeta['SRA_run']),
        genemeta=hg38_ref('genemeta.org.Hs.eg.db.RDS')
    output: sexp='saved_data/SummarizedExperiment_rnaseq_hisat2_grch38_snp_tran_knownGene.RDS'
    params:
        expected_bam_files=','.join(expand(
            'aligned/rnaseq_hisat2_grch38_snp_tran/{SRA_run}/Aligned.bam',
            SRA_run=rnaseq_samplemeta['SRA_run'])),
        bam_file_pattern='aligned/rnaseq_hisat2_grch38_snp_tran/{{SAMPLE}}/Aligned.bam',
    version: R_package_version('RSubread')
    threads: len(rnaseq_samplemeta)
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['rnaseq_count']
    shell: '''
    scripts/rnaseq-count.R \
      --samplemeta-file {input.samplemeta:q} \
      --sample-id-column SRA_run \
      --bam-file-pattern {params.bam_file_pattern:q} \
      --output-file {output.sexp:q} \
      --expected-bam-files {params.expected_bam_files:q} \
      --threads {threads:q} \
      --annotation-txdb 'TxDb.Hsapiens.UCSC.hg38.knownGene' \
      --additional-gene-info {input.genemeta:q}
    '''

rule count_rnaseq_star_ensembl:
    '''Assign & count reads reads aligned to Ensembl genes by STAR.

    https://bioconductor.org/packages/release/bioc/html/Rsubread.html

    '''
    input:
        samplemeta='saved_data/samplemeta-RNASeq.RDS',
        bam_files=expand(
            'aligned/rnaseq_star_hg38.analysisSet_ensembl.{{release}}/{SRA_run}/Aligned.sortedByCoord.out.bam',
            SRA_run=rnaseq_samplemeta['SRA_run']),
        bai_files=expand(
            'aligned/rnaseq_star_hg38.analysisSet_ensembl.{{release}}/{SRA_run}/Aligned.sortedByCoord.out.bam.bai',
            SRA_run=rnaseq_samplemeta['SRA_run']),
        txdb=hg38_ref('TxDb.Hsapiens.ensembl.hg38.v{release}.sqlite3'),
        genemeta=hg38_ref('genemeta.ensembl.{release}.RDS')
    output: sexp='saved_data/SummarizedExperiment_rnaseq_star_hg38.analysisSet_ensembl.{release,\\d+}.RDS'
    params:
        expected_bam_files=','.join(expand(
            'aligned/rnaseq_star_hg38.analysisSet_ensembl.{{release}}/{SRA_run}/Aligned.sortedByCoord.out.bam',
            SRA_run=rnaseq_samplemeta['SRA_run'])),
        bam_file_pattern=lambda wildcards: expand('aligned/rnaseq_star_hg38.analysisSet_ensembl.{release}/{{SAMPLE}}/Aligned.sortedByCoord.out.bam', **wildcards)
    version: R_package_version('RSubread')
    threads: len(rnaseq_samplemeta)
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['rnaseq_count']
    shell: '''
    scripts/rnaseq-count.R \
      --samplemeta-file {input.samplemeta:q} \
      --sample-id-column SRA_run \
      --bam-file-pattern {params.bam_file_pattern:q} \
      --output-file {output.sexp:q} \
      --expected-bam-files {params.expected_bam_files:q} \
      --threads {threads:q} \
      --annotation-txdb {input.txdb:q} \
      --additional-gene-info {input.genemeta:q}
    '''
rule count_rnaseq_star_knownGene:
    '''Assign & count reads reads aligned to UCSC known genes by STAR.

    https://bioconductor.org/packages/release/bioc/html/Rsubread.html

    '''
    input:
        samplemeta='saved_data/samplemeta-RNASeq.RDS',
        bam_files=expand(
            'aligned/rnaseq_star_hg38.analysisSet_knownGene/{SRA_run}/Aligned.sortedByCoord.out.bam',
            SRA_run=rnaseq_samplemeta['SRA_run']),
        bai_files=expand(
            'aligned/rnaseq_star_hg38.analysisSet_knownGene/{SRA_run}/Aligned.sortedByCoord.out.bam.bai',
            SRA_run=rnaseq_samplemeta['SRA_run']),
        genemeta=hg38_ref('genemeta.org.Hs.eg.db.RDS')
    output: sexp='saved_data/SummarizedExperiment_rnaseq_star_hg38.analysisSet_knownGene.RDS'
    params:
        expected_bam_files=','.join(expand(
            'aligned/rnaseq_star_hg38.analysisSet_knownGene/{SRA_run}/Aligned.sortedByCoord.out.bam',
            SRA_run=rnaseq_samplemeta['SRA_run']))
    version: R_package_version('RSubread')
    threads: len(rnaseq_samplemeta)
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['rnaseq_count']
    shell: '''
    scripts/rnaseq-count.R \
      --samplemeta-file {input.samplemeta:q} \
      --sample-id-column SRA_run \
      --bam-file-pattern aligned/rnaseq_star_hg38.analysisSet_knownGene/%s/Aligned.sortedByCoord.out.bam \
      --output-file {output.sexp:q} \
      --expected-bam-files {params.expected_bam_files:q} \
      --threads {threads:q} \
      --annotation-txdb TxDb.Hsapiens.UCSC.hg38.knownGene \
      --additional-gene-info {input.genemeta:q}
    '''

rule quant_rnaseq_with_salmon:
    '''Quantify genes from reads using Salmon.

    https://combine-lab.github.io/salmon/

    '''
    input:
        salmon_index=hg38_ref('Salmon_index_{genome_build}_{transcriptome}/sa.bin'),
        genemap_file=hg38_ref('Salmon_index_{genome_build}_{transcriptome}/genemap.txt'),
        fastq='fastq_files/{SRA_run}.fq.gz',
    output:
        list_salmon_output_files('salmon_quant/{genome_build}_{transcriptome}/{SRA_run}')
    params:
        index_dir=hg38_ref('Salmon_index_{genome_build}_{transcriptome}'),
        outdir='salmon_quant/{genome_build}_{transcriptome}/{SRA_run}',
        libtype=lambda wildcards: rnaseq_sample_libtypes[wildcards.SRA_run]
    version: SOFTWARE_VERSIONS['SALMON']
    threads: 16
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['salmon']
    shell: '''
    timeout 1h \
      salmon quant \
      --index {params.index_dir:q} \
      --libType {params.libtype:q} \
      --unmatedReads {input.fastq:q} \
      --threads {threads:q} \
      --seqBias --gcBias --useVBOpt \
      --dumpEq --dumpEqWeights \
      --geneMap {input.genemap_file:q} \
      --output {params.outdir:q} \
      --auxDir aux_info \
      --numGibbsSamples 100
    '''

rule convert_salmon_to_hdf5:
    '''Produce abundance.h5 file for Salmon using wasabi.

    https://github.com/COMBINE-lab/wasabi

    '''
    input: list_salmon_output_files('{salmon_quant_dir}')
    output: '{salmon_quant_dir}/abundance.h5'
    shell: ''' scripts/convert-salmon-to-hdf5.R {wildcards.salmon_quant_dir:q} '''

rule run_shoal:
    '''Run shoal on the output of Salmon.'''
    input:
        samples=list_salmon_output_files(
            expand('salmon_quant/{{genome}}_{{transcriptome}}/{sample}',
                   sample=rnaseq_samplemeta['SRA_run'])),
    output:
        samples=expand('shoal_quant/{{genome}}_{{transcriptome}}/{sample}_adapt.sf',
                       sample=rnaseq_samplemeta['SRA_run']),
        prior='shoal_quant/{genome}_{transcriptome}/prior/prior.tsv',
    params:
        quantdir='salmon_quant/{genome}_{transcriptome}',
        outdir='shoal_quant/{genome}_{transcriptome}',
    threads: 16
    shell: '''
    run_shoal.sh -j {threads:q} -q {params.quantdir:q} -o {params.outdir:q}
    '''

rule quant_rnaseq_with_kallisto:
    '''Quantify genes from reads using Kallisto.

    https://pachterlab.github.io/kallisto/about

    '''
    input:
        kallisto_index=hg38_ref('Kallisto_index_{genome_build}_{transcriptome}'),
        fastq='fastq_files/{SRA_run}.fq.gz',
    output:
        list_kallisto_output_files('kallisto_quant/{genome_build}_{transcriptome}/{SRA_run}')
    params:
        outdir='kallisto_quant/{genome_build}_{transcriptome}/{SRA_run}',
        libtype=lambda wildcards: rnaseq_sample_libtypes[wildcards.SRA_run]
    version: SOFTWARE_VERSIONS['KALLISTO']
    threads: 16
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['kallisto']
    run:
        libType = list(rnaseq_samplemeta['libType'][rnaseq_samplemeta['SRA_run'] == wildcards.SRA_run])[0]
        if libType == 'SF':
            lib_opt = '--fr-stranded'
        elif libType == 'SR':
            lib_opt = '--rf-stranded'
        else:
            raise ValueError('Unknown kallisto libtype: {}'.format(libType))
        shell('''
        mkdir -p {params.outdir:q}
        kallisto quant \
          --index {input.kallisto_index:q} --output-dir {params.outdir:q} \
          {lib_opt:q} --single --threads {threads:q} --bootstrap-samples 100 \
          --bias --fragment-length 200 --sd 80 {input.fastq:q}
        ''')

# TODO: Write R script to convert bootstraps into SummarizedExperiment
# RDS file, and write a rule for it.

rule align_chipseq_with_bowtie2:
    '''Align ChIP-seq reads to genome using bowtie2.

    http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

    '''
    input:
        fastq='fastq_files/{SRA_run}.fq.gz',
        index_file=hg38_ref('BT2_index_{genome_build}/index.1.bt2l')
    output:
        bam='aligned/chipseq_bowtie2_{genome_build}/{SRA_run}/Aligned.bam'
    params:
        index_basename=hg38_ref('BT2_index_{genome_build}/index')
    version: SOFTWARE_VERSIONS['BOWTIE2']
    threads: 8
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['bowtie2']
    shell: '''
    bowtie2 --threads {threads:q} --mm \
      -U {input.fastq:q} -x {params.index_basename:q} -q \
      --end-to-end --sensitive | \
    picard-tools SortSam I=/dev/stdin O={output.bam:q} \
      SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT
    '''

rule get_liftover_chain:
    '''Download LiftOver chain file from UCSC.

    http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/

    '''
    input: FTP.remote('hgdownload.cse.ucsc.edu/goldenPath/{src_genome}/liftOver/{src_genome}ToHg38.over.chain.gz', static=True)
    output: 'saved_data/{src_genome}ToHg38.over.chain'
    shell: 'zcat {input:q} > {output:q}'

rule get_motifmap:
    '''Download MotifMap data for hg19

    http://motifmap.igb.uci.edu/

    '''
    input: tar=HTTP.remote('www.igb.uci.edu/~motifmap/motifmap/HUMAN/hg19/multiz46way_placental/HUMAN.hg19_multiz46way.tar.bz2', insecure=True, static=True)
    output: bed='saved_data/MotifMap_HUMAN_hg19_BBLS_1_00_FDR_0_10.bed'
    params: file_inside_tar='HUMAN_hg19_BBLS_1_00_FDR_0_10.bed'
    shell: '''tar -O -xjf {input.tar:q} {params.file_inside_tar:q} > {output.bed:q}'''

rule liftover_motifmap:
    '''Use LiftOver to translate MotifMap BED to hg38 coordinates.'''
    input: bed='saved_data/MotifMap_HUMAN_hg19_BBLS_1_00_FDR_0_10.bed',
           chain='saved_data/hg19ToHg38.over.chain',
    output: bed='saved_data/MotifMap_HUMAN_hg38_BBLS_1_00_FDR_0_10.bed'
    params: allow_gap=2
    script: 'scripts/liftOver-MotifMap.R'

rule get_cpg:
    output: 'saved_data/UCSC_hg38_cpgIslandExtUnmasked.RDS'
    script: 'scripts/get-CpG.R'

# http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability
rule get_blacklist_regions:
    '''Download UCSC "consensus excludable regions" tracks.

    http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability

    '''
    input: FTP.remote('hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/{track_name}.bed.gz', static=True)
    output: 'saved_data/{track_name,wgEncode.*}_hg19.bed'
    shell: '''zcat {input:q} > {output:q}'''

rule liftover_blacklist_regions:
    '''Use LiftOver to translate blacklists to hg38 coordinates.

    http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/

    '''
    input: bed='saved_data/{track_name}_hg19.bed',
           chain='saved_data/hg19ToHg38.over.chain',
    output: bed='saved_data/{track_name,[^_]+}.bed'
    shell: '''liftOver {input.bed:q} {input.chain:q} {output.bed:q} /dev/null'''

rule generate_greylist:
    '''Generate greylist regions from ChIP-Seq input files.

    This uses the methodology from the GreyListChIP Bioconductor
    package, but uses a custom implementation rather than using the
    package itself.

    http://bioconductor.org/packages/release/bioc/html/GreyListChIP.html

    '''
    input:
        samplemeta='saved_data/samplemeta-ChIPSeq.RDS',
        chip_input=expand('aligned/chipseq_bowtie2_hg38.analysisSet/{SRA_run}/Aligned.bam',
                          SRA_run=dfselect(chipseq_samplemeta, 'SRA_run',
                                           chip_antibody='input')),
    output:
        nbfit_data='saved_data/ChIPSeq-input-depth-NBGLM-fits.RDS',
        counts_data='saved_data/window-counts-input-unfiltered-1kb.RDS',
        greylist_data='saved_data/ChIPSeq-input-greylist.RDS',
        greylist_bed='saved_data/ChIPSeq-input-greylist.bed',
    version: SOFTWARE_VERSIONS['BIOC']
    threads: 8
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['greylist']
    shell: '''MC_CORES={threads:q} scripts/generate-greylists.R'''

rule merge_blacklists:
    '''Merge ChIP-Seq blacklists from multiple sources into a single file.'''
    input:
        'saved_data/wgEncodeDacMapabilityConsensusExcludable.bed',
        'saved_data/wgEncodeDukeMapabilityRegionsExcludable.bed',
        'saved_data/ChIPSeq-input-greylist.bed',
    output:
        'saved_data/ChIPSeq-merged-blacklist.bed'
    shell: '''cat {input:q} > {output:q}'''

rule generate_promoter_regions_ensembl:
    '''Generate a file describing promoter regions for Ensembl genes.'''
    input:
        txdb=hg38_ref('TxDb.Hsapiens.ensembl.hg38.v{release}.sqlite3'),
        genemeta=hg38_ref('genemeta.ensembl.{release}.RDS'),
    output:
        rds="saved_data/promoter-regions_hg38.analysisSet_ensembl.{release}_{radius,[0-9.]+.*?bp}.RDS",
    shell: '''
    scripts/generate-promoters.R \
      --txdb {input.txdb:q} \
      --promoter-radius {wildcards.radius:q} \
      --additional-gene-info {input.genemeta:q} \
      --output-file {output.rds:q}
    '''

rule generate_promoter_regions_knownGene:
    '''Generate a file describing promoter regions for UCSC knownGene.'''
    input:
        genemeta=hg38_ref('genemeta.org.Hs.eg.db.RDS'),
    params:
        txdb='TxDb.Hsapiens.UCSC.hg38.knownGene'
    output:
        rds="saved_data/promoter-regions_hg38.analysisSet_knownGene_{radius,[0-9.]+.*?bp}.RDS",
    shell: '''
    scripts/generate-promoters.R \
      --txdb {params.txdb:q} \
      --promoter-radius {wildcards.radius:q} \
      --additional-gene-info {input.genemeta:q} \
      --output-file {output.rds:q}
    '''

rule macs_predictd:
    '''Determine ChIP-Seq fragment length using 'macs2 predictd'.

    https://github.com/taoliu/MACS

    '''
    input: bam_files=aligned_chipseq_bam_files,
    output:
        rfile='results/macs_predictd/predictd',
        pdf='plots/predictd_model.pdf',
        # This is actually the desired output file, despite being
        # a log file, so it is listed in the output files instead of declaring it as a log file
        logfile='results/macs_predictd/output.log'
    version: SOFTWARE_VERSIONS['MACS']
    run:
        rfile_basename = os.path.basename(output.rfile)
        rfile_dirname = os.path.dirname(output.rfile)
        shell('''
        macs2 predictd -i {input.bam_files:q} -f BAM -g hs \
          --outdir {rfile_dirname:q} --rfile {rfile_basename:q} \
          &>{output.logfile:q}
        cd plots
        Rscript ../{rfile_dirname:q}/{rfile_basename:q}
        ''')

rule callpeak_macs_all_conditions_all_donors:
    '''Call peaks using macs on all samples.

    https://github.com/taoliu/MACS

    '''
    input:
        chip_input=lambda wildcards:
        expand('aligned/chipseq_bowtie2_{genome_build}/{SRA_run}/Aligned.bam',
               SRA_run=dfselect(chipseq_samplemeta, 'SRA_run',
                                chip_antibody='input'),
               genome_build=wildcards.genome_build),
        chip_pulldown=lambda wildcards:
        expand('aligned/chipseq_bowtie2_{genome_build}/{SRA_run}/Aligned.bam',
               SRA_run=dfselect(chipseq_samplemeta, 'SRA_run',
                                chip_antibody=wildcards.chip_antibody),
               genome_build=wildcards.genome_build)
    output:
        outfiles=list_macs_callpeak_output_files('peak_calls/macs_{genome_build}/{chip_antibody}_condition.ALL_donor.ALL'),
    log: 'peak_calls/macs_{genome_build}/{chip_antibody}_condition.ALL_donor.ALL/peakcall.log'
    params:
        outdir='peak_calls/macs_{genome_build}/{chip_antibody}_condition.ALL_donor.ALL',
    version: SOFTWARE_VERSIONS['MACS']
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['macs_callpeak']
    shell: '''
    macs2 callpeak \
      --treatment {input.chip_pulldown:q} \
      --control {input.chip_input:q} \
      --format BAM \
      --gsize hs \
      --keep-dup auto \
      --outdir {params.outdir:q} \
      --name peakcall \
      --nomodel --extsize 147 \
      --pvalue=0.5 \
      2>&1 | tee {log:q};
    prename 's/peakcall_//' {params.outdir:q}/*
    '''

rule callpeak_macs_all_conditions_single_donor:
    '''Call peaks using macs on all samples from a single donor.

    https://github.com/taoliu/MACS

    '''
    input:
        chip_input=lambda wildcards:
        expand('aligned/chipseq_bowtie2_{genome_build}/{SRA_run}/Aligned.bam',
               SRA_run=dfselect(chipseq_samplemeta, 'SRA_run',
                                chip_antibody='input'),
               genome_build=wildcards.genome_build),
        chip_pulldown=lambda wildcards:
        expand('aligned/chipseq_bowtie2_{genome_build}/{SRA_run}/Aligned.bam',
               SRA_run=dfselect(chipseq_samplemeta, 'SRA_run',
                                chip_antibody=wildcards.chip_antibody,
                                donor_id=wildcards.donor),
               genome_build=wildcards.genome_build)
    output:
        outfiles=list_macs_callpeak_output_files('peak_calls/macs_{genome_build}/{chip_antibody}_condition.ALL_donor.{donor,D[0-9]+}'),
    log: 'peak_calls/macs_{genome_build}/{chip_antibody}_condition.ALL_donor.{donor,D[0-9]+}/peakcall.log'
    params:
        outdir='peak_calls/macs_{genome_build}/{chip_antibody}_condition.ALL_donor.{donor}',
    version: SOFTWARE_VERSIONS['MACS']
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['macs_callpeak']
    shell: '''
    macs2 callpeak \
      --treatment {input.chip_pulldown:q} \
      --control {input.chip_input:q} \
      --format BAM \
      --gsize hs \
      --keep-dup auto \
      --outdir {params.outdir:q} \
      --name peakcall \
      --bdg \
      --nomodel --extsize 147 \
      --pvalue=0.5 \
      2>&1 | tee {log:q}
    prename 's/peakcall_//' {params.outdir:q}/*
    '''

rule callpeak_macs_single_condition_all_donors:
    '''Call peaks using macs on all samples from a single condition.

    https://github.com/taoliu/MACS

    '''
    input:
        chip_input=lambda wildcards:
        expand('aligned/chipseq_bowtie2_{genome_build}/{SRA_run}/Aligned.bam',
               SRA_run=dfselect(chipseq_samplemeta, 'SRA_run',
                                chip_antibody='input'),
               genome_build=wildcards.genome_build),
        chip_pulldown=lambda wildcards:
        expand('aligned/chipseq_bowtie2_{genome_build}/{SRA_run}/Aligned.bam',
               SRA_run=dfselect(chipseq_samplemeta, 'SRA_run',
                                chip_antibody=wildcards.chip_antibody,
                                cell_type=wildcards.cell_type,
                                time_point=wildcards.time_point),
               genome_build=wildcards.genome_build)
    output:
        outfiles=list_macs_callpeak_output_files('peak_calls/macs_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_donor.ALL'),
    log: 'peak_calls/macs_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_donor.ALL/peakcall.log'
    params:
        outdir='peak_calls/macs_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point}_donor.ALL',
    version: SOFTWARE_VERSIONS['MACS']
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['macs_callpeak']
    shell: '''
    macs2 callpeak \
      --treatment {input.chip_pulldown:q} \
      --control {input.chip_input:q} \
      --format BAM \
      --gsize hs \
      --keep-dup auto \
      --outdir {params.outdir:q} \
      --name peakcall \
      --bdg \
      --nomodel --extsize 147 \
      --pvalue=0.5 \
      2>&1 | tee {log:q}
    prename 's/peakcall_//' {params.outdir:q}/*
    '''

rule callpeak_macs_single_condition_single_donor:
    '''Call peaks using macs on a single sample.

    https://github.com/taoliu/MACS

    '''
    input:
        chip_input=lambda wildcards:
        expand('aligned/chipseq_bowtie2_{genome_build}/{SRA_run}/Aligned.bam',
               SRA_run=dfselect(chipseq_samplemeta, 'SRA_run',
                                chip_antibody='input'),
               genome_build=wildcards.genome_build),
        chip_pulldown=lambda wildcards:
        expand('aligned/chipseq_bowtie2_{genome_build}/{SRA_run}/Aligned.bam',
               SRA_run=dfselect(chipseq_samplemeta, 'SRA_run',
                                chip_antibody=wildcards.chip_antibody,
                                donor_id=wildcards.donor,
                                cell_type=wildcards.cell_type,
                                time_point=wildcards.time_point),
               genome_build=wildcards.genome_build)
    output:
        outfiles=list_macs_callpeak_output_files('peak_calls/macs_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_donor.{donor,D[0-9]+}'),
    log: 'peak_calls/macs_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_donor.{donor,D[0-9]+}/peakcall.log'
    params:
        outdir='peak_calls/macs_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point}_donor.{donor}',
    version: SOFTWARE_VERSIONS['MACS']
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['macs_callpeak']
    shell: '''
    macs2 callpeak \
      --treatment {input.chip_pulldown:q} \
      --control {input.chip_input:q} \
      --format BAM \
      --gsize hs \
      --keep-dup auto \
      --outdir {params.outdir:q} \
      --name peakcall \
      --nomodel --extsize 147 \
      --pvalue=0.5 \
      2>&1 | tee {log:q}
    prename 's/peakcall_//' {params.outdir:q}/*
    '''

rule callpeak_epic_all_conditions_all_donors:
    '''Call peaks using epic on all samples.

    https://github.com/endrebak/epic/

    '''
    input:
        chip_input=lambda wildcards:
        expand('aligned/chipseq_bowtie2_{genome_build}/{SRA_run}/Aligned_reads_macs_filterdup.bed',
               SRA_run=dfselect(chipseq_samplemeta, 'SRA_run',
                                chip_antibody='input'),
               genome_build=wildcards.genome_build),
        chip_pulldown=lambda wildcards:
        expand('aligned/chipseq_bowtie2_{genome_build}/{SRA_run}/Aligned_reads_macs_filterdup.bed',
               SRA_run=dfselect(chipseq_samplemeta, 'SRA_run',
                                chip_antibody=wildcards.chip_antibody),
               genome_build=wildcards.genome_build)
    output:
        peaks='peak_calls/epic_{genome_build}/{chip_antibody}_condition.ALL_donor.ALL/peaks.tsv',
        chip_bw='peak_calls/epic_{genome_build}/{chip_antibody}_condition.ALL_donor.ALL/bigwig/sum_treatment.bw',
        control_bw='peak_calls/epic_{genome_build}/{chip_antibody}_condition.ALL_donor.ALL/bigwig/sum_input.bw',
    log: 'peak_calls/epic_{genome_build}/{chip_antibody}_condition.ALL_donor.ALL/peakcall.log',
    params:
        outdir='peak_calls/epic_{genome_build}/{chip_antibody}_condition.ALL_donor.ALL',
    version: SOFTWARE_VERSIONS['EPIC']
    threads: 4
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['epic_callpeak']
    run:
        with open(output.peaks, 'wb') as outfile, open(log[0], 'wb') as logfile:
            cmd = ['epic'] + \
                  ['--treatment'] + input.chip_pulldown + \
                  ['--control'] + input.chip_input + \
                  [
                      '--number-cores', threads,
                      '--genome', 'hg38',
                      '--fragment-size', 147,
                      '--keep-duplicates', 'True',
                      '--bigwig', os.path.join(params.outdir, 'bigwig'),
                  ]
            cmd = [str(x) for x in cmd]
            p = Popen(cmd, stdout=outfile, stderr=PIPE)
            for logline in p.stderr:
                logfile.write(logline)
                sys.stderr.write(logline.decode(sys.getdefaultencoding()))

rule callpeak_epic_all_conditions_single_donor:
    '''Call peaks using epic on all samples from a single donor.

    https://github.com/endrebak/epic/

    '''
    input:
        chip_input=lambda wildcards:
        expand('aligned/chipseq_bowtie2_{genome_build}/{SRA_run}/Aligned_reads_macs_filterdup.bed',
               SRA_run=dfselect(chipseq_samplemeta, 'SRA_run',
                                chip_antibody='input'),
               genome_build=wildcards.genome_build),
        chip_pulldown=lambda wildcards:
        expand('aligned/chipseq_bowtie2_{genome_build}/{SRA_run}/Aligned_reads_macs_filterdup.bed',
               SRA_run=dfselect(chipseq_samplemeta, 'SRA_run',
                                chip_antibody=wildcards.chip_antibody,
                                donor_id=wildcards.donor),
               genome_build=wildcards.genome_build)
    output:
        peaks='peak_calls/epic_{genome_build}/{chip_antibody}_condition.ALL_donor.{donor,D[0-9]+}/peaks.tsv',
        chip_bw='peak_calls/epic_{genome_build}/{chip_antibody}_condition.ALL_donor.{donor,D[0-9]+}/bigwig/sum_treatment.bw',
        control_bw='peak_calls/epic_{genome_build}/{chip_antibody}_condition.ALL_donor.{donor,D[0-9]+}/bigwig/sum_input.bw',
    log: 'peak_calls/epic_{genome_build}/{chip_antibody}_condition.ALL_donor.{donor,D[0-9]+}/peakcall.log',
    params:
        outdir='peak_calls/epic_{genome_build}/{chip_antibody}_condition.ALL_donor.{donor}',
    version: SOFTWARE_VERSIONS['EPIC']
    threads: 4
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['epic_callpeak']
    run:
        with open(output.peaks, 'wb') as outfile, open(log[0], 'wb') as logfile:
            cmd = ['epic'] + \
                  ['--treatment'] + input.chip_pulldown + \
                  ['--control'] + input.chip_input + \
                  [
                      '--number-cores', threads,
                      '--genome', 'hg38',
                      '--fragment-size', 147,
                      '--keep-duplicates', 'True',
                      '--bigwig', os.path.join(params.outdir, 'bigwig'),
                  ]
            cmd = [str(x) for x in cmd]
            p = Popen(cmd, stdout=outfile, stderr=PIPE)
            for logline in p.stderr:
                logfile.write(logline)
                sys.stderr.write(logline.decode(sys.getdefaultencoding()))

rule callpeak_epic_single_condition_all_donors:
    '''Call peaks using epic on all samples from a single condition.

    https://github.com/endrebak/epic/

    '''
    input:
        chip_input=lambda wildcards:
        expand('aligned/chipseq_bowtie2_{genome_build}/{SRA_run}/Aligned_reads_macs_filterdup.bed',
               SRA_run=dfselect(chipseq_samplemeta, 'SRA_run',
                                chip_antibody='input'),
               genome_build=wildcards.genome_build),
        chip_pulldown=lambda wildcards:
        expand('aligned/chipseq_bowtie2_{genome_build}/{SRA_run}/Aligned_reads_macs_filterdup.bed',
               SRA_run=dfselect(chipseq_samplemeta, 'SRA_run',
                                chip_antibody=wildcards.chip_antibody,
                                cell_type=wildcards.cell_type,
                                time_point=wildcards.time_point),
               genome_build=wildcards.genome_build)
    output:
        peaks='peak_calls/epic_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_donor.ALL/peaks.tsv',
        chip_bw='peak_calls/epic_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_donor.ALL/bigwig/sum_treatment.bw',
        control_bw='peak_calls/epic_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_donor.ALL/bigwig/sum_input.bw',
    log: 'peak_calls/epic_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_donor.ALL/peakcall.log'
    params:
        outdir='peak_calls/epic_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point}_donor.ALL',
    version: SOFTWARE_VERSIONS['EPIC']
    threads: 4
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['epic_callpeak']
    run:
        with open(output.peaks, 'wb') as outfile, open(log[0], 'wb') as logfile:
            cmd = ['epic'] + \
                  ['--treatment'] + input.chip_pulldown + \
                  ['--control'] + input.chip_input + \
                  [
                      '--number-cores', threads,
                      '--genome', 'hg38',
                      '--fragment-size', 147,
                      '--keep-duplicates', 'True',
                      '--bigwig', os.path.join(params.outdir, 'bigwig'),
                  ]
            cmd = [str(x) for x in cmd]
            p = Popen(cmd, stdout=outfile, stderr=PIPE)
            for logline in p.stderr:
                logfile.write(logline)
                sys.stderr.write(logline.decode(sys.getdefaultencoding()))

rule callpeak_epic_single_condition_single_donor:
    '''Call peaks using epic on a single sample.

    https://github.com/endrebak/epic/

    '''
    input:
        chip_input=lambda wildcards:
        expand('aligned/chipseq_bowtie2_{genome_build}/{SRA_run}/Aligned_reads_macs_filterdup.bed',
               SRA_run=dfselect(chipseq_samplemeta, 'SRA_run',
                                chip_antibody='input'),
               genome_build=wildcards.genome_build),
        chip_pulldown=lambda wildcards:
        expand('aligned/chipseq_bowtie2_{genome_build}/{SRA_run}/Aligned_reads_macs_filterdup.bed',
               SRA_run=dfselect(chipseq_samplemeta, 'SRA_run',
                                chip_antibody=wildcards.chip_antibody,
                                donor_id=wildcards.donor,
                                cell_type=wildcards.cell_type,
                                time_point=wildcards.time_point),
               genome_build=wildcards.genome_build)
    output:
        peaks='peak_calls/epic_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_donor.{donor,D[0-9]+}/peaks.tsv',
        chip_bw='peak_calls/epic_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_donor.{donor,D[0-9]+}/bigwig/sum_treatment.bw',
        control_bw='peak_calls/epic_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_donor.{donor,D[0-9]+}/bigwig/sum_input.bw',
    log: 'peak_calls/epic_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_donor.{donor,D[0-9]+}/peakcall.log',
    params:
        outdir='peak_calls/epic_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point}_donor.{donor}',
    version: SOFTWARE_VERSIONS['EPIC']
    threads: 4
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['epic_callpeak']
    run:
        with open(output.peaks, 'wb') as outfile, open(log[0], 'wb') as logfile:
            cmd = ['epic'] + \
                  ['--treatment'] + input.chip_pulldown + \
                  ['--control'] + input.chip_input + \
                  [
                      '--number-cores', threads,
                      '--genome', 'hg38',
                      '--fragment-size', 147,
                      '--keep-duplicates', 'True',
                      '--bigwig', os.path.join(params.outdir, 'bigwig'),
                  ]
            cmd = [str(x) for x in cmd]
            p = Popen(cmd, stdout=outfile, stderr=PIPE)
            for logline in p.stderr:
                logfile.write(logline)
                sys.stderr.write(logline.decode(sys.getdefaultencoding()))

rule convert_epic_to_narrowpeak:
    '''Convert epic results to narrowPeak format.

    This is important so that peak calls from both macs and epic are
    in the same format.

    '''
    input:
        'peak_calls/{dir}/peaks.tsv'
    output:
        'peak_calls/{dir,epic_.*}/peaks.narrowPeak'
    run:
        peaks = pd.DataFrame.from_csv(input[0], header=1, sep=' ', index_col=None)
        ndigits = ceil(log10(peaks.shape[0]+1))
        name_format = 'epic_peak_{{:0{}d}}'.format(ndigits)
        tiny_float = np.finfo(float).tiny
        narrowpeak = pd.DataFrame.from_items((
            ('chrom', peaks['Chromosome']),
            ('chromStart', peaks['Start']),
            ('chromEnd', peaks['End']),
            ('name', pd.Series(name_format.format(x) for x in range(peaks.shape[0]))),
            ('score', peaks['Score']),
            ('strand', '.'),
            ('signalValue', peaks['Fold_change']),
            ('pValue', -np.log10(peaks['P'].where(peaks['P'] > tiny_float, other=tiny_float))),
            ('qValue', -np.log10(peaks['FDR'].where(peaks['FDR'] > tiny_float, other=tiny_float))),
            # Epic doesn't call a peak, so just use the midpoint
            ('peak', np.array(np.around((peaks['End'] - peaks['Start']) / 2), dtype=np.int64)),
        ))
        narrowpeak.to_csv(output[0], sep='\t',
                          header=False, index=False,
                          quoting=csv.QUOTE_NONE,)

rule filter_blacklisted_peaks:
    '''Remove any peaks that overlap blacklisted regions.'''
    input:
        peaks='peak_calls/{dirname}/peaks.narrowPeak',
        blacklist='saved_data/ChIPSeq-merged-blacklist.bed',
    output:
        peaks='peak_calls/{dirname}/peaks_noBL.narrowPeak',
    version: SOFTWARE_VERSIONS['BEDTOOLS']
    shell: '''
    bedtools subtract -A -a {input.peaks:q} -b {input.blacklist:q} > {output.peaks:q}
    if [ ! -s {output.peaks:q} ]; then
      rm -f {output.peaks:q}
    fi
    '''

rule run_idr_macs_all_conditions:
    '''Run IDR on macs all-condition peak calls.

    https://github.com/nboley/idr

    '''
    input:
        donorA_peaks='peak_calls/macs_{genome_build}/{chip_antibody}_condition.ALL_donor.{donorA}/peaks_noBL.narrowPeak',
        donorB_peaks='peak_calls/macs_{genome_build}/{chip_antibody}_condition.ALL_donor.{donorB}/peaks_noBL.narrowPeak',
    output:
        temp_donorA_peaks=temp('idr_analysis/macs_{genome_build}/{chip_antibody}_condition.ALL_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/donorA_temp.narrowPeak'),
        temp_donorB_peaks=temp('idr_analysis/macs_{genome_build}/{chip_antibody}_condition.ALL_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/donorB_temp.narrowPeak'),
        outfile='idr_analysis/macs_{genome_build}/{chip_antibody}_condition.ALL_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/idrValues.txt',
        plotfile='idr_analysis/macs_{genome_build}/{chip_antibody}_condition.ALL_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/idrValues.png',

    log: 'idr_analysis/macs_{genome_build}/{chip_antibody}_condition.ALL_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/idr.log',
    version: SOFTWARE_VERSIONS['IDR']
    run:
        pick_top_peaks(input.donorA_peaks, output.temp_donorA_peaks, by='pValue', number=150000)
        pick_top_peaks(input.donorB_peaks, output.temp_donorB_peaks, by='pValue', number=150000)
        shell('''
        idr --samples {output.temp_donorA_peaks:q} {output.temp_donorB_peaks:q} \
          --input-file-type narrowPeak \
          --rank p.value \
          --output-file {output.outfile:q} \
          --output-file-type bed \
          --log-output-file {log:q} \
          --plot \
          --random-seed 1986
        mv {output.outfile:q}.png {output.plotfile:q}
        ''')

rule run_idr_macs_single_condition:
    '''Run IDR on macs single-condition peak calls.

    https://github.com/nboley/idr

    '''
    input:
        donorA_peaks='peak_calls/macs_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point}_donor.{donorA}/peaks_noBL.narrowPeak',
        donorB_peaks='peak_calls/macs_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point}_donor.{donorB}/peaks_noBL.narrowPeak',
    output:
        temp_donorA_peaks=temp('idr_analysis/macs_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/donorA_temp.narrowPeak'),
        temp_donorB_peaks=temp('idr_analysis/macs_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/donorB_temp.narrowPeak'),
        outfile='idr_analysis/macs_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/idrValues.txt',
        plotfile='idr_analysis/macs_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/idrValues.png',
    log: 'idr_analysis/macs_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/idr.log',
    version: SOFTWARE_VERSIONS['IDR']
    run:
        pick_top_peaks(input.donorA_peaks, output.temp_donorA_peaks, by='pValue', number=150000)
        pick_top_peaks(input.donorB_peaks, output.temp_donorB_peaks, by='pValue', number=150000)
        shell('''
        idr --samples {output.temp_donorA_peaks:q} {output.temp_donorB_peaks:q} \
          --input-file-type narrowPeak \
          --rank p.value \
          --output-file {output.outfile:q} \
          --output-file-type bed \
          --log-output-file {log:q} \
          --plot \
          --random-seed 1986
        mv {output.outfile:q}.png {output.plotfile:q}
        ''')

rule run_idr_epic_all_conditions:
    '''Run IDR on epic all-condition peak calls.

    https://github.com/nboley/idr

    '''
    input:
        donorA_peaks='peak_calls/epic_{genome_build}/{chip_antibody}_condition.ALL_donor.{donorA}/peaks_noBL.narrowPeak',
        donorB_peaks='peak_calls/epic_{genome_build}/{chip_antibody}_condition.ALL_donor.{donorB}/peaks_noBL.narrowPeak',
    output:
        temp_donorA_peaks=temp('idr_analysis/epic_{genome_build}/{chip_antibody}_condition.ALL_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/donorA_temp.narrowPeak'),
        temp_donorB_peaks=temp('idr_analysis/epic_{genome_build}/{chip_antibody}_condition.ALL_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/donorB_temp.narrowPeak'),
        outfile='idr_analysis/epic_{genome_build}/{chip_antibody}_condition.ALL_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/idrValues.txt',
        plotfile='idr_analysis/epic_{genome_build}/{chip_antibody}_condition.ALL_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/idrValues.png',
    log: 'idr_analysis/epic_{genome_build}/{chip_antibody}_condition.ALL_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/idr.log',
    version: SOFTWARE_VERSIONS['IDR']
    run:
        pick_top_peaks(input.donorA_peaks, output.temp_donorA_peaks, by='score', number=150000)
        pick_top_peaks(input.donorB_peaks, output.temp_donorB_peaks, by='score', number=150000)
        shell('''
        idr --samples {output.temp_donorA_peaks:q} {output.temp_donorB_peaks:q} \
          --input-file-type bed \
          --rank score \
          --output-file {output.outfile:q} \
          --output-file-type bed \
          --log-output-file {log:q} \
          --plot \
          --random-seed 1986
        mv {output.outfile:q}.png {output.plotfile:q}
        ''')

rule run_idr_epic_single_condition:
    '''Run IDR on epic single-condition peak calls.

    https://github.com/nboley/idr

    '''
    input:
        donorA_peaks='peak_calls/epic_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point}_donor.{donorA}/peaks_noBL.narrowPeak',
        donorB_peaks='peak_calls/epic_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point}_donor.{donorB}/peaks_noBL.narrowPeak',
    output:
        temp_donorA_peaks=temp('idr_analysis/epic_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/donorA_temp.narrowPeak'),
        temp_donorB_peaks=temp('idr_analysis/epic_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/donorB_temp.narrowPeak'),
        outfile='idr_analysis/epic_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/idrValues.txt',
        plotfile='idr_analysis/epic_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/idrValues.png',
    log: 'idr_analysis/epic_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point,Day[0-9]+}_{donorA,D[0-9]+}vs{donorB,D[0-9]+}/idr.log',
    version: SOFTWARE_VERSIONS['IDR']
    run:
        pick_top_peaks(input.donorA_peaks, output.temp_donorA_peaks, by='score', number=150000)
        pick_top_peaks(input.donorB_peaks, output.temp_donorB_peaks, by='score', number=150000)
        shell('''
        idr --samples {output.temp_donorA_peaks:q} {output.temp_donorB_peaks:q} \
          --input-file-type bed \
          --rank score \
          --output-file {output.outfile:q} \
          --output-file-type bed \
          --log-output-file {log:q} \
          --plot \
          --random-seed 1986
        mv {output.outfile:q}.png {output.plotfile:q}
        ''')

rule plot_idr:
    '''Reproduce IDR pots using ggplot2.'''
    input:
        'idr_analysis/{peak_caller}_{genome_build}/{chip_antibody}_condition.{condition}_{donorA}vs{donorB}/idrValues.txt'
    output:
        'plots/IDR/{peak_caller}_{genome_build}/{chip_antibody}/condition.{condition}/{donorA,D[0-9]+}vs{donorB,D[0-9]+}_idrplots.pdf'
    params:
        sampleA='{chip_antibody}-{condition}-{donorA}',
        sampleB='{chip_antibody}-{condition}-{donorB}',
        common_prefix='{chip_antibody}-{condition}',
    version: SOFTWARE_VERSIONS['R']
    shell: '''
    scripts/plot-idr.R -i {input:q} -o {output:q} \
      -A {params.sampleA:q} -B {params.sampleB:q} \
      -P {params.common_prefix:q}
    '''

rule idr_filter_peaks_one_condition:
    '''Annotate single-condition peak calls with IDR threshold.'''
    input:
        combined_peaks='peak_calls/{peak_caller}_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point}_donor.ALL/peaks_noBL.narrowPeak',
        idr_results_files=lambda wildcards:
        expand('idr_analysis/{peak_caller}_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point}_{donorPair}/idrValues.txt',
               **wildcards,
               donorPair=dfselect(idr_sample_pairs, chip_antibody=wildcards.chip_antibody,
                                  cell_type=wildcards.cell_type, time_point=wildcards.time_point) \
               .apply(lambda x: '%svs%s' % (x['donorA'], x['donorB']), axis=1))
    output:
        filtered_peaks='peak_calls/{peak_caller}_{genome_build}/{chip_antibody}_condition.{cell_type}.{time_point}_donor.ALL/peaks_noBL_IDR.narrowPeak',
    run:
        idr_results = ','.join(input.idr_results_files)
        shell('''scripts/filter-by-idr.R -p {input.combined_peaks:q} -o {output.filtered_peaks:q} -i {idr_results:q} -r''')

rule idr_filter_peaks_all_conditions:
    '''Annotate all-condition peak calls with IDR threshold.'''
    input:
        combined_peaks='peak_calls/{peak_caller}_{genome_build}/{chip_antibody}_condition.ALL_donor.ALL/peaks_noBL.narrowPeak',
        idr_results_files=lambda wildcards:
        expand('idr_analysis/{peak_caller}_{genome_build}/{chip_antibody}_condition.ALL_{donorPair}/idrValues.txt',
               **wildcards,
               donorPair=dfselect(idr_sample_pairs, chip_antibody=wildcards.chip_antibody) \
               .apply(lambda x: '%svs%s' % (x['donorA'], x['donorB']), axis=1).drop_duplicates())
    output:
        filtered_peaks='peak_calls/{peak_caller}_{genome_build}/{chip_antibody}_condition.ALL_donor.ALL/peaks_noBL_IDR.narrowPeak',
    run:
        idr_results = ','.join(input.idr_results_files)
        shell('''scripts/filter-by-idr.R -p {input.combined_peaks:q} -o {output.filtered_peaks:q} -i {idr_results:q} -r''')

rule csaw_compute_ccf:
    '''Compute cross-strand correlation values for ChIP-Seq samples.

    https://bioconductor.org/packages/release/bioc/html/csaw.html

    '''
    input:
        samplemeta='saved_data/samplemeta-ChIPSeq.RDS',
        bam_files=expand('aligned/chipseq_bowtie2_hg38.analysisSet/{sra_run}/Aligned.{ext}',
                         sra_run=chipseq_samplemeta['SRA_run'],
                         ext=['bam', 'bam.bai']),
        blacklist='saved_data/ChIPSeq-merged-blacklist.bed'
    output:
        'saved_data/csaw-ccf.RDS', 'saved_data/csaw-ccf-noBL.RDS'
    version: R_package_version('csaw')
    threads: 8
    shell: 'MC_CORES={threads:q} scripts/csaw-compute-ccf.R'

rule csaw_plot_ccf:
    '''Plot cross-strand correlation values for all ChIP-Seq samples.

    https://bioconductor.org/packages/release/bioc/html/csaw.html

    '''
    input:
        samplemeta='saved_data/samplemeta-ChIPSeq.RDS',
        ccf_data='saved_data/csaw-ccf.RDS',
        ccf_noBL_data='saved_data/csaw-ccf-noBL.RDS',
    output:
        expand('plots/csaw/CCF-plots{suffix}.pdf',
               suffix=('', '-relative', '-noBL', '-relative-noBL')),
        'plots/csaw/CCF-max-plot.pdf'
    version: SOFTWARE_VERSIONS['R']
    shell: 'scripts/csaw-plot-ccf.R'

rule csaw_profile_sites:
    '''Profile read coverage around local coverage maxima in ChIP-Seq samples.

    https://bioconductor.org/packages/release/bioc/html/csaw.html

    '''
    input:
        samplemeta='saved_data/samplemeta-ChIPSeq.RDS',
        bam_files=expand('aligned/chipseq_bowtie2_hg38.analysisSet/{sra_run}/Aligned.{ext}',
                         sra_run=chipseq_samplemeta['SRA_run'],
                         ext=['bam', 'bam.bai']),
        blacklist='saved_data/ChIPSeq-merged-blacklist.bed'
    output:
        'saved_data/csaw-sample-maxima.RDS',
        'saved_data/csaw-siteprof.RDS',
        'plots/csaw/site-profile-plots.pdf'
    version: R_package_version('csaw')
    threads: 8
    shell: 'MC_CORES={threads:q} scripts/csaw-profile-sites.R'

rule csaw_count_windows:
    '''Count ChIP-Seq reads overlapping windows in each sample.

    https://bioconductor.org/packages/release/bioc/html/csaw.html

    '''
    input:
        samplemeta='saved_data/samplemeta-ChIPSeq.RDS',
        bam_files=expand('aligned/chipseq_bowtie2_hg38.analysisSet/{sra_run}/Aligned.{ext}',
                         sra_run=chipseq_samplemeta['SRA_run'],
                         ext=['bam', 'bam.bai']),
        blacklist='saved_data/ChIPSeq-merged-blacklist.bed',
    output:
        'saved_data/csaw-counts_{window_size,[0-9.]+.*?bp}-windows_{read_ext,[0-9.]+.*?bp}-reads.RDS'
    version: R_package_version('csaw')
    threads: 1
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['csaw_count_windows']
    shell: '''
    scripts/csaw-count-windows.R \
      --samplemeta-file {input.samplemeta:q} \
      --sample-id-column SRA_run \
      --bam-file-pattern 'aligned/chipseq_bowtie2_hg38.analysisSet/%s/Aligned.bam' \
      --window-width {wildcards.window_size:q} \
      --read-extension {wildcards.read_ext:q} \
      --blacklist {input.blacklist:q} \
      --threads {threads:q} \
      --output-file {output:q}
    '''

rule csaw_count_bins:
    '''Count ChIP-Seq reads overlapping windows in each sample.

    https://bioconductor.org/packages/release/bioc/html/csaw.html

    '''
    input:
        samplemeta='saved_data/samplemeta-ChIPSeq.RDS',
        bam_files=expand('aligned/chipseq_bowtie2_hg38.analysisSet/{sra_run}/Aligned.{ext}',
                         sra_run=chipseq_samplemeta['SRA_run'],
                         ext=['bam', 'bam.bai']),
        blacklist='saved_data/ChIPSeq-merged-blacklist.bed',
    output:
        'saved_data/csaw-counts_{window_size,[0-9.]+.*?bp}-bigbins.RDS'
    version: R_package_version('csaw')
    threads: 8
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['csaw_count_bins']
    shell: '''
    scripts/csaw-count-windows.R \
      --samplemeta-file {input.samplemeta:q} \
      --sample-id-column SRA_run \
      --bam-file-pattern 'aligned/chipseq_bowtie2_hg38.analysisSet/%s/Aligned.bam' \
      --window-width {wildcards.window_size:q} \
      --blacklist {input.blacklist:q} \
      --bin \
      --threads {threads:q} \
      --output-file {output:q}
    '''

rule csaw_count_promoters:
    '''Count ChIP-Seq reads in promoters in each sample.

    https://bioconductor.org/packages/release/bioc/html/csaw.html

    '''
    input:
        samplemeta='saved_data/samplemeta-ChIPSeq.RDS',
        bam_files=expand('aligned/chipseq_bowtie2_{{genome}}/{sra_run}/Aligned.{ext}',
                         sra_run=chipseq_samplemeta['SRA_run'],
                         ext=['bam', 'bam.bai']),
        promoters='saved_data/promoter-regions_{genome}_{transcriptome}_{radius}.RDS',
        blacklist='saved_data/ChIPSeq-merged-blacklist.bed',
    output:
        'saved_data/promoter-counts_{genome,[^_]+}_{transcriptome,[^_]+}_{radius,[0-9.]+.*?bp}-radius_{read_ext,[0-9.]+.*?bp}-reads.RDS'
    version: R_package_version('csaw')
    threads: 16
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['csaw_count_regions']
    shell: '''
    scripts/csaw-count-regions.R \
      --samplemeta-file {input.samplemeta:q} \
      --sample-id-column SRA_run \
      --bam-file-pattern 'aligned/chipseq_bowtie2_hg38.analysisSet/{{SAMPLE}}/Aligned.bam' \
      --regions {input.promoters:q} \
      --read-extension {wildcards.read_ext:q} \
      --blacklist {input.blacklist:q} \
      --threads {threads:q} \
      --output-file {output:q}
    '''

rule csaw_count_peaks:
    '''Count ChIP-Seq reads in peak regions in each sample.

    https://bioconductor.org/packages/release/bioc/html/csaw.html

    '''
    input:
        samplemeta='saved_data/samplemeta-ChIPSeq.RDS',
        bam_files=lambda wildcards:
        expand('aligned/chipseq_bowtie2_{genome}/{sra_run}/Aligned.{ext}',
               genome = wildcards.genome,
               sra_run = dfselect(chipseq_samplemeta, 'SRA_run', chip_antibody=['input', wildcards.chip]),
               ext = [ 'bam', 'bam.bai' ]),
        peaks='peak_calls/{peak_caller}_{genome}/{chip}_condition.ALL_donor.ALL/peaks_noBL_IDR.narrowPeak',
        blacklist='saved_data/ChIPSeq-merged-blacklist.bed',
    output:
        'saved_data/peak-counts_{genome,[^_]+}_{peak_caller}_{chip}_{read_ext,[0-9.]+.*?bp}-reads.RDS'
    params:
        sample_id_list = lambda wildcards: ",".join(dfselect(chipseq_samplemeta, 'SRA_run', chip_antibody=['input', wildcards.chip]))
    version: R_package_version('csaw')
    threads: 16
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['csaw_count_regions']
    shell: '''
    scripts/csaw-count-regions.R \
      --samplemeta-file {input.samplemeta:q} \
      --sample-id-column SRA_run \
      --filter-sample-ids={params.sample_id_list:q} \
      --bam-file-pattern 'aligned/chipseq_bowtie2_hg38.analysisSet/{{SAMPLE}}/Aligned.bam' \
      --regions {input.peaks:q} \
      --read-extension {wildcards.read_ext:q} \
      --blacklist {input.blacklist:q} \
      --threads {threads:q} \
      --output-file {output:q}
    '''

rule split_csaw_window_counts:
    '''Split csaw window count data by histone mark.

    https://bioconductor.org/packages/release/bioc/html/csaw.html

    '''
    input:
        'saved_data/csaw-counts_{window_size}-windows_{read_ext}-reads.RDS',
    output:
        expand('saved_data/csaw-counts_{{window_size,[0-9.]+.*?bp}}-windows_{{read_ext,[0-9.]+.*?bp}}-reads_{chip}.RDS',
               chip=set(chipseq_samplemeta['chip_antibody'])),
    version: SOFTWARE_VERSIONS['BIOC']
    shell: '''
    scripts/split-sexp.R \
      -i {input:q} \
      -o 'saved_data/csaw-counts_{wildcards.window_size:q}-windows_{wildcards.read_ext:q}-reads_{{chip_antibody}}.RDS'
    '''

rule split_csaw_bigbin_counts:
    '''Split csaw bin count data by histone mark.

    https://bioconductor.org/packages/release/bioc/html/csaw.html

    '''
    input:
        'saved_data/csaw-counts_{window_size}-bigbins.RDS',
    output:
        expand('saved_data/csaw-counts_{{window_size,[0-9.]+.*?bp}}-bigbins_{chip}.RDS',
               chip=set(chipseq_samplemeta['chip_antibody'])),
    version: SOFTWARE_VERSIONS['BIOC']
    shell: '''
    scripts/split-sexp.R \
      -i {input:q} \
      -o 'saved_data/csaw-counts_{wildcards.window_size:q}-bigbins_{{chip_antibody}}.RDS'
    '''

rule split_csaw_promoter_counts:
    '''Split csaw bin count data by histone mark.

    https://bioconductor.org/packages/release/bioc/html/csaw.html

    '''
    input:
        'saved_data/promoter-counts_{base}_{read_ext}-reads.RDS',
    output:
        expand('saved_data/promoter-counts_{{base}}_{{read_ext,[0-9.]+.*?bp}}-reads_{chip}.RDS',
               chip=set(chipseq_samplemeta['chip_antibody'])),
    version: SOFTWARE_VERSIONS['BIOC']
    shell: '''
    scripts/split-sexp.R \
      -i {input:q} \
      -o 'saved_data/promoter-counts_{wildcards.base:q}_{wildcards.read_ext:q}-reads_{{chip_antibody}}.RDS'
    '''

rule collect_gene_abundance_ensembl:
    '''Generate a SummarizedExperiment object from kallisto's abundance.h5 format.

    This uses the tximport and sleuth R packages.'''
    input:
        samplemeta='saved_data/samplemeta-RNASeq.RDS',
        txdb=hg38_ref('TxDb.Hsapiens.ensembl.hg38.v{release}.sqlite3'),
        genemeta=hg38_ref('genemeta.ensembl.{release}.RDS'),
        samples=expand('{{quantifier}}_quant/hg38.analysisSet_ensembl.{{release}}/{sample}/abundance.h5',
                       sample=rnaseq_samplemeta['SRA_run']),
    output:
        sexp='saved_data/SummarizedExperiment_rnaseq_{quantifier}_hg38.analysisSet_ensembl.{release,\\d+}.RDS'
    version: (R_package_version('tximport'), R_package_version('sleuth'))
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['rnaseq_count']
    run:
        cmd = [
            'scripts/convert-quant-to-sexp.R',
            '--samplemeta-file', input.samplemeta,
            '--sample-id-column', 'SRA_run',
            '--abundance-file-pattern', *expand('{quantifier}_quant/hg38.analysisSet_ensembl.{release}/{{SAMPLE}}/abundance.h5', **wildcards),
            '--output-file', output.sexp,
            '--expected-abundance-files', ','.join(input.samples),
            '--aggregate-level', 'gene',
            '--annotation-txdb', input.txdb,
            '--gene-info', input.genemeta,
        ]
        check_call(cmd)

rule collect_gene_abundance_knownGene:
    '''Generate a SummarizedExperiment object from kallisto's abundance.h5 format.

    This uses the tximport and sleuth R packages.'''
    input:
        samplemeta='saved_data/samplemeta-RNASeq.RDS',
        genemeta=hg38_ref('genemeta.org.Hs.eg.db.RDS'),
        samples=expand('{{quantifier}}_quant/hg38.analysisSet_knownGene/{sample}/abundance.h5',
                       sample=rnaseq_samplemeta['SRA_run']),
    output:
        sexp='saved_data/SummarizedExperiment_rnaseq_{quantifier}_hg38.analysisSet_knownGene.RDS',
    version: (R_package_version('tximport'), R_package_version('sleuth'))
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['rnaseq_count']
    run:
        cmd = [
            'scripts/convert-quant-to-sexp.R',
            '--samplemeta-file', input.samplemeta,
            '--sample-id-column', 'SRA_run',
            '--abundance-file-pattern', *expand('{quantifier}_quant/hg38.analysisSet_knownGene/{{SAMPLE}}/abundance.h5', **wildcards),
            '--output-file', output.sexp,
            '--expected-abundance-files', ','.join(input.samples),
            '--aggregate-level', 'gene',
            '--annotation-txdb', 'TxDb.Hsapiens.UCSC.hg38.knownGene',
            '--gene-info', input.genemeta,
        ]
        check_call(cmd)

rule collect_shoal_gene_abundance_ensembl:
    '''Generate a SummarizedExperiment object from kallisto's abundance.h5 format.

    This uses the tximport and sleuth R packages.'''
    input:
        samplemeta='saved_data/samplemeta-RNASeq.RDS',
        txdb=hg38_ref('TxDb.Hsapiens.ensembl.hg38.v{release}.sqlite3'),
        genemeta=hg38_ref('genemeta.ensembl.{release}.RDS'),
        samples=expand('shoal_quant/hg38.analysisSet_ensembl.{{release}}/{sample}_adapt.sf',
                       sample=rnaseq_samplemeta['SRA_run']),
    output:
        sexp='saved_data/SummarizedExperiment_rnaseq_shoal_hg38.analysisSet_ensembl.{release,\\d+}.RDS'
    version: (R_package_version('tximport'), R_package_version('sleuth'))
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['rnaseq_count']
    run:
        cmd = [
            'scripts/convert-shoal-to-sexp.R',
            '--samplemeta-file', input.samplemeta,
            '--sample-id-column', 'SRA_run',
            '--shoal-dir', *expand('shoal_quant/hg38.analysisSet_ensembl.{release}', **wildcards),
            '--output-file', output.sexp,
            '--aggregate-level', 'gene',
            '--annotation-txdb', input.txdb,
            '--gene-info', input.genemeta,
        ]
        check_call(cmd)

rule collect_shoal_gene_abundance_knownGene:
    '''Generate a SummarizedExperiment object from kallisto's abundance.h5 format.

    This uses the tximport and sleuth R packages.'''
    input:
        samplemeta='saved_data/samplemeta-RNASeq.RDS',
        genemeta=hg38_ref('genemeta.org.Hs.eg.db.RDS'),
        samples=expand('shoal_quant/hg38.analysisSet_knownGene/{sample}_adapt.sf',
                       sample=rnaseq_samplemeta['SRA_run']),
    output:
        sexp='saved_data/SummarizedExperiment_rnaseq_shoal_hg38.analysisSet_knownGene.RDS',
    version: (R_package_version('tximport'), R_package_version('sleuth'))
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['rnaseq_count']
    run:
        cmd = [
            'scripts/convert-shoal-to-sexp.R',
            '--samplemeta-file', input.samplemeta,
            '--sample-id-column', 'SRA_run',
            '--shoal-dir', *expand('shoal_quant/hg38.analysisSet_knownGene', **wildcards),
            '--output-file', output.sexp,
            '--aggregate-level', 'gene',
            '--annotation-txdb', 'TxDb.Hsapiens.UCSC.hg38.knownGene',
            '--gene-info', input.genemeta,
        ]
        check_call(cmd)

rule collect_transcript_abundance_ensembl:
    '''Generate a SummarizedExperiment object from kallisto's abundance.h5 format.

    This uses the tximport and sleuth R packages.'''
    input:
        samplemeta='saved_data/samplemeta-RNASeq.RDS',
        txdb=hg38_ref('TxDb.Hsapiens.ensembl.hg38.v{release}.sqlite3'),
        samples=expand('{{quantifier}}_quant/hg38.analysisSet_ensembl.{{release}}/{sample}/abundance.h5',
                       sample=rnaseq_samplemeta['SRA_run']),
    output:
        sexp='saved_data/SummarizedExperiment_rnaseq_transcript_{quantifier}_hg38.analysisSet_ensembl.{release,\\d+}.RDS'
    version: (R_package_version('tximport'), R_package_version('sleuth'))
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['rnaseq_count']
    run:
        cmd = [
            'scripts/convert-quant-to-sexp.R',
            '--samplemeta-file', input.samplemeta,
            '--sample-id-column', 'SRA_run',
            '--abundance-file-pattern', *expand('{quantifier}_quant/hg38.analysisSet_ensembl.{release}/{{SAMPLE}}/abundance.h5', **wildcards),
            '--output-file', output.sexp,
            '--expected-abundance-files', ','.join(input.samples),
            '--aggregate-level', 'transcript',
            '--annotation-txdb', input.txdb,
        ]
        check_call(cmd)

rule collect_transcript_abundance_knownGene:
    '''Generate a SummarizedExperiment object from kallisto's abundance.h5 format.

    This uses the tximport and sleuth R packages.'''
    input:
        samplemeta='saved_data/samplemeta-RNASeq.RDS',
        samples=expand('{{quantifier}}_quant/hg38.analysisSet_knownGene/{sample}/abundance.h5',
                       sample=rnaseq_samplemeta['SRA_run']),
    output:
        sexp='saved_data/SummarizedExperiment_rnaseq_transcript_{quantifier}_hg38.analysisSet_knownGene.RDS',
    version: (R_package_version('tximport'), R_package_version('sleuth'))
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['rnaseq_count']
    run:
        cmd = [
            'scripts/convert-quant-to-sexp.R',
            '--samplemeta-file', input.samplemeta,
            '--sample-id-column', 'SRA_run',
            '--abundance-file-pattern', *expand('{quantifier}_quant/hg38.analysisSet_knownGene/{{SAMPLE}}/abundance.h5', **wildcards),
            '--output-file', output.sexp,
            '--expected-abundance-files', ','.join(input.samples),
            '--aggregate-level', 'transcript',
            '--annotation-txdb', 'TxDb.Hsapiens.UCSC.hg38.knownGene',
        ]
        check_call(cmd)

rule collect_shoal_transcript_abundance_ensembl:
    '''Generate a SummarizedExperiment object from kallisto's abundance.h5 format.

    This uses the tximport and sleuth R packages.'''
    input:
        samplemeta='saved_data/samplemeta-RNASeq.RDS',
        txdb=hg38_ref('TxDb.Hsapiens.ensembl.hg38.v{release}.sqlite3'),
        samples=expand('shoal_quant/hg38.analysisSet_ensembl.{{release}}/{sample}_adapt.sf',
                       sample=rnaseq_samplemeta['SRA_run']),
    output:
        sexp='saved_data/SummarizedExperiment_rnaseq_transcript_shoal_hg38.analysisSet_ensembl.{release,\\d+}.RDS'
    version: (R_package_version('tximport'), R_package_version('sleuth'))
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['rnaseq_count']
    run:
        cmd = [
            'scripts/convert-shoal-to-sexp.R',
            '--samplemeta-file', input.samplemeta,
            '--sample-id-column', 'SRA_run',
            '--shoal-dir', *expand('shoal_quant/hg38.analysisSet_ensembl.{release}', **wildcards),
            '--output-file', output.sexp,
            '--aggregate-level', 'transcript',
            '--annotation-txdb', input.txdb,
        ]
        check_call(cmd)

rule collect_shoal_transcript_abundance_knownGene:
    '''Generate a SummarizedExperiment object from kallisto's abundance.h5 format.

    This uses the tximport and sleuth R packages.'''
    input:
        samplemeta='saved_data/samplemeta-RNASeq.RDS',
        samples=expand('shoal_quant/hg38.analysisSet_knownGene/{sample}_adapt.sf',
                       sample=rnaseq_samplemeta['SRA_run']),
    output:
        sexp='saved_data/SummarizedExperiment_rnaseq_transcript_shoal_hg38.analysisSet_knownGene.RDS',
    version: (R_package_version('tximport'), R_package_version('sleuth'))
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['rnaseq_count']
    run:
        cmd = [
            'scripts/convert-shoal-to-sexp.R',
            '--samplemeta-file', input.samplemeta,
            '--sample-id-column', 'SRA_run',
            '--shoal-dir', *expand('shoal_quant/hg38.analysisSet_knownGene', **wildcards),
            '--output-file', output.sexp,
            '--aggregate-level', 'transcript',
            '--annotation-txdb', 'TxDb.Hsapiens.UCSC.hg38.knownGene',
        ]
        check_call(cmd)

rule rnaseq_explore:
    '''Perform exploratory data analysis on RNA-seq dataset'''
    input:
        rmd='scripts/rnaseq-explore.Rmd',
        sexp='saved_data/SummarizedExperiment_rnaseq_{dataset}.RDS',
    output:
        html='reports/RNA-seq/{dataset}-exploration.html',
        plots=expand('plots/RNA-seq/{{dataset}}/{plotfile}',
                     plotfile=['AveLogCPM-plots.pdf',
                               'disp-plots.pdf', 'qc-weights.pdf',
                               'rnaseq-ComBat-qc.pdf', 'rnaseq-MDSPlots.pdf',
                               'rnaseq-MDSPlots-BatchCorrect.pdf',
                               'weights-vs-covars.pdf']),
    version: R_package_version('rmarkdown')
    threads: 2
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['rnaseq_analyze']
    run:
        os.environ['MC_CORES'] = str(threads)
        rmd_render(input=input.rmd, output_file=os.path.join(os.getcwd(), output.html),
                   output_format='html_notebook',
                   params={ 'dataset': wildcards.dataset, })

rule rnaseq_compare:
    '''Perform basic comparisons between RNA-seq quantification methods'''
    input:
        rmd='scripts/rnaseq-compare.Rmd',
        sexps=expand('saved_data/SummarizedExperiment_rnaseq_{aligner_and_genome}_{transcriptome}.RDS',
                           aligner_and_genome=[
                               'star_hg38.analysisSet', 'hisat2_grch38_snp_tran', 'kallisto_hg38.analysisSet',
                               'salmon_hg38.analysisSet', 'shoal_hg38.analysisSet'
                           ],
                           transcriptome=['ensembl.85','knownGene']),
    output:
        html='reports/RNA-seq/rnaseq-compare.html',
    version: R_package_version('rmarkdown')
    threads: 10
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['rnaseq_analyze']
    run:
        os.environ['MC_CORES'] = str(threads)
        rmd_render(input=input.rmd, output_file=os.path.join(os.getcwd(), output.html),
                   output_format='html_notebook')

rule rnaseq_diffexp:
    '''Perform differential expression analysis on RNA-seq dataset'''
    input:
        rmd='scripts/rnaseq-diffexp.Rmd',
        sexp='saved_data/SummarizedExperiment_rnaseq_{quant_method}_{genome}_{transcriptome}.RDS',
    output:
        html='reports/RNA-seq/{quant_method,[^_]+}_{genome}_{transcriptome}-diffexp.html',
        table='results/RNA-seq/{quant_method,[^_]+}_{genome}_{transcriptome}-diffexp.xlsx',
        rds='saved_data/RNA-seq/{quant_method,[^_]+}_{genome}_{transcriptome}-diffexp-tables.RDS',
        rda='saved_data/RNA-seq/{quant_method,[^_]+}_{genome}_{transcriptome}-diffexp.rda',
    version: (R_package_version('rmarkdown'), R_package_version('limma'))
    threads: 2
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['rnaseq_analyze']
    run:
        os.environ['MC_CORES'] = str(threads)
        rmd_render(input=input.rmd, output_file=os.path.join(os.getcwd(), output.html),
                   output_format='html_notebook',
                   params={
                       'quant_method': wildcards.quant_method,
                       'genome': wildcards.genome,
                       'transcriptome': wildcards.transcriptome,
                   })

rule chipseq_peak_size_analysis:
    input:
        peaks=expand('peak_calls/epic_hg38.analysisSet/{chip}_condition.ALL_donor.ALL/peaks_noBL_IDR.narrowPeak',
                     chip=set(chipseq_samplemeta_noinput['chip_antibody'])),
        rmd='scripts/chipseq-peak-sizes.Rmd'
    output:
        html='reports/ChIP-seq/peak-sizes.html'
    version: R_package_version('rmarkdown')
    threads: len(set(chipseq_samplemeta_noinput['chip_antibody']))
    run:
        os.environ['MC_CORES'] = str(threads)
        rmd_render(input=input.rmd, output_file=os.path.join(os.getcwd(), output.html),
                   output_format='html_notebook')

rule chipseq_explore:
    '''Perform exploratory data analysis on ChIP-seq dataset'''
    input:
        rmd='scripts/chipseq-explore-{chip_antibody}.Rmd',
        sexp='saved_data/csaw-counts_500bp-windows_147bp-reads_{chip_antibody}.RDS',
        bigbin_sexp='saved_data/csaw-counts_10kbp-bigbins_{chip_antibody}.RDS',
        peaks='peak_calls/epic_hg38.analysisSet/{chip_antibody}_condition.ALL_donor.ALL/peaks_noBL_IDR.narrowPeak',
    output:
        html='reports/ChIP-seq/{chip_antibody}-exploration.html',
    version: R_package_version('rmarkdown')
    threads: 4
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['chipseq_analyze']
    run:
        os.environ['MC_CORES'] = str(threads)
        rmd_render(input=input.rmd, output_file=os.path.join(os.getcwd(), output.html),
                   output_format='html_notebook',
                   params={
                       'histone_mark': wildcards.chip_antibody,
                       'window_size': '500bp',
                       'fragment_length': '147bp',
                       'bigbin_size': '10kbp',
                   })

rule chipseq_diffmod:
    '''Perform differential modification analysis on ChIP-seq dataset'''
    input:
        rmd='scripts/chipseq-diffmod.Rmd',
        sexp='saved_data/csaw-counts_500bp-windows_147bp-reads_{chip_antibody}.RDS',
        peaks='peak_calls/epic_hg38.analysisSet/{chip_antibody}_condition.ALL_donor.ALL/peaks_noBL_IDR.narrowPeak',
    output:
        html='reports/ChIP-seq/{chip_antibody}-diffmod.html',
        peaks_xlsx='results/ChIP-seq/{chip_antibody}-peak-diffmod.xlsx',
        peaks_rds='saved_data/ChIP-seq/{chip_antibody}-peak-diffmod.RDS',
        windows_rds='saved_data/ChIP-seq/{chip_antibody}-window-diffmod.RDS',
        rda='saved_data/ChIP-seq/{chip_antibody}-diffmod.rda',
    version: R_package_version('rmarkdown')
    threads: 4
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['chipseq_analyze']
    run:
        os.environ['MC_CORES'] = str(threads)
        rmd_render(input=input.rmd, output_file=os.path.join(os.getcwd(), output.html),
                   output_format='html_notebook',
                   params={
                       'histone_mark': wildcards.chip_antibody,
                       'window_size': '500bp',
                       'fragment_length': '147bp',
                   })


rule chipseq_NvM_diminish_analysis:
    '''Test the hypothesis of diminishing NvM differences over time.'''
    input:
        rmd='scripts/chipseq-NvM-diminish-{chip_antibody}.Rmd',
        rda='saved_data/ChIP-seq/{chip_antibody}-diffmod.rda',
    output:
        html='reports/ChIP-seq/{chip_antibody}-NvM-diminish-analysis.html',
    version: R_package_version('rmarkdown')
    threads: 1
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['chipseq_analyze']
    run:
        rmd_render(input=input.rmd, output_file=os.path.join(os.getcwd(), output.html))

rule chipseq_promoter_explore:
    '''Perform exploratory data analysis on ChIP-seq dataset'''
    input:
        rmd='scripts/chipseq-promoter-explore-{chip_antibody}.Rmd',
        sexp='saved_data/promoter-counts_{genome}_{transcriptome}_{promoter_radius}-radius_147bp-reads_{chip_antibody}.RDS'
    output:
        html='reports/ChIP-seq/{genome,[^_]+}_{transcriptome,[^_]+}_{chip_antibody}-{promoter_radius,[0-9.]+\\w*?bp}-promoter-exploration.html',
    version: R_package_version('rmarkdown')
    threads: 4
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['chipseq_analyze']
    run:
        os.environ['MC_CORES'] = str(threads)
        rmd_render(input=input.rmd, output_file=os.path.join(os.getcwd(), output.html),
                   output_format='html_notebook',
                   params={
                       'genome': wildcards.genome,
                       'transcriptome': wildcards.transcriptome,
                       'histone_mark': wildcards.chip_antibody,
                       'promoter_radius': wildcards.promoter_radius,
                       'fragment_length': '147bp',
                       'bigbin_size': '10kbp',
                   })

rule chipseq_promoter_diffmod:
    '''Perform differential modification analysis on ChIP-seq dataset'''
    input:
        rmd='scripts/chipseq-promoter-diffmod.Rmd',
        sexp='saved_data/promoter-counts_{genome}_{transcriptome}_{promoter_radius}-radius_147bp-reads_{chip_antibody}.RDS'
    output:
        html='reports/ChIP-seq/{genome,[^_]+}_{transcriptome,[^_]+}_{chip_antibody}-{promoter_radius,[0-9.]+\\w*?bp}-promoter-diffmod.html',
        xlsx='results/ChIP-seq/{genome,[^_]+}_{transcriptome,[^_]+}_{chip_antibody}-{promoter_radius,[0-9.]+\\w*?bp}-promoter-diffmod.xlsx',
        rds='saved_data/ChIP-seq/{genome,[^_]+}_{transcriptome,[^_]+}_{chip_antibody}-{promoter_radius,[0-9.]+\\w*?bp}-promoter-diffmod.RDS',
        rda='saved_data/ChIP-seq/{genome,[^_]+}_{transcriptome,[^_]+}_{chip_antibody}-{promoter_radius,[0-9.]+\\w*?bp}-promoter-diffmod.rda',
    version: R_package_version('rmarkdown')
    threads: 4
    resources: mem_gb=MEMORY_REQUIREMENTS_GB['chipseq_analyze']
    run:
        os.environ['MC_CORES'] = str(threads)
        rmd_render(input=input.rmd,
                   output_file=os.path.join(os.getcwd(), output.html),
                   output_format='html_notebook',
                   params={
                       'genome': wildcards.genome,
                       'transcriptome': wildcards.transcriptome,
                       'histone_mark': wildcards.chip_antibody,
                       'promoter_radius': wildcards.promoter_radius,
                       'fragment_length': '147bp',
                   })

rule lamere_2016_fig7:
    '''Reproduce Lamere 2016 Figure 7'''
    input:
        rmd='scripts/promoter-cpg-figure7.Rmd',
        rda_h3k4me2='saved_data/ChIP-seq/hg38.analysisSet_ensembl.85_H3K4me2-1kbp-promoter-diffmod.rda',
        rda_h3k4me3='saved_data/ChIP-seq/hg38.analysisSet_ensembl.85_H3K4me3-1kbp-promoter-diffmod.rda',
        rda_rnaseq='saved_data/RNA-seq/hisat2_grch38_snp_tran_ensembl.85-diffexp.rda',
        cpgi='saved_data/UCSC_hg38_cpgIslandExtUnmasked.RDS',
    output:
        html='reports/lamere_2016_fig7.html'
    threads: 1
    run:
        rmd_render(input=input.rmd,
                   output_file=os.path.join(os.getcwd(), output.html),
                   output_format='html_notebook')

# Potential TFBS sources:
# JASPAR:
# ORegAnno
# http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=647218275_SuDu8UNSQBfBbQQBa7YlbaviM9r4&c=chr11&g=oreganno
# UCSC TFBS Conserved track (generated from TRANSFAC):
# http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=647218275_SuDu8UNSQBfBbQQBa7YlbaviM9r4&c=chr11&g=tfbsConsSites
rule get_jaspar_all:
    '''Download JASPAR predicted TFBS.

    http://jaspar.genereg.net/genome-tracks/
    http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=647218737_MhrHME1Kg3Iy3wRbHaaIa9F3L7GJ&c=chr1&g=hub_186875_JasparTFBS
    http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2018/

    '''
    input:
        HTTP.remote('expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2018/hg38/JASPAR2018_hg38_all_chr.bed.gz',
                    insecure=True, static=True)
    output: 'saved_data/JASPAR2018_hg38_all_chr.bed.gz'
    shell: '''cp {input:q} {output:q}'''

rule filter_jaspar:
    input: 'saved_data/JASPAR2018_hg38_all_chr.bed.gz'
    output: 'saved_data/JASPAR2018_hg38_all_chr_score{score_threshold}.bed'
    shell: '''zcat {input:q} | perl -lane 'print if $A[4] >= {wildcards.score_threshold};' > {output:q}'''

rule get_tfbs_conserved:
    '''Download UCSC TFBS Conserved track (from hg19).

    http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=647218275_SuDu8UNSQBfBbQQBa7YlbaviM9r4&g=tfbsConsSites&hgTracksConfigPage=configure

    '''
    # Note: The other inputs are UCSC tables (via rtracklayer), not files
    input: chain='saved_data/hg19ToHg38.over.chain'
    output: 'saved_data/UCSC_tfbsCons.RDS'
    script: 'scripts/get-tfbs-conserved.R'

# rule get_oreganno_tfbs:
#     pass

rule generate_tfbs_overlap:
    input:
        rmd='scripts/tfbs-overlap.Rmd',
        tfbs_data='saved_data/UCSC_tfbsCons.RDS',
        promoter_sexps=set(
            expand('saved_data/promoter-counts_hg38.analysisSet_ensembl.85_{promoter_radius}-radius_147bp-reads_{chip_antibody}.RDS',
                   zip, **dict(chipseq_samplemeta_noinput))),
        peak_sexps=set(
            expand('saved_data/peak-counts_hg38.analysisSet_epic_{chip_antibody}_147bp-reads.RDS',
                   zip, **dict(chipseq_samplemeta_noinput))),
    output:
        html='reports/tfbs-overlap.html',
        promoter_tfbs_overlap='saved_data/promoter-tfbs-overlap_hg38.analysisSet_ensembl.85.RDS',
        peak_tfbs_overlap='saved_data/peak-tfbs-overlap_hg38.analysisSet.RDS'
    threads: 3
    run:
        os.environ['MC_CORES'] = str(threads)
        rmd_render(input=input.rmd,
                   output_file=os.path.join(os.getcwd(), output.html),
                   output_format='html_notebook')

rule run_mofa_promoter:
    '''Run MOFA using RNA-seqs and ChIP-Seq promoter data.'''
    input:
        rmd='scripts/promoter-mofa-run.Rmd',
        promoter_sexps=set(
            expand('saved_data/promoter-counts_hg38.analysisSet_ensembl.85_{promoter_radius}-radius_147bp-reads_{chip_antibody}.RDS',
                   zip, **dict(chipseq_samplemeta_noinput))),
        rna_sexp='saved_data/SummarizedExperiment_rnaseq_shoal_hg38.analysisSet_ensembl.85.RDS',
    output:
        html='reports/promoter-mofa-run.html',
        final_model=expand('saved_data/mofa/mofa-model_hg38.analysisSet_ensembl.85_rna+promoter.{ext}',
                           ext={'hdf5', 'RDS'}),
        test_models=expand('saved_data/mofa/mofa-model_hg38.analysisSet_ensembl.85_rna+promoter_explore{i}.{ext}',
                           i=range(1,5),
                           ext={'hdf5', 'RDS'}),
    threads: 8
    run:
        os.environ['MC_CORES'] = str(threads)
        rmd_render(input=input.rmd,
                   output_file=os.path.join(os.getcwd(), output.html),
                   output_format='html_notebook')

rule run_mofa_peak:
    '''Run MOFA using RNA-seqs and ChIP-Seq peak data.'''
    input:
        rmd='scripts/peak-mofa-run.Rmd',
        peak_sexps=set(
            expand('saved_data/peak-counts_hg38.analysisSet_epic_{chip_antibody}_147bp-reads.RDS',
                   zip, **dict(chipseq_samplemeta_noinput))),
        rna_sexp='saved_data/SummarizedExperiment_rnaseq_shoal_hg38.analysisSet_ensembl.85.RDS',
    output:
        html='reports/peak-mofa-run.html',
        final_model=expand('saved_data/mofa/mofa-model_hg38.analysisSet_ensembl.85_rna+peak.{ext}',
                           ext={'hdf5', 'RDS'}),
        test_models=expand('saved_data/mofa/mofa-model_hg38.analysisSet_ensembl.85_rna+peak_explore{i}.{ext}',
                           i=range(1,5),
                           ext={'hdf5', 'RDS'}),
    threads: 8
    run:
        os.environ['MC_CORES'] = str(threads)
        rmd_render(input=input.rmd,
                   output_file=os.path.join(os.getcwd(), output.html),
                   output_format='html_notebook')

rule analyze_mofa_promoter:
    '''Analyze MOFA results.'''
    input:
        rmd = 'scripts/promoter-mofa-analyze.Rmd',
        mofa_model = 'saved_data/mofa/mofa-model_hg38.analysisSet_ensembl.85_rna+promoter.RDS',
        tfbs_overlap_sets = 'saved_data/promoter-tfbs-overlap_hg38.analysisSet_ensembl.85.RDS',
    output:
        html='reports/promoter-mofa-analyze.html'
    threads: 1
    run:
        os.environ['MC_CORES'] = str(threads)
        rmd_render(input=input.rmd,
                   output_file=os.path.join(os.getcwd(), output.html),
                   output_format='html_notebook')

rule analyze_mofa_peak:
    '''Analyze MOFA results.'''
    input:
        rmd = 'scripts/peak-mofa-analyze.Rmd',
        mofa_model = 'saved_data/mofa/mofa-model_hg38.analysisSet_ensembl.85_rna+peak.RDS',
        tfbs_overlap_sets = 'saved_data/peak-tfbs-overlap_hg38.analysisSet.RDS',
    output:
        html='reports/peak-mofa-analyze.html'
    threads: 1
    run:
        os.environ['MC_CORES'] = str(threads)
        rmd_render(input=input.rmd,
                   output_file=os.path.join(os.getcwd(), output.html),
                   output_format='html_notebook')

rule prepare_msigdb:
    '''Prepare R data files from the MSigDB XML file.

    Note that this XML file must be manually downloaded, since it
    requires a login.'''
    input:
        xml='saved_data/msigdb_v6.1.xml'
    output:
        xml_cache='saved_data/msigdb_v6.1.RDS',
        sets=expand("saved_data/msigdb-{identifier}.RDS",
                    identifier=['symbol', 'entrez', 'ensembl']),
        meta='saved_data/msigdb-meta.RDS',
    shell: '''Rscript scripts/prepare-msigdb.R'''

rule prepare_graphite:
    '''Prepare graphite-provided pathways for analysis.'''
    # No input, since the graphite package already includes the data
    output:
        DB=expand('saved_data/graphite-{idtype}.RDS', idtype=['symbol', 'entrez', 'ensembl']),
        SPIA=expand('saved_data/SPIA/graphite-{idtype}-{dbname}ExSPIA.RData',
                    idtype=['symbol', 'entrez', 'ensembl'],
                    dbname = list(r('graphite:::.dbs$hsapiens')))
    shell: '''Rscript scripts/prepare-graphite.R'''

rule rnaseq_gst:
    '''Run gene set tests on RNA-seq results.'''
    input:
        rmd='scripts/rnaseq-gene-set-tests.Rmd',
        diffexp='saved_data/RNA-seq/shoal_hg38.analysisSet_ensembl.85-diffexp.rda',
        msigdb='saved_data/msigdb-ensembl.RDS',
        graphite='saved_data/graphite-ensembl.RDS',
        tfbs_overlap='saved_data/promoter-tfbs-overlap_hg38.analysisSet_ensembl.85.RDS'
    output:
        rds='saved_data/CAMERA-results-RNA.RDS',
        xlsx=expand('results/RNA-seq/CAMERA-results-{collection}.xlsx',
                    collection=[
                        "MSigDB.h", "MSigDB.c2.CP", "MSigDB.c2.CGP", "MSigDB.c3.MIR",
                        "MSigDB.c3.TFT", "MSigDB.c5", "MSigDB.c7", 'TFBS_overlap', 'graphite',
                    ]),
    threads: 9
    run:
        os.environ['MC_CORES'] = str(threads)
        rmd_run_without_rendering(input.rmd)

rule promoter_gst:
    '''Run gene set tests on RNA-seq results.'''
    input:
        rmd='scripts/promoter-gene-set-tests.Rmd',
        diffmod='saved_data/ChIP-seq/hg38.analysisSet_ensembl.85_{histone_mark}-{promoter_radius}-promoter-diffmod.rda',
        msigdb='saved_data/msigdb-ensembl.RDS',
        graphite='saved_data/graphite-ensembl.RDS',
        tfbs_overlap='saved_data/promoter-tfbs-overlap_hg38.analysisSet_ensembl.85.RDS'
    output:
        rds='saved_data/CAMERA-results-{histone_mark}-{promoter_radius}-promoter.RDS',
        xlsx=expand('results/ChIP-seq/CAMERA-results-{{histone_mark}}-{{promoter_radius}}-promoter-{collection}.xlsx',
                    collection=[
                        "MSigDB.h", "MSigDB.c2.CP", "MSigDB.c2.CGP", "MSigDB.c3.MIR",
                        "MSigDB.c3.TFT", "MSigDB.c5", "MSigDB.c7", 'TFBS_overlap', 'graphite',
                    ]),
    threads: 9
    run:
        os.environ['MC_CORES'] = str(threads)
        rmd_run_without_rendering(
            input.rmd,
            params={
                'genome': 'hg38.analysisSet',
                'transcriptome': 'ensembl.85',
                'histone_mark': wildcards.histone_mark,
                'promoter_radius': wildcards.promoter_radius,
                'fragment_length': '147bp',
            })

rule select_abundant_tss_ensembl:
    '''Select the most abundant TSS for each gene.'''
    input:
        txdb=hg38_ref('TxDb.Hsapiens.ensembl.hg38.v{release}.sqlite3'),
        sexp='saved_data/SummarizedExperiment_rnaseq_transcript_{quant_method}_{genome}_ensembl.{release,\\d+}.RDS',
        genemeta=hg38_ref('genemeta.ensembl.{release}.RDS'),
    output:
        tss='saved_data/tss_{quant_method,[^_]+}_{genome,[^_]+}_ensembl.{release,\\d+}.RDS'
    shell: '''
    Rscript scripts/select-abundant-tss.R \
      --transcript-quant {input.sexp:q} \
      --annotation-txdb {input.txdb:q} \
      --additional-gene-info {input.genemeta:q} \
      --output-file {output.tss:q}
    '''

rule select_abundant_tss_knownGene:
    '''Select the most abundant TSS for each gene.'''
    input:
        sexp='saved_data/SummarizedExperiment_rnaseq_transcript_{quant_method}_{genome}_knownGene.RDS',
        genemeta=hg38_ref('genemeta.org.Hs.eg.db.RDS'),
    output:
        tss='saved_data/tss_{quant_method,[^_]+}_{genome,[^_]+}_knownGene.RDS'
    shell: '''
    Rscript scripts/select-abundant-tss.R \
      --transcript-quant {input.sexp:q} \
      --annotation-txdb 'TxDb.Hsapiens.UCSC.hg38.knownGene' \
      --additional-gene-info {input.genemeta:q} \
      --output-file {output.tss:q}
    '''

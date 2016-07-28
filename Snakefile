import contextlib
import os
import os.path
import shutil
import subprocess

from subprocess import check_call, Popen, PIPE, CalledProcessError, list2cmdline
from atomicwrites import atomic_write, AtomicWriter
from rpy2 import robjects
from rpy2.robjects import pandas2ri

pandas2ri.activate()

# Run a separate Snakemake workflow to fetch the sample metadata,
# which must be avilable before evaluating the rules below
from processify import processify
from snakemake import snakemake
snakemake = processify(snakemake)
snakemake('pre.Snakefile', targets=expand(os.path.join('saved_data', 'samplemeta-{dataset}.RDS'), dataset=('RNASeq', 'ChIPSeq')))

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

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
def ensure_empty_dir(path):
    if os.path.exists(path):
        shutil.rmtree(path)
    ensure_dir(path)

def Popen_pipeline(cmds, stdin=None, stdout=None, *args, **kwargs):
    '''Popen a pipeline of several commands.

    Returns a list with all the process objects returned by Popen.

    Each command's stdout becomes the next command's stdin. The stdin
    argument becomes the stdin of the first command, while the stdout
    argument becomes the stdout of the last command. All other
    arguments are passed to every invocation of Popen, so ensure that
    they make sense in that context.

    '''
    cmds = list(cmds)
    if len(cmds) == 0:
        raise ValueError("Cannot run a pipeline with zero commands")
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

def read_R_dataframe(rdsfile):
    '''Read an R data frame stored in an RDS file.

    Returns the result as a Pandas DataFrame.

    '''
    readRDS = robjects.r('readRDS')
    df = readRDS((robjects.StrVector([rdsfile])))
    return(pandas2ri.ri2py(df))

subworkflow hg38_ref:
    workdir: os.path.expanduser("~/references/hg38")

rule all:
    output: "temp"
    shell: 'false'

rule fetch_sra_run:
    '''Script to fetch the .sra file for an SRA run

    (An SRA run identifier starts with SRR.)

    '''
    output: 'sra_files/{sra_run}.sra'
    shell: 'scripts/get-sra-run-files.R {wildcards.sra_run:q}'
    resources: concurrent_downloads=1

rule extract_fastq:
    '''Extract FASTQ from SRA files.'''
    input: 'sra_files/{sra_run}.sra'
    output: 'fastq_files/{sra_run,SRR\\d+}.{fqext,fq(|\\.gz|\\.bz2|\\.qp)}'
    run:
        cmds = [
            ['fastq-dump', '--stdout', input[0]],
            ['scripts/fill-in-empty-fastq-qual.py'],
            fastq_compression_cmds[wildcards.fqext]['compress'],
        ]
        with atomic_write(output[0], mode="wb", overwrite=True) as outfile:
            pipeline = Popen_pipeline(cmds, stdout=outfile)
            retcodes = [ p.wait() for p in pipeline ]
            for retcode, cmd in zip(retcodes, cmds):
                if retcode != 0:
                    raise CalledProcessError(rercode, cmd)

# TODO: Transcriptome-space BAM and gene counts
rule align_rnaseq_with_star_single_end:
    '''Align fastq file with star'''
    input: fastq='fastq_files/{samplename}.fq.gz',
           index_sa=hg38_ref('STAR_index_hg38.analysisSet_knownGene/SA'),
           transcriptome_gff=hg38_ref('knownGene.gff'),
    output: bam='aligned/rnaseq_star_hg38.analysisSet_knownGene/{samplename}/Aligned.sortedByCoord.out.bam',
            sj='aligned/rnaseq_star_hg38.analysisSet_knownGene/{samplename}/SJ.out.tab',
    threads: 8
    run:
        index_genomedir = os.path.dirname(input.index_sa)
        outdir = os.path.dirname(output.bam) + os.path.sep
        ensure_dir(outdir)
        read_cmd = list2cmdline(fastq_compression_cmds['fq.gz']['decompress'])
        star_cmd = [
            'STAR',
            '--runThreadN', threads,
            '--runMode', 'alignReads',
            '--genomeDir', index_genomedir,
            '--sjdbGTFfile', input.transcriptome_gff,
            '--sjdbGTFfeatureExon', 'CDS',
            '--sjdbGTFtagExonParentTranscript', 'Parent',
            '--sjdbOverhang', '100',
            '--readFilesIn', input.fastq,
            '--readFilesCommand', read_cmd,
            '--outSAMattributes', 'Standard',
            '--outSAMunmapped', 'Within',
            '--outFileNamePrefix', outdir,
            '--outSAMtype', 'BAM', 'SortedByCoordinate',
        ]
        # Run STAR
        shell(list2cmdline(map(str, star_cmd)))

rule align_rnaseq_with_hisat2_single_end:
    '''Align fastq file with HISAT2'''
    input: fastq='fastq_files/{samplename}.fq.gz', index_f1=hg38_ref('HISAT2_index_grch38_snp_tran/index.1.ht2'), transcriptome_gff=hg38_ref('knownGene.gff'),
    output: bam='aligned/rnaseq_hisat2_grch38_snp_tran/{samplename}/Aligned.bam', log='aligned/rnaseq_hisat2_grch38_snp_tran/{samplename}/hisat2.log'
    threads: 8
    run:
        index_basename = re.sub('\\.1\\.ht2', "", input.index_f1)
        outdir = os.path.dirname(output.bam)
        ensure_dir(outdir)
        hisat2_cmd = [
            'hisat2',
            '--threads', str(threads),
            '-q', '--phred33',
            '--very-sensitive',
            '--dta-cufflinks',
            '-x', index_basename,
            '-U', input.fastq,
            '-k', '20',
            '--time',
        ]
        bam_sort_cmd = [
            'picard-tools', 'SortSam', 'I=/dev/stdin', 'O=/dev/stdout',
            'SORT_ORDER=coordinate', 'VALIDATION_STRINGENCY=LENIENT',
        ]
        with atomic_write(output.bam, mode="wb", overwrite=True) as outfile, \
             open(output.log, mode="wb") as logfile:
            pipeline = Popen_pipeline([hisat2_cmd, bam_sort_cmd], stdout=outfile, stderr=logfile)
            retcodes = [ p.wait() for p in pipeline ]
            for retcode, cmd in zip(retcodes, cmds):
                if retcode != 0:
                    raise CalledProcessError(rercode, cmd)

rule extract_bam:
    '''Extract SRA file to BAM.'''
    input: 'sra_files/{sra_run}.sra'
    output: 'bam_files/{sra_run}.bam'
    # The call to samtools is theoretically redundant, but needed
    # because it munges the empty quality strings from sam-dump into
    # "*", which picard will accept.
    shell: '''
    sam-dump -u -c {input:q} | \
        samtools view -bS - | \
        picard-tools SortSam I=/dev/stdin O={output:q} \
            SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT
    '''

rule index_bam:
    '''Create .bai file for a bam file.'''
    input: '{basename}.bam'
    output: '{basename}.bam.bai'
    shell: '''
    picard-tools BuildBamIndex I={input:q} O={output:q} \
        VALIDATION_STRINGENCY=LENIENT
    '''

rnaseq_samplemeta = read_R_dataframe("saved_data/samplemeta-RNASeq.RDS")
rnaseq_bam_files = [ "bam_files/{}.bam".format(basename) for basename in rnaseq_samplemeta["SRA_run"] ]
chipseq_samplemeta = read_R_dataframe("saved_data/samplemeta-ChIPSeq.RDS")
chipseq_bam_files = [ "bam_files/{}.bam".format(basename) for basename in chipseq_samplemeta["SRA_run"] ]

rule count_rnaseq_reads:
    input: samplemeta="saved_data/samplemeta-RNASeq.RDS", bam_files=rnaseq_bam_files
    output: "saved_data/SummarizedExperiment-RNASeq.RDS"
    threads: 8
    run:
        cmd = [
            "scripts/rnaseq-count.R",
            "SAMPLEMETA_FILE={}".format(input.samplemeta),
            "BAM_FILES={}".format(",".join(input.bam_files)),
            "SUMEXP_OUTPUT_FILE={}".format(output[0]),
            "THREADS={}".format(threads)
        ]
        check_call(cmd)

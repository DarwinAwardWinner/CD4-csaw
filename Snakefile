import contextlib
import os
import os.path
import rpy2.rinterface
import shutil
import subprocess

from atomicwrites import atomic_write, AtomicWriter
from subprocess import check_call, Popen, PIPE, CalledProcessError, list2cmdline
from rpy2 import robjects
from rpy2.robjects import pandas2ri

from snakemake.io import expand
from snakemake.utils import min_version
min_version("3.7.1")

pandas2ri.activate()
rpy2.rinterface.set_writeconsole_warnerror(lambda x: sys.stderr.write(x))

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

# http://www.gencodegenes.org/releases/; this is all the hg38 releases
# so far
gencode_to_ens_version = {
    20: 76,
    21: 77,
    22: 79,
    23: 81,
    24: 83,
    25: 85
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

# Run a separate Snakemake workflow to fetch the sample metadata,
# which must be avilable before evaluating the rules below. Without
# this two-step workflow, the below rules would involve quite complex
# use of multiple dynamic() inputs and outputs.
from processify import processify
from snakemake import snakemake
snakemake = processify(snakemake)
result = snakemake(
    'pre.Snakefile',
    targets=expand(os.path.join('saved_data', 'samplemeta-{dataset}.RDS'),
                   dataset=('RNASeq', 'ChIPSeq')),
    quiet=True)
if not result:
    raise Exception("Could not retrieve experiment metadata from GEO")

rnaseq_samplemeta = read_R_dataframe("saved_data/samplemeta-RNASeq.RDS")
chipseq_samplemeta = read_R_dataframe("saved_data/samplemeta-ChIPSeq.RDS")

rnaseq_star_outdirs = [
    'rnaseq_star_hg38.analysisSet_knownGene',
    'rnaseq_star_hg38.analysisSet_gencode.v25',
]
rnaseq_hisat_outdir = 'rnaseq_hisat2_grch38_snp_tran'

aligned_rnaseq_star_bam_files = []
for d in rnaseq_star_outdirs:
    for samp in rnaseq_samplemeta["SRA_run"]:
        aligned_rnaseq_star_bam_files.append(
            'aligned/{dirname}/{samp}/Aligned.sortedByCoord.out.bam'.format(
                dirname=d,
                samp=samp))

aligned_rnaseq_hisat_bam_files = []
for samp in rnaseq_samplemeta["SRA_run"]:
    aligned_rnaseq_hisat_bam_files.append(
        'aligned/{dirname}/{samp}/Aligned.bam'.format(
            dirname=rnaseq_hisat_outdir,
            samp=samp))

aligned_rnaseq_bam_files = aligned_rnaseq_star_bam_files + aligned_rnaseq_hisat_bam_files
aligned_rnaseq_bai_files = [ bam + '.bai' for bam in aligned_rnaseq_bam_files ]

subworkflow hg38_ref:
    workdir: os.path.expanduser("~/references/hg38")

include: 'rulegraph.Snakefile'

rule all:
    input:
        rnaseq_counts=[
            'saved_data/rnaseq_star_hg38.analysisSet_gencode.v25_counts.RDS',
            'saved_data/rnaseq_star_hg38.analysisSet_knownGene_counts.RDS',
            'saved_data/rnaseq_hisat2_grch38_snp_tran_counts.RDS',
        ],

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
            wait_for_subprocs(pipeline)

# TODO: Transcriptome-space BAM and gene counts
rule align_rnaseq_with_star_single_end:
    '''Align fastq file with star'''
    input: fastq='fastq_files/{samplename}.fq.gz',
           index_sa=hg38_ref('STAR_index_{genome_build}_{transcriptome}/SA'),
           transcriptome_gff=hg38_ref('{transcriptome}.gff3'),
    output: bam='aligned/rnaseq_star_{genome_build}_{transcriptome}/{samplename}/Aligned.sortedByCoord.out.bam',
            sj='aligned/rnaseq_star_{genome_build}_{transcriptome}/{samplename}/SJ.out.tab',
            tx_bam='aligned/rnaseq_star_{genome_build}_{transcriptome}/{samplename}/Aligned.toTranscriptome.out.bam',
            gene_counts='aligned/rnaseq_star_{genome_build}_{transcriptome}/{samplename}/ReadsPerGene.out.tab',
    params: temp_sam='aligned/rnaseq_star_{genome_build}_{transcriptome}/{samplename}/Aligned.out.sam',
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
            '--sjdbGTFtagExonParentGene', 'gene_id',
            '--sjdbOverhang', '100',
            '--readFilesIn', input.fastq,
            '--readFilesCommand', read_cmd,
            '--outSAMattributes', 'Standard',
            '--outSAMunmapped', 'Within',
            '--outFileNamePrefix', outdir,
            '--outSAMtype', 'SAM',
            '--quantMode', 'TranscriptomeSAM', 'GeneCounts',
        ]
        # Run STAR
        shell(list2cmdline(map(str, star_cmd)))
        # Sort BAM
        picard_cmd = 'picard-tools SortSam I={params.temp_sam:q} O={output.bam:q} SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT'
        shell(picard_cmd)
        # Delete sam file
        os.remove(params.temp_sam)

rule align_rnaseq_with_hisat2_single_end:
    '''Align fastq file with HISAT2'''
    input: fastq='fastq_files/{samplename}.fq.gz',
           index_f1=hg38_ref('HISAT2_index_grch38_snp_tran/index.1.ht2'),
           transcriptome_gff=hg38_ref('knownGene.gff3'),
           chrom_mapping=hg38_ref('chrom_mapping_GRCh38_ensembl2UCSC.txt'),
    output: bam='aligned/rnaseq_hisat2_grch38_snp_tran/{samplename}/Aligned.bam',
            log='aligned/rnaseq_hisat2_grch38_snp_tran/{samplename}/hisat2.log'
    threads: 8
    run:
        index_basename = re.sub('\\.1\\.ht2', "", input.index_f1)
        outdir = os.path.dirname(output.bam)
        ensure_dir(outdir)
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
        with atomic_write(output.bam, mode="wb", overwrite=True) as outfile, \
             open(output.log, mode="wb") as logfile:
            pipeline = Popen_pipeline(cmds, stdout=outfile, stderr=logfile)
            wait_for_subprocs(pipeline)

rule index_bam:
    '''Create .bai file for a bam file.'''
    input: '{basename}.bam'
    output: '{basename}.bam.bai'
    shell: '''
    picard-tools BuildBamIndex I={input:q} O={output:q} \
        VALIDATION_STRINGENCY=LENIENT
    '''

# The hisat2 documentation doesn't specify which version of Ensembl
# they used to build the prebuilt index, so I guess we'll just use the
# latest.
rule count_rnaseq_hiseq2:
    input: samplemeta='saved_data/samplemeta-RNASeq.RDS',
           bam_files=expand('aligned/rnaseq_hisat2_grch38_snp_tran/{SRA_run}/Aligned.bam',
                            SRA_run=rnaseq_samplemeta['SRA_run']),
           bai_files=expand('aligned/rnaseq_hisat2_grch38_snp_tran/{SRA_run}/Aligned.bam.bai',
                            SRA_run=rnaseq_samplemeta['SRA_run']),
           txdb=hg38_ref('TxDb.Hsapiens.Ensembl.hg38.v85.sqlite3'),
           genemeta=hg38_ref('genemeta.ensembl.v85.RDS')
    output: sexp='saved_data/rnaseq_hisat2_grch38_snp_tran_counts.RDS'
    threads: 4
    run:
        cmd = [
            'scripts/rnaseq-count.R',
            '--samplemeta-file', input.samplemeta,
            '--sample-id-column', 'SRA_run',
            '--bam-file-pattern', 'aligned/rnaseq_hisat2_grch38_snp_tran/%s/Aligned.bam',
            '--output-file', output.sexp,
            '--expected-bam-files', ','.join(input.bam_files),
            '--threads', str(threads),
            '--annotation-txdb', input.txdb,
            '--additional-gene-info', input.genemeta,
        ]
        check_call(cmd)

rule count_rnaseq_star_gencode:
    input: samplemeta='saved_data/samplemeta-RNASeq.RDS',
           bam_files=expand(
               'aligned/rnaseq_star_hg38.analysisSet_gencode.v{{gencode_ver}}/{SRA_run}/Aligned.sortedByCoord.out.bam',
               SRA_run=rnaseq_samplemeta['SRA_run']),
           bai_files=expand(
               'aligned/rnaseq_star_hg38.analysisSet_gencode.v{{gencode_ver}}/{SRA_run}/Aligned.sortedByCoord.out.bam.bai',
               SRA_run=rnaseq_samplemeta['SRA_run']),
           txdb=lambda wildcards: hg38_ref(
               'TxDb.Hsapiens.Ensembl.hg38.v{ens_ver}.sqlite3'.format(
                   ens_ver=gencode_to_ens_version[int(wildcards.gencode_ver)])),
           genemeta=lambda wildcards: hg38_ref(
               'genemeta.ensembl.v{ens_ver}.RDS'.format(
               ens_ver=gencode_to_ens_version[int(wildcards.gencode_ver)]))
    output: sexp='saved_data/rnaseq_star_hg38.analysisSet_gencode.v{gencode_ver,\\d+}_counts.RDS'
    threads: 4
    run:
        cmd = [
            'scripts/rnaseq-count.R',
            '--samplemeta-file', input.samplemeta,
            '--sample-id-column', 'SRA_run',
            '--bam-file-pattern',
            'aligned/rnaseq_star_hg38.analysisSet_gencode.v{gencode_ver}/%s/Aligned.sortedByCoord.out.bam'.format(
                   gencode_ver=wildcards.gencode_ver,
            ),
            '--output-file', output.sexp,
            '--expected-bam-files', ','.join(input.bam_files),
            '--threads', str(threads),
            '--annotation-txdb', input.txdb,
            '--additional-gene-info', input.genemeta,
        ]
        print(cmd)
        check_call(cmd)

rule count_rnaseq_star_knownGene:
    input: samplemeta='saved_data/samplemeta-RNASeq.RDS',
           bam_files=expand(
               'aligned/rnaseq_star_hg38.analysisSet_knownGene/{SRA_run}/Aligned.sortedByCoord.out.bam',
               SRA_run=rnaseq_samplemeta['SRA_run']),
           bai_files=expand(
               'aligned/rnaseq_star_hg38.analysisSet_knownGene/{SRA_run}/Aligned.sortedByCoord.out.bam.bai',
               SRA_run=rnaseq_samplemeta['SRA_run']),
           genemeta=hg38_ref('genemeta.org.Hs.eg.db.RDS')
    output: sexp='saved_data/rnaseq_star_hg38.analysisSet_knownGene_counts.RDS'
    threads: 4
    run:
        cmd = [
            'scripts/rnaseq-count.R',
            '--samplemeta-file', input.samplemeta,
            '--sample-id-column', 'SRA_run',
            '--bam-file-pattern', 'aligned/rnaseq_star_hg38.analysisSet_knownGene/%s/Aligned.sortedByCoord.out.bam',
            '--output-file', output.sexp,
            '--expected-bam-files', ','.join(input.bam_files),
            '--threads', str(threads),
            '--annotation-txdb', 'TxDb.Hsapiens.UCSC.hg38.knownGene',
            '--additional-gene-info', input.genemeta,
        ]
        check_call(cmd)

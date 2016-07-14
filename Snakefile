import os.path

from rpy2 import robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()

# Run a separate Snakemake workflow to fetch the sample metadata,
# which must be avilable before evaluating the rules below
from processify import processify
from snakemake import snakemake
snakemake = processify(snakemake)
snakemake('pre.Snakefile', targets=expand(os.path.join('saved_data', 'samplemeta-{dataset}.RDS'), dataset=('RNASeq', 'ChIPSeq')))

def read_R_dataframe(rdsfile):
    '''Read an R data frame stored in an RDS file.

    Returns the result as a Pandas DataFrame.

    '''
    readRDS = robjects.r('readRDS')
    df = readRDS((robjects.StrVector([rdsfile])))
    return(pandas2ri.ri2py(df))

rule fetch_sra_run:
    output: 'sra_files/{sra_run}.sra'
    shell: 'scripts/get-sra-run-files.R {wildcards.sra_run:q}'

rule extract_bam:
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
    input: '{basename}.bam'
    output: '{basename}.bam.bai'
    shell: '''
    picard-tools BuildBamIndex I={input:q} O={output:q} \
        VALIDATION_STRINGENCY=LENIENT'''

rnaseq_samplemeta = read_R_dataframe("saved_data/samplemeta-RNASeq.RDS")
rnaseq_bam_files = [ "bam_files/{}.bam".format(basename) for basename in rnaseq_samplemeta["SRA_run"] ]
chipseq_samplemeta = read_R_dataframe("saved_data/samplemeta-ChIPSeq.RDS")
chipseq_bam_files = [ "bam_files/{}.bam".format(basename) for basename in chipseq_samplemeta["SRA_run"] ]

rule count_rnaseq_reads:
    input: samplemeta="saved_data/samplemeta-RNASeq.RDS", bam_files=rnaseq_bam_files
    output: "saved_data/SummarizedExperiment-RNASeq.RDS"

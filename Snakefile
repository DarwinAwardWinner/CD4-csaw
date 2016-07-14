import os.path

# Run a separate Snakemake workflow to fetch the sample metadata,
# which must be avilable before evaluating the rules below
from processify import processify
from snakemake import snakemake
snakemake = processify(snakemake)
snakemake('pre.Snakefile', targets=expand(os.path.join('saved_data', 'samplemeta-{dataset}.RDS'), dataset=('RNASeq', 'ChIPSeq')))

rule fetch_sra_run:
    output: 'sra_files/{sra_run}.sra'
    shell: 'scripts/get-sra-run-files.R {wildcards.sra_run:q}'

rule extract_bam:
    input: 'sra_files/{sra_run}.sra'
    output: 'bam_files/{sra_run}.bam'
    shell: '''
    sam-dump -u -c {input:q} | \
        picard-tools SortSam I=/dev/stdin O={output:q} \
            SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT
    '''

rule index_bam:
    input: '{basename}.bam'
    output: '{basename}.bam.bai'
    shell: '''
    picard-tools BuildBamIndex I={input}, O={output} \
        VALIDATION_STRINGENCY=LENIENT'''

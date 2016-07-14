import os.path

# Run a separate Snakemake workflow to
from processify import processify
from snakemake import snakemake
snakemake = processify(snakemake)
snakemake("pre.Snakefile", targets=expand(os.path.join("saved_data", "samplemeta-{dataset}.RDS"), dataset=("RNASeq", "ChIPSeq")))

rule fetch_sra_run:
    output: "sra_files/{sra_run}.sra"
    shell: "scripts/get-sra-run-files.R {wildcards.sra_run:q}"

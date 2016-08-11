import os.path

rule get_sra_metadb:
    output: os.path.join("saved_data", "SRAmetadb.sqlite")
    shell: 'scripts/get-sra-metadb.R'

rule get_geo_meta:
    input: os.path.join("saved_data", "SRAmetadb.sqlite")
    output: expand(os.path.join("saved_data", "samplemeta-{dataset}.RDS"), dataset=("RNASeq", "ChIPSeq"))
    shell: "scripts/get-geo-metadata.R"

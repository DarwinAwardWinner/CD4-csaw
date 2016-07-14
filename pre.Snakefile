import os.path

rule get_geo_meta:
    output: expand(os.path.join("saved_data", "samplemeta-{dataset}.RDS"), dataset=("RNASeq", "ChIPSeq"))
    shell: "scripts/get-geo-metadata.R"

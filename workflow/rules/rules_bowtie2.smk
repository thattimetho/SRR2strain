## ----- Bowtie2 mapping and indexing rules ------
from snakemake.io import expand

rule bowtie2_index:
## bowtie2_index                                : Index fasta file for mapping later
    input:
        genome_db_in=expand("data/{genome_fasta_dir}/{genome_ID}.fa",
                            genome_fasta_dir=config["data_dir_settings"]["genome_dir"],
                            genome_ID=config["metadata_settings"]["genome_ID"])
    output:
        genome_index_out=expand("data/{genome_index_dir}/{genome_ID}_index.1.bt2",
                                genome_index_dir=config["data_dir_settings"]["index_dir"])
    params:
        index_name=config["metadata_settings"]["genome_ID"]
    conda:
        "../envs/manual-tools.yml"
    shell:
        "bowtie2-build {input.genome_db_in} {params.index_name}_index"

rule bowtie2_map:
## bowtie2_map                                  : Map reads unto reference genome database
    input:
        
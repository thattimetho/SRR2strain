## ----- Samtools helper rules ------
from snakemake.io import expand

rule samtools_convert_to_bam:
## samtools_convert_to_bam                      : Index fasta file for mapping later
    input:
        genome_db_in=expand("data/{genome_fasta_dir}/{genome_ID}.fa",
                            genome_fasta_dir=config["data_dir_settings"]["genome_dir"],
                            genome_ID=config["metadata_settings"]["genome_ID"])
    output:
        genome_index_out=expand("data/{genome_index_dir}/{genome_ID}_index.1.bt2",
                                genome_index_dir=config["data_dir_settings"]["index_dir"],
                                genome_ID=config["metadata_settings"]["genome_ID"])
    params:
        index_name=config["metadata_settings"]["genome_ID"]
    conda:
        "../envs/manual-samtools.yml"
    shell:
        "samtools view -b {input} > {output}"


rule samtools_sort:
## samtools_sort                                : Sort bam file
    input:
        unsorted_bam_in="data/{run_ID}/mappings/{mapping_ID}.bam"
    output:
        sorted_bam_out="data/{run_ID}/mappings/{mapping_ID}.sorted.bam"
    conda:
        "../envs/manual-samtools.yml"
    shell:
        "samtools sort -o {output} {input}"
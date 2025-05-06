## ----- Bowtie2 mapping and indexing rules ------
from snakemake.io import expand, temp

rule bowtie2_index_db:
## bowtie2_index_db                             : Index fasta file for mapping later
    input:
        genome_db_in = expand("data/{genome_fasta_dir}/{genome_ID}.fa",
                              genome_fasta_dir=config["data_dir_settings"]["genome_dir"],
                              genome_ID=config["metadata_settings"]["genome_ID"])
    output:
        genome_index_out = expand("data/{genome_index_dir}/{genome_ID}_index.1.bt2",
                                  genome_index_dir=config["data_dir_settings"]["index_dir"],
                                  genome_ID=config["metadata_settings"]["genome_ID"])
    params:
        index_name = config["metadata_settings"]["genome_ID"]
    log:
        stdout = expand("data/{run_ID}/logs/bowtie2_build_{genome_ID}.log",
                        run_ID = config["metadata_settings"]["dataset_ID"],
                        genome_ID=config["metadata_settings"]["genome_ID"]),
        stderr = expand("data/{run_ID}/logs/bowtie2_build_{genome_ID}.err.log",
                        run_ID = config["metadata_settings"]["dataset_ID"],
                        genome_ID=config["metadata_settings"]["genome_ID"])
    threads:
        max(workflow.cores, 4)
    conda:
        "manual-tools"
    shell:
        "bowtie2-build --threads {threads} {input.genome_db_in} {params.index_name}_index > {log.stdout} 2> {log.stderr}"

rule bowtie2_map:
## bowtie2_map                                  : Map reads unto reference genome database
    input:
        genome_index_in = expand("data/{genome_index_dir}/{genome_ID}_index.1.bt2",
                                 genome_index_dir=config["data_dir_settings"]["index_dir"],
                                 genome_ID=config["metadata_settings"]["genome_ID"]),
        sra_dataset_reads_in_1 = "data/{run_ID}/raw_reads/{SRR_ID}_1.fastq",
        sra_dataset_reads_in_2 = "data/{run_ID}/raw_reads/{SRR_ID}_2.fastq"
    output:
        mapped_reads_sam_unsorted_out = temp("data/{run_ID}/mappings/{SRR_ID}_{genome_ID}.sam"),
    params:
        genome_index_in = expand("data/{genome_index_dir}/{genome_ID}_index",
                                 genome_index_dir=config["data_dir_settings"]["index_dir"],
                                 genome_ID=config["metadata_settings"]["genome_ID"])
    log:
        stdout = "data/{run_ID}/logs/bowtie2_{SRR_ID}_{genome_ID}.log",
        stderr = "data/{run_ID}/logs/bowtie2_{SRR_ID}_{genome_ID}.err.log"
    threads:
        max(workflow.cores/2, 4)
    conda:
        "manual-tools"
    shell:
        """
        bowtie2 --threads {threads} --very-sensitive-local --no-unal -x {params.genome_index_in} \
        -1 {input.sra_dataset_reads_in_1} -2 {input.sra_dataset_reads_in_2} -S {output.mapped_reads_sam_unsorted_out} > \
        {log.stdout} 2> {log.stderr}
        """

##
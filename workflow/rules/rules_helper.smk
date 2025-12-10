## ----- Rules related to data or additional analyses ------

from snakemake.io import expand

rule checkv_end_to_end:
## checkv_end_to_end                            : Run CheckV genome analysis (end_to_end pipeline)
    input:
        genome_db_in = expand("data/{genome_fasta_dir}/{genome_ID}.fa",
                              genome_fasta_dir = config["data_dir_settings"]["genome_dir"],
                              genome_ID = config["metadata_settings"]["genome_ID"])
    output:
        genome_quality_out = expand("data/{genome_index_dir}/checkv_{genome_ID}/quality_summary.tsv",
                                    genome_index_dir = config["data_dir_settings"]["index_dir"],
                                    genome_ID = config["metadata_settings"]["genome_ID"])
    params:
        index_name = config["metadata_settings"]["genome_ID"],
        index_dir = config["data_dir_settings"]["index_dir"],
        checkv_db = config["database_dir_settings"]["checkv_db_dir"]
    log:
        stdout = expand("data/{run_ID}/logs/checkv-end_to_end-{genome_ID}.log",
                        run_ID = config["metadata_settings"]["dataset_ID"],
                        genome_ID = config["metadata_settings"]["genome_ID"]),
        stderr = expand("data/{run_ID}/logs/checkv-end_to_end-{genome_ID}.err.log",
                        run_ID = config["metadata_settings"]["dataset_ID"],
                        genome_ID = config["metadata_settings"]["genome_ID"])
    threads:
        max(workflow.cores, 10)
    conda:
        "manual-checkv"
    shell:
        "checkv end_to_end {input.genome_db_in} checkv_{params.index_name} -d {params.checkv_db} -t {threads} > {log.stdout} 2> {log.stderr}"

rule coverm_coverage:
## coverm_coverage                              : Run CoverM to calculate coverage per contig, per SRR_ID
    input:
        genome_db_in = expand("data/{genome_fasta_dir}/{genome_ID}.fa",
                              genome_fasta_dir = config["data_dir_settings"]["genome_dir"],
                              genome_ID = config["metadata_settings"]["genome_ID"]),
        mapped_reads_bam_sorted_in = "data/{run_ID}/mappings/{SRR_ID}-{genome_ID}.sorted.bam"
    output:
        genome_coverage_out = "data/{run_ID}/coverm-{SRR_ID}-{genome_ID}.tsv"
    log:
        stdout = "data/{run_ID}/logs/coverm-{SRR_ID}-{genome_ID}.log",
        stderr = "data/{run_ID}/logs/coverm-{SRR_ID}-{genome_ID}.err.log"
    threads:
        1
    conda:
        "manual-coverm"
    shell:
        "coverm contig -b {input.mapped_reads_bam_sorted_in} -r {input.genome_db_in} -o {output} > {log.stdout} 2> {log.stderr}"


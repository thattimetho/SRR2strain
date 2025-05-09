## ----- Rules related to gene prediction (Prodigal) ------
from snakemake.io import expand

rule prodigal_run:
## prodigal_run                                 : Run prodigal protein prediction tool
    input:
        genome_db_in = expand("data/{genome_fasta_dir}/{genome_ID}.fa",
                              genome_fasta_dir=config["data_dir_settings"]["genome_dir"],
                              genome_ID=config["metadata_settings"]["genome_ID"])
    output:
        genes_db_out = expand("data/{genes_dir}/{genome_ID}.genes",
                               genes_dir=config["data_dir_settings"]["genes_dir"],
                               genome_ID=config["metadata_settings"]["genome_ID"]),
        genes_fna_out = expand("data/{genes_dir}/{genome_ID}.genes.fna",
                               genes_dir=config["data_dir_settings"]["genes_dir"],
                               genome_ID=config["metadata_settings"]["genome_ID"]),
        genes_faa_out = expand("data/{genes_dir}/{genome_ID}.genes.faa",
                               genes_dir=config["data_dir_settings"]["genes_dir"],
                               genome_ID=config["metadata_settings"]["genome_ID"])
    log:
        stdout=expand("data/{run_ID}/logs/prodigal_{genome_ID}.log",
                      run_ID = config["metadata_settings"]["dataset_ID"],
                      genome_ID=config["metadata_settings"]["genome_ID"]),
        stderr=expand("data/{run_ID}/logs/prodigal_{genome_ID}.err.log",
                      run_ID = config["metadata_settings"]["dataset_ID"],
                      genome_ID=config["metadata_settings"]["genome_ID"])
    threads:
        4
    conda:
        "manual-tools"
    shell:
        """
        prodigal -i {input.genome_db_in} -o {output.genes_db_out} -a {output.genes_faa_out} \
        -d {output.genes_fna_out} -p meta > {log.stdout} 2> {log.stderr}
        """
##
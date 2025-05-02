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
    conda:
        "../envs/manual-tools.yml"
    shell:
        """
        prodigal -i {input.genome_db_in} -o {output.genes_db_out} -a {output.genes_faa_out} 
        -d {output.genes_fna_out} -p anon 
        """
##
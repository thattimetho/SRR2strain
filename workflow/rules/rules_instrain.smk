## ----- Rules related to microdiversity analysis (inStrain) ------
from snakemake.io import expand

rule instrain_run:
## instrain_run                                 : Run inStrain to calculate sample diversity statistics
    input:
        srr_mapped_in = "data/{run_ID}/mappings/{SRR_ID}_{genome_ID}.sorted.bam",
        genome_db_in = expand("data/{genome_fasta_dir}/{genome_ID}.fa",
                              genome_fasta_dir=config["data_dir_settings"]["genome_dir"],
                              genome_ID=config["metadata_settings"]["genome_ID"]),
        genes_fna_in = expand("data/{genes_dir}/{genome_ID}.genes.fna",
                               genes_dir=config["data_dir_settings"]["genes_dir"],
                               genome_ID=config["metadata_settings"]["genome_ID"])
    output:
        instrain_scaffold_out = "data/{run_ID}/instrain/instrain_{run_ID}_{SRR_ID}_{genome_ID}/output/instrain_{run_ID}_{SRR_ID}_scaffold_info.tsv"
    params:
        instrain_dir_out = "data/{run_ID}/instrain/instrain_{run_ID}_{SRR_ID}_{genome_ID}"
    log:
        stdout="data/{run_ID}/logs/instrain_{SRR_ID}_{genome_ID}.log",
        stderr="data/{run_ID}/logs/instrain_{SRR_ID}_{genome_ID}.err.log"
    threads:
        max(workflow.cores/2, 4)
    conda:
        "manual-instrain"
    shell:
        """
        inStrain profile --skip_plot_generation --skip_genome_wide -p {threads} {input.srr_mapped_in} {input.genome_db_in} 
        -o {params.instrain_dir_out} -g {input.genes_fna_in} > {log.stdout} 2> {log.stderr}
        """

rule phabox_lifestyle_inference:
## phabox_lifestyle_inference                   : Run phage lifestyle prediction analysis
    input:
        genome_db_in = expand("data/{genome_fasta_dir}/{genome_ID}.fa",
                              genome_fasta_dir=config["data_dir_settings"]["genome_dir"],
                              genome_ID=config["metadata_settings"]["genome_ID"])
    params:
        phabox_dir_out = config["data_dir_settings"]["genome_dir"]
    log:
        stdout="data/{run_ID}/logs/phabox2_{genome_ID}.log",
        stderr="data/{run_ID}/logs/phabox2_{genome_ID}.err.log"
    threads:
        max(workflow.cores/2, 4)
    conda:
        "manual-phabox"
    shell:
        """
        phabox2 --task phatyp --threads {threads} --dbdir --outpth {params.phabox_dir_out} 
        --contigs {input.genome_db_in} > {log.stdout} 2> {log.stderr}
        """
##
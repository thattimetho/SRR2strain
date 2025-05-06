## ----- SRA-related helper rules ------
from snakemake.io import temp, expand

rule sra_prefetch:
## sra_prefetch                                 : Prefetches SRA datasets from NCBI SRA repo
    output:
        sra_dataset_out = "data/{run_ID}/sra_temp/{SRR_ID}.sra"
    params:
        sra_ID_prefetch = lambda wildcards: wildcards.SRR_ID,
        sra_dataset_out_dir="data/{run_ID}/sra_temp/"
    log:
        stdout="data/{run_ID}/logs/sra_prefetch_{SRR_ID}.log",
        stderr="data/{run_ID}/logs/sra_prefetch_{SRR_ID}.err.log"
    threads:
        1
    shell:
        """
        prefetch {params.sra_ID_prefetch} -O {params.sra_dataset_out_dir} -f yes > {log.stdout} 2> {log.stderr}
        """


rule sra_fasterq_dump:
## sra_fasterq_dump                             : Converts (or dumps) SRA datasets into .fastq files
    input:
        sra_dataset_in = "data/{run_ID}/sra_temp/{SRR_ID}.sra"
    output:
        sra_dataset_reads_out_1 = temp("data/{run_ID}/raw_reads/{SRR_ID}_1.fastq"),
        sra_dataset_reads_out_2 = temp("data/{run_ID}/raw_reads/{SRR_ID}_2.fastq")
    params:
        sra_dataset_reads_dir = "data/{run_ID}/raw_reads/"
    log:
        stdout="data/{run_ID}/logs/sra_fasterq_dump_{SRR_ID}.log",
        stderr="data/{run_ID}/logs/sra_fasterq_dump_{SRR_ID}.err.log"
    threads:
        1
    shell:
        """
        fasterq-dump --threads {threads} --out-dir {params.sra_dataset_reads_dir} {input.sra_dataset_in} > {log.stdout}
        2> {log.stderr}
        """

##
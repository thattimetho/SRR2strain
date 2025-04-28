## ----- SRA-related helper rules ------

rule sra_prefetch:
## sra_prefetch                                 : Prefetches SRA datasets from NCBI SRA repo
    input:
        sra_prefetch_list_in="data/{run_ID}/metadata/SRR_Acc_List_{bioproject_ID}.txt"
    output:
        sra_dataset_out="data/{run_ID}/sra_temp_{bioproject_ID}/{SRR_ID}.sra"
    params:
        sra_dataset_out_dir="data/{run_ID}/sra_temp_{bioproject_ID}/"
    threads:
        max(int(workflow.cores/2), 4)
    shell:
        "parallel --verbose -j {threads} -a {input.sra_prefetch_list_in} prefetch -O {params.sra_dataset_out_dir}"


rule sra_fasterq_dump:
## sra_fasterq_dump                             : Converts (or dumps) SRA datasets into .fastq files
    input:
        sra_dataset_in = "data/{run_ID}/sra_temp/{SRR_ID}.sra"
    output:
        sra_dataset_reads_out_1 = "data/{run_ID}/raw_reads/{SRR_ID}_1.fastq",
        sra_dataset_reads_out_2 = "data/{run_ID}/raw_reads/{SRR_ID}_2.fastq"
    params:
        sra_dataset_reads_dir = "data/{run_ID}/raw_reads/"
    threads:
        max(int(workflow.cores/10), 1)
    shell:
        "fasterq-dump --threads {threads} --out-dir {params.sra_dataset_reads_dir} {input.sra_dataset_in}"

##
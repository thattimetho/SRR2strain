## ----- Samtools helper rules ------
from snakemake.io import temp

rule samtools_convert_to_bam:
## samtools_convert_to_bam                      : Index fasta file for mapping later
    input:
        mapped_reads_sam_unsorted_in = "data/{run_ID}/mappings/{SRR_ID}_{genome_ID}.sam"
    output:
        mapped_reads_bam_unsorted_out = temp("data/{run_ID}/mappings/{SRR_ID}_{genome_ID}.bam")
    conda:
        "../envs/manual-samtools.yml"
    shell:
        "samtools view -b {input} > {output}"


rule samtools_sort:
## samtools_sort                                : Sort bam file
    input:
        mapped_reads_bam_unsorted_in = "data/{run_ID}/mappings/{SRR_ID}_{genome_ID}.bam"
    output:
        mapped_reads_bam_sorted_out = "data/{run_ID}/mappings/{SRR_ID}_{genome_ID}.sorted.bam"
    conda:
        "../envs/manual-samtools.yml"
    shell:
        "samtools sort -o {output} {input}"

##
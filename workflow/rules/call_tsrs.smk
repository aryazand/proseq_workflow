rule call_tsrs:
    input:
        "results/deeptools/5prime_coverage/{sample}_{strand}.bw",
    output:
        "results/tsrs/{sample}_{strand}.bed",
    conda:
        "../envs/tsr_detectr.yml"
    log:
        "results/tsrs/{sample}_{strand}.log",
    script:
        "../scripts/call_tsrs.R"


rule consolidate_tsrs:
    input:
        tsrs=lambda wildcards: expand(
            "results/tsrs/{sample}_{strand}.bed",
            sample=fastq_process_align.samples.index,
            strand=["forward", "reverse"],
        ),
    output:
        "results/tsrs/TSRs.bed",
    params:
        extra="-s -d 10 -c 4,5,6 -o distinct,count,distinct",
    log:
        "results/tsrs/TSRs.log",
    wrapper:
        "v8.1.1/bio/bedtools/merge"

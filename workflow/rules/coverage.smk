rule deeptools_coverage:
    input:
        bam=fastq_process_align.get_bam_2,
        bai=fastq_process_align.get_bai,
    output:
        "results/deeptools/coverage/{sample}_{direction}.bw",
    threads: 4
    params:
        effective_genome_size=config["mapping_stats"]["deeptools_coverage"][
            "genome_size"
        ],
        extra=config["mapping_stats"]["deeptools_coverage"]["extra"]
        + " --filterRNAstrand {direction}",
    log:
        "results/deeptools/coverage/{sample}_{direction}.log",
    message:
        "generate normalized coverage using deeptools"
    wrapper:
        "v7.0.0/bio/deeptools/bamcoverage"

rule deeptools_5prime_coverage:
    input:
        bam=fastq_process_align.get_bam_2,
        bai=fastq_process_align.get_bai,
    output:
        "results/deeptools/5prime_coverage/{sample}_{direction}.bw",
    threads: 4
    params:
        effective_genome_size=config["mapping_stats"]["deeptools_coverage"][
            "genome_size"
        ],
        extra=config["mapping_stats"]["deeptools_coverage"]["extra"]
        + " --Offset 1 --filterRNAstrand {direction}",
    log:
        "results/deeptools/5prime_coverage/{sample}_{direction}.log",
    message:
        "generate normalized 5prime coverage using deeptools"
    wrapper:
        "v7.0.0/bio/deeptools/bamcoverage"

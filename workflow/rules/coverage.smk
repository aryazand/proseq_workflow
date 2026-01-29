rule deeptools_coverage:
    input:
        bam=fastq_process_align.get_cram,
        bai=fastq_process_align.get_crai,
    output:
        "results/deeptools/coverage/{sample}_{direction}.bw",
    wildcard_constraints:
        direction="for|forward|plus|rev|reverse|minus"
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
        bam=fastq_process_align.get_cram,
        bai=fastq_process_align.get_crai,
    output:
        "results/deeptools/5prime_coverage/{sample}_{direction}.bw",
    wildcard_constraints:
        direction="for|forward|plus|rev|reverse|minus"
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


rule deeptools_3prime_coverage:
    input:
        bam=fastq_process_align.get_cram,
        bai=fastq_process_align.get_crai,
    output:
        "results/deeptools/3prime_coverage/{sample}_{direction}.bw",
    wildcard_constraints:
        direction="for|forward|plus|rev|reverse|minus"
    threads: 4
    params:
        effective_genome_size=config["mapping_stats"]["deeptools_coverage"][
            "genome_size"
        ],
        extra=config["mapping_stats"]["deeptools_coverage"]["extra"]
        + " --Offset -1 --filterRNAstrand {direction}",
    log:
        "results/deeptools/3prime_coverage/{sample}_{direction}.log",
    message:
        "generate normalized 3prime coverage using deeptools"
    wrapper:
        "v7.0.0/bio/deeptools/bamcoverage"

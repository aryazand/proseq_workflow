rule deeptools_coverage_forward:
    input:
        bam=fastq_process_align.get_bam_2,
        bai=fastq_process_align.get_bai,
    output:
        "results/deeptools/coverage/{sample}_forward.bw",
    threads: 4
    params:
        effective_genome_size=config["mapping_stats"]["deeptools_coverage"][
            "genome_size"
        ],
        extra=config["mapping_stats"]["deeptools_coverage"]["extra"]
        + " --filterRNAstrand forward",
    log:
        "results/deeptools/coverage/{sample}_forward.log",
    message:
        "generate normalized forward coverage using deeptools"
    wrapper:
        "v7.0.0/bio/deeptools/bamcoverage"


rule deeptools_coverage_reverse:
    input:
        bam=fastq_process_align.get_bam_2,
        bai=fastq_process_align.get_bai,
    output:
        "results/deeptools/coverage/{sample}_reverse.bw",
    threads: 4
    params:
        effective_genome_size=config["mapping_stats"]["deeptools_coverage"][
            "genome_size"
        ],
        extra=config["mapping_stats"]["deeptools_coverage"]["extra"]
        + " --filterRNAstrand reverse",
    log:
        "results/deeptools/coverage/{sample}_reverse.log",
    message:
        "generate normalized reverse coverage using deeptools"
    wrapper:
        "v7.0.0/bio/deeptools/bamcoverage"


rule deeptools_5prime_coverage_forward:
    input:
        bam=fastq_process_align.get_bam_2,
        bai=fastq_process_align.get_bai,
    output:
        "results/deeptools/5prime_coverage/{sample}_forward.bw",
    threads: 4
    params:
        effective_genome_size=config["mapping_stats"]["deeptools_coverage"][
            "genome_size"
        ],
        extra=config["mapping_stats"]["deeptools_coverage"]["extra"]
        + "--Offset 1 --filterRNAstrand forward",
    log:
        "results/deeptools/5prime_coverage/{sample}_forward.log",
    message:
        "generate normalized forward 5prime coverage using deeptools"
    wrapper:
        "v7.0.0/bio/deeptools/bamcoverage"


rule deeptools_5prime_coverage_reverse:
    input:
        bam=fastq_process_align.get_bam_2,
        bai=fastq_process_align.get_bai,
    output:
        "results/deeptools/5prime_coverage/{sample}_reverse.bw",
    threads: 4
    params:
        effective_genome_size=config["mapping_stats"]["deeptools_coverage"][
            "genome_size"
        ],
        extra=config["mapping_stats"]["deeptools_coverage"]["extra"]
        + "--Offset 1 --filterRNAstrand reverse",
    log:
        "results/deeptools/5prime_coverage/{sample}_reverse.log",
    message:
        "generate normalized reverse 5prime coverage using deeptools"
    wrapper:
        "v7.0.0/bio/deeptools/bamcoverage"

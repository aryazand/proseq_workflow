rule get_chr_sizes:
    input:
        "results/get_genome/genome.fasta.fai",
    output:
        "results/get_genome/genome.chrom.sizes",
    log:
        "results/get_genome/chrom_sizes.log",
    shell:
        """
        cut -f1,2 {input} > {output}
        """


rule bed_to_bigbed:
    input:
        bed="results/tsrs/TSRs.bed",
        chrom_sizes="results/get_genome/genome.chrom.sizes",
    output:
        bigbed="results/tsrs/TSRs.bb",
    conda:
        "../envs/ucsc_tools.yaml"
    log:
        "results/tsrs/bigbed.log",
    shell:
        """
        bedToBigBed {input.bed} {input.chrom_sizes} {output.bigbed}
        """


rule ucsc_trackhub:
    input:
        tsrs="results/tsrs/TSRs.bb",
        bw=lambda wildcards: expand(
            "results/deeptools/5prime_coverage/{sample}_{strand}.bw",
            sample=fastq_process_align.samples.index,
            strand=["forward", "reverse"],
        ),
    output:
        dir=directory(trackhub_dir),
        trackdb=os.path.join(
            trackhub_dir, config["ucsc_trackhub"]["genome"], "trackDb.txt"
        ),
        hubtxt=os.path.join(
            trackhub_dir, config["ucsc_trackhub"]["hub_name"] + ".hub.txt"
        ),
        genomestxt=os.path.join(
            trackhub_dir, config["ucsc_trackhub"]["hub_name"] + ".genomes.txt"
        ),
    conda:
        "../envs/trackhub.yaml"
    log:
        os.path.join(trackhub_dir, "trackhub.log"),
    params:
        hub_name=config["ucsc_trackhub"]["hub_name"],
        defaultPos=config["ucsc_trackhub"]["defaultPos"],
        short_label=config["ucsc_trackhub"]["short_label"],
        long_label=config["ucsc_trackhub"]["long_label"],
        genome=config["ucsc_trackhub"]["genome"],
        email=config["ucsc_trackhub"]["email"],
    script:
        "../scripts/trackhub.py"

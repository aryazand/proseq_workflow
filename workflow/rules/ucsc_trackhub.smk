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


rule bed_to_bigBed:
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


rule gff3_to_GenePred:
    input:
        gff="results/get_genome/genome.gff",
    output:
        genePred="results/get_genome/genome.GenePred",
    conda:
        "../envs/ucsc_tools.yaml"
    log:
        "results/get_genome/gff3_to_GenePred.log",
    shell:
        """
        gff3ToGenePred {input.gff} {output.genePred}
        """


rule GenePred_to_bigGenePred:
    input:
        GenePred="results/get_genome/genome.GenePred",
    output:
        bigGenePred="results/get_genome/genome.bb",
    conda:
        "../envs/ucsc_tools.yaml"
    log:
        "results/get_genome/GenePred_to_bigGenePred.log",
    shell:
        """
        genePredToBigGenePred {input.GenePred} {output.bigGenePred}
        """


rule faToTwoBit:
    input:
        "results/get_genome/genome.fasta",
    output:
        "results/get_genome/genome.2bit",
    log:
        "results/get_genome/fa_to_2bit.log",
    wrapper:
        "v7.1.0/bio/ucsc/faToTwoBit"


rule ucsc_trackhub:
    input:
        genome_2bit="results/get_genome/genome.2bit",
        genome_genePred="results/get_genome/genome.bb",
        tsrs="results/tsrs/TSRs.bb",
        plus_bw=lambda wildcards: expand(
            "results/deeptools/5prime_coverage/{sample}_forward.bw",
            sample=fastq_process_align.samples.index,
        ),
        minus_bw=lambda wildcards: expand(
            "results/deeptools/5prime_coverage/{sample}_reverse.bw",
            sample=fastq_process_align.samples.index,
        ),
    output:
        dir=directory(trackhub_dir),
        trackdb=expand(
            os.path.join(trackhub_dir, "{org}", "trackDb.txt"),
            org=config["ucsc_trackhub"]["genomes"],
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
    script:
        "../scripts/trackhub.py"

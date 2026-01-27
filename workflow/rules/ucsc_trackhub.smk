rule get_chr_sizes:
    input:
        fastq_process_align.get_fasta_index,
    output:
        "results/get_genome/genome.chrom.sizes",
    log:
        "results/get_genome/chrom_sizes.log",
    shell:
        """
        cut -f1,2 {input} > {output}
        """


rule gff3ToGenePred:
    input:
        gff="{file}.gff",
    output:
        genePred="{file}.genePred",
    params:
        extra=lambda wc: config["ucsc_trackhub"]["gff3ToGenePred"].get(
            get_config_key(wc.file), ""
        ),
    container:
        "docker://quay.io/biocontainers/ucsc-gff3togenepred:482--h0b57e2e_0"
    log:
        "{file}.gff3ToGenePred.log",
    shell:
        """
        gff3ToGenePred {params.extra} {input.gff} {output.genePred}            
        """


rule gtfToGenePred:
    input:
        gtf="{file}.gtf",
    output:
        genePred="{file}.genePred",
    params:
        extra=lambda wc: config["ucsc_trackhub"]["gtfToGenePred"].get(
            get_config_key(wc.file), ""
        ),
    container:
        "docker://quay.io/biocontainers/ucsc-gtftogenepred:482--h0b57e2e_0"
    log:
        "{file}.gtfToGenePred.log",
    shell:
        """
        gtfToGenePred {params.extra} {input.gtf} {output.genePred}            
        """


rule genePredToBigGenePred:
    input:
        genePred="{file}.genePred",
    output:
        bgInput="{file}.bgInput",
    params:
        extra=lambda wc: config["ucsc_trackhub"]["genePredToBigGenePred"].get(
            get_config_key(wc.file), ""
        ),
    container:
        "docker://quay.io/biocontainers/ucsc-genepredtobiggenepred:482--h0b57e2e_0"
    log:
        "{file}.genePredToBigGenePred.log",
    shell:
        """
        genePredToBigGenePred {params.extra} {input.genePred}  stdout | sort -k1,1 -k2,2n > {output.bgInput}            
        """


rule bedToBigBed:
    input:
        bed=branch(
            lambda wc: "tsrs" in wc.file, then="{file}.bed", otherwise="{file}.bgInput"
        ),
        chrom_sizes="results/get_genome/genome.chrom.sizes",
    output:
        bigBed="{file}.bb",
    params:
        extra=lambda wc: config["ucsc_trackhub"]["bedToBigBed"].get(
            get_config_key(wc.file), ""
        ),
    container:
        "docker://quay.io/biocontainers/ucsc-bedtobigbed:482--hdc0a859_0"
    log:
        "{file}.bedToBigBed.log",
    shell:
        """
        bedToBigBed {params.extra} {input.bed} {input.chrom_sizes} {output.bigBed}          
        """


rule faToTwoBit:
    input:
        "results/get_genome/genome.fasta",
    output:
        "results/get_genome/genome.2bit",
    params:
        extra=config["ucsc_trackhub"]["faToTwoBit"]["genome"],
    log:
        "results/get_genome/fa_to_2bit.log",
    wrapper:
        "v7.1.0/bio/ucsc/faToTwoBit"


rule ucsc_trackhub:
    input:
        genome_2bit="results/get_genome/genome.2bit",
        genome_genePred="results/get_genome/genome.bb",
        tsrs=["results/tsrs/TSRs.bb", "results/tsrs/TSRs_merged.bb"],
        stringtie=lambda wildcards: expand(
            "results/stringtie/{sample}.bb",
            sample=fastq_process_align.samples.index,
        ),
        fiveprime_plus_bw=lambda wildcards: expand(
            "results/deeptools/5prime_coverage/{sample}_forward.bw",
            sample=fastq_process_align.samples.index,
        ),
        fiveprime_minus_bw=lambda wildcards: expand(
            "results/deeptools/5prime_coverage/{sample}_reverse.bw",
            sample=fastq_process_align.samples.index,
        ),
        plus_bw=lambda wildcards: expand(
            "results/deeptools/coverage/{sample}_forward.bw",
            sample=fastq_process_align.samples.index,
        ),
        minus_bw=lambda wildcards: expand(
            "results/deeptools/coverage/{sample}_reverse.bw",
            sample=fastq_process_align.samples.index,
        ),
    output:
        dir=directory(trackhub_dir),
        trackdb=expand(
            os.path.join(trackhub_dir, "{org}", "trackDb.txt"),
            org=config["ucsc_trackhub"]["genomes"],
        ),
        hubtxt=os.path.join(
            trackhub_dir, config["ucsc_trackhub"]["hub_file"]["hub_name"] + ".hub.txt"
        ),
        genomestxt=os.path.join(
            trackhub_dir,
            config["ucsc_trackhub"]["hub_file"]["hub_name"] + ".genomes.txt",
        ),
    container:
        "docker://quay.io/biocontainers/trackhub:1.0--pyh7cba7a3_0"
    log:
        os.path.join(trackhub_dir, "trackhub.log"),
    script:
        "../scripts/trackhub.py"

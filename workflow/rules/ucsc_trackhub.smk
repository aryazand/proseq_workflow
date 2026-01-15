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


rule tsr_bed_to_bigBed:
    input:
        bed="results/tsrs/TSRs.bed",
        chrom_sizes="results/get_genome/genome.chrom.sizes",
    output:
        bigbed="results/tsrs/TSRs.bb",
    params:
        extra=config["ucsc_trackhub"]["process_tsrs"]["tsr_bed_to_bigBed"],
    conda:
        "../envs/ucsc_tools.yaml"
    log:
        "results/tsrs/bigbed.log",
    shell:
        """
        bedToBigBed {params.extra} {input.bed} {input.chrom_sizes} {output.bigbed}
        """


rule ref_gff_to_bigGenePred:
    input:
        gff="results/get_genome/genome.gff",
        chrom_sizes="results/get_genome/genome.chrom.sizes",
    output:
        bigGenePred="results/get_genome/genome.bb",
    params:
        intermediate_files=lambda wildcards, input: multiext(
            subpath(input.gff, strip_suffix=".gff"), ".genePred", ".bgpInput"
        ),
        gff3ToGenePred_extra=config["ucsc_trackhub"]["process_genome_annotation"][
            "gff3_to_GenePred"
        ],
        genePredToBigGenePred_extra=config["ucsc_trackhub"][
            "process_genome_annotation"
        ]["GenePred_to_bgpInput"],
        bedToBigBed_extra=config["ucsc_trackhub"]["process_genome_annotation"][
            "bgpInput_to_bigGenePred"
        ],
    conda:
        "../envs/ucsc_tools.yaml"
    log:
        "results/get_genome/bigbed.log",
    shell:
        """
        gff3ToGenePred {params.gff3ToGenePred_extra} {input.gff} {params.intermediate_files[0]}
        genePredToBigGenePred {params.genePredToBigGenePred_extra} {params.intermediate_files[0]} stdout | sort -k1,1 -k2,2n > {params.intermediate_files[1]}
        bedToBigBed {params.bedToBigBed_extra} {params.intermediate_files[1]} {input.chrom_sizes} {output.bigGenePred}
        """


rule stringtie_gtf_to_bigGenePred:
    input:
        gtf="results/stringtie/stringtie_merged.gtf",
        chrom_sizes="results/get_genome/genome.chrom.sizes",
    output:
        bigGenePred="results/stringtie/stringtie_merged.bb",
    params:
        intermediate_files=lambda wildcards, input: multiext(
            subpath(input.gtf, strip_suffix=".gtf"), ".genePred", ".bgpInput"
        ),
        gtfToGenePred_extra=config["ucsc_trackhub"]["process_stringtie"][
            "gtf_to_genePred"
        ],
        genePredToBigGenePred_extra=config["ucsc_trackhub"]["process_stringtie"][
            "GenePred_to_bgpInput"
        ],
        bedToBigBed_extra=config["ucsc_trackhub"]["process_stringtie"][
            "bgpInput_to_bigGenePred"
        ],
    conda:
        "../envs/ucsc_tools.yaml"
    log:
        "results/stringtie/bigbed.log",
    shell:
        """
        gtfToGenePred {params.gtfToGenePred_extra} {input.gtf} {params.intermediate_files[0]}
        genePredToBigGenePred {params.genePredToBigGenePred_extra} {params.intermediate_files[0]} stdout | sort -k1,1 -k2,2n > {params.intermediate_files[1]}
        bedToBigBed {params.bedToBigBed_extra} {params.intermediate_files[1]} {input.chrom_sizes} {output.bigGenePred}
        """


rule gffcompare_gtf_to_bigGenePred:
    input:
        gtf="results/gffcompare/gffcmp.annotated.gtf",
        chrom_sizes="results/get_genome/genome.chrom.sizes",
    output:
        bigGenePred="results/gffcompare/gffcmp.annotated.bb",
    params:
        intermediate_files=lambda wildcards, input: multiext(
            subpath(input.gtf, strip_suffix=".gtf"), ".genePred", ".bgpInput"
        ),
        gtfToGenePred_extra=config["ucsc_trackhub"]["process_gffcompare"][
            "gtf_to_genePred"
        ],
        genePredToBigGenePred_extra=config["ucsc_trackhub"]["process_gffcompare"][
            "GenePred_to_bgpInput"
        ],
        bedToBigBed_extra=config["ucsc_trackhub"]["process_gffcompare"][
            "bgpInput_to_bigGenePred"
        ],
    conda:
        "../envs/ucsc_tools.yaml"
    log:
        "results/gffcompare/bigbed.log",
    shell:
        """
        gtfToGenePred {params.gtfToGenePred_extra} {input.gtf} {params.intermediate_files[0]}
        genePredToBigGenePred {params.genePredToBigGenePred_extra} {params.intermediate_files[0]} stdout | sort -k1,1 -k2,2n > {params.intermediate_files[1]}
        bedToBigBed {params.bedToBigBed_extra} {params.intermediate_files[1]} {input.chrom_sizes} {output.bigGenePred}
        """


rule faToTwoBit:
    input:
        "results/get_genome/genome.fasta",
    output:
        "results/get_genome/genome.2bit",
    params:
        extra=config["ucsc_trackhub"]["process_fasta"]["faToTwoBit"],
    log:
        "results/get_genome/fa_to_2bit.log",
    wrapper:
        "v7.1.0/bio/ucsc/faToTwoBit"


rule ucsc_trackhub:
    input:
        genome_2bit="results/get_genome/genome.2bit",
        genome_genePred="results/get_genome/genome.bb",
        tsrs="results/tsrs/TSRs.bb",
        stringtie="results/stringtie/stringtie_merged.bb",
        gffcompare="results/gffcompare/gffcmp.annotated.bb",
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

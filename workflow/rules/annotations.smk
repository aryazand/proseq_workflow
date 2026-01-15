rule stringtie:
    input:
        bam=fastq_process_align.get_bam_2,
    output:
        "results/stringtie/{sample}.gtf",
    conda:
        "../envs/stringtie.yml"
    log:
        "results/stringtie/{sample}.log",
    params:
        ref_gff=config["annotation"]["stringtie"]["ref_anno"],
        strandedness=lambda wildcards: (
            "--rf"
            if config["annotation"]["stringtie"]["strandedness"] == "reverse"
            else (
                "--fr"
                if config["annotation"]["stringtie"]["strandedness"] == "forward"
                else ""
            )
        ),
    shell:
        "stringtie {input.bam} -G {params.ref_gff} {params.strandedness} -o {output}"


rule stringtie_merge:
    input:
        gtfs=expand(
            "results/stringtie/{sample}.gtf", sample=fastq_process_align.samples.index
        ),
    output:
        gtf="results/stringtie/stringtie_merged.gtf",
    conda:
        "../envs/stringtie.yml"
    log:
        "results/stringtie/stringtie_merge.log",
    params:
        ref_gff=config["annotation"]["stringtie"]["ref_anno"],
    shell:
        """
        stringtie --merge -G {params.ref_gff} -o {output.gtf} {input.gtfs}
        """


rule gff_compare:
    input:
        gtfs=expand(
            "results/stringtie/{sample}.gtf", sample=fastq_process_align.samples.index
        ),
        ref_gff="results/get_genome/genome.gff",
    output:
        gtf="results/gffcompare/gffcmp.annotated.gtf",
    params:
        extra=config["annotation"]["gffcompare"]["extra"],
        output_prefix=subpath(output.gtf, strip_suffix=".annotated.gtf"),
    conda:
        "../envs/gffcompare.yml"
    log:
        "results/gffcompare/gffcompare.log",
    shell:
        """
        gffcompare {params.extra} -o {params.output_prefix} -r {input.ref_gff} {input.gtfs}
        """

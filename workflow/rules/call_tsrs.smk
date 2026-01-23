rule call_tsrs:
    input:
        "results/deeptools/5prime_coverage/{sample}_{strand}.bw",
    output:
        "results/tsrs/{sample}_{strand}.bed",
    container:
        "http://pricenas.biochem.uiowa.edu/arya/container_images/tsrDetectR.sif"
    log:
        "results/tsrs/{sample}_{strand}.log",
    script:
        "../scripts/call_tsrs.R"


rule consolidate_sample_tsrs:
    input:
        tsrs_plus="results/tsrs/{sample}_forward.bed",
        tsrs_minus="results/tsrs/{sample}_reverse.bed",
    output:
        "results/tsrs/{sample}.bed",
    log:
        "results/tsrs/{sample}.log",
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n | awk -F'\\t' -v sample={wildcards.sample} 'BEGIN{{OFS="\\t"}} {{$4=sample; print}}' > {output} 
        """


rule consolidate_all_tsrs:
    input:
        tsrs=expand(
            "results/tsrs/{sample}.bed", sample=fastq_process_align.samples.index
        ),
    output:
        "results/tsrs/TSRs.bed",
    log:
        "results/tsrs/TSRs.log",
    shell:
        """
        TEMP_FILE="results/tsrs/temp.bed"
        cat {input} | sort -k1,1 -k2,2n | awk -F'\\t' 'BEGIN{{OFS="\\t"}} {{$5=log($5); print}}' > $TEMP_FILE
        MAX_VAL=$(awk 'BEGIN {{max=0}} {{if ($5>max) max=$5}} END {{print max}}' $TEMP_FILE)
        awk -F'\\t' -v max_val=$MAX_VAL 'BEGIN{{OFS="\\t"}} {{$5=1000*$5/max_val; print}}' $TEMP_FILE > {output}
        rm $TEMP_FILE
        """


rule merge_tsrs:
    input:
        "results/tsrs/TSRs.bed",
    output:
        "results/tsrs/TSRs_merged.bed",
    params:
        extra=config["tsrs"]["bedtools_merge"]["extra"],
    log:
        "results/tsrs/TSRs.log",
    wrapper:
        "v8.1.1/bio/bedtools/merge"

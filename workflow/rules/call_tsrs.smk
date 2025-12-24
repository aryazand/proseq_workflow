rule call_tsrs:
    input:
        "results/deeptools/5prime_coverage/{sample}_{strand}.bw"
    output:
        "results/tsrs/{sample}_{strand}.bed"
    params:
        method = config["tsrs"]["tsrDetectR"]["method"],
        window = config["tsrs"]["tsrDetectR"]["window"],
        background = config["tsrs"]["tsrDetectR"]["background"],
        threshold = config["tsrs"]["tsrDetectR"]["threshold"],
        basedir = workflow.basedir
    conda:
        "../envs/tsr_detectr.yml"
    log:
        "results/tsrs/{sample}_{strand}.log"
    shell:
        """
        Rscript -e '
        options(timeout = 600)
        if (!requireNamespace("tsrDetectR", quietly=TRUE)) {{
          message("Installing tsrDetectR...")
          Sys.setenv(R_COMPILE_AND_INSTALL_PACKAGES="never")
          remotes::install_github("aryazand/tsrDetectR", upgrade="never", dependencies=TRUE)
        }}
        '
        mkdir -p results/tsrs
        Rscript {params.basedir}/scripts/call_tsrs.R -m {params.method} -w {params.window} -b {params.background} -t {params.threshold} -s {wildcards.strand} -i {input} -o {output}
        """
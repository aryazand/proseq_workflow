###############################
# This R script calls tsrs based on commandline arguments
################################

# Install required packages if not already installed
if (!requireNamespace("tsrDetectR", quietly = TRUE)) {
  message("Installing tsrDetectR...")
  remotes::install_github("aryazand/tsrDetectR", upgrade="never", dependencies=TRUE, INSTALL_opts =c("--no-lock"))
}

##################################
# PROCESS COMMANDLINE ARGUMENTS
##################################

in_file = snakemake@input[[1]]
out_file = snakemake@output[[1]]
strand = snakemake@wildcards[["strand"]]
param_method = snakemake@config[["tsrs"]][["tsrDetectR"]][["method"]]
param_window = snakemake@config[["tsrs"]][["tsrDetectR"]][["window"]]
param_background = snakemake@config[["tsrs"]][["tsrDetectR"]][["background"]]
param_threshold = snakemake@config[["tsrs"]][["tsrDetectR"]][["threshold"]]

################################
# CALL TSRs
################################

input_bw <- rtracklayer::import.bw(in_file, as = "RleList")

if(param_method == "maxtss"){

  tsr_func <- function(i){
    tsrs <- tsrDetectR::findtsr_maxtss(x = input_bw[[i]],
                          w = param_window,
                          background = param_background,
                          rel_threshold = param_threshold)  
    
    GenomicRanges::GRanges(
      seqnames = i,
      ranges = tsrs,
      strand = ifelse(strand %in% c("plus", "forward", "for"), "+", "-")
    )
  }

  tsrs <- purrr::map(names(input_bw), tsr_func) 

} 

# Collapse list of GRanges into one GRanges
tsrs <- do.call(c, tsrs)

# Export to BED file
rtracklayer::export.bed(tsrs, con = out_file)
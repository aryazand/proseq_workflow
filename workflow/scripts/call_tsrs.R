###############################
# This R script calls tsrs based on commandline arguments
################################

# Install required packages if not already installed
if (!requireNamespace("tsrDetectR", quietly = TRUE)) {
  message("Installing tsrDetectR...")
  remotes::install_github("aryazand/tsrDetectR", upgrade="never", dependencies=TRUE, INSTALL_opts =c("--no-lock"))
}

##################################
# PROCESS INPUT ARGUMENTS FROM SNAKEMAKE
##################################

in_file = snakemake@input[[1]]
out_file = snakemake@output[[1]]
strand = snakemake@wildcards[["strand"]]
param_method = snakemake@config[["tsrs"]][["tsrDetectR"]][["method"]]
param_window = snakemake@config[["tsrs"]][["tsrDetectR"]][["window"]]
param_background = snakemake@config[["tsrs"]][["tsrDetectR"]][["background"]]
param_threshold = snakemake@config[["tsrs"]][["tsrDetectR"]][["threshold"]]

if (exists(snakemake@config[["tsrs"]][["tsrDetectR"]][["regex"]])) {
  param_regex = snakemake@config[["tsrs"]][["tsrDetectR"]][["regex"]]
  if (is.null(param_regex) || param_regex == "") {
    param_regex = ".*"
  }
} else {
  param_regex = ".*"
}

################################
# LOAD COVERAGE DATA
################################

input_bw <- rtracklayer::import.bw(in_file, as = "RleList")

# Filter bw 
names_match <- grepl(param_regex, names(input_bw))
input_bw <- input_bw[names_match]
input_bw <- input_bw[sapply(input_bw, sum) > 0]


################################
# CALL TSRs
################################

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
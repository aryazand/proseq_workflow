###############################
# This R script calls tsrs based on commandline arguments
################################

##################################
# PROCESS COMMANDLINE ARGUMENTS
##################################

library("optparse")

parser <- OptionParser()
parser <- add_option(parser, c("-m", "--method"), action="store", type="character",
                    default="maxtss", help="method for identifying tsrs [default %maxtss]",
                    metavar="character") 
parser <- add_option(parser, c("-w", "--window"), action="store", type="integer",
                    default=11, help="window size [default %11]",
                    metavar="number") 
parser <- add_option(parser, c("-b", "--background"), action="store", type="integer",
                    default=101, help="window size for calculating background [default %101]",
                    metavar="number")
parser <- add_option(parser, c("-t", "--rel_threshold"), action="store", type="double",
                    help="threshold value for calling TSR [no default]",
                    metavar="number")
parser <- add_option(parser, c("-i", "--in_file"), action="store", type="character",
                    help="BigWig file to parse [no default]",
                    metavar="character")
parser <- add_option(parser, c("-o", "--out_file"), action="store", type="character",
                    help="output BED file [no default]",
                    metavar="character")
parser <- add_option(parser, c("-s", "--strand"), action="store", type="character",
                    help="strand [no default]",
                    metavar="character")

parsed_args <- parse_args(parser)

################################
# CALL TSRs
################################

input_bw <- rtracklayer::import.bw(parsed_args$in_file, as = "RleList")

if(parsed_args$method == "maxtss"){

  tsr_func <- function(i){
    tsrs <- tsrDetectR::findtsr_maxtss(x = input_bw[[i]],
                          w = parsed_args$window,
                          background = parsed_args$background,
                          rel_threshold = parsed_args$rel_threshold)  
    
    GenomicRanges::GRanges(
      seqnames = i,
      ranges = tsrs,
      strand = ifelse(parsed_args$strand == "plus", "+", "-")
    )
  }

  tsrs <- purrr::map(names(input_bw), tsr_func) 

} 

# Collapse list of GRanges into one GRanges
tsrs <- do.call(c, tsrs)

# Export to BED file
rtracklayer::export.bed(tsrs, con = parsed_args$out_file)
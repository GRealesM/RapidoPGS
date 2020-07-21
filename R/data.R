#' LD block architecture for European populations (hg19).
#'
#' A GRanges object containing the LD block for European ancestry, in hg19 build
#' This dataset was obtained from Belisa and Pickrell (2016), in bed format,
#' then converted to GRanges. See manuscript for more details.
#'
#' @format A GRanges object containing 1703 ranges
#' \describe{
#'   \item{seqnames}{chromosome}
#'   \item{ranges}{start and stop positions for the block}
#'   \item{strand}{genomic strand, irrelevant here
#' }
#' @source \url{https://bitbucket.org/nygcresearch/ldetect-data/src} 
"EUR_ld.blocks"


#' LD block architecture for European populations (hg38).
#'
#' A GRanges object containing the LD block for European ancestry, in hg38 build
#' This dataset was obtained from Belisa and Pickrell (2016), in bed format,
#' then liftovered to hg38 using UCSC liftOver tool, then converted to GRanges. 
#' See manuscript for more details.
#'
#' @format A GRanges object containing 1625 ranges
#' \describe{
#'   \item{seqnames}{chromosome}
#'   \item{ranges}{start and stop positions for the block}
#'   \item{strand}{genomic strand, irrelevant here
#' }
#' @source \url{https://bitbucket.org/nygcresearch/ldetect-data/src} 
#' 
"EUR_ld.blocks38"


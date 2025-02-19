#' LD block architecture for European populations (hg19).
#'
#' A GRanges object containing the LD block for European ancestry, in hg19 build.
#' This dataset was obtained from \doi{10.1093/bioinformatics/btv546}, in bed format,
#' then converted to GRanges. See manuscript for more details.
#'
#' @format A GRanges object containing 1703 ranges
#' \describe{
#'   \item{seqnames}{chromosome}
#'   \item{ranges}{start and stop positions for the block}
#'   \item{strand}{genomic strand, irrelevant here}
#' }
#' @source \url{https://bitbucket.org/nygcresearch/ldetect-data/src} 
"EUR_ld.blocks19"


#' LD block architecture for European populations (hg38).
#'
#' A GRanges object containing the LD block for European ancestry, in hg38 build.
#' This dataset was obtained from \doi{10.1101/2022.03.04.483057}, in bed format,
#'  then transformed to GRanges.
#' See manuscript for more details.
#'
#' @format A GRanges object containing 1361 ranges
#' \describe{
#'   \item{seqnames}{chromosome}
#'   \item{ranges}{start and stop positions for the block}
#'   \item{strand}{genomic strand, irrelevant here}
#' }
#' @source \url{https://github.com/jmacdon/LDblocks_GRCh38/tree/master/data} 
#' 
"EUR_ld.blocks38"

#' Subset of Michailidou BRCA GWAS sumstat dataset.
#'
#' A data.table containing a subset of \doi{10.1038/nature24284} breast cancer summary statistic dataset, in hg38 build.
#' This dataset is freely available in GWAS catalog (see link below). We removed unnecessary and all-missing columns, and rows
#' with missing data at hm_beta and hm_effect_allele_frequency, and took a random sample of 100,000 SNPs without replacement.
#'
#' @format A data.table object containing 100,000 SNPs
#' \describe{
#'   \item{hm_rsid}{rsids, or SNP ids}
#'   \item{hm_chrom}{chromosome}
#'   \item{hm_pos}{base position, in hg38}
#'   \item{hm_other_allele}{reference, or non-effect allele}
#'   \item{hm_effect_allele}{alternative, or effect allele}
#'   \item{hm_beta}{beta, log(OR), or effect size}
#'   \item{hm_effect_allele_frequency}{effect allele frequency}
#'   \item{standard_error}{standard error of beta}
#'   \item{p_value}{p-value}
#' }
#' @source \url{http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004988/harmonised/29059683-GCST004988-EFO_0000305.h.tsv.gz} 
#' 
"michailidou38"

#' Subset of Michailidou BRCA GWAS sumstat dataset.
#'
#' A data.table containing a subset of \doi{10.1038/nature24284} breast cancer summary statistic dataset, in hg19 build.
#' This dataset is freely available in GWAS catalog (see link below). We used "chromosome", "base_pair_location columns", removed unnecessary and all-missing columns, and took a random sample of 100,000 SNPs without replacement.
#'
#' @format A data.table object containing 100,000 SNPs
#' \describe{
#' SNPID, CHR, BP, REF, ALT, ALT_FREQ, BETA, SE, P
#'   \item{SNPID}{rsids, or SNP ids}
#'   \item{CHR}{chromosome}
#'   \item{BP}{base position, in hg38}
#'   \item{REF}{reference, or non-effect allele}
#'   \item{ALT}{alternative, or effect allele}
#'   \item{ALT_FREQ}{effect allele frequency}
#'   \item{BETA}{beta, log(OR), or effect size}
#'   \item{SE}{standard error of beta}
#'   \item{P}{p-value}
#' }
#' @source \url{http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004988/harmonised/29059683-GCST004988-EFO_0000305-build37.f.tsv.gz} 
#' 
"michailidou19"

#'@importFrom susieR susie_rss
NULL


##' Posterior inclusion probabilities under a multiple causal variant model
##'
##' A wrapper around susie_RSS from the susieR package. See
##' \url{https://stephenslab.github.io/susieR} for the package and Wang et al,
##' JRSSB 2020 for the paper describing the model
##' \url{https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssb.12388}
##'
##' NB: intended as an internal function, called from other functions in the
##' RapidoPGS package.
##' @param d data.table with columns SNP, BETA, SE (other columns allowed, and
##'   will be ignored)
##' @param LD LD matrix with dimnames matching "snp" column
##' @param nref number of individuals used to create the LD matrix
##' @param max_it maximum number of iterations for susie. The default of 100
##'   should be plenty, but increase if there is not convergence
##' @inheritParams computePGS
##' @return returns pip for each SNP, which we will use as a snp weight in
##'   generating the PGS (=pip * BETA)
##' @examples
##' data(michailidou) # load example data
##' d=michailidou[hm_chrom==3 & abs(hm_pos-27303612) < 1e+5] # focus on a window of association
##' setnames(d, old = c("hm_rsid", "hm_chrom", "hm_pos", "hm_other_allele",
##'   "hm_effect_allele", "hm_beta", "hm_effect_allele_frequency",
##'   "standard_error", "p_value"), new=c("SNPID","CHR", "BP","REF", "ALT",
##'   "BETA", "ALT_FREQ", "SE", "P")) # rename
##' LD=(1 - abs(outer(d$ALT_FREQ,d$ALT_FREQ, "-"))) *
##'    outer(sign(d$BETA), sign(d$BETA), "*")
##' dimnames(LD)=list(d$SNPID,d$SNPID)
##' susie_pip(d, LD)
##' @author Chris Wallace
susie_pip=function(d,LD,nref=503,pi_i=1e-4,max_it=100) {
  ## checks
  nd <- names(d)
  if(!all(c("SNPID","BETA","SE") %in% nd))
    stop("require columns SNP, BETA, SE")
   ## snps should be unique
   if("SNPID" %in% nd && is.factor(d$SNPID))
     stop("dataset ",suffix,": SNPID should be a character vector but is a factor")
   if("SNPID" %in% nd && any(duplicated(d$SNPID)))
     stop("dataset ",suffix,": duplicated SNPIDs found")
  ## SNPIDs should be in LD matrix
  if(!is.matrix(LD) ||
     !all(d$SNPID %in% colnames(LD)) ||
     !identical(rownames(LD),colnames(LD)) ||
     any(LD < -1) || any(LD > 1))
    stop("LD should be a correlation matrix containing all SNPIDs listed in d$SNPID (match by rownames, colnames)")

  z=d$BETA/d$SE
  res=susie_rss(z, LD[d$SNPID,d$SNPID], z_ld_weight = 1/nref,
               ,null_weight=1 - length(d$SNPID)*pi_i
               ,estimate_prior_method="EM"
                ## ,verbose=TRUE
               ,max_iter=max_it )
  pip=res$pip[ -length(res$pip) ]
  names(pip)=d$SNPID
  pip
}

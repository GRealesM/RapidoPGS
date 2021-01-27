##' @importFrom stats coef lm na.omit pnorm
NULL

##' Helper function to sum logs without loss of precision
##'
##' Sums logs without loss of precision
##' This function is verbatim of its namesake in cupcake package (github.com/ollyburren/cupcake/)
##'
##' @param x a vector of logs to sum
##' @return a scalar
##' @author Chris Wallace
logsum <- function(x) {
    my.max <- max(x) ##take out the maximum value in log form)
    my.res <- my.max + log(sum(exp(x - my.max )))
    return(my.res)
}


##' Estimate trait standard deviation given vectors of variance of coefficients,  MAF and sample size
##'
##' Estimate is based on var(beta-hat) = var(Y) / (n * var(X))
##' var(X) = 2* maf*(1-maf)
##' so we can estimate var(Y) by regressing n*var(X) against 1/var(beta)
##' This function is verbatim from its namesake in coloc package (github.com/chr1swallace/coloc/), by Chris Wallace
##' 
##' @title Estimate trait variance, internal function
##' @param vbeta vector of variance of coefficients
##' @param maf vector of MAF (same length as vbeta)
##' @param n sample size
##' @return estimated standard deviation of Y
##' @author Chris Wallace
sdY.est <- function(vbeta, maf, n) {
  warning("estimating sdY from maf and varbeta, please directly supply sdY if known")
  oneover <- 1/vbeta
  nvx <- 2 * n * maf * (1-maf)
  m <- lm(nvx ~ oneover - 1)
  cf <- coef(m)[['oneover']]
  if(cf < 0)
    stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
  return(sqrt(cf))
}


##' Compute posterior probabilities using Wakefield's approximate Bayes Factors for quantitative traits
##'
##' \code{wakefield_pp_quant} computes posterior probabilities for a given SNP to be causal for a given SNP under the assumption of a single causal variant.
##'
##' This function was adapted from \code{wakefield_pp} in cupcake package (github.com/ollyburren/cupcake/)
##'
##' @param beta a vector of effect sizes (\eqn{\beta}) from a quantitative trait GWAS
##' @param se vector of standard errors of effect sizes (\eqn{\beta})
##' @param sdY a scalar of the standard deviation given vectors of variance of coefficients,  MAF and sample size. Can be calculated using \code{sdY.est}
##' @param sd.prior a scalar representing our prior expectation of \eqn{\beta} (DEFAULT 0.15).
##' @param pi_i a scalar representing the prior probability (DEFAULT \eqn{1 \times 10^{-4}})
##' The method assumes a normal prior on the population log relative risk centred at 0 and the DEFAULT
##' value sets the variance of this distribution to 0.04, equivalent to a 95\%  belief that the true relative risk
##' is in the range of 0.66-1.5 at any causal variant.
##' @return a vector of posterior probabilities.
##' @author Guillermo Reales, Chris Wallace
wakefield_pp_quant <- function(beta, se, sdY, sd.prior=0.15, pi_i=1e-4) { 
  # compute V
  V <- se^2
  # Compute z too
  z <- beta/se
  # Multiply prior by sdY
  prior <- sd.prior*sdY
  ## Shrinkage factor: ratio of the prior variance to the total variance
  r <- prior^2 / (prior^2 + V)
  ## Approximate BF
  lABF = 0.5 * (log(1-r) + (r * z^2))
  ## tABF - to add one we create another element at the end of zero for which pi_i is 1
  tABF <- c(lABF,0)
  vpi_i<-c(rep(pi_i,length(lABF)),1)
  sBF <- logsum(tABF + log(vpi_i))
  exp(lABF+log(pi_i)-sBF)
}


##' compute posterior probabilities using Wakefield's approximate Bayes Factors
##' \code{wakefield_pp} computes posterior probabilities for a given SNP to be causal for a given SNP under the assumption of a single causal variant.
##'
##' This function was adapted from its namesake in cupcake package (github.com/ollyburren/cupcake/) to no longer require allele frequencies.
##'
##' @param beta a vector of effect sizes (\eqn{\beta}) from a quantitative trait GWAS
##' @param se vector of standard errors of effect sizes (\eqn{\beta})
##' @param pi_i a scalar representing the prior probability (DEFAULT \eqn{1 \times 10^{-4}})
##' @param sd.prior a scalar representing our prior expectation of \eqn{\beta} (DEFAULT 0.2).
##' The method assumes a normal prior on the population log relative risk centred at 0 and the DEFAULT
##' value sets the variance of this distribution to 0.04, equivalent to a 95\%  belief that the true relative risk
##' is in the range of 0.66-1.5 at any causal variant.
##' @return a vector of posterior probabilities.
##' @author Olly Burren, Chris Wallace, Guillermo Reales
wakefield_pp <- function(beta, se, pi_i=1e-4,sd.prior=0.2) {
  if(length(beta) != length(se))
    stop("beta and se must be vectors of the same size")
  # compute V
  V <- se^2
  # compute z
  z <- beta/se
  ## Shrinkage factor: ratio of the prior variance to the total variance
  r <- sd.prior^2 / (sd.prior^2 + V)
  ## Approximate BF
  lABF = 0.5 * (log(1-r) + (r * z^2))
  ## tABF - to add one we create another element at the end of zero for which pi_i is 1
  tABF <- c(lABF,0)
  vpi_i<-c(rep(pi_i,length(lABF)),1)
  sBF <- logsum(tABF + log(vpi_i))
  exp(lABF+log(pi_i)-sBF)
}



##' Retrieve GWAS summary datasets from GWAS catalog
##' '\code{gwascat.download} takes a PMID from the user and downloads the associated summary statistics datasets published in GWAS catalog
##' 
##' This function, takes PUBMED ids as an input, searches at the GWAS catalog
##' for harmonised datasets associated to that, interactively asking the  
##' user to choose if there are more than one, and fetches the dataset. 
##'
##' If multiple files are available for the same study, R will prompt an interactive
##' dialogue to select a specific file, by number. If you know the number and
##' prefer to select it automatically, you can provide it using file argument.
##'
##' @param ID a numeric. A PubMed ID (PMID) reference number from a GWAS paper.
##' @param filenum a numeric. If multiple files are available, which one to choose? If NULL (DEFAULT), R will prompt an interactive prompt, asking for the number.
##' @param hm_only a logical. Should GWAS catalog harmonised columns be retained?
##' @return a data.table containing the dataset.
##' @import data.table curl RCurl
##' @export
##' @author Guillermo Reales
##' @examples 
##' 
##' \dontrun{
##'	ds <- gwascat.download(29059683, hm_only = FALSE) # This should work: Michailidou dataset
##'	wrongds <- gwascat.download(01223247236) # This shouldn't work: The Empress pub phone number
##'}
##'

gwascat.download <- function(ID, filenum = NULL, hm_only=TRUE){
		gwc.manifest <- fread("https://www.ebi.ac.uk/gwas/api/search/downloads/studies_alternative")
		
		study.manifest <- gwc.manifest[ID == PUBMEDID,] 
		
		if(nrow(study.manifest) == 0) stop("Please provide a valid PUBMED ID")
		if(nrow(study.manifest) > 1){
			message("There are ",nrow(study.manifest)," datasets associated with this PUBMED ID")
			if(!is.null(filenum)){
				if(!is.numeric(filenum) || length(filenum) > 1 || filenum > nrow(study.manifest)){
					stop("Please provide a valid file number!")
				} else{
					study_acc  <- filenum
				}
			} else{ 
				print(study.manifest[,"STUDY ACCESSION"])
				study_acc  <- readline(prompt=paste("Please select accession (from 1 to ", nrow(study.manifest),"):", sep=""))
				study_acc  <- as.integer(study_acc)
				while(!study_acc %in% 1:nrow(study.manifest)){
					message("Oops! Seems that your number is not in the options. Please try again.")
					print(study.manifest[,"STUDY ACCESSION"])
					study_acc  <- readline(prompt=paste("Please select accession (from 1 to ", nrow(study.manifest),"):", sep=""))
					study_acc  <- as.integer(study_acc)
					}

			}
		study.manifest <- study.manifest[study_acc,]
		}
		url <- paste("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/", gsub(" ","", study.manifest$`FIRST AUTHOR`), "_", ID,"_",study.manifest$`STUDY ACCESSION`, "/harmonised/",ID,"-",study.manifest$`STUDY ACCESSION`,"-",sapply(strsplit(study.manifest$MAPPED_TRAIT_URI, "/"), `[`,5),".h.tsv.gz", sep="")
		if(!RCurl::url.exists(url)) stop("The file you requested is unavailable. This may be due to the fact that public and harmonised summary statistics do not exist. Please check at GWAS catalog website.")
		message("Retrieving dataset for ",study.manifest$`DISEASE/TRAIT`,", by ", study.manifest$`FIRST AUTHOR`,", from ",study.manifest$DATE, ", published at ", study.manifest$JOURNAL,", with accession ", study.manifest$`STUDY ACCESSION`,".")
		
		ds <- fread(url)
		if(hm_only){
		hmcols <- grep("hm_",names(ds), value=TRUE)
		ds  <- ds[,..hmcols]
		}
		return(ds)
}



##' Compute PGS from GWAS summary statistics using posteriors from Wakefield's approximate Bayes Factors
##' 
##' '\code{rapidopgs_single} computes PGS from a from GWAS summary statistics using posteriors from Wakefield's approximate Bayes Factors
##' 
##' This function will take a GWAS summary statistic dataset as an input,
##' will assign align it to a reference panel file (if provided), then it will assign 
##' SNPs to LD blocks and compute Wakefield's ppi by LD block, then will use it 
##' to generate PGS weights by multiplying those posteriors by effect sizes (\eqn{\beta}). 
##' Optionally, it will filter SNPs by a custom filter on ppi and then recalculate weights, to improve accuracy.
##' 
##' Alternatively, if filt_threshold is larger than one, RápidoPGS will select the top
##' \code{filt_threshold} SNPs by absolute weights (note, not ppi but weights).
##' 
##' The GWAS summary statistics file to compute PGS using our method must contain the following minimum columns, with these exact column names:
##' \describe{
##'   \item{CHR}{Chromosome}
##'   \item{BP}{Base position (in GRCh37/hg19 or GRCh38/hg38). If using hg38, use build = "hg38" in parameters}
##'   \item{REF}{Reference, or non-effect allele}
##'   \item{ALT}{Alternative, or effect allele, the one \eqn{\beta} refers to}
##'   \item{ALT_FREQ}{Minor/ALT allele frequency in the tested population, or in a close population from a reference panel. Required for Quantitative traits only}
##'   \item{BETA}{\eqn{\beta} (or log(OR)), or effect sizes}
##'   \item{SE}{standard error of \eqn{\beta}}
##' }
##'
##' If a reference is provided,  it should have 5 columns: CHR, BP,
##' SNPID, REF, and ALT. Also, it should be in the same build as 
##' the summary statistics. In both files, column order does not matter.
##' @param data a data.table containing GWAS summary statistic dataset
##'   with all required information.
##' @param N a scalar representing the sample in the study, or a string indicating 
##'   the column name containing it. Required for quantitative traits only.
##' @param trait a string specifying if the dataset corresponds to a case-control
##'   ("cc") or a quantitative trait ("quant") GWAS. If trait = "quant", an 
##'   ALT_FREQ column is required.
##' @param build a string containing the genome build of the dataset,
##'   either "hg19" (for hg19/GRCh37) or "hg38" (hg38/GRCh38). DEFAULT
##'   "hg19".
##' @param pi_i a scalar representing the prior probability (DEFAULT:
##'   \eqn{1 \times 10^{-4}}).
##' @param sd.prior the prior specifies that BETA at causal SNPs
##'   follows a centred normal distribution with standard deviation
##'   sd.prior. Sensible and widely used DEFAULTs are 0.2 for case
##'   control traits, and 0.15 * var(trait) for quantitative (selected
##'   if trait == "quant").
##' @param filt_threshold a scalar indicating the ppi threshold (if
##'   \code{filt_threshold} < 1) or the number of top SNPs by absolute
##'   weights (if \code{filt_threshold} >= 1) to filter the dataset
##'   after PGS computation. If NULL (DEFAULT), no thresholding will
##'   be applied.
##' @param recalc a logical indicating if weights should be
##'   recalculated after thresholding. Only relevant if \code{filt_threshold}
##'   is defined.
##' @param reference a string indicating the path of the reference file 
##'   SNPs should be filtered and aligned to, see Details.
##' @return a data.table containing the formatted sumstats dataset with
##'   computed PGS weights.
##' @import data.table 
##' @importFrom bigsnpr snp_match
##' @importFrom GenomicRanges GRanges findOverlaps
##' @importFrom IRanges IRanges
##' @export
##' @author Guillermo Reales, Chris Wallace
##' @examples
##' sumstats <- data.table(SNPID=c("rs139096444","rs3843766","rs61977545", "rs544733737",
##'			"rs2177641", "rs183491817", "rs72995775","rs78598863", "rs1411315"), 
##'			CHR=c("4","20","14","2","4","6","6","21","13"), 
##'			BP=c(1479959, 13000913, 29107209, 203573414, 57331393, 11003529, 149256398, 
##'					25630085, 79166661), 
##'			REF=c("C","C","C","T","G","C","C","G","T"), 
##'			ALT=c("A","T","T","A","A","A","T","A","C"), 
##'			BETA=c(0.012,0.0079,0.0224,0.0033,0.0153,0.058,0.0742,0.001,-0.0131),
##'			SE=c(0.0099,0.0066,0.0203,0.0171,0.0063,0.0255,0.043,0.0188,0.0074))
##'
##' PGS  <- rapidopgs_single(sumstats,  trait = "cc")
rapidopgs_single <- function(data,
                             N=NULL,
                             trait=c("cc","quant"),
                             build = "hg19",
                             pi_i= 1e-04,
                             sd.prior=if(trait == "quant") {0.15} else {0.2},
                             filt_threshold = NULL,
                             recalc=TRUE,
                             reference=NULL
){
  
  if(!"data.table" %in% class(data))
    data <- as.data.table(data)
  if(!trait %in% c("cc", "quant")) stop("Please, specify your study type, choose case-control ('cc') or quantitative ('quant').")
  if(length(trait) != 1) stop("Please select only one study type")
  if(trait == "quant" && is.null(N)) stop("N (sample size) is required for quantitative traits, please provide them, either as an integer or as column name containing it.")
  if(trait == "quant"){
    mincol <- c("CHR","BP", "REF","ALT","BETA", "SE", "ALT_FREQ")
    if(!all(mincol %in% names(data)))
      stop("All minimum columns should be present in the summary statistics dataset. Please check, missing columns: ",
           paste(setdiff(mincol,names(data)), collapse=", "))  
    ds <- copy(data) # avoid modifying input data.table
    ds  <- na.omit(ds, cols=c("CHR","BP", "BETA", "SE", "ALT_FREQ")) # Remove NA in relevant columns
  } else{
    ## Here's a list of columns that the dataset must have, otherwise it will fail
    mincol <- c("CHR","BP", "REF","ALT","BETA", "SE")
    if(!all(mincol %in% names(data)))
      stop("All minimum columns should be present in the summary statistics dataset. Please check, missing columns: ",
           paste(setdiff(mincol,names(data)), collapse=", "))
    ds <- copy(data) # avoid modifying input data.table
    ds  <- na.omit(ds, cols=c("CHR","BP", "BETA", "SE")) # Remove NA in relevant columns
    ds <- ds[SE != 0,] # Remove SNPs with exactly zero SE (sign of problems in the data)
  }	
  if(!"SNPID" %in% names(ds)){
    ds[,SNPID:=paste(CHR,BP, sep=":")] # SNPID is not strictly required to be provided. If it's not, we create it using CHR:BP
  }
  ## First step, align to reference, if reference is provided
  if(!is.null(reference)){
    ## Next step is to filter and align our alleles and their effects to the hapmap3 reference, which I have already formatted for our purposes.
    message("Filtering SNPs...")
    refset <- fread(reference)
    mincolref <- c("CHR","BP", "SNPID", "REF","ALT")
    if(!all(mincolref %in% names(refset))) stop("All minimum columns should be present in the reference file. Please check, missing columns: ", paste(setdiff(mincol,names(refset)), collapse=", "))
    
    setnames(refset, old=c("CHR","BP", "SNPID", "REF","ALT"), new=c("chr","pos","id","a0","a1"))
    setnames(ds, old=c("CHR","BP", "BETA", "REF","ALT"), new=c("chr","pos","beta","a0","a1"))
    
    message("Matching and aligning SNPs to the reference")
    info_snp <- snp_match(ds,refset)
    
    message("Original sumstat had ", nrow(sumstats.chr)," for chr",chr,". After matching ", nrow(info_snp.chr)," remained, and will be used for further steps.")
    ds <- data.table(info_snp)
    ds[,SNPID:=id][, c("_NUM_ID_.ss","_NUM_ID_", "id"):=NULL] # Replaces SNPID in the sumstat by SNPID in the reference, and removes snp_match cols.
    setnames(ds, old=c("chr","pos","beta","a0","a1"), new=c("CHR","BP", "BETA", "REF","ALT"))
  }
  
  # Assign ld.blocks, in case they werent computed yet
  if(!"ld.block" %in% names(ds)){
    if(build == "hg19"){ 
      blranges <- RapidoPGS::EUR_ld.blocks
    }else if(build == "hg38"){ 
      blranges <- EUR_ld.blocks38
    }else{ 
      stop("RápidoPGS only accepts hg19 or hg38 at the moment, please check.")
    }		
    message("Assigning LD blocks...")
    snpranges <- GRanges(seqnames=paste("chr",ds$CHR, sep=""), ranges=IRanges(start=ds$BP, end=ds$BP, names=ds$SNPID), strand="*")
    ds[,ld.block:=findOverlaps(snpranges, blranges, select='last')]
    message("Done!")
  }
  ds <- ds[!is.na(ld.block),] # Remove SNPs without assigned block
  
  if(trait == "quant"){
    if(is.numeric(N) && length(N) == 1) { # In case N is supplied as a number
      message("Computing a RápidoPGS-single model for a quantitative trait with ", N, " individuals...")
      ds[,sdY:=sdY.est(vbeta=SE^2, maf=ALT_FREQ, n=N)][,ppi:=wakefield_pp_quant(BETA,SE,sdY,sd.prior), by="ld.block"][,weight:=ppi*BETA]
      if(!is.null(filt_threshold)){
        ds  <- ds[ds$ppi > filt_threshold,]
        if(recalc){
          ds[,sdY:=sdY.est(vbeta=SE^2, maf=ALT_FREQ, n=N)][,ppi:=wakefield_pp_quant(BETA,SE,sdY,sd.prior), by="ld.block"][,weight:=ppi*BETA]		
        }
      }
    } else{ # In case column name is supplied
      if(is.character(N) && length(N) == 1){
        Nco <- N
        message("Computing a RápidoPGS-single model for a quantitative trait multiple with multiple N, supplied by column ", Nco,". Mmax N: ", max(ds[,get(Nco)]), ", min N = ", min(ds[,get(Nco)]), ", and mean N: ",mean(ds[,get(Nco)]), "...")
        ds[,sdY:=sdY.est(vbeta=SE^2, maf=ALT_FREQ, n=get(Nco))][,ppi:=wakefield_pp_quant(BETA,SE,sdY,sd.prior), by="ld.block"][,weight:=ppi*BETA]
        if(!is.null(filt_threshold)){
          ds  <- ds[ds$ppi > filt_threshold,]
          if(recalc){
            ds[,sdY:=sdY.est(vbeta=SE^2, maf=ALT_FREQ, n=get(Nco))][,ppi:=wakefield_pp_quant(BETA,SE,sdY,sd.prior), by="ld.block"][,weight:=ppi*BETA]		
          }
        }
      }
    }
  }
  if(trait == "cc"){
    message("Computing a RápidoPGS-single model for a case-control dataset...")
    ds[,ppi:=wakefield_pp(beta=BETA, se=SE,pi_i, sd.prior), by = "ld.block"][, weight:=ppi*BETA]
    if(!is.null(filt_threshold)){
      if(filt_threshold < 1){
        ds  <- ds[ds$ppi > filt_threshold,]
      } else {
        ds <- ds[order(-rank(abs(weight))),][1:filt_threshold,] 
      }
      if(recalc){
        ds[,ppi:=wakefield_pp(beta=BETA, se=SE,pi_i, sd.prior), by = "ld.block"][, weight:=ppi*BETA]
      }
    }
  }
  return(ds)
}	


##' Posterior inclusion probabilities under a multiple causal variant model
##'
##' A wrapper around susie_RSS from the susieR package. See
##' \url{https://stephenslab.github.io/susieR} for the package and Wang et al.
##' JRSSB 2020 for the paper describing the model
##' \url{https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssb.12388}
##'
##' NB: intended as an internal function, called from other functions in the
##' RapidoPGS package.
##' @param d data.table with columns SNPID, BETA, SE (other columns allowed, d
##'   will be ignored)
##' @param LD LD matrix with dimnames matching "SNPID" column
##' @param nref number of individuals used to create the LD matrix
##' @param pi_i prior probability that a given variant is causal. Default 1e-4
##' @param sd.prior Standard deviation prior of the trait, if NULL (default), it  will be estimated
##' @param max_it maximum number of iterations for susie. The default of 100
##'   should be plenty, but increase if there is not convergence
##' @importFrom susieR susie_rss
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
##' @author Chris Wallace, Guillermo Reales

susie_pip <- function(d,LD,nref=503,pi_i=1e-4, sd.prior=NULL, max_it=100) {
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
     !all(LD >= -1 && LD <= 1))
    stop("LD should be a correlation matrix containing all SNPIDs listed in d$SNPID (match by rownames, colnames)")
  
  z <- d$BETA/d$SE
  if(is.null(sd.prior)){
    res=susie_rss(z, LD[d$SNPID,d$SNPID], z_ld_weight = 1/nref,
                  ,null_weight=1 - length(d$SNPID)*pi_i
                  ,estimate_prior_method="EM"
                  ## ,verbose=TRUE
                  ,max_iter=max_it )
  } else{
    res=susie_rss(z, LD[d$SNPID,d$SNPID], z_ld_weight = 1/nref,
                  ,null_weight=1 - length(d$SNPID)*pi_i
                  ,estimate_prior_method="EM"
                  ,estimate_prior_variance=FALSE  ### <- stop internal estimation
                  ,prior_variance = sd.prior^2 ### <- use custom sd.prior
                  ## ,verbose=TRUE
                  ,max_iter=max_it )
  }
  pip <- res$pip[ -length(res$pip) ]
  names(pip) <- d$SNPID
  pip
}

##' Compute PGS from GWAS summary statistics using Bayesian sum of single-effect 
##' (SuSiE) linear regression using z scores
##' 
##' '\code{rapidopgs_multi} computes PGS from a from GWAS summary statistics 
##' using Bayesian sum of single-effect (SuSiE) linear regression using z scores
##' 
##' This function will take a GWAS summary statistic dataset as an input,
##' will assign LD blocks to it, then use a preset reference panel in Plink format 
##' to compute LD matrices for each block. Then SuSiE method will be used to 
##' compute posterior probabilities of variants to be causal and generate PGS
##' weights by multiplying those posteriors by effect sizes (\eqn{\beta}). 
##' Unlike \code{rapidopgs_single}, this approach will assume one or more causal variants
##' Also, \code{rapidopgs_multi} does not require allele frequencies nor sample sizes
##' 
##' The GWAS summary statistics file to compute PGS using our method must contain
##' the following minimum columns, with these exact column names:
##' \describe{
##'   \item{CHR}{Chromosome}
##'   \item{BP}{Base position (in GRCh37/hg19).
##'   \item{REF}{Reference, or non-effect allele}
##'   \item{ALT}{Alternative, or effect allele, the one \eqn{\beta} refers to}
##'   \item{BETA}{\eqn{\beta} (or log(OR)), or effect sizes}
##'   \item{SE}{standard error of \eqn{\beta}}
##'   \item{P}{P-value for the association test}
##' }
##' Other columns are allowed, and will be ignored.
##' 
##' Reference panel should be divided by chromosome, in Plink format.
##' Both reference panel and summary statistic dataset should be in GRCh37/hg19. 
##' 
##' @param data a data.table containing GWAS summary statistic dataset
##'   with all required information.
##' @param reference a string representing the path to the directory containing 
##'   the reference panel (eg. "../ref-data/").
##' @param ancestry a string indicating the ancestral population (DEFAULT: "EUR")
##' @param pi_i a scalar representing the prior probability (DEFAULT:
##'   \eqn{1 \times 10^{-4}}).
##' @param iterations number of maximum iterations at SuSie step. DEFAULT:100, 
##'    which should be plenty to ensure convergence. 
##' @param ncores a numeric specifying the number of cores to be used.
##'    If not set, it won't use parallelisation.
##' @param alpha.block a numeric threshold for minimum P-value in LD blocks.
##'    Blocks with minimum P above \code{alpha.block} will be skipped. Default: 1e-4.
##' @param alpha.snp a numeric threshold for P-value pruning within LD block.
##'    SNPs with P above \{alpha.snp} will be removed. Default: 0.01.
##' @return a data.table containing the sumstats dataset with
##'   computed PGS weights.
##' @import data.table 
##' @importFrom bigsnpr snp_match snp_cor snp_readBed snp_attach
##' @importFrom GenomicRanges GRanges findOverlaps
##' @importFrom IRanges IRanges
##' @export
##' @author Guillermo Reales, Chris Wallace
##' @examples
##' \dontrun{
##' sumstats <- data.table(SNPID=c("rs139096444","rs3843766","rs61977545", "rs544733737",
##'			"rs2177641", "rs183491817", "rs72995775","rs78598863", "rs1411315"), 
##'			CHR=c("4","20","14","2","4","6","6","21","13"), 
##'			BP=c(1479959, 13000913, 29107209, 203573414, 57331393, 11003529, 149256398, 
##'					25630085, 79166661), 
##'			REF=c("C","C","C","T","G","C","C","G","T"), 
##'			ALT=c("A","T","T","A","A","A","T","A","C"), 
##'			ALT_FREQ=c(0.2611,0.4482,0.0321,0.0538,0.574,0.0174,0.0084,0.0304,0.7528),
##'			BETA=c(0.012,0.0079,0.0224,0.0033,0.0153,0.058,0.0742,0.001,-0.0131),
##'			SE=c(0.0099,0.0066,0.0203,0.0171,0.0063,0.0255,0.043,0.0188,0.0074),
##'			P=c(0.2237,0.2316,0.2682,0.8477,0.01473,0.02298,0.08472,0.9573,0.07535))
##' PGS  <- rapidopgs_multi(sumstats, reference = "ref-data/", ncores=2)
##'}

rapidopgs_multi <- function(data, reference, ancestry="EUR", pi_i = 1e-04, iterations = 100, ncores=1, alpha.block=1e-4, alpha.snp=0.01, sd.prior=NULL){
  
  mincol <- c("CHR","BP", "REF","ALT","BETA", "SE","P")
  if(!all(mincol %in% names(data)))
    stop("All minimum columns should be present in the summary statistics dataset. Please check, missing columns: ",
         paste(setdiff(mincol,names(data)), collapse=", "))
  if(!file.exists(paste(reference,"chr1.bed", sep=""))) # Temp fix to detect panel		
    stop("No reference panel detected. Please check.")
  
  message("Running RápidoPGS with multiple causal variant assumption.")
  
  ds <- copy(data) # avoid modifying input data.table
  ds[,SNPID:=paste(CHR,BP, sep = ":")]
  ds  <- na.omit(ds, cols=c("CHR","BP", "BETA", "SE", "P"))
  # Assign ld.blocks, in case they werent computed yet.
  # Note that in this case only hg19 is admitted.
  if(!"ld.block" %in% names(ds)){
    blranges <- EUR_ld.blocks	
    message("Assigning LD blocks...")
    snpranges <- GRanges(seqnames=paste("chr",ds$CHR, sep=""), ranges=IRanges(start=ds$BP, end=ds$BP, names=ds$SNPID), strand="*")
    ds[,ld.block:=findOverlaps(snpranges, blranges, select='last')]
    message("Done!")
  }
  ## We ensure that rows without assigned block are removed
  if(any(is.na(ds$ld.block))){
    # warning(length(ds$ld.block[is.na(ds$ld.block)]), " SNPs didn't have any LD block assigned and will be removed")
    ds <- na.omit(ds, cols="ld.block")
  }
  
  # We don't need all columns, so we can filter them by name
  ds <- ds[,c("CHR","BP","SNPID","REF","ALT","BETA","SE","P","n_eff", "ld.block")]
  # Dear snp_match require specific names, so let's abide
  names(ds)[c(1:2,4:7)]  <- c("chr","pos","a0","a1","beta", "beta_se")
  results <- data.table()
  
  for(chrs in 1:22){
    message("Working on chromosome ", chrs,".")
    # First thing is to check if we already have .rds files for our chr (ie. if we have read the bed files before). If not, we'll read it. This will create a .bk file, but it won't be a problem since we won't call this function again.
    if(!file.exists(paste0(reference,"chr",chrs,".rds"))){
      snp_readBed(paste0(reference,"chr",chrs,".bed"))
    }
    # Attach the "bigSNP" object in R session
    # This object contain SNP data from a single chromosome from 1000 genomes phase 1.
    # This object contains a matrix of genotypes (analogous to bed), as well as fam and map slots, PLINK format.
    obj.bigSNP <- snp_attach(paste(reference,"chr",chrs, ".rds", sep=""))
    # In this case, we'll focus only on individuals of european ancestry
    euridx  <- grep(ancestry, obj.bigSNP$fam$family.ID)
    # Get aliases for useful slots
    G   <- obj.bigSNP$genotypes
    # Filter sumstats and panel by the SNPs that we're going to use
    ds.chr <- as.data.frame(ds[ds$chr == chrs,])
    map.chr <- obj.bigSNP$map[-3]
    names(map.chr) <- c("chr", "SNPID", "pos", "a0", "a1")
    map.chr$SNPID <- paste(map.chr$chr,map.chr$pos, sep=":") # We're using CHR:BP to match, rather than rsIDs
    
    map.chr <- map.chr[map.chr$SNPID %in% ds.chr$SNPID,]	
    message("Matching and aligning SNPs in chr",chrs," to the reference")
    # Align sumstats to panel
    snp.chr <- snp_match(ds.chr, map.chr)
    
    pb <- txtProgressBar(min = 0, max = length(unique(snp.chr$ld.block)), style = 3) # Show a nice progress bar
    
    for(block in unique(snp.chr$ld.block)){
      
      idxbar <- which(unique(snp.chr$ld.block) == block)
      
      snp.block <- snp.chr[snp.chr$ld.block == block,]
      # Skip block if min P is above a certain threshold
      if(min(snp.block$P) > alpha.block){
        setTxtProgressBar(pb,idxbar)
        next
      }
      # Remove SNPs with P above a certain threshold
      snp.block <- snp.block[snp.block$P < alpha.snp,]
      
      # Skip block if has only one SNP after filtering
      if(nrow(snp.block) < 2){
        #				message ("\nWarning: Block ", block," has only one SNP. Skipping...")
        setTxtProgressBar(pb,idxbar)
        next
      }
      
      # Recover SNP indices	
      snp.idx <- which(paste(obj.bigSNP$map$chromosome,obj.bigSNP$map$physical.pos, sep=":") %in% snp.block$SNPID) 
      # Remove duplicates, which sometimes appear
      if(length(snp.idx) != length(snp.block$SNPID)){
        du <- paste(obj.bigSNP$map$chromosome,obj.bigSNP$map$physical.pos, sep=":") 
        dup <- du[du %in% snp.block$SNPID]
        dup <- dup[duplicated(dup)]
        snp.block <- snp.block[!snp.block$SNPID %in% dup,]
        snp.idx <- which(du %in% snp.block$SNPID)
      }
      
      # Compute LD matrix for that block
      LDmat.block <- snp_cor(G, ind.col = snp.idx, ind.row= euridx, ncores = ncores, size = length(snp.idx))
      dimnames(LDmat.block) <- list(snp.block$SNPID,snp.block$SNPID)
      LDmat.block <- as.matrix(LDmat.block)	
      
      snp.block <- snp.block[,c("chr","pos", "a0","a1","SNPID","beta","beta_se")]
      names(snp.block) <- c("CHR","BP","REF","ALT","SNPID","BETA","SE")
      # Compute Susie
      snp.block$ppi_susie <- susie_pip(snp.block,LDmat.block,nref=length(euridx), pi_i=pi_i, max_it = iterations, sd.prior = sd.prior)
      # Append to results
      results <- rbind(results, snp.block)
      
      # Progress bar
      
      setTxtProgressBar(pb, idxbar)
    }
    close(pb)
  }
  results[,weight:=BETA*ppi_susie]
  return(results)
  
}

##' Download 1000 Genomes Phase III panel  
##' 
##' \code{create_1000G} downloads and gets 1000 Genomes Phase III panel in 
##' PLINK format, and apply quality control for being used to compute PGS using 
##' \code{rapidopgs_multi}.
##' Given the size of the files, running this function can take long, depending
##' on broadband speed and server status. We also recommend to ensure that there
##' is at least 60GB free space available in disk.
##'
##' @param directory a string indicating the directory to download the panel
##' @param remove.related a logical stating if related individuals should be removed.
##' Default TRUE. 
##' @param qc.maf a numeric to set the MAF threshold for variants to be removed. DEFAULT 0.05
##' @param qc.hwe a numeric indicating the threshold for Hardy-Weinberg exact test 
##' p-value, below which variants will be removed. DEFAULT 1e-10.
##' @param qc.geno a numeric to set maximum missing call rates for variants. DEFAULT = 0.
##' @param autosomes.only If FALSE, it will include X and Y chromosomes, too.
##' @return bed, fam and bim files for each chromosome in the chosen directory.
##' @import data.table bigsnpr bigreadr
##' @export
##' @author Guillermo Reales
##' @examples
##' \dontrun{
##' create_1000G()
##'}

create_1000G <- function(directory = "ref-data", remove.related=TRUE, qc.maf = 0.05, qc.hwe=1e-10, qc.geno=0, autosomes.only=TRUE){
  
  # Remove annoying timeout limit in download.file
  timeout <- getOption('timeout')
  options(timeout=10000)
  
  # Max chromosomes to be considered, if autosomes.only is set to FALSE it will
  # be modified below
  max.chr=22
  
  dir.create(directory)
  plink <- download_plink(directory)
  message("Downloading 1000 Genomes Phase III files in VCF format.")
  message("These are (very) big files, so they'll take a while to download. Time for a cuppa.")
  for(chr in c(1:max.chr)){
    message("Downloading chromosome ",chr,"...")
    download.file(paste0("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"),
                  destfile = paste0(directory,"/chr",chr,".vcf.gz"), mode = "wb")
    message("Done!")
  }
  if(!autosomes.only){
  # X and Y chromosomes have a bit different link
  message("Downloading chromosome X. Note that current RápidoPGS-multi version uses autosomes only.")
  download.file("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz", destfile = paste0(directory,"/chr23.vcf.gz"), mode = "wb")
  message("Done!")
  message("Downloading chromosome Y. Note that current RápidoPGS-multi version uses autosomes only.")
  download.file("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz", destfile = paste0(directory,"/chr24.vcf.gz"), mode = "wb")
  message("Done!")
  max.chr=24
  }
  message("Transforming files into Plink format...")
  for(chr in 1:max.chr){
    system(paste0(plink, " --vcf ",directory, "/chr", chr, ".vcf.gz --make-bed --out ", directory, "/chr", chr), ignore.stdout = FALSE, ignore.stderr = FALSE)
  }
  # We don't need our hard-earned vcf files anymore, so we can delete them
  unlink(paste0(directory,"/*vcf.gz"))
  message("Done!")
  message("Starting QC procedure...")
  ped <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20200731.ALL.ped")
  pop <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv")
  
  for(chr in 1:max.chr){
    chrno <- ifelse(chr <23, chr, ifelse(chr==23, "X", "Y"))
    message("Processing Chr ", chrno, "...")
    bed <- snp_plinkQC(plink, prefix.in = paste0(directory,"/chr",chr),
                       prefix.out = paste0(directory,"/chr",chr,"_QC"),
                       geno = qc.geno, maf = qc.maf, hwe = qc.hwe)
    rds <- snp_readBed(bed)
    snp <- snp_attach(paste0(directory,"/chr",chr,"_QC.rds"))  
    G <- snp$genotypes
    fam <- snp$fam
    fam <- dplyr::left_join(fam, ped, by = c("sample.ID" = "Individual ID"))
    fam <- dplyr::left_join(fam, pop, by = c("Population" = "Population Code"))
    
    snp$fam$family.ID <- paste(fam$`Super Population`, fam$Population, sep = "_")
    snp$fam$paternal.ID <- snp$fam$maternal.ID <- 0L
    if(remove.related){
      unrelated_tags <- c("unrel","unrels","father","mother")
      ind_norel <- which(fam$Relationship %in% unrelated_tags & fam$Siblings == 0  & fam$`Second Order` == 0 & fam$`Third Order` == 0)
      
      maf <- snp_MAF(G, ind.row = ind_norel)
      
      bed <- snp_writeBed(snp, tempfile(fileext = ".bed"),
                          ind.row = ind_norel, ind.col = which(maf > 0.05))
    } else{
      
      maf <- snp_MAF(G)
      
      bed <- snp_writeBed(snp, paste0(directory, "/chr",chr,"_tmp.bed"),
                          ind.col = which(maf > qc.maf))
    }
    rds <- snp_readBed(bed)
    snp <- snp_attach(rds)
    unlink(paste0(directory, "/chr", chr,"\\.*"))
    snp_writeBed(snp, paste0(directory, "/chr",chr,".bed"))
    message("Done!")
  }
  unlink(paste0(directory, "/*_QC.*"))
  unlink(paste0(directory, "/*_tmp.*"))
  unlink(paste0(directory, "/plink"))
  message("Done! Now you can find your reference panel ready at ", directory,"/.")
  on.exit(options(timeout=timeout))
}


##' Compute Standard deviation prior (SD prior) for quantitative traits
##'  using pre-computed heritability.
##' 
##' \code{sd.prior.est} function will take the dataset as an input, a \eqn{h^2} 
##' value obtained from a public repository such as LDhub, 
##' (http://ldsc.broadinstitute.org/ldhub/), sample size and number of variants,
##' and will provide a sd.prior estimate that can be used to improve prediction
##' performance of RápidoPGS functions on quantitative traits.
##' 
##' @param data a data.table containing the GWAS summary statistic input dataset. 
##' Must contain SNPID and SE columns.
##' @param h2 a numeric. Heritability estimate or h^2 (See details).
##' @param pi_i a numeric. Prior that a given variant is causal. DEFAULT = 1e-4.
##' @export
##' @author Guillermo Reales, Elena Vigorito, Chris Wallace
##' @examples 
##' sumstats <- data.table(SNPID=c("4:1479959","20:13000913","14:29107209","2:203573414",
##' "4:57331393","6:11003529","6:149256398","21:25630085","13:79166661"), 
##'			REF=c("C","C","C","T","G","C","C","G","T"), 
##'			ALT=c("A","T","T","A","A","A","T","A","C"), 
##'			ALT_FREQ=c(0.2611,0.4482,0.0321,0.0538,0.574,0.0174,0.0084,0.0304,0.7528),
##'			BETA=c(0.012,0.0079,0.0224,0.0033,0.0153,0.058,0.0742,0.001,-0.0131),
##'			SE=c(0.0099,0.0066,0.0203,0.0171,0.0063,0.0255,0.043,0.0188,0.0074),
##'			P=c(0.2237,0.2316,0.2682,0.8477,0.01473,0.02298,0.08472,0.9573,0.07535)) 
##' sd.prior <- sd.prior.est(sumstats, h2 = 0.2456, N = 45658, pi_i=1e-4)
##' 
##' 
sd.prior.est <- function(data, h2, N, pi_i=1e-4){
  data <- copy(data)
  data <- unique(data, by="SNPID")
  if(is.character(N)){
    sample.size <- N
    var_b <- h2 / (pi_i * sum(1/(data[,get(sample.size)] * data$SE^2)))
  } else{
    var_b <- h2 / (pi_i * sum(1/(N * data$SE^2)))
  }
  return(sqrt(var_b))
  
}






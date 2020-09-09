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
##' var(X) = 2*maf*(1-maf)
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
##' @param se.beta vector of standard errors of effect sizes (\eqn{\beta})
##' @param sdY a scalar of the standard deviation given vectors of variance of coefficients,  MAF and sample size. Can be calculated using \code{sdY.est}
##' @param sd.prior a scalar representing our prior expectation of \eqn{\beta} (DEFAULT 0.15).
##' @param pi_i a scalar representing the prior probability (DEFAULT \eqn{1 \times 10^{-4}})
##' The method assumes a normal prior on the population log relative risk centred at 0 and the DEFAULT
##' value sets the variance of this distribution to 0.04, equivalent to a 95\%  belief that the true relative risk
##' is in the range of 0.66-1.5 at any causal variant.
##' @return a vector of posterior probabilities.
##' @author Guillermo Reales, Chris Wallace
wakefield_pp_quant <- function(beta, se.beta, sdY, sd.prior=0.15, pi_i=1e-4) { 
  # compute V
  V <- se.beta^2
  # Compute z too
  z <- beta/se.beta
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
##' This function is verbatim of its namesake in cupcake package (github.com/ollyburren/cupcake/)
##'
##' @param p a vector of univariate pvalues from a GWAS
##' @param f a vector of minor allele frequencies taken from some reference population.
##' @param N a scalar or vector for total sample size of GWAS
##' @param s a scalar representing the proportion of cases (n.cases/N)
##' @param pi_i a scalar representing the prior probability (DEFAULT \eqn{1 \times 10^{-4}})
##' @param sd.prior a scalar representing our prior expectation of \eqn{\beta} (DEFAULT 0.2).
##' The method assumes a normal prior on the population log relative risk centred at 0 and the DEFAULT
##' value sets the variance of this distribution to 0.04, equivalent to a 95\%  belief that the true relative risk
##' is in the range of 0.66-1.5 at any causal variant.
##' @param log.p if FALSE (DEFAULT), p is a p value. If TRUE, p is a log(p) value.  Use this if your dataset holds p values too small to be accurately stored without using logs
##' @return a vector of posterior probabilities.
##' @author Olly Burren, Chris Wallace

wakefield_pp <- function(p,f, N, s,pi_i=1e-4,sd.prior=0.2,log.p=FALSE) {
    if(length(p) != length(f))
      stop("p and f must be vectors of the same size")
    # compute V
    V <- 1 / (2 * N * f * (1 - f) * s * (1 - s))
    # convert p vals to z
    z <- stats::qnorm(0.5 * p, lower.tail = FALSE,log.p=log.p)
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


##' Compute PGS from GWAS summary statistics using posteriors from Wakefield's approximate Bayes Factors
##' 
##' '\code{computePGS} computes PGS from a from GWAS summary statistics using posteriors from Wakefield's approximate Bayes Factors
##' 
##' Main RápidoPGS function. This function will take a GWAS summary statistic dataset as an input,
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
##'   \item{SNPID}{rsids, or SNP identifiers. If not available, they can be anything (eg. CHR_BP)}
##'   \item{REF}{Reference, or non-effect allele}
##'   \item{ALT}{Alternative, or effect allele, the one \eqn{\beta} refers to}
##'   \item{ALT_FREQ}{Minor/ALT allele frequency in the tested population, or in a close population from a reference panel}
##'   \item{BETA}{\eqn{\beta} (or log(OR)), or effect sizes}
##'   \item{SE}{standard error of \eqn{\beta}}
##'   \item{P}{P-value for the association test}
##' }
##'
##' If a reference is provided. It should have 5 columns: CHR, BP,
##' SNPID, REF, and ALT. Also, it should be in the same build as 
##' the summary statistics. In both files, column order does not matter.
##' @param data a data.table containing GWAS summary statistic dataset
##'   with all required information.
##' @param N0 a scalar representing the number of controls in the
##'   study (or the number of subjects in quantitative trait GWAS),
##'   or a string indicating the column name containing it.
##' @param N1 a scalar representing the number of cases in the
##'   case-control study, or a string indicating the column name containing it. 
##'   If NULL (DEFAULT), quantitative trait will be assumed.
##' @param build a string containing the genome build of the dataset,
##'   either "hg19" (for hg19/GRCh37) or "hg38" (hg38/GRCh38). DEFAULT
##'   "hg19".
##' @param pi_i a scalar representing the prior probability (DEFAULT:
##'   \eqn{1 \times 10^{-4}}).
##' @param sd.prior the prior specifies that BETA at causal SNPs
##'   follows a centred normal distribution with standard deviation
##'   sd.prior. Sensible and widely used DEFAULTs are 0.2 for case
##'   control traits, and 0.15 * var(trait) for quantitative (selected
##'   if N1 is NULL).
##' @param log.p if FALSE (DEFAULT), p is a p value. If TRUE, p is a
##'   log(p) value. Use this if your dataset holds p values too small
##'   to be accurately stored without using logs.
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
##' @param forsAUC a logical indicating if output should be in sAUC
##'   evaluation format as we used it for the paper.
##' @param altformat a logical indicating if output should be in a
##'   format containing pid (chr:pos), ALT, and weights only. DEFAULT
##'   FALSE
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
##'			ALT_FREQ=c(0.2611,0.4482,0.0321,0.0538,0.574,0.0174,0.0084,0.0304,0.7528),
##'			BETA=c(0.012,0.0079,0.0224,0.0033,0.0153,0.058,0.0742,0.001,-0.0131),
##'			SE=c(0.0099,0.0066,0.0203,0.0171,0.0063,0.0255,0.043,0.0188,0.0074),
##'			P=c(0.2237,0.2316,0.2682,0.8477,0.01473,0.02298,0.08472,0.9573,0.07535))
##'
##' PGS  <- computePGS(sumstats,  N0= 119078 ,N1=137045, build = "hg38")
##'
	
computePGS <- function(data,
                       N0,
                       N1=NULL,
                       build = "hg19",
                       pi_i= 1e-04,
                       sd.prior=if(is.null(N1)) {0.15} else {0.2},
                       log.p=FALSE,
                       filt_threshold = NULL,
                       recalc=TRUE,
                       reference=NULL,
                       forsAUC=FALSE,
                       altformat=FALSE){
	
  ## Here's a list of columns that the dataset must have, otherwise it will fail
  mincol <- c("CHR","BP", "REF","ALT","BETA", "SE", "P", "ALT_FREQ", "SNPID")
  if(!all(mincol %in% names(data)))
    stop("All minimum columns should be present in the summary statistics dataset. Please check, missing columns: ",
         paste(setdiff(mincol,names(data)), collapse=", "))

  ds <- copy(data) # avoid modifying input data.table
  ds  <- na.omit(ds, cols=c("CHR","BP", "BETA", "SE", "P","ALT_FREQ")) # Remove NA in relevant columns
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
			blranges <- EUR_ld.blocks
		}else if(build == "hg38"){ 
			blranges <- EUR_ld.blocks38
		}else{ 
			stop("RapidoPGS only accepts hg19 or hg38 at the moment, please check.")
		}		
	message("Assigning LD blocks...")
	snpranges <- GRanges(seqnames=paste("chr",ds$CHR, sep=""), ranges=IRanges(start=ds$BP, end=ds$BP, names=ds$SNPID), strand="*")
        ds[,ld.block:=findOverlaps(snpranges, blranges, select='last')]
	message("Done!")
	}

  if(is.null(N1)){
    if(is.numeric(N0) && length(N0) == 1) { # In case N0 is supplied
      message("N1 not supplied. Assuming quantitative trait with ", N0, " individuals. Computing PGS.")
      ds[,sdY:=sdY.est(vbeta=SE^2, maf=ALT_FREQ, n=N0)][,ppi:=wakefield_pp_quant(BETA,SE,sdY,sd.prior), by="ld.block"][,weight:=ppi*BETA]
      if(!is.null(filt_threshold)){
        ds  <- ds[ds$ppi > filt_threshold,]
        if(recalc){
          ds[,sdY:=sdY.est(vbeta=SE^2, maf=ALT_FREQ, n=N0)][,ppi:=wakefield_pp_quant(BETA,SE,sdY,sd.prior), by="ld.block"][,weight:=ppi*BETA]		
        }
      }
    } else{ # In case column name is supplied
      if(is.character(N0) && length(N0) == 1){
        Nco <- N0
        message("N1 not supplied.  Assuming quantitative trait with possibly multiple N, supplied by column ", Nco,". Mmax N: ", max(ds[,get(Nco)]), ", min N = ", min(ds[,get(Nco)]), ", and mean N: ",mean(ds[,get(Nco)]), ". Computing PGS.")
        ds[,sdY:=sdY.est(vbeta=SE^2, maf=ALT_FREQ, n=get(Nco))][,ppi:=wakefield_pp_quant(BETA,SE,sdY,sd.prior), by="ld.block"][,weight:=ppi*BETA]
	if(!is.null(filt_threshold)){
          ds  <- ds[ds$ppi > filt_threshold,]
          if(recalc){
            ds[,sdY:=sdY.est(vbeta=SE^2, maf=ALT_FREQ, n=get(Nco))][,ppi:=wakefield_pp_quant(BETA,SE,sdY,sd.prior), by="ld.block"][,weight:=ppi*BETA]		
          }
	}
      }
    }
  } else  { # If both are supplied
    if(is.character(N0) && is.character(N1)){
      message("Computing PGS for a case-control dataset. Both N0 and N1 columns provided.")
      Nco  <- N0
      Nca <- N1
      ds <- na.omit(ds, cols=c(N0,N1))
      if(log.p){
        ds[,P:=pnorm(-abs(BETA/SE),log.p=TRUE)*2]	
      }
      ds[,ppi:=wakefield_pp(p = P, f = ALT_FREQ, N=get(Nco)+get(Nca), s = get(Nca)/(get(Nco)+get(Nca)), pi_i, sd.prior, log.p), by = "ld.block"][, weight:=ppi*BETA]	
      if(!is.null(filt_threshold)){
        if(filt_threshold < 1){
          ds  <- ds[ds$ppi > filt_threshold,]
        } else {
          ds <- ds[order(-rank(abs(weight))),][1:filt_threshold,] 
        }
        
        if(recalc){
          ds[,ppi:=wakefield_pp(p = P, f = ALT_FREQ, N= get(Nco)+get(Nca) , s = get(Nca)/(get(Nca)+get(Nca)), pi_i, sd.prior, log.p), by = "ld.block"][, weight:=ppi*BETA]
        }
      }
    }else{
      message("Computing PGS for a case-control dataset, with ", N0," controls, and ", N1, " cases.")
      if(log.p){
        ds[,P:=pnorm(-abs(BETA/SE),log.p=TRUE)*2]	
      }
      ds[,ppi:=wakefield_pp(p = P, f = ALT_FREQ, N= N0+N1, s = N1/(N0+N1), pi_i, sd.prior, log.p), by = "ld.block"][, weight:=ppi*BETA]	
      if(!is.null(filt_threshold)){
        if(filt_threshold < 1){
          ds  <- ds[ds$ppi > filt_threshold,]
        } else {
          ds <- ds[order(-rank(abs(weight))),][1:filt_threshold,] 
        }	
        if(recalc){
          ds[,ppi:=wakefield_pp(p = P, f = ALT_FREQ, N= N0+N1 , s = N1/(N0+N1), pi_i, sd.prior, log.p), by = "ld.block"][, weight:=ppi*BETA]
        }
      }
    }	
  }
  if(forsAUC){
    ds[,pid:=paste(CHR,BP,sep=":")]
    ds <- ds[,c("pid", "ALT", "weight")]
  }
  if(altformat){
    ds  <- ds[,c("SNPID", "ALT", "weight")]
  }
  return(ds)
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




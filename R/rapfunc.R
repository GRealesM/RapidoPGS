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

##' Finding a file in an FTP directory
##' This is an internal function to help \code{gwascat.download} find the right file
##' 
##'
##' @param ftp_address a string. An FTP address provided by \code{gwascat.download}.
##' @param acc a string containing the accession for the desired study.
##' @param hm a logical. Should it look in the harmonised directory?
##' @return a data.table containing the dataset.
##' @import data.table curl RCurl magrittr
##' @author Guillermo Reales
find_file_in_ftp <- function(ftp_address, acc, hm){
		
		flist <- getURL(ftp_address, dirlistonly = TRUE) %>% strsplit(., "\n")  %>% unlist 
		flist <- flist[grepl("GCST",flist)]

		accn <- gsub("GCST", "", acc)  %>% as.numeric
		ac1 <- floor(accn/1000)*1000 + 1

		sfl <- flist[grepl(paste0(ac1, "-"), flist)][1] # Choose the smallest one, if multiple

		tadd <- paste0(ftp_address, sfl, "/", acc, "/")	

		if(!RCurl::url.exists(tadd)) stop("There seem to be no summary statistics available for this accession in GWAS catalog. Please check.\n")
		tlist <- getURL(tadd, dirlistonly = TRUE) %>% strsplit(., "\n")  %>% unlist
	
		ziplist <- tlist[grepl(".csv$|.tsv$|.gz$", perl = TRUE ,tlist)]
		fileurls <- paste0(tadd, ziplist)

		if(hm && "harmonised" %in% tlist) {
			message("You have selected 'harmonised' dataset.")
			qadd <- paste0(ftp_address, sfl, "/", acc, "/harmonised/")
			qlist <- getURL(qadd, dirlistonly = TRUE) %>% strsplit(., "\n")  %>% unlist
			ziplist <- qlist[grepl(".csv$|.tsv$|.gz$",  perl = TRUE , qlist)]
			fileurls <- paste0(qadd, ziplist)
		}else if (!"harmonised" %in% tlist) {
		   message("There were no harmonised files available for this accession. Trying the not harmonised ones.")
		}

		return(fileurls)
}

##' Retrieve GWAS summary datasets from GWAS catalog
##' '\code{gwascat.download} takes a PMID from the user and downloads the associated summary statistics datasets published in GWAS catalog
##' 
##' This function, takes PUBMED ids as an input, searches at the GWAS catalog
##' for harmonised datasets associated to that, interactively asking the  
##' user to choose if there are more than one, and fetches the dataset. 
##'
##' If multiple files are available for the same study, R will prompt an interactive
##' dialogue to select a specific file, by number. 
##'
##' @param ID a numeric. A PubMed ID (PMID) reference number from a GWAS paper.
##' @param harmonised a logical. Should GWAS catalog harmonised files be pursued? 
##'     If not available, the function will fall back to non-harmonised
##' @return a character vector containing the url(s) to the dataset(s).
##' @import data.table curl RCurl
##' @export
##' @author Guillermo Reales
##'
gwascat.download <- function(ID, harmonised = TRUE){
		gwc.manifest <- fread("https://www.ebi.ac.uk/gwas/api/search/downloads/studies_alternative", quote = "")
		
		study.manifest <- gwc.manifest[ PUBMEDID == ID ,] 
		study_acc  <- 1 # For cases when the study has one accession only. Will be updated below otherwise
		
		if(nrow(study.manifest) == 0) stop("Please provide a valid PUBMED ID from a study published on GWAS catalog.")
		message("There are ",nrow(study.manifest)," datasets associated with this PUBMED ID (Note that summary statistics might not be available.)")
		if(nrow(study.manifest) > 1){	
			print(study.manifest[,c("STUDY ACCESSION", "DISEASE/TRAIT", "INITIAL SAMPLE SIZE")])
				study_acc  <- readline(prompt=paste("Please select accession (from 1 to ", nrow(study.manifest),"):", sep=""))
				study_acc  <- as.integer(study_acc)
				while(!study_acc %in% seq_len(nrow(study.manifest))){
					message("Oops! Seems that your number is not in the options. Please try again.")
					print(study.manifest[,c("STUDY ACCESSION", "DISEASE/TRAIT", "INITIAL SAMPLE SIZE")])
					study_acc  <- readline(prompt=paste("Please select accession (from 1 to ", nrow(study.manifest),"):", sep=""))
					study_acc  <- as.integer(study_acc)
					message("You selected option ", study_acc, ".")
					}

			}
		
		study.manifest <- study.manifest[study_acc,]
		# The filesystem changed in GWAS catalog, and now files are organised by the study accession.
		# We need to build the url accordingly

		ftp_address <- "ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
		access  <- study.manifest$`STUDY ACCESSION`
		message("Accession ", access, " selected.")

		ff <- find_file_in_ftp(ftp_address, acc = access, hm = harmonised)
		if(length(ff) == 0) stop("Oh! I couldn't retrieve any file from this accession. Are you sure the data is available at GWAS catalog?")
		
		furl <- ff # In case there's only one url

		if(length(ff) > 1){
			message("There are ", length(ff)," datasets associated with this Accession number.")
				fdf <- data.table(files = ff)
				print(fdf)
				ds_url  <- readline(prompt=paste0("\nNote that 'build37' means hg19/GRCh37, while 'h' means hg38.\nPlease select the file you want to download (from 1 to ", nrow(fdf),"):"))
				ds_url  <- as.integer(ds_url)
				while(!ds_url %in% seq_len(nrow(fdf))) {
					message("Oops! Seems that your number is not in the options. Please try again.")
					print(fdf)
					ds_url  <- readline(prompt=paste0("\nNote that 'build37' means hg19/GRCh37, while 'h' means hg38.\nPlease select the file you want to download (from 1 to ", nrow(fdf),"):"))
					ds_url  <- as.integer(ds_url)
					}
					furl <- ff[ds_url]
			}

		if(!RCurl::url.exists(furl)) stop("The file you requested is unavailable. This may be due to the fact that public and harmonised summary statistics do not exist. Please check at GWAS catalog website.")
		message("Retrieving dataset for ",study.manifest$`DISEASE/TRAIT`,", by ", study.manifest$`FIRST AUTHOR`,", from ",study.manifest$DATE, ", published at ", study.manifest$JOURNAL,", with accession ", study.manifest$`STUDY ACCESSION`," and url:", furl, ".\n")
		message("Note that some files indicate the genome build in the file name. For example 'GRCh37' or 'build37' refer to hg19/GRCh37, whereas 'h' and 'GRCh38' indicate hg38/GRCh38 build.")
		
		ds <- fread(furl)
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
##' Alternatively, if filt_threshold is larger than one, RapidoPGS will select the top
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
##'			CHR=c(4,20,14,2,4,6,6,21,13), 
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
    info_snp <- snp_match(ds,refset, match.min.prop=0)
    
    message("Original sumstat had ", nrow(sumstats.chr)," for chr",chr,". After matching ", nrow(info_snp.chr)," remained, and will be used for further steps.")
    ds <- data.table(info_snp)
    ds[,SNPID:=id][, c("_NUM_ID_.ss","_NUM_ID_", "id"):=NULL] # Replaces SNPID in the sumstat by SNPID in the reference, and removes snp_match cols.
    setnames(ds, old=c("chr","pos","beta","a0","a1"), new=c("CHR","BP", "BETA", "REF","ALT"))
  }
  
  # Assign ld.blocks, in case they werent computed yet
  if(!"ld.block" %in% names(ds)){
    if(build == "hg19"){ 
      blranges <- EUR_ld.blocks19
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
  ds <- ds[!is.na(ld.block),] # Remove SNPs without assigned block
  
  if(trait == "quant"){
    if(is.numeric(N) && length(N) == 1) { # In case N is supplied as a number
      message("Computing a RapidoPGS-single model for a quantitative trait with N = ", N, ".")
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
        message("Computing a RapidoPGS-single model for a quantitative trait with N supplied by column ", Nco,". Mmax N: ", max(ds[,get(Nco)]), ", min N = ", min(ds[,get(Nco)]), ", and mean N: ",mean(ds[,get(Nco)]), "...")
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
    message("Computing a RapidoPGS-single model for a case-control dataset...")
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

##' Compute PGS from GWAS summary statistics using Bayesian sum of single-effect 
##' (SuSiE) linear regression using z scores
##' 
##' '\code{rapidopgs_multi} computes PGS from a from GWAS summary statistics 
##' using Bayesian sum of single-effect (SuSiE) linear regression using z scores
##' 
##' This function will take a GWAS summary statistic dataset as an input,
##' will assign LD blocks to it, then use user-provided LD matrices or a preset 
##' reference panel in Plink format to compute LD matrices for each block. 
##' Then SuSiE method will be used to compute posterior probabilities of variants to be causal 
##' and generate PGS weights by multiplying those posteriors by effect sizes (\eqn{\beta}). 
##' Unlike \code{rapidopgs_single}, this approach will assume one or more causal variants.
##' 
##' 
##' The GWAS summary statistics file to compute PGS using our method must contain
##' the following minimum columns, with these exact column names:
##' \describe{
##'   \item{CHR}{Chromosome}
##'   \item{BP}{Base position (in GRCh37/hg19).}
##'   \item{REF}{Reference, or non-effect allele}
##'   \item{ALT}{Alternative, or effect allele, the one \eqn{\beta} refers to}
##'   \item{BETA}{\eqn{\beta} (or log(OR)), or effect sizes}
##'   \item{SE}{standard error of \eqn{\beta}}
##'   \item{P}{P-value for the association test}
##' }
##' In addition, quantitative traits must have the following extra column:
##' \describe{
##'   \item{ALT_FREQ}{Minor allele frequency.}
##' }
##' Also, for quantitative traits, sample size must be supplied, either as a number,
##' or indicating the column name, for per-SNP sample size datasets (see below).
##' Other columns are allowed, and will be ignored.
##' 
##' Reference panel should be divided by chromosome, in Plink format.
##' Both reference panel and summary statistic dataset should be in GRCh37/hg19.
##' For 1000 Genomes panel, you can use \code{create_1000G} function to set it up
##' automatically.
##' 
##' If prefer to use LD matrices, you must indicate the path to the directory 
##' where they are stored. They must be in RDS format, named LD_chrZ.rds (where
##' Z is the 1-22 chromosome number). If you don't have LD matrices already,
##' we recommend downloading those gently provided by Prive et al., at 
##' \url{https://figshare.com/articles/dataset/European_LD_reference/13034123}.
##' These matrices were computed using for 1,054,330 HapMap3 variants based on 
##' 362,320 European individuals of the UK biobank.
##'  
##' 
##' @param data a data.table containing GWAS summary statistic dataset
##'   with all required information.
##' @param reference a string representing the path to the directory containing 
##'   the reference panel (eg. "../ref-data").
##' @param LDmatrices a string representing the path to the directory containing 
##'   the pre-computed LD matrices.
##' @param N a numeric indicating the number of individuals used to generate input
##'  GWAS dataset, or a string indicating the column name containing per-SNP sample size.
##' @param build a string indicating the genome build. 'hg19' and 'hg38' are supported. Note that your LD matrices or reference panel should match the build.
##' @param trait a string indicating if trait is a case-control ("cc") or quantitative ("quant").
##' @param ncores a numeric specifying the number of cores (CPUs) to be used.
##'    If using pre-computed LD matrices, one core is enough for best performance.
##' @param alpha.block a numeric threshold for minimum P-value in LD blocks.
##'    Blocks with minimum P above \code{alpha.block} will be skipped. Default: 1e-4.
##' @param alpha.snp a numeric threshold for P-value pruning within LD block.
##'    SNPs with P above \code{alpha.snp} will be removed. Default: 0.01.
##' @param sd.prior the prior specifies that BETA at causal SNPs
##'   follows a centred normal distribution with standard deviation
##'   sd.prior.
##'   If NULL (default) it will be automatically estimated (recommended).
##' @param ancestry a string indicating the ancestral population (DEFAULT: "EUR", European). 
##'   If using an alternative population, bear in mind that your LD matrices or reference must be from the same population.
##'   You'll also need to provide matching LD.blocks via the LDblocks argument.
##' @param LDblocks a string indicating the path to an alternative LD block file in .RData format. Only required for non-European PGS.
##' @return a data.table containing the sumstats dataset with
##'   computed PGS weights.
##' @import data.table
##' @importFrom bigsnpr snp_match snp_cor snp_readBed snp_attach
##' @importFrom GenomicRanges GRanges findOverlaps
##' @importFrom IRanges IRanges 
##' @importFrom utils download.file setTxtProgressBar txtProgressBar
##' @importFrom coloc runsusie
##' @export
##' @author Guillermo Reales, Chris Wallace
##' @examples
##' \dontrun{
##' ss <- data.table(
##'			CHR=c(4,20,14,2,4,6,6,21,13), 
##'			BP=c(1479959, 13000913, 29107209, 203573414, 57331393, 11003529, 149256398, 
##'					25630085, 79166661), 
##'			REF=c("C","C","C","T","G","C","C","G","T"), 
##'			ALT=c("A","T","T","A","A","A","T","A","C"), 
##'			BETA=c(0.012,0.0079,0.0224,0.0033,0.0153,0.058,0.0742,0.001,-0.0131),
##'			SE=c(0.0099,0.0066,0.0203,0.0171,0.0063,0.0255,0.043,0.0188,0.0074),
##'			P=c(0.2237,0.2316,0.2682,0.8477,0.01473,0.02298,0.08472,0.9573,0.07535))
##' PGS <- rapidopgs_multi(ss, reference = "ref-data/", N = 20000, build = "hg19", trait="cc", ncores=5)
##'}
rapidopgs_multi <- function(data, reference=NULL, LDmatrices=NULL, N=NULL, build=c("hg19", "hg38"), trait=c("cc", "quant"), ncores=1, alpha.block=1e-4, alpha.snp=0.01, sd.prior=NULL, ancestry="EUR", LDblocks = NULL){
 
  # Sanity checks
  if(!"data.table" %in% class(data)) data <- as.data.table(data)
  if(is.null(N)) stop("N (sample size) is now required for all GWAS summary statistics. Please provide it.")
  if(!is.null(N) && length(N) != 1) stop("Please provide a single N value.")
  if((is.null(reference) && is.null(LDmatrices)) | (!is.null(reference) && !is.null(LDmatrices))) stop("Please provide either a reference panel or LD matrices.")
  if(ancestry != "EUR" && is.null(LDblocks)) stop("You selected non-European ancestry but provided no alternative LDblock file.\nAt the moment, RapidoPGS multi only provides LD blocks for Europeans.\nIf you want to create your own, check EUR_ld.blocks.RData in the package files for reference.")
  if(!is.null(reference) && ncores == 1) message("You selected 1 core. When using a reference you can increase the number of cores")
  if(!is.null(LDmatrices) && ncores > 1) message("Note: You selected more than 1 core, but RapidoPGS-multi does not use parallelisation when using LD matrices.")
  if(length(build) != 1 || !build %in% c("hg19", "hg38")) stop("Please provide a valid genomic build for your data: either 'hg19' or 'hg38' are supported.")
  if(length(trait) != 1 || !trait %in% c("cc", "quant")) stop("Please provide a valid specify your trait type: either case-control ('cc') or quantitative ('quant') are supported.") 
  
  mincol <- c("CHR","BP", "REF","ALT","BETA", "SE","P")
  
  if(!all(mincol %in% names(data)))
    stop("All minimum columns should be present in the summary statistics dataset. Please check, missing columns: ",
         paste(setdiff(mincol,names(data)), collapse=", "))
  
  if(!is.null(reference)){ # If reference is provided, make safety check
    if(!dir.exists(reference)){
      stop("Reference directory doesn't exists. Please select a valid directory. Need to create a reference panel first? Check `create_1000G()`.")
    }
    if(!file.exists(paste(reference,"/chr1.bed", sep=""))) # Temp fix to detect panel		
      stop("No reference panel detected. Please check.")
  }
    
  if(!is.null(LDmatrices)){
    if(!dir.exists(LDmatrices)){
      stop("LDmatices directory doesn't exists. Please select a valid directory. Need to download the LDmatrices?")
    }
    if(!all(paste0("LD_chr", 1:22, ".rds") %in% dir(path=LDmatrices))){
      stop("Not all matrices LD_chr...rds are present in the directory. Please check.")
    }
  }
  
  message("You selected ", build, " build. Please note that your reference or LDmatrices must be in the same build as your dataset.")

  ds <- copy(data) # avoid modifying input data.table
  ds[,SNPID:=paste(CHR,BP, sep = ":")]
  ds  <- na.omit(ds, cols=c("CHR","BP", "BETA", "SE", "P"))
  
  # Compute LD blocks
  # Assign ld.blocks, in case they werent computed yet.
  # Note that in this case only hg19 is admitted.
  if(!"ld.block" %in% names(ds)){

    if(build == "hg19" && ancestry == "EUR"){
        blranges <- EUR_ld.blocks19	
    } else if (build == "hg38" && ancestry == "EUR") {
        blranges <- EUR_ld.blocks38	
    }

    # Special case for non-European LD blocks
    if(ancestry != "EUR"){
      if(!file.exists(LDblocks) || !grepl(".rda$", LDblocks, ignore.case = TRUE)) stop("Please provide a valid LDblock file as ldblock.rda")
      load(LDblocks)
      blranges <- ldblock
      message("You selected non-European ancestry. Plase ensure that either your reference panel or your LD matrices correspond to this ancestry.")
    } 
  
    ###########
    message("Assigning LD blocks using pre-computed ", build, " LD blocks.")
    snpranges <- GRanges(seqnames=paste("chr",ds$CHR, sep=""), ranges=IRanges(start=ds$BP, end=ds$BP, names=ds$SNPID), strand="*")
    ds[,ld.block:=findOverlaps(snpranges, blranges, select='last')]
    message("Done!")
  }

  ## We ensure that rows without assigned block are removed
  if(any(is.na(ds$ld.block))){
    # warning(length(ds$ld.block[is.na(ds$ld.block)]), " SNPs didn't have any LD block assigned and will be removed")
    ds <- na.omit(ds, cols="ld.block")
  }
  
  # Dear snp_match require specific names, so let's abide
  setnames(ds, c("CHR", "BP","REF", "ALT", "BETA","SE"), c("chr","pos","a0","a1","beta", "beta_se"))

  #######
  
  results <- data.table()
  
  message("Running RapidoPGS-multi model with multiple causal variant assumption for a GWAS summary statistics dataset in ", build,", of ", ancestry," ancestry, and with N = ", N, ".")

  if(!is.null(LDmatrices)){
      # In case LDmatrices are provided
      # Import HapMap3 manifest
      map <- as.data.table(readRDS(paste(LDmatrices, "map.rds", sep="/")))
      map[,SNPID:=paste(chr,pos, sep=":")] 
      # Check if ds is aligned to map, and align if not
      ds.snp <- paste(ds$chr, ds$pos, ds$a0, ds$a1, sep=":")
      map.snp <- paste(map$chr, map$pos, map$a0, map$a1, sep=":") 

      if(!all(ds.snp %in% map.snp)){
        ds <- as.data.table(snp_match(ds, map, match.min.prop = 0))
      }
  }
    
  ### Here the loop by chromosome starts
  for(chrs in 1:22){

    if(!is.null(reference)){
      message("Working on chromosome ", chrs,".")
    # First thing is to check if we already have .rds files for our chr (ie. if we have read the bed files before). If not, we'll read it. This will create a .bk file, but it won't be a problem since we won't call this function again.
      if(!file.exists(paste0(reference,"/chr",chrs,".rds"))){
        snp_readBed(paste0(reference,"/chr",chrs,".bed"))
      }
      # Attach the "bigSNP" object in R session
      # This object contain SNP data from a single chromosome from 1000 genomes phase III.
      # This object contains a matrix of genotypes (analogous to bed), as well as fam and map slots, PLINK format.
      obj.bigSNP <- snp_attach(paste(reference,"/chr",chrs, ".rds", sep=""))
      # By default EUR, but this can be changed, just ensure you supply the right ancestry coding to the function
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
      ds.chr <- snp_match(ds.chr, map.chr, match.min.prop=0)
      ds.chr <- data.table(ds.chr)
      
      pb <- txtProgressBar(min = 0, max = length(unique(ds.chr$ld.block)), style = 3) # Show a nice progress bar

    } else if (!is.null(LDmatrices)) {

      message("Working on chromosome ", chrs,".")
      LD.chr <- readRDS(paste0(LDmatrices, "/LD_chr",chrs, ".rds"))
      ds.chr  <-  ds[chr == chrs]
      map.chr <- map[chr == chrs]
    
      # Check LDmatrix has same dimensions as manifest
      if(nrow(LD.chr) != nrow(map.chr)){
      stop("Something is wrong. LD matrix doesn't have the same rows as the manifest for this chromosome, please check.")
      }   
      pb <- txtProgressBar(min = 0, max = length(unique(ds.chr$ld.block)), style = 3) # Show a nice progress bar
    }

    # Here the loop by block starts
    for(block in unique(ds.chr$ld.block)){
        #  message("This is block ",block)
      idxbar <- which(unique(ds.chr$ld.block) == block)
      snp.block <- ds.chr[ds.chr$ld.block == block]
      
      # Skip block if min P is above a certain threshold
      if(min(snp.block$P) > alpha.block){
        #				message ("\nBlock ", block," has no SNP with P-value below ",alpha.block," threshold. Skipping block.")
        setTxtProgressBar(pb,idxbar)
        next
      }
      # Remove SNPs with P above a certain threshold
      snp.block <- snp.block[snp.block$P < alpha.snp]
      
      # Skip block if has only one SNP after filtering
      if(nrow(snp.block) < 2){
        #				message ("\nWarning: Block ", block," has only one SNP. Skipping...")
        setTxtProgressBar(pb,idxbar)
        next
      }
      
      if(!is.null(reference)){

            # Recover SNP indices. To avoid duplicates, match by alleles too
            map.idx <- paste(obj.bigSNP$map$chromosome,obj.bigSNP$map$physical.pos, obj.bigSNP$map$allele1, obj.bigSNP$map$allele2, sep=":")
            block.idx <- paste(snp.block$chr, snp.block$pos, snp.block$a0, snp.block$a1, sep=":")
            snp.idx <- which(map.idx %in% block.idx) 

            # Compute LD matrix for that block
            LD.block <- snp_cor(G, ind.col = snp.idx, ind.row= euridx, ncores = ncores, size = length(snp.idx))
            dimnames(LD.block) <- list(snp.block$SNPID,snp.block$SNPID)
            LD.block <- as.matrix(LD.block)	

            snp.block <- copy(snp.block)
            setnames(snp.block ,c("chr","pos", "a0","a1","SNPID","beta","beta_se") ,c("CHR","BP","REF","ALT","SNPID","BETA","SE"))

      } else if (!is.null(LDmatrices)) {

            snp.idx <- match(snp.block$SNPID, map.chr$SNPID) 
            LD.block <- as.matrix(LD.chr[snp.idx,snp.idx])
            dimnames(LD.block) <- list(snp.block$SNPID, snp.block$SNPID)
            snp.block <- snp.block[,c("chr","pos", "a0","a1","SNPID","beta","beta_se")]
            snp.block <- copy(snp.block)
            # Prepare dataset for runsusie
            setnames(snp.block ,c("chr","pos", "a0","a1","SNPID","beta","beta_se") ,c("CHR","BP","REF","ALT","SNPID","BETA","SE"))

      }

      # Prepare snp.block for susie
      # Match ids of resulting SNPs with those in map manifest (and hence, LD matrix)
      # Remove duplicates
      if(all(!duplicated(snp.block$SNPID))){
        snp.block <- snp.block[!duplicated(snp.block$SNPID)]
      }
      
      # susie.ds <- list(snp=snp.block$SNPID, beta=snp.block$BETA, varbeta=snp.block$SE^2, LD=LD.block, type=trait)
      # susie.ds <- list(snp=snp.block$SNPID, beta=snp.block$BETA, varbeta=snp.block$SE^2, LD=LD.block, N = as.numeric(N), type = trait)
      # check_dataset requires sdY if type == "quant", even though susie doesn't need it. We added sdY = 1.
      if (trait == "cc") {
        susie.ds <- list(snp = snp.block$SNPID, beta = snp.block$BETA, varbeta = snp.block$SE^2, LD = LD.block, N = as.numeric(N), type = trait)
      } else if (trait == "quant") {
        susie.ds <- list(snp = snp.block$SNPID, beta = snp.block$BETA, varbeta = snp.block$SE^2, LD = LD.block, N = as.numeric(N), type = trait, sdY = 1)
      }

      if(is.null(sd.prior)){
        prior_est <-  TRUE
        prior_var <-  50 # Susie default.
      } else{
        prior_est  <-  FALSE
        prior_var  <-  sd.prior^2
      }
      
      ppi_susie <- suppressMessages(runsusie(susie.ds, prior_variance=prior_var, estimate_prior_variance=prior_est))
      snp.block$ppi_susie <- ppi_susie$pip
      
      # Append to results
      results <- rbind(results, snp.block)
      
      # Progress bar
      setTxtProgressBar(pb, idxbar)
      
    } # End of block loop
    close(pb)
  } # End of loop 
  results[,weight:=BETA*ppi_susie]
  return(results)
}



##' Download 1000 Genomes Phase III panel  
##' 
##' \code{create_1000G} downloads and gets 1000 Genomes Phase III panel (hg19) in 
##' PLINK format, and apply quality control for being used to compute PGS using 
##' \code{rapidopgs_multi}.
##' Given the size of the files, running this function can take long, depending
##' on broadband speed and server status. We also recommend to ensure that there
##' is at least 60GB free space available in disk.
##'
##' @param directory a string indicating the directory to download the panel
##' @param remove.related a logical stating if related individuals should be removed.
##' Default TRUE. 
##' @param qc.maf a numeric to set the MAF threshold for variants to be removed. DEFAULT 0.01
##' @param qc.hwe a numeric indicating the threshold for Hardy-Weinberg exact test 
##' p-value, below which variants will be removed. DEFAULT 1e-10.
##' @param qc.geno a numeric to set maximum missing call rates for variants. DEFAULT = 0.
##' @param autosomes.only If FALSE, it will include X and Y chromosomes, too.
##' @return bed, fam and bim files for each chromosome in the chosen directory.
##' @import data.table bigsnpr bigreadr 
##' @importFrom utils download.file
##' @export
##' @author Guillermo Reales
##' @examples
##' \dontrun{
##' create_1000G()
##'}
create_1000G <- function(directory = "ref-data", remove.related=TRUE, qc.maf = 0.01, qc.hwe=1e-10, qc.geno=0, autosomes.only=TRUE){
  
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
    download.file(paste0("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"),
                  destfile = paste0(directory,"/chr",chr,".vcf.gz"), mode = "wb")
    message("Done!")
  }
  if(!autosomes.only){
  # X and Y chromosomes have a bit different link
  message("Downloading chromosome X. Note that current RapidoPGS-multi version uses autosomes only.")
  download.file("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz", destfile = paste0(directory,"/chr23.vcf.gz"), mode = "wb")
  message("Done!")
  message("Downloading chromosome Y. Note that current RapidoPGS-multi version uses autosomes only.")
  download.file("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz", destfile = paste0(directory,"/chr24.vcf.gz"), mode = "wb")
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
##' performance of RapidoPGS functions on quantitative traits.
##' 
##' @param data a data.table containing the GWAS summary statistic input dataset. 
##' Must contain SNPID and SE columns.
##' @param h2 a numeric. Heritability estimate or h^2 (See details).
##' @param pi_i a numeric. Prior that a given variant is causal. DEFAULT = 1e-4.
##' @param N a numeric. Sample size of the GWAS input dataset.
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






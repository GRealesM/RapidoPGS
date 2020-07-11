##' Complement genotypes
##'
##' ie A/G -> T/C
##' This function was copycat from annotSnpStats package (github.com/chr1swallace/annotSnpStats/), by Chris Wallace
##'
##' @title g.complement
##' @param x character vector of genotypes
##' @export
##' @return character vector of genotypes on the alternative strand
##' @examples
##' g.complement(c("A/G","A/T"))
##' @author Chris Wallace
g.complement <- function (x) {
  x <- toupper(x)
  switches <- c(A = "t", T = "a", C = "g", G = "c")
  for (i in seq_along(switches)) x <- sub(names(switches)[i], switches[i], x)
  toupper(x)
}


##' Reverse alleles in a genotype
##'
##' ie A/G -> G/A
##' This function is verbatim of its namesake from annotSnpStats package (github.com/chr1swallace/annotSnpStats/), by Chris Wallace
##'
##' @title g.rev
##' @param x character vector of genotypes
##' @param sep character with which to separate alleles. Default is "/".
##' @export
##' @return character vector of reversed genotypes 
##' @examples
##' g.rev(c("A/G","A/T"))
##' @author Chris Wallace
g.rev <- function (x, sep = "/") {
  sapply(strsplit(x, sep), function(g) paste(rev(g), collapse = "/"))
}


##' Define possible allele switching classes
##'
##' This function is verbatim of its namesake in annotSnpStats package (github.com/chr1swallace/annotSnpStats/), by Chris Wallace
##'
##' @title g.class
##' @param x vector of allele codes from dataset X
##' @param y vector of allele codes from dataset Y, same length as x
##' @return character vector of allele switching classes
##' @export
##' @examples
##' alleles.X <- c(snp1="A/G",snp2="A/G",snp3="A/G",snp4="A/G",snp5="A/T",snp6="A/T")
##' alleles.Y <- c(snp1="A/G",snp2="G/A",snp3="T/C",snp4="C/T",snp5="A/T",snp6="T/A")
##' classes <- g.class(x=alleles.X,y=alleles.Y)
##' cbind(alleles.X,alleles.Y,classes)
##' @author Chris Wallace
g.class <- function (x, y) {
  if (!identical(names(x), names(y))) 
    stop("x and y must relate to same SNPs")
  mat <- matrix(FALSE, length(x), 4, dimnames = list(names(x), c("nochange", "rev", "comp", "revcomp")))
  mat[, "nochange"] <- x == y
  mat[, "rev"] <- x == g.rev(y)
  mat[, "comp"] <- x == g.complement(y)
  mat[, "revcomp"] <- x == g.rev(g.complement(y))
  indels <- x %in% c("I/D", "D/I")
  if (any(indels)) 
    mat[indels, c("comp", "revcomp")] <- FALSE
  ret <- character(nrow(mat))
  rs <- rowSums(mat)
  if (length(wh <- which(rs > 1))) 
    ret[wh] <- "ambig"
  if (length(wh <- which(rs == 0))) 
    ret[wh] <- "impossible"
  if (length(wh <- which(rs == 1))) 
    ret[wh] <- colnames(mat)[apply(mat[wh, , drop = FALSE], 1, which)]
  return(ret)
}

##' Helper function to sum logs without loss of precision
##'
##' Sums logs without loss of precision
##' This function is verbatim of its namesake in cupcake package (github.com/ollyburren/cupcake/)
##'
##' @param x a vector of logs to sum
##' @return a scalar
##' @export
##' @author Olly Burren, Chris Wallace
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
##' @export
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
##' @export
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
##' @return a vector of posterior probabilities.
##' @export
##' @author Olly Burren, Chris Wallace

wakefield_pp <- function(p,f, N, s,pi_i=1e-4,sd.prior=0.2) {
    if(length(p) != length(f))
      stop("p and f must be vectors of the same size")
    # compute V
    V <- 1 / (2 * N * f * (1 - f) * s * (1 - s))
    # convert p vals to z
    z <- stats::qnorm(0.5 * p, lower.tail = FALSE)
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


##' Preprocess GWAS summary statistics datasets for PGS computation
##' 
##' \code{pgs.file.preprocess} imports the dataset, checks for minimum information, assign LD blocks to SNPs from an externally calculated LD architecture,
##' and optionally filters SNPs and aligns alleles to the reference HapMap3 panel.
##' 
##' The GWAS summary statistics file to compute PGS using our method must contain the following minimum columns, with these exact column names:
##' CHR19: Chromosome 
##' BP19: Base position (in GRCh37/hg19)
##' REF: Reference, or non-effect allele
##' ALT: Alternative, or effect allele, the one \eqn{\beta} refers to
##' BETA: Beta (or log(OR)), or effect sizes.
##' SE: Standard error of \eqn{\beta}
##' P: P-values for the association test
##' ALT_FREQ: Minor/ALT allele frequency in the tested population, or in a close population from a reference panel. No problem it flipped.
##'
##' If a ref is provided. It should have 4 columns: CHR19, BP19, SNPID, REF, and ALT
##' @param dataset a string with the name or path to the GWAS summary statistic dataset, or the path to it
##' @param blockfile a string with the name or path to a RDS file containing the genomic coordinates of LD blocks, computed for EUR populations.
##' @param ref a string indicating the reference file SNPs should be filtered and aligned to.
##' @return a data.table to be used for PGS computation.
##' @export
##' @author Guillermo Reales, Chris Wallace

pgs.file.preprocess  <- function(dataset, blockfile="ld.blocks.RDS", ref=NULL){

	message("Loading dataset...")
	ds <- fread(dataset)
	# Here's a list of columns that the dataset must have, otherwise it will fail
	mincol <- c("CHR19","BP19", "REF","ALT","BETA", "SE", "P", "ALT_FREQ")
	if(!all(mincol %in% names(ds))) stop("All minimum columns should be present. Please check, missing columns: ", setdiff(mincol,names(ds)))
	
	if(!is.null(ref)){
	# Next step is to filter and align our alleles and their effects to the hapmap3 reference, which I have already formatted for our purposes.
	message("Filtering SNPs...")
	refset <- fread(ref)
	refset[, alleles:=paste(REF,ALT, sep="/")][,pid:=paste(CHR19, BP19, sep=":")]
	ds[,alleles:=paste(REF,ALT, sep="/")][,pid:=paste(CHR19, BP19, sep=":")]
	ds <- ds[pid %in% refset$pid,]
	ds  <- merge(ds, refset[,.(SNPID, pid, alleles)], by ='pid', suffixes=c("", ".reference"))
	ds[, alleles:=toupper(alleles)][, c("REF", "ALT"):=list(toupper(REF), toupper(ALT))][, SNPID:=SNPID.reference]

	message("Aligning alleles...")
	# Check if everything is alright  
	if(!all(g.class(ds$alleles.reference, ds$alleles)== "nochange")){
	    allele_diagnostics <- g.class(ds$alleles.reference, ds$alleles)
	    alleles_to_flip <-  allele_diagnostics == "rev"
	    alleles_to_comp <- allele_diagnostics == "comp"
	    alleles_to_revcomp <- allele_diagnostics == "revcomp"
	    cat("Some SNPs have to be flipped. ", sum(alleles_to_flip), " to flip, ", sum(alleles_to_comp), " to find their complement, and ", sum(alleles_to_revcomp), " to find their reverse complement.\n")
	    ds$alleles[alleles_to_flip] <- unlist(g.rev(ds$alleles[alleles_to_flip]))
	    ds$alleles[alleles_to_comp] <- g.complement(ds$alleles[alleles_to_comp])
	    ds$alleles[alleles_to_revcomp] <- unlist(g.rev(g.complement(ds$alleles[alleles_to_revcomp])))
	    ds$REF <- sapply(strsplit(ds$alleles, "/"), `[`, 1)
	    ds$ALT <- sapply(strsplit(ds$alleles, "/"), `[`, 2)
	    ds$BETA[alleles_to_flip] <- ds$BETA[alleles_to_flip]*-1
	    ds$BETA[alleles_to_revcomp] <- ds$BETA[alleles_to_revcomp]*-1
	   # NOTE: I introduced the following bit from milcytokine_basis on to guarantee that we have no ambiguity nor duplicated SNPs
	#if(!all(g.class(ds$alleles.reference, ds$alleles)== "nochange")){
	#	ds  <- ds[g.class(ds$alleles.reference, ds$alleles)== "nochange",]
	#	}
	
	    rm(alleles_to_flip, alleles_to_comp, alleles_to_revcomp)
	  }
	ds[, c("alleles", "alleles.reference"):=NULL]
	}

	ds <- unique(ds)

	# I made ld.blocks file from fourier_ls-all.bed, so now I can simply load the object
	# blocks <- fread("fourier_ls-all.bed")
	# blranges <- GRanges(seqnames=blocks$chr, ranges=IRanges(blocks$start, end=blocks$stop, names=1:nrow(blocks)), strand="*") 
	# saveRDS(blranges, "ld.blocks.RDS")
	message("Assigning LD blocks...")
	blranges <- readRDS(blockfile)
	ds  <- ds[!is.na(CHR19) & !is.na(BP19),] 
	snpranges <- GRanges(seqnames=paste("chr",ds$CHR19, sep=""), ranges=IRanges(ds$BP19, end=ds$BP19, names=ds$SNPID), strand="*")
        ds[,ld.block:=findOverlaps(snpranges, blranges, select='last')]
	ds  <- na.omit(ds)
	message("Done!")
	return(ds)
}




##' Compute PGS from GWAS summary statistics using posteriors from Wakefield's approximate Bayes Factors
##' '\code{computePGS} computes PGS from a from GWAS summary statistics using posteriors from Wakefield's approximate Bayes Factors
##' 
##' This function will take a GWAS summary statistic dataset with SNPs assigned to LD blocks and compute Wakefield's ppi by LD block, then will use it 
##' to generate PGS weights by multiplying those posteriors by effect sizes (\eqn{\beta}). 
##' Optionally, it will filter SNPs by a custom filter on ppi and then recalculate weights, to improve accuracy.
##'
##' @param ds a data.table containing GWAS summary statistic dataset with all information, including ld.blocks. Usually an output from \code{pgs.file.preprocess}
##' @param N0 a scalar representing the number of controls in the study (or the number of subjects in quantitative trait GWAS)
##' @param N1 a scalar representing the number of cases in the case-control study. If NULL (default), quantitative trait will be assumed
##' @param pi_i a scalar representing the prior probability (DEFAULT \eqn{1 \times 10^{-4}})
##' @param sd.prior a scalar representing our prior expectation of \eqn{\beta} (DEFAULT 0.2)
##' @param filt_threshold a scalar indicating the ppi threshold to filter the dataset after PGS computation. If NULL (Default), nothresholding will be applied
##' @param recalc a logical indicating if weights should be recalculated after thresholding. If TRUE, \code{filt_threshold} should be defined
##' @param forsAUC a logical indicating if output should be in sAUC evaluation format as we used it for the paper.
##' @param altformat a logical indicating if output should be in a format containing pid (chr:pos), ALT, and weights only. Default FALSE
##' @return a data.table containing the formatted sumstat dataset with computed PGS weights.
##' @export
##' @author Guillermo Reales, Chris Wallace

computePGS <- function(ds, N0,N1=NULL,pi_i= 1e-04, sd.prior=0.2, filt_threshold = NULL, recalc=FALSE, forsAUC=FALSE, altformat=FALSE){

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
		message("N1 not supplied.  Assuming quantitative trait with multiple N, supplied by column ", Nco,". Mmax N: ", max(ds[,get(Nco)]), ", min N = ", min(ds[,get(Nco)]), ", and mean N: ",mean(ds[,get(Nco)]), ". Computing PGS.")
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
			ds[,ppi:=wakefield_pp(p = P, f = ALT_FREQ, N=get(Nco)+get(Nca), s = get(Nca)/(get(Nco)+get(Nca)), pi_i, sd.prior), by = "ld.block"][, weight:=ppi*BETA]	
				if(!is.null(filt_threshold)){
					ds  <- ds[ds$ppi > filt_threshold,]
				if(recalc){
					ds[,ppi:=wakefield_pp(p = P, f = ALT_FREQ, N= get(Nco)+get(Nca) , s = get(Nca)/(get(Nca)+get(Nca)), pi_i = 1e-04, sd.prior=0.2), by = "ld.block"][, weight:=ppi*BETA]
				}
				}
		}else{
			message("Computing PGS for a case-control dataset, with ", N0," controls, and ", N1, " cases.")
			ds[,ppi:=wakefield_pp(p = P, f = ALT_FREQ, N= N0+N1, s = N1/(N0+N1), pi_i, sd.prior), by = "ld.block"][, weight:=ppi*BETA]	
				if(!is.null(filt_threshold)){
					ds  <- ds[ds$ppi > filt_threshold,]
				if(recalc){
					ds[,ppi:=wakefield_pp(p = P, f = ALT_FREQ, N= N0+N1 , s = N1/(N0+N1), pi_i = 1e-04, sd.prior=0.2), by = "ld.block"][, weight:=ppi*BETA]
				}
				}
		  }	
		}
	if(forsAUC){
		ds[,pid:=paste(CHR19,BP19,sep=":")]
		ds <- ds[,c("pid", "ALT", "weight")]
	}
	if(altformat){
		ds  <- ds[,c("SNPID", "ALT", "weight")]
	}

	return(ds)
}	


##' Get SNP correlation summation
##' '\code{get.correlation.adj} computes a summation of the genotypic correlation between PGS model and a reference panel
##' 
##' This function, taken from SummaryAUC (https://github.com/lsncibb/SummaryAUC)
##' source code, will use a reference panel from 1000 Genomes Project (379 individuals 
##' of European ancestry) to calculate the genotypic correlation between SNPs in the PGS
##' model and the panel, outputting a summation that will be later use at the denominator
##' for delta calculation in \code{auc} computation. 
##' For more details, please see original publication (Song et al., 2019; https://academic.oup.com/bioinformatics/article/35/20/4038/5419857). 
##'
##' @param prs.model a data.table containing the PGS model, ideally an output from \code{computePGS}
##' @param tau a numeric vector representing the genotypic variance in the GWAS validation dataset (sqrt(2 * MAF * (1 - MAF) * INFO)
##' @param KG.plink.pre a string representing the prefix of the reference plink files, provided by SummaryAUC
##' @param soFile a string of the path to the compiled .so file for computation, provided by SummaryAUC
##' @param pos_thr a scalar representing the threshold for distance of SNPs, correlation will be caculated for any SNP pairs within this distance
##' @return a scalar. A summation of beta * tau * cor for SNPs
##' @export

get.correlation.adj <- function(prs.model, tau, KG.plink.pre = 'KG.all.chr', soFile = 'scripts/getAdjCorrelation.so', pos_thr = 5e8){
	kg.fam = read.table(file = paste0(KG.plink.pre, '.fam'), header = F, stringsAsFactors = F)
	NSample = nrow(kg.fam)
	
	kg.bim = read.table(file = paste0(KG.plink.pre, '.bim'), header = F, stringsAsFactors = F)
	# Let's create pids!
	kg.bim$pid  = paste(kg.bim$V1, kg.bim$V4, sep = ":")
	pids = intersect(prs.model$pid, kg.bim$pid)
	Npids = length(pids)
	idx = match(pids, prs.model$pid)
	prs.model = prs.model[idx, ]
	tau = tau[idx]
	
	idx = match(prs.model$pid, kg.bim$pid)
	chrs = kg.bim[idx,1]
	pos = kg.bim[idx, 4]
	flip.idx = which(prs.model$ALT != kg.bim[idx, 'V5']) # Here V5 is bim A1 (usually minor allele)
	idx[flip.idx] = -idx[flip.idx]
	
	beta_tau = prs.model$weight * tau
	kg.bed.file = paste0(KG.plink.pre, '.bed')
	
	if(!file.exists(soFile)){
		cFile = gsub('.so$', '.c', soFile)
		system(paste('R CMD SHLIB', cFile))
	}
	dyn.load(soFile)
	
	adj.cor.results = .C("getAdjCorrelation", plinkBed = kg.bed.file, NumSample = as.integer(NSample), NumSNP = as.integer(Npids), idx = as.integer(idx), chrs = as.integer(chrs), pos = as.integer(pos), pos_thr = as.integer(pos_thr), beta_tau = as.double(beta_tau), adj_cor = as.double(0.1))
	
	return(adj.cor.results$adj_cor)
}


##' Compute AUC using a GWAS summary statistic dataset for validation
##' '\code{get.correlation.adj} computes a summation of the genotypic correlation between PGS model and a reference panel
##' 
##' This function, adapted from SummaryAUC (https://github.com/lsncibb/SummaryAUC)
##' source code, will use a PGS model of a binary trait and a GWAS summary statistics dataset for validation. 
##' Ideally, PGS model will have a specific format, containing 3 columns: pid (chr:pos in hg19 build), ALT, and weights.
##' GWAS validation dataset will contain: pid, ALT, MAF, BETA, and P. 
##' 
##' For more details on the method, please see original publication (Song et al., 2019; https://academic.oup.com/bioinformatics/article/35/20/4038/5419857). 
##'
##' @param prs.model.file a data.table containing the PGS model, ideally an output from \code{computePGS}
##' @param gwas.summary.stats.file a vector representing the name of the GWAS sumstats validation dataset
##' @param N0 a scalar. Number of controls in the GWAS validation dataset
##' @param N1 a scalar. Number of cases in the GWAS validation dataset
##' @param soFile a string of the path to the compiled .so file for computation, provided by SummaryAUC
##' @param pos_thr a scalar representing the threshold for distance of SNPs, correlation will be caculated for any SNP pairs within this distance
##' @param flag.correction.ajd.inputted.data a logical. Are there inputted SNPs in the GWAS validation dataset? (Not tested in our implementation)
##' @return a list containing AUC, AUC variance, adj.cor, and number of SNPs in common in the PRS model and the validaton dataset
##' @export


auc <- function(prs.model.file, gwas.summary.stats.file, N0,N1, soFile = 'getAdjCorrelation.so', KG.plink.pre = 'KG.all.chr', pos_thr = 5e8, flag.correlation.adj.inputted.data = FALSE){
	#browser()
	#library("mvtnorm") ## As included in the package, no need to load it for this function
	#work.dir = dirname(prs.model.file)
	#setwd(work.dir)
	
	prs.model = read.delim(file = prs.model.file, header = T, sep = '\t', stringsAsFactors = F)
	
	gwas.summary.stats = read.table(file = gwas.summary.stats.file, header = T, stringsAsFactors = F, check.names = F, fill=TRUE)
	# CHR	SNP	A1	MAF	BETA	P	INFO
	# 14	rs7160549	C	0.472562	0.206895	2.36804e-07	0.997
	pids = intersect(prs.model[, 'pid'], gwas.summary.stats[, 'pid'])

	idx = match(pids, gwas.summary.stats[, 'pid'])
	gwas.summary.stats = gwas.summary.stats[idx, ]
	idx = match(pids, prs.model[, 'pid'])
	prs.model = prs.model[idx, ]

	if("OR" %in% colnames(gwas.summary.stats)){
		beta_gwas = log(gwas.summary.stats[, 'OR'])
	}else{
		beta_gwas = gwas.summary.stats[, 'BETA']
	}
	idx = which(gwas.summary.stats[, 'ALT'] != prs.model[, 'ALT'])  
	beta_gwas[idx] = -beta_gwas[idx]
	P = gwas.summary.stats[, 'P']
	MAF = gwas.summary.stats[, 'MAF']
	Z = sign(beta_gwas) * qnorm(P/2, lower.tail=FALSE)

	if("INFO" %in% colnames(gwas.summary.stats)){
		INFO = as.numeric(gwas.summary.stats[, 'INFO'])
	}else{
		INFO = rep(1, length(MAF))
	}
	tau = sqrt(2 * MAF * (1 - MAF) * INFO)
	beta_prs = prs.model[,'weight']
	nSNPs = nrow(prs.model)	
	
	if(flag.correlation.adj.inputted.data){
		
		if(!file.exists(soFile)){
			cFile = gsub('.so$', '.c', soFile)
			system(paste('R CMD SHLIB', cFile))
		}
		dyn.load(soFile)
		
		adj.cor = 0
		for(chr in 1:22){
			cur.chr.dosage.file = paste0("dosage.chr.", chr, ".txt.gz")
			if(!file.exists(cur.chr.dosage.file))
				next
			
			# n.cur.chr.snp = as.integer(system(paste0("zcat dosage.chr.", chr, ".txt.gz | wc -l"), intern = TRUE))
			# if(n.cur.chr.snp < 1)
			# 	next
			
			cat("Chromosome ", chr, ' Reading dosage data:\n', sep = "")
			dosage = read.delim(file = paste0("dosage.chr.", chr, ".txt.gz"), header = F, stringsAsFactors = F, sep = '\t')
			# 1     rs61769350      693731  A       G       1.944 ......
			# https://www.cog-genomics.org/plink/1.9/assoc#dosage
			# 'format=1' normally indicates a single 0..2 A1 expected count, here A1 is the 4th column with skip0=1 skip1=1 format=1
			
			cur.snps = intersect(snps, dosage[, 2])
			if(length(cur.snps) < 1)
				next
			
			idx = match(cur.snps, dosage[, 2])
			dosage = dosage[idx, ]
			dosage.chr = dosage[, 1]
			dosage.pos = dosage[, 3]
			dosage.datamatrix = t(as.matrix(dosage[, -(1:5), drop = F]))
			storage.mode(dosage.datamatrix) = 'double'
			
			idx = dosage[, 4] != prs.model[match(cur.snps, snps), 2]
			dosage.datamatrix[, idx] = 2 - dosage.datamatrix[, idx]
			
			cur.tau = rep(NA, ncol(dosage.datamatrix))
			for(i in 1:ncol(dosage.datamatrix))
				cur.tau[i] = sqrt(var(dosage.datamatrix[, i], na.rm = TRUE))
			
			idx = match(cur.snps, snps)
			cur.beta_prs = beta_prs[idx]
			cur.beta_tau = cur.beta_prs * cur.tau
			tau[idx] = cur.tau
			
			if(ncol(dosage.datamatrix) < 2){
				rm(dosage, dosage.datamatrix, dosage.chr, dosage.pos, cur.tau, cur.beta_tau)
				next
			}
			
			idx = is.na(dosage.datamatrix)
			dosage.datamatrix[idx] = -9

			adj.cor.results = .C("getAdjCorrelationDosage", D = as.double(dosage.datamatrix), NumSample = as.integer(nrow(dosage.datamatrix)), NumSNP = as.integer(ncol(dosage.datamatrix)), chrs = as.integer(dosage.chr), pos = as.integer(dosage.pos), pos_thr = as.integer(pos_thr), beta_tau = as.double(cur.beta_tau), adj_cor = as.double(0.1))
			adj.cor = adj.cor + adj.cor.results$adj_cor
			
			rm(dosage, dosage.datamatrix, dosage.chr, dosage.pos, cur.tau, cur.beta_tau, adj.cor.results)
		}
	} else{
		
		adj.cor = get.correlation.adj(prs.model, tau, KG.plink.pre, soFile, pos_thr)
		# cat(paste0('adj.cor:', adj.cor, '\n'))
		# delta = sqrt(1/N1 + 1/N0) * sum(beta_prs * Z * tau) / sqrt(2 * sum((beta_prs * tau)^2) + 4 * adj.cor)
		# cat(paste0('delta:', delta, '\n'))
	}

		# cat(paste0('adj.cor:', adj.cor, '\n'))

		delta = sqrt(1/N1 + 1/N0) * sum(beta_prs * Z * tau) / sqrt(2 * sum((beta_prs * tau)^2) + 4 * adj.cor)
		# cat(paste0('delta:', delta, '\n'))

		auc0 = pnorm(delta)
		pxy = pmvnorm(lower=c(-delta,-delta),upper=Inf,mean=c(0,0),sigma=matrix(c(1,0.5,0.5,1),2,2))[1]
		var0 = (pxy - auc0^2) * (N1+N0) / (N1*N0)
		return(list(AUC=auc0,AUC.Var=var0, Adj.Cor=adj.cor, SNPs=length(pids)))
}


##' Retrieve GWAS summary datasets from GWAS catalog
##' '\code{gwascat.download} takes a PMID from the user and downloads the associated summary statistics datasets published in GWAS catalog
##' 
##' This function, takes PUBMED ids as an input, searches at the GWAS catalog
##' for harmonised datasets associated to that, interactively asking the  
##' user to choose if there are more than one, and fetches the dataset. 
##'
##' @param ID a numeric. A PubMed ID (PMID) reference number from a GWAS paper.
##' @param hm_only a logical. Should GWAS catalog harmonised columns be retained?
##' @return a data.table containing the dataset, if available.
##' @examples
##' test1 <- gwascat.download(27863252) # Astle, should work
##' negcontrol <- gwascat.download(01223247236) # The Empress Pub phone number. Shouldn't work!
##' @export

gwascat.download <- function(ID, hm_only=TRUE){
		gwc.manifest <- fread("https://www.ebi.ac.uk/gwas/api/search/downloads/studies_alternative")
		
		
		study.manifest <- gwc.manifest[ID == PUBMEDID,] 
		
		if(nrow(study.manifest) == 0) stop("Please provide a valid PUBMED ID")
		if(nrow(study.manifest) > 1){
			message("There are multiple datasets associated to this PUBMED ID")
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
		url <- paste("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/", gsub(" ","", study.manifest$`FIRST AUTHOR`), "_", ID,"_",study.manifest$`STUDY ACCESSION`, "/harmonised/",ID,"-",study.manifest$`STUDY ACCESSION`,"-",sapply(strsplit(study.manifest$MAPPED_TRAIT_URI, "/"), `[`,5),".h.tsv.gz", sep="")
		if(!RCurl::url.exists(url)) stop("The file you requested is unavailable. This may be due to the fact that public and harmonised summary statistics do not exist. Please check at GWAS catalog website.")
		message("Retrieving dataset for ",study.manifest$`DISEASE/TRAIT`,", by ", study.manifest$`FIRST AUTHOR`,", from ",study.manifest$DATE, ", published at ", study.manifest$JOURNAL,", with accession ", study.manifest$`STUDY ACCESSION`,".")
		
		ds <- fread(url)
		if(hm_only){
		hmcols <- grep("hm_",names(ds), value=TRUE)
		ds  <- ds[,..hmcols]
		}
}






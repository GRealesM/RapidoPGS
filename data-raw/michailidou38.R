## code to prepare `michailidou38` dataset goes here


library(data.table)
set.seed(1)
options(timeout=10000)

ds <- fread("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004988/harmonised/29059683-GCST004988-EFO_0000305.h.tsv.gz")

# If running interactively, you can also run the next line, selecting option 2 -- hg38
# ds <- gwascat.download(29059683)
setnames(ds, old = c("hm_rsid","hm_chrom","hm_pos", "hm_other_allele", "hm_effect_allele", "hm_effect_allele_frequency", "hm_beta", "standard_error", "p_value"), new = c("SNPID","CHR", "BP", "REF","ALT","ALT_FREQ", "BETA", "SE", "P"))
ds <- ds[,.(SNPID, CHR, BP, REF, ALT, ALT_FREQ, BETA, SE, P)]
ds <- ds[CHR !="X"]
ds$CHR <- as.numeric(ds$CHR)
ds <- ds[order(CHR, BP)]
ds <- na.omit(ds, cols = c("BETA", "ALT_FREQ"))
tokeep <- sort(sample(1:nrow(ds), 100000, replace = FALSE))
michailidou38 <- ds[tokeep]
michailidou38[, N:= 137045 + 119078] # Cases and controls for Michailidou, required for RápidoPGS_multi
usethis::use_data(michailidou38, overwrite = TRUE)

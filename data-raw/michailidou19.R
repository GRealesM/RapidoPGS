## code to prepare `michailidou19` dataset goes here

library(data.table)
set.seed(1)
options(timeout=10000)

ds <- fread("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004988/harmonised/29059683-GCST004988-EFO_0000305-build37.f.tsv.gz")

# If running interactively, you can also run the next line, selecting option 1 -- hg19
# ds <- gwascat.download(29059683)
setnames(ds, old = c("variant_id","chromosome","base_pair_location", "other_allele", "effect_allele", "effect_allele_frequency", "beta", "standard_error", "p_value"), new = c("SNPID","CHR", "BP", "REF","ALT","ALT_FREQ", "BETA", "SE", "P"))
ds <- ds[,.(SNPID, CHR, BP, REF, ALT, ALT_FREQ, BETA, SE, P)]
ds <- ds[CHR !="X"]
ds$CHR <- as.numeric(ds$CHR)
ds <- ds[order(CHR, BP)]
ds <- na.omit(ds, cols = c("BETA", "ALT_FREQ"))
tokeep <- sort(sample(1:nrow(ds), 100000, replace = FALSE))
michailidou19 <- ds[tokeep]
michailidou19[, N:= 137045 + 119078] # Cases and controls for Michailidou, required for RápidoPGS_multi
usethis::use_data(michailidou19, overwrite = TRUE)

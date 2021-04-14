## code to prepare `michailidou19` dataset goes here

setDTthreads(8)
set.seed(1)
options(timeout=10000)
ds <- fread("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/MichailidouK_29059683_GCST004988/harmonised/29059683-GCST004988-EFO_0000305-build37.f.tsv.gz")
setnames(ds, old = c("variant_id","chromosome","base_pair_location", "other_allele", "effect_allele", "effect_allele_frequency", "beta", "standard_error", "p_value"), new = c("SNPID","CHR", "BP", "REF","ALT","ALT_FREQ", "BETA", "SE", "P"))
ds <- ds[,.(SNPID, CHR, BP, REF, ALT, ALT_FREQ, BETA, SE, P)]
ds <- ds[CHR !="X"]
ds$CHR <- as.numeric(ds$CHR)
tokeep <- sort(sample(1:nrow(ds), 100000))
michailidou19 <- ds[tokeep]
usethis::use_data(michailidou19, overwrite = TRUE)

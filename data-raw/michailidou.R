## code to prepare `michailidou` dataset goes here


library(data.table)
set.seed(1)
ds <- gwascat.download(29059683, hm_only = FALSE)
ds <- ds[,c(2:7,11,20:21)]
ds <- na.omit(ds, cols = c("hm_beta", "hm_effect_allele_frequency"))
snpselect <- sample(1:nrow(ds), 100000, replace = FALSE)
michailidou <- ds[snpselect,]
usethis::use_data(michailidou, overwrite = TRUE)

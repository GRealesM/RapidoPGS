# Code to prepare MacDonald's hg38 LD block data

library(data.table)
library(GenomicRanges)
bl <- fread("https://raw.githubusercontent.com/jmacdon/LDblocks_GRCh38/master/data/EUR_LD_blocks.bed")
EUR_ld.blocks38 <- GRanges(seqnames=bl$chr, ranges=IRanges(start=bl$start, end=bl$stop), strand="*")
usethis::use_data(EUR_ld.blocks38, overwrite = TRUE)
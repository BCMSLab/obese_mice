# load required libraries
library(tidyverse)
library(DESeq2)

# load data
ob_counts <- read_rds('data/ob_counts.rds')
info <- as_tibble(colData(ob_counts))
se2 <- ob_counts[rowRanges(ob_counts)$biotype == 'protein_coding',]

# remove batch effects

# create DESeq2 object from gene-level counts and metadata
dds <- DESeqDataSetFromMatrix(countData = assay(se2),
                              colData = info,
                              design = ~ 1)

# Estimate library size correction scaling factors
dds <- estimateSizeFactors(dds)

# write data
write_rds(dds, 'data/ob_counts_batched.rds')

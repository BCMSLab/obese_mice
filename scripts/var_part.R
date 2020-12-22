# load required libraries
library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)
library(cowplot)
library(variancePartition)

# load data
dds <- read_rds('data/ob_counts_batched.rds')
info <- as_tibble(colData(dds))

# identify genes that pass expression cutoff
isexpr <- rowSums(fpm(dds)>1) >= 0.5*ncol(dds)

# compute log2 Fragments Per Million
dds_transform <- vst(dds)[isexpr,]

# Define formula
form <- ~ (1|group) + (1|diet) + (1|tissue) + (1|mouse) + (1|batch)

# Run variancePartition analysis
varPart <- fitExtractVarPartModel(assay(dds_transform), form, info)
write_rds(varPart, 'data/variancePartition.rds')

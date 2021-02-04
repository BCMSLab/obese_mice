# load required libraries
library(tidyverse)
library(SummarizedExperiment)
library(variancePartition)

# load data
ob_counts <- read_rds('data/ob_counts_batched.rds')
mat <- assay(ob_counts, 'transformed')
info <- as.data.frame(colData(ob_counts))

# remove low counts genes
mat2 <- mat[round(rowVars(mat), 4) > 0,]
dim(mat)
dim(mat2)

# Define formula
form <- ~ (1|group) + (1|diet) + (1|tissue) + (1|mouse)

# Run variancePartition analysis
varPart <- fitExtractVarPartModel(mat2, form, info)

# write object
write_rds(varPart, 'data/variance_partitioned.rds')

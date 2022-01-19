# load required libraries
library(tidyverse)
library(DESeq2)
library(sva)
library(limma)

# load data
ob_counts <- read_rds('data/ob_counts.rds')
info <- as_tibble(colData(ob_counts))
se2 <- ob_counts[rowRanges(ob_counts)$biotype == 'protein_coding',]

# remove batch effects

# create DESeq2 object from gene-level counts and metadata
dds <- DESeqDataSetFromMatrix(countData = assay(se2),
                              colData = info,
                              design = ~ 1)
dds <- estimateSizeFactors(dds)

# estimate library size correction scaling factors
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]

mod  <- model.matrix(~ group + diet+ tissue, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))

#n = num.sv(dat,mod,method="leek")

svseq <- svaseq(dat, mod, mod0, n.sv = 2)
sv <- as.data.frame(svseq$sv) %>%
    setNames(paste0('SV', 1:ncol(svseq$sv)))

# add variables to dds object
colData(dds) <- cbind(colData(dds), sv)

# transformed counts
assay(dds, 'transformed') <- assay(vst(dds))
assay(dds, 'removed_batch') <- removeBatchEffect(assay(dds, 'transformed'), dds$batch)
assay(dds, 'removed_sv') <- removeBatchEffect(assay(dds, 'transformed'), dds$SV1)
assay(dds, 'removed_tissue') <- removeBatchEffect(assay(dds, 'removed_batch'), dds$tissue)

# write data
write_rds(dds, 'data/ob_counts_batched.rds')

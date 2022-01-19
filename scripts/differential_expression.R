# load libraries
library(tidyverse)
library(DESeq2)
library(openxlsx)

# load and format data
ob_counts <- read_rds('data/ob_counts.rds')

# extract unique tissue names
tissues <- unique(ob_counts$tissue)

# apply differential expression to each tissue separately

tissue_res <- map(tissues, function(x) {
    # subset by tissue
    ob_sub <- ob_counts[,ob_counts$tissue == x]
    
    # relevel factors
    ob_sub$diet <- factor(ob_sub$diet, levels = c("ND", "HFD"))
    ob_sub$group <- factor(ob_sub$group, levels = c("WT", "ob/ob"))
    
    # create a model matrix
    mod <- model.matrix(~ group * diet, colData(ob_sub))
    
    # make deseq dataset
    dds <- DESeqDataSet(ob_sub, mod)
    
    # remove low counts
    keep <- rowSums(counts(dds)) >= 50
    dds <- dds[keep,]
    
    # run deseq
    dds <- DESeq(dds)
    
    # extract results for groupob.ob.dietHFD
    res_names <- resultsNames(dds)
    res <- results(dds, name = res_names[4], tidy = TRUE)
    res <- select(res, geneID = row, everything()) %>% as_tibble()
    
    # return
    res
})

names(tissue_res) <- tissues

openxlsx::write.xlsx(tissue_res, 'data/differential_expression.xlsx')

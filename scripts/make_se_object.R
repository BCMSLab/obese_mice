# load required libraries
library(tidyverse)
library(SummarizedExperiment)
library(GenomicFeatures)
library(biomaRt)

# load data: raw data
raw_data <- list.files('data/1st RNAanalysis/1.Rawdata_information',
                       full.names = TRUE) %>%
    map_df(read_tsv)

# load data: alignment rate
align_rate <- list.files('data/1st RNAanalysis/2.ReadAlignmentRate/',
                         full.names = TRUE) %>%
    map_df(read_tsv)

# load data: counts
counts_summary <- read_tsv('data/1st RNAanalysis/3.FeatureCount_result/MouseRNA-seq.txt.summary') %>%
    gather(file, value, -Status)

counts <- read_tsv('data/1st RNAanalysis/3.FeatureCount_result/MouseRNA-seq.txt',
                   skip = 1)

# make se object
# count matrix
mat <- as.matrix(dplyr::select(counts, ends_with('.bam')))
rownames(mat) <- counts$Geneid

# pheno type data
pd <- tibble(id = colnames(mat)) %>%
    mutate(id = str_replace(id, '_ob', ''),
           id = str_replace(id, 'HDF', 'HFD'),
           group = str_split(id, '_', simplify = TRUE)[, 1],
           diet = str_split(id, '_', simplify = TRUE)[, 2],
           replicate = str_split(id, '_', simplify = TRUE)[, 3],
           tissue = str_split(id, '_', simplify = TRUE)[, 4],
           batch = str_split(id, '_|\\.', simplify = TRUE)[, 5],
           mouse = paste(group, diet, replicate, sep = '_')) %>%
    as.data.frame()

rownames(pd) <- colnames(mat)

# gene info
ens_mart <- useMart('ensembl')
ens_mm10 <- useMart(biomart = 'ensembl',
                    dataset = 'mmusculus_gene_ensembl')

ens_df <- getBM(
    attributes = c('ensembl_gene_id',
                   'chromosome_name',
                   'start_position',
                   'end_position',
                   'strand',
                   'entrezgene_id',
                   'mgi_symbol',
                   'gene_biotype'),
    mart = ens_mm10
)

# subset to available ids
mat2 <- mat[intersect(rownames(mat), ens_df$ensembl_gene_id),]

# make gr object
gr <- ens_df[match(rownames(mat2), ens_df$ensembl_gene_id),]

gr$strand <- ifelse(gr$strand == 1, '+', '-')
names(gr) <- c('gene_id', 'seqnames', 'start', 'end', 'strand',
               'entrez_id', 'symbol', 'biotype')

row_ranges <- makeGRangesFromDataFrame(gr,
                               keep.extra.columns = TRUE)

names(row_ranges) <- row_ranges$gene_id

all(rownames(mat2) == names(row_ranges))

# make object
se <- SummarizedExperiment(assays = list(counts = mat2),
                           colData = pd,
                           rowRanges = row_ranges)

# write object
write_rds(se, 'data/ob_counts.rds')

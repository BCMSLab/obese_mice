# load libraries
library(tidyverse)
library(reshape2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GO.db)

# load and format data
ob_counts <- read_rds('data/ob_counts.rds')

var_part <- read_rds('data/variance_partitioned.rds')
var_part <- var_part@.Data %>%
    bind_cols() %>%
    setNames(var_part@names) %>%
    mutate(gene_id = var_part@row.names) %>%
    as_tibble() %>%
    left_join(as_tibble(mcols(ob_counts))) %>%
    filter(!is.na(symbol))

# prepare input
vars <- c('tissue', 'group', 'diet')
stats <- map(vars, function(x) {
    vec <- pull(var_part, x)
    names(vec) <- var_part$symbol
    vec[order(vec, decreasing = TRUE)][1:1000]
}) %>%
    set_names(vars)
str(stats)

pathways <- select(org.Mm.eg.db,
                   na.omit(unique(var_part$symbol)),
                   'GO',
                   'SYMBOL')
pathways <- split(pathways$SYMBOL, pathways$GO)
#pathways <- pathways[lengths(pathways) > 10 & lengths(pathways) < 1000]
pathways <- as_tibble(melt(pathways)) %>%
    setNames(c('gene', 'term')) %>%
    dplyr::select(term, gene)

pathway_names <- tibble(
    term = unique(pathways$term),
    name = map(term, ~GOTERM[[.x]]@Term)
)

ov <- map_df(stats, function(x) {
    enricher(names(x),
             TERM2GENE = pathways,
             TERM2NAME = pathway_names,
             pAdjustMethod = 'fdr')@result
}, .id = 'variable') %>%
    as_tibble()

write_rds(ov, 'data/over_representation.rds')

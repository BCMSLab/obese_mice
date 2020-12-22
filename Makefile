#!/bin/bash

all: directories \
	data/ob_counts.rds \
	data/ob_counts_batched.rds \
	data/variance_partitioned.rds \
	data/over_representation.rds \
	manuscript.pdf
	
directories:
	mkdir -p manuscript manuscript/tables manuscript/figures
	mkdir -p log
	
data/ob_counts.rds: scripts/make_se_object.R 
	R CMD BATCH --vanilla $< log/ob_counts.txt
	
data/ob_counts_batched.rds: scripts/remove_batch_effects.R \
	data/ob_counts.rds 
	R CMD BATCH --vanilla $< log/ob_counts_batched.txt
	
data/variance_partitioned.rds: scripts/var_part.R \
	data/ob_counts_batched.rds
	R CMD BATCH --vanilla $< log/variance_partitioned.txt
	
data/over_representation.rds: scripts/over_representation.R \
	data/ob_counts.rds \
	data/variancePartition.rds
	R CMD BATCH --vanilla $< log/over_representation.txt
	
manuscript.pdf: manuscript.Rmd \
	$(bash find data/*) \
	$(bash find scripts/*)
	Rscript -e 'rmarkdown::render("manuscript.Rmd")' &> log/manuscript.txt

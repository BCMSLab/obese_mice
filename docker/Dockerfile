FROM rocker/verse

RUN sudo apt-get update -y
RUN sudo apt-get install -y git-all

RUN Rscript -e "update.packages(ask = FALSE)"
RUN Rscript -e "install.packages('lme4')"
RUN Rscript -e "install.packages('cowplot')"
RUN Rscript -e "install.packages('BiocManager')"

RUN Rscript -e "BiocManager::install('GenomicFeatures')"
RUN Rscript -e "BiocManager::install('SummarizedExperiment')"
RUN Rscript -e "BiocManager::install('biomaRt')"
RUN Rscript -e "BiocManager::install('ComplexHeatmap')"
RUN Rscript -e "BiocManager::install('circlize')"
RUN Rscript -e "BiocManager::install('variancePartition')"
RUN Rscript -e "BiocManager::install('clusterProfiler')"
RUN Rscript -e "BiocManager::install('fgsea')"
RUN Rscript -e "BiocManager::install('DESeq2')"
RUN Rscript -e "BiocManager::install('org.Mm.eg.db')"
RUN Rscript -e "BiocManager::install('GO.db')"

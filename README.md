![Build and push Docker image](https://github.com/MahShaaban/obese_mice/workflows/Build%20and%20push%20Docker%20image/badge.svg)

# obese_mice

Code for the data descriptor paper (Tissue-specific gene expression in obese
hyperglycemic mice)

## Setting up the docker environment

The analysis was run on a [docker](https://hub.docker.com/r/bcmslab/obsesmice/)
image based on the the latest **rocker/verse**. Other R packages were added to
the image and were made available as an image that can be obtained and launched
on any local machine running
[docker](https://hub.docker.com/r/bcmslab/obsesmice/).

```bash
docker pull bcmslab/obsesmice:latest
docker run -it bcmslab/obsesmice:latest bash
```

Alternatively, the image can be built locally form the provided Docerfile.
Pass `docker/Dockerfile` to the `docker build` command. 

## Obtaining the source code

The source code is hosted publicly on this repository in the form of a research
compendium. This includes the scripts to reproduce the figures and tables in 
this manuscript.

```bash
git clone https://github.com/MahShaaban/obese_mice
```

## Runing the analysis

In the directory `obese_mice`, run `make`

```bash
cd curated_adipo_descriptor
make all
```

This file contains the recipes to generate several `rds` objects in `data` and
tables and figures in `manuscript/`.

FROM bioconductor/bioconductor_docker:RELEASE_3_21

LABEL maintainer="Jimmy Lucas"
LABEL description="Docker image for Nextflow scRNA-seq pipeline with R/Bioconductor/Seurat"

ENV DEBIAN_FRONTEND=noninteractive
ENV R_LIBS_USER=/usr/local/lib/R/site-library
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    curl \
    wget \
    git \
    ca-certificates \
    procps \
    make \
    cmake \
    build-essential \
    pkg-config \
    pandoc \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libpng-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libxt-dev \
    libx11-dev \
    libglu1-mesa-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses-dev \
    libfftw3-dev \
    libtiff-dev \
    libjpeg-dev \
    libhdf5-dev \
    libglpk-dev \
    libgmp3-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages(c( \
    'optparse', \
    'yaml', \
    'tidyverse', \
    'data.table', \
    'ggplot2', \
    'pheatmap', \
    'patchwork', \
    'cowplot', \
    'Matrix', \
    'hdf5r', \
    'Seurat', \
    'SeuratObject', \
    'remotes', \
    'bookdown', \
    'rmarkdown', \
    'knitr' \
    ), repos='https://cloud.r-project.org')"

RUN Rscript -e "BiocManager::install(c( \
    'SingleCellExperiment', \
    'SummarizedExperiment', \
    'scater', \
    'scran', \
    'DropletUtils', \
    'BiocParallel', \
    'DESeq2', \
    'edgeR', \
    'limma', \
    'ComplexHeatmap', \
    'clusterProfiler', \
    'enrichplot', \
    'DOSE', \
    'GO.db', \
    'org.Hs.eg.db', \
    'org.Mm.eg.db', \
    'AnnotationDbi', \
    'AnnotationHub', \
    'biomaRt' \
    ), ask=FALSE, update=FALSE)"

RUN Rscript -e "options(timeout = 600); remotes::install_git('https://github.com/chris-mcginnis-ucsf/DoubletFinder.git', dependencies = TRUE, upgrade = 'never')"

RUN Rscript -e "library(optparse); library(yaml); library(Seurat); library(SingleCellExperiment); library(DoubletFinder); sessionInfo()"

WORKDIR /workspace

CMD ["/bin/bash"]
FROM continuumio/miniconda3

LABEL maintainer="Jimmy Lucas" \
      description="Conda environment for the Nextflow scRNA-seq pipeline (R, Bioconductor 3.21, Seurat)"

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

USER root

# Git is used to install DoubletFinder; process tools help diagnose pipeline runs.
RUN apt-get update && apt-get install -y --no-install-recommends \
    git ca-certificates procps build-essential pkg-config cmake \
    libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
    libpng-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
    libxt-dev libx11-dev libglu1-mesa-dev zlib1g-dev libbz2-dev \
    liblzma-dev libncurses-dev libfftw3-dev libtiff-dev libjpeg-dev \
    libhdf5-dev libglpk-dev libgmp3-dev \
    && rm -rf /var/lib/apt/lists/*

COPY envs/scrnaseq.yml /tmp/scrnaseq.yml

RUN conda env create -f /tmp/scrnaseq.yml \
    && conda clean --all --yes \
    && rm /tmp/scrnaseq.yml

ENV PATH=/opt/conda/envs/scrnaseq/bin:$PATH \
    CONDA_PREFIX=/opt/conda/envs/scrnaseq

# Some Conda dependencies include packages from a newer Bioconductor
# release. Switch the whole R library to 3.21 before adding study packages.
RUN Rscript -e "options(timeout=600, Ncpus=2); BiocManager::install(version='3.21', ask=FALSE, update=TRUE); stopifnot(as.character(BiocManager::version()) == '3.21')"

RUN Rscript -e "options(timeout=600, Ncpus=2); pkgs <- c('SingleCellExperiment','SummarizedExperiment','scater','scran','DropletUtils','BiocParallel','DESeq2','edgeR','limma','ComplexHeatmap','clusterProfiler','enrichplot','DOSE','GO.db','org.Hs.eg.db','org.Mm.eg.db','AnnotationDbi','AnnotationHub','biomaRt'); BiocManager::install(pkgs, version='3.21', ask=FALSE, update=FALSE); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing)) stop('Bioconductor packages missing after installation: ', paste(missing, collapse=', '))"

# Install DoubletFinder into the same R library.
RUN Rscript -e "options(timeout=600); remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', dependencies=NA, upgrade='never', lib=.libPaths()[1])"

# Fail the build if any requested package is missing or fails to load.
RUN Rscript -e "pkgs <- c('optparse','yaml','tidyverse','data.table','ggplot2','pheatmap','patchwork','cowplot','Matrix','hdf5r','Seurat','SeuratObject','remotes','bookdown','rmarkdown','knitr','SingleCellExperiment','SummarizedExperiment','scater','scran','DropletUtils','BiocParallel','DESeq2','edgeR','limma','ComplexHeatmap','clusterProfiler','enrichplot','DOSE','GO.db','org.Hs.eg.db','org.Mm.eg.db','AnnotationDbi','AnnotationHub','biomaRt','DoubletFinder'); invisible(lapply(pkgs, function(pkg) { library(pkg, character.only=TRUE); message('OK: ',pkg) })); stopifnot(as.character(BiocManager::version()) == '3.21'); sessionInfo()"

# Record the actual versions after every installation step has completed.
COPY scripts/export_environment_manifest.R /usr/local/bin/export_environment_manifest.R
RUN Rscript /usr/local/bin/export_environment_manifest.R /opt/environment-manifest.yml

RUN echo 'conda activate scrnaseq' >> /root/.bashrc

WORKDIR /workspace

CMD ["/bin/bash"]

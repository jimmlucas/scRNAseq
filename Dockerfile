FROM continuumio/miniconda3

LABEL maintainer="Jimmy Lucas" \
      description="Conda environment for the Nextflow scRNA-seq pipeline (R, Bioconductor 3.21, Seurat)"

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

USER root

# Git is used to install DoubletFinder; process tools help diagnose pipeline runs.
RUN apt-get update && apt-get install -y --no-install-recommends \
    git ca-certificates procps \
    && rm -rf /var/lib/apt/lists/*

COPY envs/scrnaseq.yml /tmp/scrnaseq.yml

RUN conda env create -f /tmp/scrnaseq.yml \
    && conda clean --all --yes \
    && rm /tmp/scrnaseq.yml

ENV PATH=/opt/conda/envs/scrnaseq/bin:$PATH \
    CONDA_PREFIX=/opt/conda/envs/scrnaseq

# Install in the same R library as the rest of the environment. Avoid a
# second dependency solver changing the pinned Bioconductor/Seurat stack.
RUN Rscript -e "options(timeout=600); remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', dependencies=NA, upgrade='never', lib=.libPaths()[1])"

# Fail the build if any requested package is missing or fails to load.
RUN Rscript -e "pkgs <- c('optparse','yaml','tidyverse','data.table','ggplot2','pheatmap','patchwork','cowplot','Matrix','hdf5r','Seurat','SeuratObject','remotes','bookdown','rmarkdown','knitr','SingleCellExperiment','SummarizedExperiment','scater','scran','DropletUtils','BiocParallel','DESeq2','edgeR','limma','ComplexHeatmap','clusterProfiler','enrichplot','DOSE','GO.db','org.Hs.eg.db','org.Mm.eg.db','AnnotationDbi','AnnotationHub','biomaRt','DoubletFinder'); invisible(lapply(pkgs, function(pkg) { library(pkg, character.only=TRUE); message('OK: ',pkg) })); stopifnot(as.character(BiocManager::version()) == '3.21'); sessionInfo()"

RUN echo 'conda activate scrnaseq' >> /root/.bashrc

WORKDIR /workspace

CMD ["/bin/bash"]

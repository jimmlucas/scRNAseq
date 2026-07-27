# Analysis Patterns From Legacy Rmds

This note summarizes recurrent patterns observed in previous scRNA-seq,
snRNA-seq-adjacent and CITE-seq analysis scripts. It is not executable code.
It is a design reference for converting past project experience into
structured pipeline causistics, policy gates and report artifacts.

## Scope

Reviewed material includes mouse scRNA-seq Rmds and human PBMC CITE-seq/BCR
Rmds. The studies share a single-cell analysis backbone, but differ in
modality, organism, platform, metadata and downstream scientific questions.

The pipeline should therefore reuse common modules while allowing modality-
specific policies.

## Shared Backbone

Common phases across projects:

```text
load data
standardize metadata
compute QC metrics
filter cells/features using reviewed thresholds
score or detect doublets
normalize data
reduce dimensionality
integrate batches/samples/runs when needed
cluster cells
annotate clusters
run downstream differential or pathway analyses
generate reports and detailed result files
```

## Mouse scRNA Patterns

Observed in the L13712 Rmd family:

```text
RNA-only Seurat objects
organism-specific mitochondrial/ribosomal/globin patterns
QC on nCount_RNA, nFeature_RNA, percent.mito, percent.ribo, percent.globin
mild filtering discussed with researchers
DoubletFinder per sample
doublets labelled before filtering
SCTransform with regressors
Harmony integration by sample
SNN clustering over Harmony reduction
SingleR annotation with mouse references
custom reference annotation for cardiac/fibroblast contexts
pathway analysis with clusterProfiler/org.Mm.eg.db
```

Important causistics:

```text
high mitochondrial signal can be biological in some tissues
ribosomal content can be biological or technical
doublet filtering should be inspected before application
custom references may be useful but require cautious interpretation
```

## Human PBMC CITE-seq Patterns

Observed in the M15214 Rmd family and report:

```text
BD Rhapsody DBEC input matrices
multiple runs
SampleTag multiplexing
RNA assay plus AB assay
BCR metadata and clonotype analysis
metadata enrichment from targets and VDJ files
Multiplet and Undetermined SampleTag handling
possible project-specific SampleTag relabeling
RNA QC plus AB QC
filtering on nCount_RNA, nFeature_RNA, percent.mito, percent.ribo,
percent.globin, nCount_AB and nFeature_AB
protein-level doublet evidence using incompatible marker co-expression
sample and run exclusion after QC review
SCTransform for RNA
CLR normalization for AB
Harmony integration by Sample_Name and Run
WNN for RNA plus AB joint clustering
cluster resolution exploration with clustree
cluster annotation with automatic tools plus researcher review
pseudobulk differential expression with edgeR/limma-voom
differential abundance with limma
BCR clonotype analysis with scRepertoire and AIRR-related tooling
```

Important causistics:

```text
repeated samples across runs require explicit sample resolution
bad-quality repeated runs may be excluded
AB markers can be removed if species/reactivity is inappropriate
AB low feature or high count outliers need modality-aware QC
SampleTag multiplets and undetermined cells should be handled explicitly
WNN should be available only when RNA and AB reductions share cells
BCR analysis should be optional and gated by VDJ metadata availability
```

## Implications For Pipeline Design

The pipeline should not become a large collection of free-form flags. Instead:

```text
metrics produce evidence
causistics summarize recurring situations
policy resolvers propose options
humans review scientific choices
validators enforce approved contracts
application modules transform data reproducibly
```

This keeps standalone Nextflow execution possible while giving external
supervisors enough structure to explain and recommend options.

## Recommended Capability Families

```text
core_input
metadata_standardization
rna_qc
ab_qc
sampletag_resolution
sample_resolution
doublet_scoring
marker_based_doublet_detection
normalization
ab_normalization
dimensionality
integration
multimodal_wnn
clustering
annotation
pseudobulk_de
differential_abundance
pathway_analysis
clonotype_analysis
report_generation
```

## Analysis Type Profiles

scRNA-seq:

```text
RNA QC
RNA doublet detection
RNA normalization
PCA/UMAP
integration when multiple samples/runs require it
clustering
annotation
DE/pathways when requested
```

snRNA-seq:

```text
RNA QC with nuclear-aware thresholds
ambient/intronic considerations
doublet detection
normalization
dimensionality
integration
clustering
annotation
```

CITE-seq:

```text
RNA QC
AB QC
SampleTag resolution when multiplexed
RNA normalization
AB CLR normalization
Harmony or equivalent integration when required
WNN joint analysis
protein marker reporting
BCR/TCR clonotype analysis when VDJ metadata exist
```

## Reporting Implications

Reports should contain both narrative interpretation and pointers to detailed
files. A report-ready phase should publish:

```text
summary tables
before/after QC figures
per-sample/per-run metrics
policy decisions and approval status
warnings and causistics
sessionInfo files
links to detailed artifacts
```

Final reports should not hide decisions. They should say what was detected,
what was proposed, who reviewed it, what was applied, and where the detailed
artifact can be inspected.


# Standalone Real Data Validation

This document freezes the current milestone and defines the next validation
step before adding more scientific modules or building BBUDY in depth.

## Milestone

```text
tag:    pre-integration-realdata-v0.1
commit: d7e6c25
scope:  standalone pipeline through integration policy resolution
```

The goal is to verify that `scRNAseq` works as a standalone scientific engine:

```text
Nextflow + Docker + analysis plan + reviewed policy files + -resume
```

An external supervisor such as BBUDY must not be required to complete this
validation.

## Why Pause Here

The pipeline has reached a meaningful scientific boundary:

```text
INPUT / LOAD
  -> QC
  -> DOUBLETS
  -> POST-DOUBLET READINESS
  -> NORMALIZATION
  -> DIMENSIONALITY / PCA
  -> PC SELECTION POLICY
  -> INTEGRATION POLICY
```

Fixture tests validate contracts, but real datasets are needed to test:

```text
large sample counts
40k-100k cells
unbalanced samples
real batch structure
snRNA-specific QC distributions
CITE-seq modalities
repeated samples across runs
Seurat v5 multi-layer objects
memory release between samples
```

## Validation Questions

The validation should answer:

```text
Can a bioinformatician run the available pipeline without knowing BBUDY exists?
Are policy gates usable manually with reviewed YAML files and -resume?
Do artifact contracts remain stable on real data?
Do resource assumptions hold before integration increases complexity?
```

## Primary Dataset

Start with the complete real L13712 dataset, not the reduced fixture.

Run phase by phase until:

```text
INTEGRATION_POLICY_RESOLUTION
```

At each human gate:

```text
1. inspect proposed policy YAML
2. edit approval.status / reviewer / timestamp manually
3. validate the reviewed policy
4. resume the next phase with nextflow -resume
```

This simulates the future BBUDY control loop without introducing BBUDY yet.

## Additional Datasets

Use at least one dataset with different causistics before continuing module
development:

```text
mouse scRNA-seq
human PBMC / CITE-seq
snRNA-seq
```

The objective is not complete biological interpretation yet. The objective is
to reach integration policy resolution and identify contract gaps.

## Resource Checks

Record per-process:

```text
CPU
duration
peak RAM
input size
output size
container image
Nextflow work hash
```

Pay special attention to:

```text
LOAD_DATA
QC_METRICS
DOUBLET_DETECTION_SCORING
NORMALIZATION
DIMENSIONALITY
```

The desired behavior is:

```text
sample/process chunk
  -> persist artifacts
  -> release memory
  -> resume safely
```

Avoid discovering after downstream expansion that an early step forces a large
global object into memory unnecessarily.

## Expected Outputs To Review

Review these contracts after each real run:

```text
control/execution_modes.yaml
control/module_catalog.yaml
control/causistics_catalog.yaml
control/*policy*/*.yaml
03_post_doublet_qc/03_analysis_readiness.yaml
04_normalization/applied_normalization_policy.yaml
05_dimensionality/05_pca_elbow_metrics.csv
control/pc_selection_policy/*.yaml
control/pc_selection_policy_approval/*.yaml
control/integration_policy/*.yaml
results/pipeline_info/execution_trace.txt
```

## Pass Criteria

The milestone is validated when:

```text
real data reaches integration_policy_resolution
manual policy review works with -resume
cells/features remain consistent across contracts
large phases complete within expected resources
published artifacts are sufficient for human review
no BBUDY-specific logic is required inside Nextflow
```

## Continue After Validation

After real-data validation, continue the scientific engine:

```text
VALIDATE_INTEGRATION_APPROVAL
  -> INTEGRATION
  -> POST_INTEGRATION_DIAGNOSTICS
  -> NEIGHBOR / UMAP
  -> CLUSTERING POLICY
  -> CLUSTERING
```

The first deeper BBUDY MVP should wait until the standalone engine reaches:

```text
INPUT -> QC -> DOUBLETS -> NORMALIZATION -> PCA -> INTEGRATION -> UMAP -> CLUSTERING
```

At that point BBUDY can safely focus on orchestration, interpretation,
recommendations, human review and reproducible `-resume` relaunches.

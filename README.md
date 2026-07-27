# scRNAseq

Modular Nextflow/R pipeline for scRNA-seq, snRNA-seq and CITE-seq analysis.

This repository is the scientific execution engine. It owns reproducible
analysis, containers, checkpoints, validation, published artifacts and
Nextflow `-resume`.

BBUDY can supervise this pipeline from outside, but the pipeline does not
depend on BBUDY.

## Architecture

The project is designed around one pipeline and multiple ways to orchestrate
it:

```text
Human-reviewed policies
        ↓
External supervisor, optional
        ↓
Nextflow scientific engine
        ↓
R / Python scientific and control scripts
```

The separation is intentional:

```text
scRNAseq
  - runs Nextflow processes
  - manages Docker containers
  - validates inputs and policies
  - computes metrics and artifacts
  - applies only approved policies
  - preserves cache and resume behavior

External supervisor
  - interprets published artifacts
  - explains decisions to users
  - requests human review
  - prepares reviewed policy files
  - relaunches Nextflow with -resume
```

The pipeline should also remain usable without any supervisor:

```text
Nextflow + Docker + configs + approved policies
```

## Current Flow

The current plan-driven path reaches normalization:

```text
VALIDATE_ANALYSIS_PLAN
        ↓
RESOLVE_ANALYSIS_PLAN
        ↓
INPUT_DISCOVERY
        ↓
RESOLVE_LOAD_CONFIG
        ↓
LOAD_DATA
        ↓
QC_METRICS
        ↓
QC_THRESHOLD_RESOLUTION
        ↓
VALIDATE_QC_THRESHOLD_APPROVAL
        ↓
QC_FILTER_POLICY
        ↓
DOUBLET_DETECTION_SCORING
        ↓
DOUBLET_POLICY_RESOLUTION
        ↓
VALIDATE_DOUBLET_APPROVAL
        ↓
DOUBLET_FILTER
        ↓
POST_DOUBLET_QC_SUMMARY
        ↓
NORMALIZATION_POLICY_RESOLUTION
        ↓
VALIDATE_NORMALIZATION_APPROVAL
        ↓
NORMALIZATION
```

The next scientific blocks are scaling/PCA and integration policy.

## Policy Pattern

Scientific decisions are separated from scientific execution:

```text
calculate metrics
        ↓
propose policy
        ↓
human review
        ↓
validate reviewed policy
        ↓
apply approved policy
        ↓
publish audit artifacts
```

This pattern is already used for:

```text
QC thresholds
Doublet filtering
Normalization strategy
```

Application phases must not recalculate decision proposals. They consume only
validated policy files.

## Execution Modes

The canonical interface is phase-level execution. This keeps debugging,
checkpoint recovery and external orchestration simple.

Examples:

```bash
nextflow run main.nf \
  --mode analysis_plan_qc_process \
  --analysis_plan configs/analysis_plan_test.yaml \
  --project_id M15214_TEST \
  --project_name CART_scRNAseq_TEST \
  --biological_run_id L13712_TEST \
  -profile docker,local_low_memory
```

```bash
nextflow run main.nf \
  --mode normalization_approval_phase \
  --analysis_plan configs/analysis_plan_test.yaml \
  --approved_normalization_policy path/to/reviewed_policy.yaml \
  --project_id M15214_TEST \
  --project_name CART_scRNAseq_TEST \
  --biological_run_id L13712_TEST \
  -profile docker,local_low_memory \
  -resume
```

Available modern control modes include:

```text
analysis_plan_qc_process
qc_threshold_approval_phase
qc_filter_policy_phase
qc_visualization_policy_phase
doublet_detection_scoring_phase
doublet_policy_resolution_phase
doublet_approval_phase
doublet_filter_phase
post_doublet_qc_summary_phase
normalization_policy_resolution_phase
normalization_approval_phase
normalization_phase
```

Legacy standalone modes are still present while the plan-driven contract
becomes canonical:

```text
qc_process
load_phase
qc_filter_phase
qc_phase_1
qc_phase_2
qc_visualization_phase
doublet_phase
qc_scores_phase
```

The formal execution contract is documented in:

```text
control/execution_modes.yaml
control/module_catalog.yaml
```

## Containers

R analysis uses:

```text
ghcr.io/jimmlucas/bioc-rnaseq-test01new2:v0.3.0
```

Python control-plane modules use:

```text
ghcr.io/jimmlucas/scrnaseq-control-python@sha256:aa6427fba278d9b83cc9f12cc4152a87c3829344f76dd339c57c31a9c70e81d4
```

The control image contains `PyYAML` and is pinned in `nextflow.config`.

## Tests

Architecture contract:

```bash
tests/architecture/run_tests.sh
```

Normalization approval contract:

```bash
tests/normalization_approval/run_tests.sh
```

Doublet and post-doublet regression tests:

```bash
tests/doublet_scoring/run_test.sh
tests/doublet_approval/run_tests.sh
tests/doublet_filter/run_test.sh
tests/post_doublet_qc/run_test.sh
tests/normalization_policy/run_test.sh
tests/normalization_approval/run_tests.sh
tests/normalization/run_test.sh
```

Several tests require Docker and the pinned R/Python images.

## Development Rule

Keep the scientific engine neutral:

```text
Do not put supervisor-specific logic inside Nextflow modules.
Expose artifacts, policies and modes.
Let external tools decide what to run next.
```

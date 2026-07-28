#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

PYTHON_BIN="${PYTHON_BIN:-python3}"
CONDA_PYTHON="/Users/jimmlucas/miniconda3/bin/python3"
CODEX_PYTHON="/Users/jimmlucas/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3"

if [[ -x "${CONDA_PYTHON}" ]]; then
  PYTHON_BIN="${CONDA_PYTHON}"
elif [[ -x "${CODEX_PYTHON}" ]]; then
  PYTHON_BIN="${CODEX_PYTHON}"
fi

cat > "${TMP_DIR}/analysis_plan.yaml" <<'YAML'
plan:
  schema_version: '1.0'
  plan_id: adaptive_alias_contract
experiment:
  project_id: TEST
  biological_run_id: TEST_RUN
  experiment_type: scRNAseq
  material: cells
  organism: human
  modalities:
    - RNA
inputs:
  matrix_type: seurat_object
  load_config: configs/test/load_config.yaml
  sample_column: sample
  biological_sample_column: sample
  group_column: group
routing:
  input_discovery: false
  load_data: false
  qc_metrics: true
  qc_visualization: false
  qc_threshold_resolution: true
  qc_filter: false
  doublet_detection: false
  sample_resolution: false
  merge_samples: false
  qc_scores: false
  normalization: false
qc:
  threshold_mode: adaptive
YAML

cat > "${TMP_DIR}/summary.csv" <<'CSV'
metric,value
cells,4
rna_features,1000
sample_count,2
group_count,1
CSV

cat > "${TMP_DIR}/by_sample.csv" <<'CSV'
sample,Group,cells,median_nCount_RNA,median_nFeature_RNA,median_percent_mito,median_percent_ribo,median_percent_globin,median_log10GenesPerUMI
s1,A,2,1000,500,5,10,0,0.90
s2,A,2,1200,600,6,11,0,0.91
CSV

cat > "${TMP_DIR}/cell_metadata.csv" <<'CSV'
cell,sample,Group,nCount_RNA,nFeature_RNA,percent.mito,percent.ribo,percent.globin,log10GenesPerUMI
c1,s1,A,1000,500,5,10,0,0.90
c2,s1,A,1100,550,7,12,0,0.91
c3,s2,A,1200,600,6,11,0,0.92
c4,s2,A,1300,650,8,13,0,0.93
CSV

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/validate_analysis_plan.py" \
  --plan "${TMP_DIR}/analysis_plan.yaml" \
  --validated-plan "${TMP_DIR}/analysis_plan.validated.yaml" \
  --summary "${TMP_DIR}/validation.tsv" \
  > "${TMP_DIR}/validate.out" \
  2> "${TMP_DIR}/validate.err"

grep -F "qc.threshold_mode='adaptive' is deprecated" \
  "${TMP_DIR}/validate.err" >/dev/null

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/resolve_qc_thresholds.py" \
  --plan "${TMP_DIR}/analysis_plan.validated.yaml" \
  --global-metrics "${TMP_DIR}/summary.csv" \
  --sample-metrics "${TMP_DIR}/by_sample.csv" \
  --cell-metadata "${TMP_DIR}/cell_metadata.csv" \
  --policy-yaml "${TMP_DIR}/qc_threshold_policy.yaml" \
  --policy-tsv "${TMP_DIR}/qc_threshold_policy.tsv" \
  --summary "${TMP_DIR}/resolution.tsv" \
  --session-info "${TMP_DIR}/session.txt"

"${PYTHON_BIN}" - "${TMP_DIR}/qc_threshold_policy.yaml" <<'PY'
from pathlib import Path
import sys
import yaml

policy = yaml.safe_load(Path(sys.argv[1]).read_text(encoding="utf-8"))
assert policy["qc_threshold_policy"]["threshold_mode"] == "data_driven"
assert policy["qc_threshold_policy"]["approval"]["status"] == "proposed"
print("adaptive alias resolves to data_driven")
PY

echo "QC threshold mode contract: VALID"

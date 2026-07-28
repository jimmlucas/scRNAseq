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

RESULTS_DIR="${TMP_DIR}/results"
POLICY_DIR="${RESULTS_DIR}/control/qc_thresholds"
IMPACT_DIR="${RESULTS_DIR}/control/qc_threshold_impact"

mkdir -p "${POLICY_DIR}" "${IMPACT_DIR}"

cat > "${TMP_DIR}/analysis_plan.yaml" <<'YAML'
plan:
  schema_version: '1.0'
  plan_id: review_contract
experiment:
  project_id: REVIEW_TEST
  biological_run_id: RUN_001
  experiment_type: scRNAseq
  material: cells
  organism: human
  modalities:
    - RNA
inputs:
  matrix_type: seurat_object
  load_config: load_config.yaml
  sample_column: SampleName
  biological_sample_column: SampleName
  group_column: Group
routing:
  input_discovery: true
  load_data: true
  qc_metrics: true
  qc_threshold_resolution: true
  qc_filter: false
qc:
  threshold_mode: data_driven
YAML

cat > "${POLICY_DIR}/qc_threshold_policy.yaml" <<'YAML'
qc_threshold_policy:
  schema_version: '1.0'
  plan_id: review_contract
  project_id: REVIEW_TEST
  biological_run_id: RUN_001
  threshold_mode: data_driven
  material: cells
  approval:
    required: true
    status: proposed
    approved_by: null
    approved_at: null
  application_semantics:
    minimum_operator: '>='
    maximum_operator: <=
    combine_metrics_with: AND
    missing_metric_action: error
  thresholds:
    nFeature_RNA:
      enabled: true
      min: 0
      max: 5000
    nCount_RNA:
      enabled: true
      min: 0
      max: 10000
    percent.mito:
      enabled: true
      min: null
      max: 20
    percent.ribo:
      enabled: false
      min: null
      max: null
YAML

cat > "${IMPACT_DIR}/qc_threshold_impact.summary.tsv" <<'TSV'
plan_id	project_id	biological_run_id	threshold_mode	approval_status	cell_count	flagged_any	retained	flagged_fraction	enabled_metric_count
review_contract	REVIEW_TEST	RUN_001	data_driven	proposed	100	7	93	0.07	3
TSV

cat > "${IMPACT_DIR}/qc_threshold_impact_by_sample.csv" <<'CSV'
sample,Group,cells_before,flagged_any,retained,flagged_fraction
S1,A,50,4,46,0.08
S2,B,50,3,47,0.06
CSV

cat > "${IMPACT_DIR}/qc_threshold_impact_by_metric.csv" <<'CSV'
sample,Group,metric,cells_before,flagged_cells,flagged_fraction,unique_flagged_cells,unique_flagged_fraction,overlapped_flagged_cells
ALL,ALL,nCount_RNA,100,5,0.05,3,0.6,2
ALL,ALL,percent.mito,100,2,0.02,2,1.0,0
CSV

cat > "${IMPACT_DIR}/qc_lower_tail_summary.csv" <<'CSV'
metric,threshold,cells_evaluated,cells_below_threshold,fraction_below_threshold
nFeature_RNA,500.0,100,1,0.01
CSV

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/review_qc_policy.py" \
  --results-dir "${RESULTS_DIR}" \
  --output-dir "${TMP_DIR}/approve" \
  --decision approve \
  --reviewer jimmy \
  --approved-at "2026-07-28T10:00:00+00:00" \
  --reason "Impact is balanced across samples."

"${PYTHON_BIN}" - "${TMP_DIR}/approve" <<'PY'
from pathlib import Path
import csv
import sys
import yaml

base = Path(sys.argv[1])
policy = yaml.safe_load(
    (base / "qc_threshold_policy.reviewed.yaml").read_text(
        encoding="utf-8",
    )
)["qc_threshold_policy"]

assert policy["approval"]["status"] == "approved"
assert policy["approval"]["reviewed_by"] == "jimmy"
assert policy["approval"]["approved_by"] == "jimmy"
assert policy["approval"]["approved_at"] == "2026-07-28T10:00:00+00:00"
assert policy["application"]["status"] == "not_applied"
assert policy["review"]["changes"] == []
assert policy["review"]["impact_snapshot"]["flagged_any"] == 7

with (base / "qc_threshold_policy.review.summary.tsv").open(
    encoding="utf-8",
    newline="",
) as handle:
    row = next(csv.DictReader(handle, delimiter="\t"))

assert row["review_status"] == "approved"
assert row["flagged_any"] == "7"
print("QC policy approve review contract: VALID")
PY

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/validate_qc_threshold_policy.py" \
  --plan "${TMP_DIR}/analysis_plan.yaml" \
  --policy "${TMP_DIR}/approve/qc_threshold_policy.reviewed.yaml" \
  --validated-policy "${TMP_DIR}/approve/approved_qc_threshold_policy.yaml" \
  --summary "${TMP_DIR}/approve/qc_threshold_approval.summary.tsv" \
  --session-info "${TMP_DIR}/approve/qc_threshold_approval.session.txt" \
  > "${TMP_DIR}/approve/validation.out"

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/review_qc_policy.py" \
  --results-dir "${RESULTS_DIR}" \
  --output-dir "${TMP_DIR}/modified" \
  --decision modify-and-approve \
  --reviewer jimmy \
  --approved-at "2026-07-28T10:05:00+00:00" \
  --set percent.mito.max=15 \
  --set thresholds.nFeature_RNA.min=200 \
  --reason "Tightened bounds after impact review."

"${PYTHON_BIN}" - "${TMP_DIR}/modified" <<'PY'
from pathlib import Path
import sys
import yaml

policy = yaml.safe_load(
    (
        Path(sys.argv[1])
        / "qc_threshold_policy.reviewed.yaml"
    ).read_text(encoding="utf-8")
)["qc_threshold_policy"]

assert policy["approval"]["status"] == "modified_and_approved"
assert policy["thresholds"]["percent.mito"]["max"] == 15.0
assert policy["thresholds"]["nFeature_RNA"]["min"] == 200.0
assert len(policy["review"]["changes"]) == 2
assert policy["application"]["status"] == "not_applied"
print("QC policy modified review contract: VALID")
PY

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/validate_qc_threshold_policy.py" \
  --plan "${TMP_DIR}/analysis_plan.yaml" \
  --policy "${TMP_DIR}/modified/qc_threshold_policy.reviewed.yaml" \
  --validated-policy "${TMP_DIR}/modified/approved_qc_threshold_policy.yaml" \
  --summary "${TMP_DIR}/modified/qc_threshold_approval.summary.tsv" \
  --session-info "${TMP_DIR}/modified/qc_threshold_approval.session.txt" \
  > "${TMP_DIR}/modified/validation.out"

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/review_qc_policy.py" \
  --results-dir "${RESULTS_DIR}" \
  --output-dir "${TMP_DIR}/rejected" \
  --decision reject \
  --reviewer jimmy \
  --approved-at "2026-07-28T10:10:00+00:00" \
  --reason "Policy requires additional review."

"${PYTHON_BIN}" - "${TMP_DIR}/rejected" <<'PY'
from pathlib import Path
import sys
import yaml

policy = yaml.safe_load(
    (
        Path(sys.argv[1])
        / "qc_threshold_policy.reviewed.yaml"
    ).read_text(encoding="utf-8")
)["qc_threshold_policy"]

assert policy["approval"]["status"] == "rejected"
assert policy["approval"]["approved_at"] is None
assert policy["approval"]["approved_by"] is None
assert policy["application"]["status"] == "not_applied"
print("QC policy reject review contract: VALID")
PY

if "${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/validate_qc_threshold_policy.py" \
  --plan "${TMP_DIR}/analysis_plan.yaml" \
  --policy "${TMP_DIR}/rejected/qc_threshold_policy.reviewed.yaml" \
  --validated-policy "${TMP_DIR}/rejected/approved_qc_threshold_policy.yaml" \
  --summary "${TMP_DIR}/rejected/qc_threshold_approval.summary.tsv" \
  --session-info "${TMP_DIR}/rejected/qc_threshold_approval.session.txt" \
  > "${TMP_DIR}/rejected/validation.out" \
  2> "${TMP_DIR}/rejected/validation.err"; then
  echo "rejected policy unexpectedly validated" >&2
  exit 1
fi

grep -q -- "not approved" "${TMP_DIR}/rejected/validation.err"

if "${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/review_qc_policy.py" \
  --results-dir "${RESULTS_DIR}" \
  --output-dir "${TMP_DIR}/bad_approve" \
  --decision approve \
  --reviewer jimmy \
  --set percent.mito.max=15 \
  2> "${TMP_DIR}/bad_approve.err"; then
  echo "approve unexpectedly accepted --set" >&2
  exit 1
fi

grep -q -- "--set is only accepted" "${TMP_DIR}/bad_approve.err"

if "${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/review_qc_policy.py" \
  --results-dir "${RESULTS_DIR}" \
  --output-dir "${TMP_DIR}/bad_modify" \
  --decision modify-and-approve \
  --reviewer jimmy \
  2> "${TMP_DIR}/bad_modify.err"; then
  echo "modify-and-approve unexpectedly accepted no changes" >&2
  exit 1
fi

grep -q -- "requires at least one --set" "${TMP_DIR}/bad_modify.err"

echo "QC policy review tests: VALID"

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

cat > "${TMP_DIR}/qc_threshold_policy.yaml" <<'YAML'
qc_threshold_policy:
  schema_version: '1.0'
  plan_id: impact_contract
  project_id: TEST
  biological_run_id: RUN
  threshold_mode: data_driven
  approval:
    required: true
    status: proposed
  application_semantics:
    minimum_operator: '>='
    maximum_operator: <=
    combine_metrics_with: AND
    missing_metric_action: error
  thresholds:
    nFeature_RNA:
      enabled: true
      min: 200
      max: 5000
    nCount_RNA:
      enabled: true
      min: 500
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

cat > "${TMP_DIR}/cell_metadata.csv" <<'CSV'
cell,sample,Group,nFeature_RNA,nCount_RNA,percent.mito,percent.ribo
c1,S1,A,1000,2000,5,10
c2,S1,A,100,2000,5,10
c3,S1,A,6000,12000,5,10
c4,S2,B,1000,2000,30,10
c5,S2,B,7000,15000,30,10
CSV

cat > "${TMP_DIR}/by_sample.csv" <<'CSV'
sample,Group,cells,median_nCount_RNA,median_nFeature_RNA,median_percent_mito
S1,A,3,2000,1000,5
S2,B,2,8500,4000,30
CSV

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/summarize_qc_threshold_impact.py" \
  --policy-yaml "${TMP_DIR}/qc_threshold_policy.yaml" \
  --cell-metadata "${TMP_DIR}/cell_metadata.csv" \
  --by-sample "${TMP_DIR}/by_sample.csv" \
  --impact-yaml "${TMP_DIR}/qc_threshold_impact.yaml" \
  --by-sample-output "${TMP_DIR}/qc_threshold_impact_by_sample.csv" \
  --by-metric-output "${TMP_DIR}/qc_threshold_impact_by_metric.csv" \
  --overlap-output "${TMP_DIR}/qc_threshold_overlap_summary.csv" \
  --summary "${TMP_DIR}/qc_threshold_impact.summary.tsv" \
  --session-info "${TMP_DIR}/qc_threshold_impact.session.txt" \
  --project-id TEST \
  --biological-run-id RUN

"${PYTHON_BIN}" - "${TMP_DIR}" <<'PY'
from pathlib import Path
import csv
import sys
import yaml

base = Path(sys.argv[1])

impact = yaml.safe_load(
    (base / "qc_threshold_impact.yaml").read_text(encoding="utf-8")
)["qc_threshold_impact"]

assert impact["application_status"] == "not_applied"
assert impact["cell_count"] == 5
assert impact["flagged_any"] == 4
assert impact["retained"] == 1
assert impact["enabled_metrics"] == [
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mito",
]

with (base / "qc_threshold_impact_by_sample.csv").open(
    encoding="utf-8",
    newline="",
) as handle:
    by_sample = {
        row["sample"]: row
        for row in csv.DictReader(handle)
    }

assert int(by_sample["S1"]["cells_before"]) == 3
assert int(by_sample["S1"]["flagged_any"]) == 2
assert int(by_sample["S1"]["retained"]) == 1
assert int(by_sample["S2"]["cells_before"]) == 2
assert int(by_sample["S2"]["flagged_any"]) == 2
assert int(by_sample["S2"]["retained"]) == 0

with (base / "qc_threshold_impact_by_metric.csv").open(
    encoding="utf-8",
    newline="",
) as handle:
    by_metric = {
        (row["sample"], row["metric"]): row
        for row in csv.DictReader(handle)
    }

assert int(by_metric[("S1", "nFeature_RNA")]["flagged_cells"]) == 2
assert int(by_metric[("S1", "nFeature_RNA")]["unique_flagged_cells"]) == 1
assert int(by_metric[("S1", "nFeature_RNA")]["overlapped_flagged_cells"]) == 1
assert int(by_metric[("S1", "nCount_RNA")]["flagged_cells"]) == 1
assert int(by_metric[("S1", "nCount_RNA")]["unique_flagged_cells"]) == 0
assert int(by_metric[("S1", "nCount_RNA")]["overlapped_flagged_cells"]) == 1
assert int(by_metric[("S1", "percent.mito")]["flagged_cells"]) == 0
assert int(by_metric[("S2", "percent.mito")]["flagged_cells"]) == 2
assert int(by_metric[("S2", "percent.mito")]["unique_flagged_cells"]) == 1

assert int(by_metric[("ALL", "nFeature_RNA")]["flagged_cells"]) == 3
assert int(by_metric[("ALL", "nFeature_RNA")]["unique_flagged_cells"]) == 1
assert int(by_metric[("ALL", "nCount_RNA")]["flagged_cells"]) == 2
assert int(by_metric[("ALL", "nCount_RNA")]["unique_flagged_cells"]) == 0
assert int(by_metric[("ALL", "percent.mito")]["flagged_cells"]) == 2
assert int(by_metric[("ALL", "percent.mito")]["unique_flagged_cells"]) == 1

with (base / "qc_threshold_overlap_summary.csv").open(
    encoding="utf-8",
    newline="",
) as handle:
    overlap = {
        (row["metric_a"], row["metric_b"]): row
        for row in csv.DictReader(handle)
    }

key = ("nCount_RNA", "nFeature_RNA")
assert int(overlap[key]["flagged_both"]) == 2

print("QC threshold impact contract: VALID")
PY

if "${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/summarize_qc_threshold_impact.py" \
  --policy-yaml "${TMP_DIR}/qc_threshold_policy.yaml" \
  --cell-metadata "${TMP_DIR}/cell_metadata.csv" \
  --by-sample "${TMP_DIR}/by_sample.csv" \
  --impact-yaml "${TMP_DIR}/identity.yaml" \
  --by-sample-output "${TMP_DIR}/identity_by_sample.csv" \
  --by-metric-output "${TMP_DIR}/identity_by_metric.csv" \
  --overlap-output "${TMP_DIR}/identity_overlap.csv" \
  --summary "${TMP_DIR}/identity.summary.tsv" \
  --session-info "${TMP_DIR}/identity.session.txt" \
  --project-id WRONG \
  2> "${TMP_DIR}/identity.err"; then
  echo "Expected identity mismatch failure" >&2
  exit 1
fi

grep -F "Policy project_id mismatch" \
  "${TMP_DIR}/identity.err" >/dev/null

cat > "${TMP_DIR}/cell_metadata_missing_metric.csv" <<'CSV'
cell,sample,Group,nFeature_RNA,nCount_RNA,percent.ribo
c1,S1,A,1000,2000,10
CSV

if "${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/summarize_qc_threshold_impact.py" \
  --policy-yaml "${TMP_DIR}/qc_threshold_policy.yaml" \
  --cell-metadata "${TMP_DIR}/cell_metadata_missing_metric.csv" \
  --by-sample "${TMP_DIR}/by_sample.csv" \
  --impact-yaml "${TMP_DIR}/missing_metric.yaml" \
  --by-sample-output "${TMP_DIR}/missing_metric_by_sample.csv" \
  --by-metric-output "${TMP_DIR}/missing_metric_by_metric.csv" \
  --overlap-output "${TMP_DIR}/missing_metric_overlap.csv" \
  --summary "${TMP_DIR}/missing_metric.summary.tsv" \
  --session-info "${TMP_DIR}/missing_metric.session.txt" \
  2> "${TMP_DIR}/missing_metric.err"; then
  echo "Expected missing metric failure" >&2
  exit 1
fi

grep -F "Enabled QC metric is missing from metadata: percent.mito" \
  "${TMP_DIR}/missing_metric.err" >/dev/null

cat > "${TMP_DIR}/cell_metadata_bad_value.csv" <<'CSV'
cell,sample,Group,nFeature_RNA,nCount_RNA,percent.mito,percent.ribo
c1,S1,A,1000,not_numeric,5,10
CSV

if "${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/summarize_qc_threshold_impact.py" \
  --policy-yaml "${TMP_DIR}/qc_threshold_policy.yaml" \
  --cell-metadata "${TMP_DIR}/cell_metadata_bad_value.csv" \
  --by-sample "${TMP_DIR}/by_sample.csv" \
  --impact-yaml "${TMP_DIR}/bad_value.yaml" \
  --by-sample-output "${TMP_DIR}/bad_value_by_sample.csv" \
  --by-metric-output "${TMP_DIR}/bad_value_by_metric.csv" \
  --overlap-output "${TMP_DIR}/bad_value_overlap.csv" \
  --summary "${TMP_DIR}/bad_value.summary.tsv" \
  --session-info "${TMP_DIR}/bad_value.session.txt" \
  2> "${TMP_DIR}/bad_value.err"; then
  echo "Expected non-numeric metric failure" >&2
  exit 1
fi

grep -F "Enabled QC metric contains missing or non-numeric values: nCount_RNA" \
  "${TMP_DIR}/bad_value.err" >/dev/null

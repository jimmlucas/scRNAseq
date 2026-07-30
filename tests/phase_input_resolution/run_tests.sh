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

OUTDIR="${TMP_DIR}/results/PHASE_TEST"
CONFIG_DIR="${TMP_DIR}/repo/configs/real_validation/PHASE_TEST"

mkdir -p \
  "${OUTDIR}/01a_qc_metrics" \
  "${OUTDIR}/control/qc_threshold_approval" \
  "${OUTDIR}/01e_qc_filter_policy" \
  "${OUTDIR}/01f_split_seurat_by_sample/sample_objects" \
  "${CONFIG_DIR}"

touch \
  "${OUTDIR}/01a_qc_metrics/seu_qc_metrics.rds" \
  "${OUTDIR}/control/qc_threshold_approval/approved_qc_threshold_policy.yaml" \
  "${OUTDIR}/01e_qc_filter_policy/seu_qc_filtered.rds" \
  "${OUTDIR}/01f_split_seurat_by_sample/01f_split_by_sample_manifest.tsv" \
  "${OUTDIR}/01f_split_seurat_by_sample/01f_split_by_sample_summary.csv" \
  "${CONFIG_DIR}/doublet_config.yaml"

pushd "${TMP_DIR}/repo" >/dev/null

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/resolve_phase_inputs.py" \
  --outdir "${OUTDIR}" \
  --phase qc_filter_policy_phase \
  --format shell \
  --require-existing \
  > "${TMP_DIR}/qc_filter.env"

grep -q "INPUT_OBJECT=.*seu_qc_metrics.rds" "${TMP_DIR}/qc_filter.env"
grep -q "APPROVED_QC_THRESHOLD_POLICY=.*approved_qc_threshold_policy.yaml" \
  "${TMP_DIR}/qc_filter.env"

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/resolve_phase_inputs.py" \
  --outdir "${OUTDIR}" \
  --phase split_seurat_by_sample_phase \
  --format shell \
  --require-existing \
  > "${TMP_DIR}/split.env"

grep -q "INPUT_OBJECT=.*seu_qc_filtered.rds" "${TMP_DIR}/split.env"

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/resolve_phase_inputs.py" \
  --outdir "${OUTDIR}" \
  --phase doublet_detection_scoring_phase \
  --project-id PHASE_TEST \
  --format json \
  --require-existing \
  > "${TMP_DIR}/doublet.json"

"${PYTHON_BIN}" - "${TMP_DIR}/doublet.json" <<'PY'
from pathlib import Path
import json
import sys

document = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
paths = document["paths"]

assert paths["input_object"]["artifact"] == "qc_filtered_object"
assert paths["input_object"]["exists"] is True
assert paths["input_object"]["path"].endswith(
    "01e_qc_filter_policy/seu_qc_filtered.rds"
)
assert paths["doublet_config"]["artifact"] == "doublet_config"
assert paths["doublet_config"]["exists"] is True
assert paths["doublet_config"]["path"].endswith(
    "configs/real_validation/PHASE_TEST/doublet_config.yaml"
)
print("Phase input JSON contract: VALID")
PY

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/resolve_phase_inputs.py" \
  --outdir "${OUTDIR}" \
  --artifact qc_filtered_object \
  --artifact doublet_config \
  --project-id PHASE_TEST \
  --format tsv \
  --require-existing \
  > "${TMP_DIR}/artifacts.tsv"

grep -q $'qc_filtered_object\tqc_filtered_object\t' "${TMP_DIR}/artifacts.tsv"
grep -q $'doublet_config\tdoublet_config\t' "${TMP_DIR}/artifacts.tsv"

if "${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/resolve_phase_inputs.py" \
  --outdir "${OUTDIR}" \
  --phase doublet_filter_phase \
  --format shell \
  --require-existing \
  > "${TMP_DIR}/missing.out" \
  2> "${TMP_DIR}/missing.err"; then
  echo "Missing doublet filter inputs unexpectedly resolved" >&2
  exit 1
fi

grep -q "Missing required resolved path" "${TMP_DIR}/missing.err"
grep -q "doublet_scored_object" "${TMP_DIR}/missing.err"
grep -q "approved_doublet_policy" "${TMP_DIR}/missing.err"

popd >/dev/null

echo "Phase input resolution tests: VALID"

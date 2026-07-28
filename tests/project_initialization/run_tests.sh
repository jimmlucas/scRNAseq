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

TARGETS="${TMP_DIR}/targets.csv"
COUNT_DIR="${TMP_DIR}/counts"
OUT_DIR="${TMP_DIR}/generated"

mkdir -p "${COUNT_DIR}/S1" "${COUNT_DIR}/S2"
touch "${COUNT_DIR}/S1/S1.rds" "${COUNT_DIR}/S2/S2.rds"

cat > "${TARGETS}" <<'CSV'
SampleName,ShortName,Group
sample_1,S1,A
sample_2,S2,B
CSV

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/init_project.py" \
  --project-id INIT_TEST \
  --biological-run-id RUN_001 \
  --experiment-type scRNAseq \
  --organism human \
  --material cells \
  --targets "${TARGETS}" \
  --count-dir "${COUNT_DIR}" \
  --output-dir "${OUT_DIR}"

test -s "${OUT_DIR}/load_config.yaml"
test -s "${OUT_DIR}/analysis_plan.yaml"

"${PYTHON_BIN}" - "${OUT_DIR}/analysis_plan.yaml" <<'PY'
from pathlib import Path
import sys
import yaml

plan = yaml.safe_load(Path(sys.argv[1]).read_text(encoding="utf-8"))
assert plan["qc"]["threshold_mode"] == "data_driven"
assert plan["experiment"]["project_id"] == "INIT_TEST"
assert plan["routing"]["qc_threshold_resolution"] is True
print("Project initializer default contract: VALID")
PY

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/validate_analysis_plan.py" \
  --plan "${OUT_DIR}/analysis_plan.yaml" \
  --validated-plan "${TMP_DIR}/validated.yaml" \
  --summary "${TMP_DIR}/validation.tsv" \
  > "${TMP_DIR}/validation.out"

cat > "${TMP_DIR}/targets_missing_shortname.csv" <<'CSV'
SampleName,Group
sample_1,A
CSV

if "${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/init_project.py" \
  --project-id INIT_TEST \
  --biological-run-id RUN_001 \
  --experiment-type scRNAseq \
  --organism human \
  --material cells \
  --targets "${TMP_DIR}/targets_missing_shortname.csv" \
  --count-dir "${COUNT_DIR}" \
  --output-dir "${TMP_DIR}/missing_shortname" \
  2> "${TMP_DIR}/missing_shortname.err"; then
  echo "Expected missing ShortName failure" >&2
  exit 1
fi

grep -F "Missing required columns in targets: ShortName" \
  "${TMP_DIR}/missing_shortname.err" >/dev/null

if "${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/init_project.py" \
  --project-id INIT_TEST \
  --biological-run-id RUN_001 \
  --experiment-type scRNAseq \
  --organism human \
  --material cells \
  --targets "${TARGETS}" \
  --count-dir "${TMP_DIR}/missing_counts" \
  --output-dir "${TMP_DIR}/missing_counts_out" \
  2> "${TMP_DIR}/missing_counts.err"; then
  echo "Expected missing count directory failure" >&2
  exit 1
fi

grep -F "Count directory not found" \
  "${TMP_DIR}/missing_counts.err" >/dev/null

touch "${COUNT_DIR}/S1/S1_extra.rds"

if "${PYTHON_BIN}" "${ROOT_DIR}/scripts/control/init_project.py" \
  --project-id INIT_TEST \
  --biological-run-id RUN_001 \
  --experiment-type scRNAseq \
  --organism human \
  --material cells \
  --targets "${TARGETS}" \
  --count-dir "${COUNT_DIR}" \
  --output-dir "${TMP_DIR}/multiple_rds" \
  2> "${TMP_DIR}/multiple_rds.err"; then
  echo "Expected multiple RDS failure" >&2
  exit 1
fi

grep -F "Multiple RDS files found for sample S1" \
  "${TMP_DIR}/multiple_rds.err" >/dev/null

echo "Project initialization contract: VALID"

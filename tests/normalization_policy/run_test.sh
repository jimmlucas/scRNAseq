#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(
  cd "$(dirname "${BASH_SOURCE[0]}")/../.." &&
  pwd
)"

cd "$PROJECT_ROOT"

READINESS="${READINESS:-results_test_repro/03_post_doublet_qc/03_analysis_readiness.yaml}"
SUMMARY="${SUMMARY:-results_test_repro/03_post_doublet_qc/03_post_doublet_qc_summary.csv}"
BY_SAMPLE="${BY_SAMPLE:-results_test_repro/03_post_doublet_qc/03_post_doublet_qc_by_sample.csv}"
WARNINGS="${WARNINGS:-results_test_repro/03_post_doublet_qc/03_post_doublet_qc_warnings.tsv}"
OUTDIR="${OUTDIR:-results_test_repro}"
EXECUTION_ID="${EXECUTION_ID:-NORMALIZATION_POLICY_REPRO_001}"

echo "=========================================="
echo "Normalization policy resolution test"
echo "=========================================="
echo "Project root:    $PROJECT_ROOT"
echo "Readiness:       $READINESS"
echo "Summary:         $SUMMARY"
echo "By sample:       $BY_SAMPLE"
echo "Warnings:        $WARNINGS"
echo "Output dir:      $OUTDIR"
echo "Execution ID:    $EXECUTION_ID"
echo

command -v nextflow >/dev/null 2>&1 || {
  echo "ERROR: nextflow is not installed or not in PATH." >&2
  exit 1
}

if [[ ! -f "$READINESS" || ! -f "$SUMMARY" || ! -f "$BY_SAMPLE" || ! -f "$WARNINGS" ]]; then
  echo "Post-doublet readiness inputs not found; running readiness test first."
  tests/post_doublet_qc/run_test.sh
fi

for required_file in "$READINESS" "$SUMMARY" "$BY_SAMPLE" "$WARNINGS"; do
  if [[ ! -f "$required_file" ]]; then
    echo "ERROR: required input not found: $required_file" >&2
    exit 1
  fi
done

echo "Running Nextflow normalization policy resolution phase..."

nextflow run main.nf \
  --mode normalization_policy_resolution_phase \
  --analysis_readiness "$READINESS" \
  --post_doublet_qc_summary "$SUMMARY" \
  --post_doublet_qc_by_sample "$BY_SAMPLE" \
  --post_doublet_qc_warnings "$WARNINGS" \
  --project_id M15214_TEST \
  --project_name CART_scRNAseq_TEST \
  --biological_run_id L13712_TEST \
  --execution_id "$EXECUTION_ID" \
  --outdir "$OUTDIR" \
  -profile docker,local_low_memory \
  -resume

OUTPUT_DIR="$OUTDIR/control/normalization_policy"

echo
echo "Validating generated normalization policy..."

python3 tests/normalization_policy/validate_normalization_policy.py \
  "$PROJECT_ROOT" \
  "$OUTPUT_DIR"

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

python3 - "$PROJECT_ROOT" "$READINESS" "$TMP_DIR" <<'PY'
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

root = Path(sys.argv[1])
readiness_path = Path(sys.argv[2])
tmp = Path(sys.argv[3])
resolver_path = root / "scripts/control/resolve_normalization_policy.py"

spec = importlib.util.spec_from_file_location(
    "resolve_normalization_policy",
    resolver_path,
)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

readiness = module.read_yaml(readiness_path, "Readiness")
readiness["analysis_readiness"]["status"] = "blocked"
module.write_yaml(tmp / "blocked_readiness.yaml", readiness)
PY

if python3 scripts/control/resolve_normalization_policy.py \
  --readiness "$TMP_DIR/blocked_readiness.yaml" \
  --summary "$SUMMARY" \
  --by-sample "$BY_SAMPLE" \
  --warnings "$WARNINGS" \
  --policy-yaml "$TMP_DIR/normalization_policy.yaml" \
  --policy-tsv "$TMP_DIR/normalization_policy.tsv" \
  --resolution-summary "$TMP_DIR/normalization_policy_resolution.summary.tsv" \
  --session-info "$TMP_DIR/normalization_policy_resolution.session.txt" \
  --project-id M15214_TEST \
  --project-name CART_scRNAseq_TEST \
  --biological-run-id L13712_TEST >/tmp/normalization_policy_blocked.stdout 2>/tmp/normalization_policy_blocked.stderr
then
  echo "Expected blocked readiness to fail" >&2
  cat /tmp/normalization_policy_blocked.stdout >&2
  exit 1
fi

echo
echo "=========================================="
echo "Normalization policy test completed successfully"
echo "=========================================="

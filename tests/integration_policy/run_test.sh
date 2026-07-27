#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(
  cd "$(dirname "${BASH_SOURCE[0]}")/../.." &&
  pwd
)"

cd "$PROJECT_ROOT"

PCA_OBJECT="${PCA_OBJECT:-results_test_repro/05_dimensionality/seu_pca.rds}"
DIMENSIONALITY_SUMMARY="${DIMENSIONALITY_SUMMARY:-results_test_repro/05_dimensionality/05_dimensionality_summary.csv}"
PC_POLICY="${PC_POLICY:-results_test_repro/control/pc_selection_policy/pc_selection_policy.yaml}"
READINESS="${READINESS:-results_test_repro/03_post_doublet_qc/03_analysis_readiness.yaml}"
NORMALIZATION_POLICY="${NORMALIZATION_POLICY:-results_test_repro/04_normalization/applied_normalization_policy.yaml}"
OUTDIR="${OUTDIR:-results_test_repro}"
EXECUTION_ID="${EXECUTION_ID:-INTEGRATION_POLICY_REPRO_001}"
TMP_DIR="$(mktemp -d)"

cleanup() {
  rm -rf "$TMP_DIR"
}

trap cleanup EXIT

echo "=========================================="
echo "Integration policy resolution test"
echo "=========================================="
echo "Project root:             $PROJECT_ROOT"
echo "PCA object:               $PCA_OBJECT"
echo "Dimensionality summary:   $DIMENSIONALITY_SUMMARY"
echo "PC policy:                $PC_POLICY"
echo "Readiness:                $READINESS"
echo "Normalization policy:     $NORMALIZATION_POLICY"
echo "Output dir:               $OUTDIR"
echo "Execution ID:             $EXECUTION_ID"
echo

command -v nextflow >/dev/null 2>&1 || {
  echo "ERROR: nextflow is not installed or not in PATH." >&2
  exit 1
}

if [[ ! -f "$PCA_OBJECT" || ! -f "$DIMENSIONALITY_SUMMARY" ]]; then
  echo "Dimensionality inputs not found; running dimensionality test first."
  tests/dimensionality/run_test.sh
fi

if [[ ! -f "$PC_POLICY" ]]; then
  echo "PC-selection policy not found; running PC-selection resolver test first."
  tests/pc_selection_policy/run_test.sh
fi

for required_file in "$PCA_OBJECT" "$DIMENSIONALITY_SUMMARY" "$PC_POLICY" "$READINESS" "$NORMALIZATION_POLICY"; do
  if [[ ! -f "$required_file" ]]; then
    echo "ERROR: required input not found: $required_file" >&2
    exit 1
  fi
done

REVIEWED_PC_POLICY="$TMP_DIR/reviewed_pc_selection_policy.yaml"
APPROVED_PC_POLICY="$TMP_DIR/approved_pc_selection_policy.yaml"

python3 - "$PROJECT_ROOT" "$PC_POLICY" "$REVIEWED_PC_POLICY" <<'PY'
from __future__ import annotations

import copy
import importlib.util
import sys
from pathlib import Path

root = Path(sys.argv[1])
policy_path = Path(sys.argv[2])
reviewed_path = Path(sys.argv[3])
validator_path = root / "scripts/control/validate_pc_selection_approval.py"

spec = importlib.util.spec_from_file_location(
    "validate_pc_selection_approval",
    validator_path,
)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

document = module.read_yaml(policy_path, "PC-selection policy")
document = copy.deepcopy(document)
approval = document["pc_selection_policy"]["approval"]
approval["status"] = "approved"
approval["reviewed_by"] = "jimmy"
approval["approved_at"] = "2026-07-27T12:00:00+02:00"
module.write_yaml(reviewed_path, document)
PY

python3 scripts/control/validate_pc_selection_approval.py \
  --plan configs/analysis_plan_test.yaml \
  --policy "$REVIEWED_PC_POLICY" \
  --validated-policy "$APPROVED_PC_POLICY" \
  --summary "$TMP_DIR/pc_selection_policy_approval.summary.tsv" \
  --session-info "$TMP_DIR/pc_selection_policy_approval.session.txt"

echo
echo "Running Nextflow integration policy resolution phase..."

nextflow run main.nf \
  --mode integration_policy_resolution_phase \
  --input_object "$PCA_OBJECT" \
  --dimensionality_summary "$DIMENSIONALITY_SUMMARY" \
  --approved_pc_selection_policy "$APPROVED_PC_POLICY" \
  --analysis_readiness "$READINESS" \
  --applied_normalization_policy "$NORMALIZATION_POLICY" \
  --project_id M15214_TEST \
  --project_name CART_scRNAseq_TEST \
  --biological_run_id L13712_TEST \
  --execution_id "$EXECUTION_ID" \
  --outdir "$OUTDIR" \
  -profile docker,local_low_memory \
  -resume

OUTPUT_DIR="$OUTDIR/control/integration_policy"

echo
echo "Validating generated integration policy..."

python3 tests/integration_policy/validate_integration_policy.py \
  "$PROJECT_ROOT" \
  "$OUTPUT_DIR"

ruby -ryaml - "$APPROVED_PC_POLICY" "$TMP_DIR/bad_pc_policy.yaml" <<'RB'
source = ARGV.fetch(0)
target = ARGV.fetch(1)

document = YAML.load_file(source)
document["pc_selection_policy"]["validation"]["status"] = "invalid"
File.write(target, document.to_yaml)
RB

if python3 scripts/control/resolve_integration_policy.py \
  --pca-object "$PCA_OBJECT" \
  --dimensionality-summary "$DIMENSIONALITY_SUMMARY" \
  --approved-pc-selection-policy "$TMP_DIR/bad_pc_policy.yaml" \
  --analysis-readiness "$READINESS" \
  --applied-normalization-policy "$NORMALIZATION_POLICY" \
  --policy-yaml "$TMP_DIR/integration_policy.yaml" \
  --policy-tsv "$TMP_DIR/integration_policy.tsv" \
  --resolution-summary "$TMP_DIR/integration_policy_resolution.summary.tsv" \
  --session-info "$TMP_DIR/integration_policy_resolution.session.txt" \
  --project-id M15214_TEST \
  --project-name CART_scRNAseq_TEST \
  --biological-run-id L13712_TEST >/tmp/integration_policy_bad.stdout 2>/tmp/integration_policy_bad.stderr
then
  echo "Expected invalid PC-selection validation to fail" >&2
  cat /tmp/integration_policy_bad.stdout >&2
  exit 1
fi

echo "Integration policy negative check: VALID"

ruby -ryaml - "$READINESS" "$TMP_DIR/single_sample_readiness.yaml" <<'RB'
source = ARGV.fetch(0)
target = ARGV.fetch(1)

document = YAML.load_file(source)
dataset = document["analysis_readiness"]["dataset"]
dataset["samples"] = 1
File.write(target, document.to_yaml)
RB

ruby -ryaml - "$NORMALIZATION_POLICY" "$TMP_DIR/single_sample_normalization_policy.yaml" <<'RB'
source = ARGV.fetch(0)
target = ARGV.fetch(1)

document = YAML.load_file(source)
document["normalization_policy"]["input"]["samples"] = 1
document["normalization_policy"]["integration"]["required"] = false
document["normalization_policy"]["integration"]["reason"] = "single_sample_dataset"
File.write(target, document.to_yaml)
RB

python3 scripts/control/resolve_integration_policy.py \
  --pca-object "$PCA_OBJECT" \
  --dimensionality-summary "$DIMENSIONALITY_SUMMARY" \
  --approved-pc-selection-policy "$APPROVED_PC_POLICY" \
  --analysis-readiness "$TMP_DIR/single_sample_readiness.yaml" \
  --applied-normalization-policy "$TMP_DIR/single_sample_normalization_policy.yaml" \
  --policy-yaml "$TMP_DIR/single_sample_integration_policy.yaml" \
  --policy-tsv "$TMP_DIR/single_sample_integration_policy.tsv" \
  --resolution-summary "$TMP_DIR/single_sample_integration_policy_resolution.summary.tsv" \
  --session-info "$TMP_DIR/single_sample_integration_policy_resolution.session.txt" \
  --project-id M15214_TEST \
  --project-name CART_scRNAseq_TEST \
  --biological-run-id L13712_TEST

python3 - "$PROJECT_ROOT" "$TMP_DIR/single_sample_integration_policy.yaml" <<'PY'
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

root = Path(sys.argv[1])
policy_path = Path(sys.argv[2])
resolver_path = root / "scripts/control/resolve_integration_policy.py"

spec = importlib.util.spec_from_file_location(
    "resolve_integration_policy",
    resolver_path,
)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

policy = module.read_yaml(policy_path, "Integration policy")[
    "integration_policy"
]
strategy = policy["strategy"]

assert strategy["integration_required"] is False
assert strategy["method"] == "none"
assert strategy["input_reduction"] == "pca"
assert strategy["output_reduction"] == "pca"
assert strategy["batch_variable"] is None
assert policy["approval"]["status"] == "proposed"
assert policy["application"]["status"] == "not_applied"

print("Integration policy no-integration contract: VALID")
PY

echo
echo "=========================================="
echo "Integration policy test completed successfully"
echo "=========================================="

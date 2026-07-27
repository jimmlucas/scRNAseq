#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(
  cd "$(dirname "${BASH_SOURCE[0]}")/../.." &&
  pwd
)"

cd "$PROJECT_ROOT"

DIMENSIONALITY_SUMMARY="${DIMENSIONALITY_SUMMARY:-results_test_repro/05_dimensionality/05_dimensionality_summary.csv}"
VARIABLE_FEATURES="${VARIABLE_FEATURES:-results_test_repro/05_dimensionality/05_variable_features_summary.csv}"
ELBOW_METRICS="${ELBOW_METRICS:-results_test_repro/05_dimensionality/05_pca_elbow_metrics.csv}"
OUTDIR="${OUTDIR:-results_test_repro}"
EXECUTION_ID="${EXECUTION_ID:-PC_SELECTION_POLICY_REPRO_001}"
TMP_DIR="$(mktemp -d)"

cleanup() {
  rm -rf "$TMP_DIR"
}

trap cleanup EXIT

echo "=========================================="
echo "PC-selection policy resolution test"
echo "=========================================="
echo "Project root:          $PROJECT_ROOT"
echo "Dimensionality summary: $DIMENSIONALITY_SUMMARY"
echo "Variable features:      $VARIABLE_FEATURES"
echo "Elbow metrics:          $ELBOW_METRICS"
echo "Output dir:             $OUTDIR"
echo "Execution ID:           $EXECUTION_ID"
echo

command -v nextflow >/dev/null 2>&1 || {
  echo "ERROR: nextflow is not installed or not in PATH." >&2
  exit 1
}

if [[ ! -f "$DIMENSIONALITY_SUMMARY" || ! -f "$VARIABLE_FEATURES" || ! -f "$ELBOW_METRICS" ]]; then
  echo "Dimensionality inputs not found; running dimensionality test first."
  tests/dimensionality/run_test.sh
fi

for required_file in "$DIMENSIONALITY_SUMMARY" "$VARIABLE_FEATURES" "$ELBOW_METRICS"; do
  if [[ ! -f "$required_file" ]]; then
    echo "ERROR: required input not found: $required_file" >&2
    exit 1
  fi
done

echo "Running Nextflow PC-selection policy resolution phase..."

nextflow run main.nf \
  --mode pc_selection_policy_resolution_phase \
  --dimensionality_summary "$DIMENSIONALITY_SUMMARY" \
  --variable_features_summary "$VARIABLE_FEATURES" \
  --pca_elbow_metrics "$ELBOW_METRICS" \
  --project_id M15214_TEST \
  --project_name CART_scRNAseq_TEST \
  --biological_run_id L13712_TEST \
  --execution_id "$EXECUTION_ID" \
  --outdir "$OUTDIR" \
  -profile docker,local_low_memory \
  -resume

OUTPUT_DIR="$OUTDIR/control/pc_selection_policy"

echo
echo "Validating generated PC-selection policy..."

python3 tests/pc_selection_policy/validate_pc_selection_policy.py \
  "$PROJECT_ROOT" \
  "$OUTPUT_DIR"

ruby -rcsv - "$DIMENSIONALITY_SUMMARY" "$TMP_DIR/bad_dimensionality_summary.csv" <<'RB'
source = ARGV.fetch(0)
target = ARGV.fetch(1)

rows = CSV.read(source, headers: true)
rows.each do |row|
  if row["metric"] == "pca_applied"
    row["value"] = "false"
  end
end

CSV.open(target, "w") do |csv|
  csv << rows.headers
  rows.each { |row| csv << row.fields }
end
RB

if python3 scripts/control/resolve_pc_selection_policy.py \
  --dimensionality-summary "$TMP_DIR/bad_dimensionality_summary.csv" \
  --variable-features "$VARIABLE_FEATURES" \
  --elbow-metrics "$ELBOW_METRICS" \
  --policy-yaml "$TMP_DIR/pc_selection_policy.yaml" \
  --policy-tsv "$TMP_DIR/pc_selection_policy.tsv" \
  --resolution-summary "$TMP_DIR/pc_selection_policy_resolution.summary.tsv" \
  --session-info "$TMP_DIR/pc_selection_policy_resolution.session.txt" \
  --project-id M15214_TEST \
  --project-name CART_scRNAseq_TEST \
  --biological-run-id L13712_TEST >/tmp/pc_selection_policy_bad.stdout 2>/tmp/pc_selection_policy_bad.stderr
then
  echo "Expected pca_applied=false to fail" >&2
  cat /tmp/pc_selection_policy_bad.stdout >&2
  exit 1
fi

echo "PC-selection negative policy check: VALID"

echo
echo "=========================================="
echo "PC-selection policy test completed successfully"
echo "=========================================="

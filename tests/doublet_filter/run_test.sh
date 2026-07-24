#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(
  cd "$(dirname "${BASH_SOURCE[0]}")/../.." &&
  pwd
)"

cd "$PROJECT_ROOT"

INPUT_OBJECT="${INPUT_OBJECT:-results_test_repro/01g_doublet_detection_scoring/seu_doublet_scored.rds}"
OUTDIR="${OUTDIR:-results_test_repro}"
EXECUTION_ID="${EXECUTION_ID:-DOUBLET_FILTER_REPRO_001}"
CONTAINER_IMAGE="${CONTAINER_IMAGE:-ghcr.io/jimmlucas/bioc-rnaseq-test01new2@sha256:b515735a14cfeb287ae628cfb3e9591492486afa521c9f0b2baba8fee2707b94}"
TMP_DIR="$(mktemp -d)"

export DOCKER_DEFAULT_PLATFORM="${DOCKER_DEFAULT_PLATFORM:-linux/amd64}"

cleanup() {
  rm -rf "$TMP_DIR"
}

trap cleanup EXIT

echo "=========================================="
echo "Doublet filter reproducibility test"
echo "=========================================="
echo "Project root:    $PROJECT_ROOT"
echo "Input object:    $INPUT_OBJECT"
echo "Output dir:      $OUTDIR"
echo "Execution ID:    $EXECUTION_ID"
echo "Container image: $CONTAINER_IMAGE"
echo

command -v docker >/dev/null 2>&1 || {
  echo "ERROR: docker is not installed or not in PATH." >&2
  exit 1
}

command -v nextflow >/dev/null 2>&1 || {
  echo "ERROR: nextflow is not installed or not in PATH." >&2
  exit 1
}

if [[ ! -f "$INPUT_OBJECT" ]]; then
  echo "Scored doublet object not found; running scoring fixture first."
  tests/doublet_scoring/run_test.sh
fi

if [[ ! -f "$INPUT_OBJECT" ]]; then
  echo "ERROR: input object not found after scoring: $INPUT_OBJECT" >&2
  exit 1
fi

APPROVED_POLICY="$TMP_DIR/approved_doublet_policy.yaml"

python3 scripts/control/validate_doublet_approval.py \
  --plan configs/analysis_plan_test.yaml \
  --policy tests/doublet_approval/fixtures/approved_policy.yaml \
  --validated-policy "$APPROVED_POLICY" \
  --summary "$TMP_DIR/doublet_policy_approval.summary.tsv" \
  --session-info "$TMP_DIR/doublet_policy_approval.session.txt"

echo
echo "Pulling the pinned R container image..."

docker pull "$CONTAINER_IMAGE"

echo
echo "Running Nextflow doublet filter phase..."

nextflow run main.nf \
  --mode doublet_filter_phase \
  --input_object "$INPUT_OBJECT" \
  --approved_doublet_policy "$APPROVED_POLICY" \
  --project_id M15214_TEST \
  --project_name CART_scRNAseq_TEST \
  --biological_run_id L13712_TEST \
  --execution_id "$EXECUTION_ID" \
  --outdir "$OUTDIR" \
  --r_image "$CONTAINER_IMAGE" \
  -profile docker,local_low_memory \
  -resume

OUTPUT_DIR="$OUTDIR/02_doublet_filter"

echo
echo "Validating generated outputs..."

docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v "$PROJECT_ROOT:$PROJECT_ROOT" \
  -w "$PROJECT_ROOT" \
  "$CONTAINER_IMAGE" \
  Rscript \
  tests/doublet_filter/validate_doublet_filter.R \
  "$INPUT_OBJECT" \
  "$OUTPUT_DIR"

python3 - "$PROJECT_ROOT" "$TMP_DIR" "$APPROVED_POLICY" <<'PY'
from __future__ import annotations

import copy
import importlib.util
import sys
from pathlib import Path

root = Path(sys.argv[1])
tmp = Path(sys.argv[2])
approved_policy = Path(sys.argv[3])
validator_path = root / "scripts/control/validate_doublet_approval.py"

spec = importlib.util.spec_from_file_location(
    "validate_doublet_approval",
    validator_path,
)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

approved = module.read_yaml(
    approved_policy,
    "Approved doublet policy",
)


def write_case(name: str, document: dict) -> None:
    module.write_yaml(tmp / f"{name}.yaml", document)


write_case("approved", approved)

proposed = copy.deepcopy(approved)
proposed["doublet_policy"]["approval"]["status"] = "proposed"
write_case("proposed", proposed)

invalid_validation = copy.deepcopy(approved)
invalid_validation["doublet_policy"]["validation"]["status"] = "invalid"
write_case("invalid_validation", invalid_validation)

already_applied = copy.deepcopy(approved)
already_applied["doublet_policy"]["application"]["status"] = "applied"
already_applied["doublet_policy"]["application"]["filtering_applied"] = True
write_case("already_applied", already_applied)

sample_missing = copy.deepcopy(approved)
sample_missing["doublet_policy"]["samples"].pop("wta-5")
write_case("sample_missing", sample_missing)

count_mismatch = copy.deepcopy(approved)
count_mismatch["doublet_policy"]["samples"]["wta-1"][
    "proposed_doublet_count"
] = 16
write_case("count_mismatch", count_mismatch)

invalid_action = copy.deepcopy(approved)
invalid_action["doublet_policy"]["samples"]["wta-1"][
    "application_action"
] = "free_text_action"
write_case("invalid_action", invalid_action)

keep_all = copy.deepcopy(approved)
for sample_policy in keep_all["doublet_policy"]["samples"].values():
    sample_policy["application_action"] = "keep_all"
write_case("keep_all", keep_all)

exclude = copy.deepcopy(approved)
for sample_policy in exclude["doublet_policy"]["samples"].values():
    sample_policy["application_action"] = "exclude_sample_from_filtering"
write_case("exclude_sample_from_filtering", exclude)
PY

run_filter_direct() {
  local case_name="$1"
  local case_out="$TMP_DIR/direct_${case_name}"
  mkdir -p "$case_out"

  docker run --rm \
    -u "$(id -u):$(id -g)" \
    -v "$PROJECT_ROOT:$PROJECT_ROOT" \
    -v "$TMP_DIR:$TMP_DIR" \
    -w "$case_out" \
    "$CONTAINER_IMAGE" \
    Rscript \
    "$PROJECT_ROOT/scripts/02_doublet_filter.R" \
    --input "$PROJECT_ROOT/$INPUT_OBJECT" \
    --policy "$TMP_DIR/${case_name}.yaml" \
    --output "seu_doublet_filtered.rds" \
    --project-id M15214_TEST \
    --biological-run-id L13712_TEST
}

expect_fail() {
  local case_name="$1"

  if run_filter_direct "$case_name" >"$TMP_DIR/${case_name}.stdout" 2>"$TMP_DIR/${case_name}.stderr"; then
    echo "Expected ${case_name} to fail" >&2
    cat "$TMP_DIR/${case_name}.stdout" >&2
    exit 1
  fi
}

expect_pass_with_counts() {
  local case_name="$1"
  local expected_removed="$2"

  run_filter_direct "$case_name" >"$TMP_DIR/${case_name}.stdout" 2>"$TMP_DIR/${case_name}.stderr"

  python3 - "$TMP_DIR/direct_${case_name}/02_doublet_filter_summary.csv" "$expected_removed" <<'PY'
import csv
import sys

summary = {}
with open(sys.argv[1], newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle):
        summary[row["metric"]] = row["value"]

expected_removed = int(sys.argv[2])
assert int(summary["cells_removed"]) == expected_removed
assert int(summary["cells_output"]) == 392 - expected_removed
assert int(summary["features_output"]) == 27522
PY
}

for case_name in \
  proposed \
  invalid_validation \
  already_applied \
  sample_missing \
  count_mismatch \
  invalid_action
do
  expect_fail "$case_name"
done

expect_pass_with_counts "keep_all" 0
expect_pass_with_counts "exclude_sample_from_filtering" 0

echo
echo "=========================================="
echo "Doublet filter test completed successfully"
echo "=========================================="

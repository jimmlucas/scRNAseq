#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(
  cd "$(dirname "${BASH_SOURCE[0]}")/../.." &&
  pwd
)"

cd "$PROJECT_ROOT"

INPUT_OBJECT="${INPUT_OBJECT:-results_test_repro/04_normalization/seu_normalized.rds}"
APPLIED_POLICY="${APPLIED_POLICY:-results_test_repro/04_normalization/applied_normalization_policy.yaml}"
OUTDIR="${OUTDIR:-results_test_repro}"
EXECUTION_ID="${EXECUTION_ID:-DIMENSIONALITY_REPRO_001}"
CONTAINER_IMAGE="${CONTAINER_IMAGE:-ghcr.io/jimmlucas/bioc-rnaseq-test01new2@sha256:b515735a14cfeb287ae628cfb3e9591492486afa521c9f0b2baba8fee2707b94}"
TMP_DIR="$(mktemp -d)"

export DOCKER_DEFAULT_PLATFORM="${DOCKER_DEFAULT_PLATFORM:-linux/amd64}"

cleanup() {
  rm -rf "$TMP_DIR"
}

trap cleanup EXIT

echo "=========================================="
echo "Dimensionality reproducibility test"
echo "=========================================="
echo "Project root:    $PROJECT_ROOT"
echo "Input object:    $INPUT_OBJECT"
echo "Applied policy:  $APPLIED_POLICY"
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

if [[ ! -f "$INPUT_OBJECT" || ! -f "$APPLIED_POLICY" ]]; then
  echo "Normalized inputs not found; running normalization test first."
  tests/normalization/run_test.sh
fi

if [[ ! -f "$INPUT_OBJECT" ]]; then
  echo "ERROR: input object not found: $INPUT_OBJECT" >&2
  exit 1
fi

if [[ ! -f "$APPLIED_POLICY" ]]; then
  echo "ERROR: applied normalization policy not found: $APPLIED_POLICY" >&2
  exit 1
fi

echo "Pulling the pinned R container image..."

docker pull "$CONTAINER_IMAGE"

echo
echo "Running Nextflow dimensionality phase..."

nextflow run main.nf \
  --mode dimensionality_phase \
  --input_object "$INPUT_OBJECT" \
  --applied_normalization_policy "$APPLIED_POLICY" \
  --project_id M15214_TEST \
  --project_name CART_scRNAseq_TEST \
  --biological_run_id L13712_TEST \
  --execution_id "$EXECUTION_ID" \
  --outdir "$OUTDIR" \
  --r_image "$CONTAINER_IMAGE" \
  --dimensionality_npcs 30 \
  -profile docker,local_low_memory \
  -resume

OUTPUT_DIR="$OUTDIR/05_dimensionality"

echo
echo "Validating generated outputs..."

docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v "$PROJECT_ROOT:$PROJECT_ROOT" \
  -w "$PROJECT_ROOT" \
  "$CONTAINER_IMAGE" \
  Rscript \
  tests/dimensionality/validate_dimensionality.R \
  "$INPUT_OBJECT" \
  "$OUTPUT_DIR"

ruby -ryaml - "$APPLIED_POLICY" "$TMP_DIR/not_applied_policy.yaml" <<'RB'
source = ARGV.fetch(0)
target = ARGV.fetch(1)

document = YAML.load_file(source)

application = document.fetch("normalization_policy").fetch("application")
application["status"] = "not_applied"
application["normalization_applied"] = false

File.write(target, document.to_yaml)
RB

echo
echo "Validating negative policy check..."

if docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v "$PROJECT_ROOT:$PROJECT_ROOT" \
  -v "$TMP_DIR:$TMP_DIR" \
  -w "$TMP_DIR" \
  "$CONTAINER_IMAGE" \
  Rscript \
  "$PROJECT_ROOT/scripts/05_dimensionality.R" \
  --input "$PROJECT_ROOT/$INPUT_OBJECT" \
  --policy "$TMP_DIR/not_applied_policy.yaml" \
  --output "$TMP_DIR/should_not_exist.rds" \
  --project-id M15214_TEST \
  --biological-run-id L13712_TEST; then
  echo "ERROR: dimensionality accepted a not_applied normalization policy" >&2
  exit 1
fi

echo "Dimensionality negative policy check: VALID"

echo
echo "=========================================="
echo "Dimensionality reproducibility test completed successfully"
echo "=========================================="

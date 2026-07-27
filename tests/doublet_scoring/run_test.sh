#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(
  cd "$(dirname "${BASH_SOURCE[0]}")/../.." &&
  pwd
)"

cd "$PROJECT_ROOT"

INPUT_OBJECT="${INPUT_OBJECT:-tests/doublet_scoring/fixtures/seu_qc_filtered.rds}"
OUTDIR="${OUTDIR:-results_test_repro}"
EXECUTION_ID="${EXECUTION_ID:-DOUBLET_SCORING_REPRO_001}"

CONTAINER_IMAGE="${CONTAINER_IMAGE:-ghcr.io/jimmlucas/bioc-rnaseq-test01new2@sha256:b515735a14cfeb287ae628cfb3e9591492486afa521c9f0b2baba8fee2707b94}"

export DOCKER_DEFAULT_PLATFORM="${DOCKER_DEFAULT_PLATFORM:-linux/amd64}"

echo "=========================================="
echo "Doublet scoring reproducibility test"
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
  echo "ERROR: input object not found: $INPUT_OBJECT" >&2
  exit 1
fi

echo "Verifying fixture checksum..."

sha256sum -c \
  tests/doublet_scoring/expected/checksums.sha256

echo
echo "Pulling the pinned container image..."

docker pull "$CONTAINER_IMAGE"

echo
echo "Running Nextflow doublet scoring phase..."

nextflow run main.nf \
  --mode doublet_detection_scoring_phase \
  --input_object "$INPUT_OBJECT" \
  --doublet_config configs/project_config_doublet_test.yaml \
  --project_id M15214_TEST \
  --project_name CART_scRNAseq_TEST \
  --biological_run_id L13712_TEST \
  --execution_id "$EXECUTION_ID" \
  --outdir "$OUTDIR" \
  --r_image "$CONTAINER_IMAGE" \
  -profile docker,local_low_memory \
  -resume

OUTPUT_DIR="$OUTDIR/01g_doublet_detection_scoring"

echo
echo "Validating generated outputs..."

docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v "$PROJECT_ROOT:$PROJECT_ROOT" \
  -w "$PROJECT_ROOT" \
  "$CONTAINER_IMAGE" \
  Rscript \
  tests/doublet_scoring/validate_doublet_scoring.R \
  "$INPUT_OBJECT" \
  "$OUTPUT_DIR"

echo
echo "=========================================="
echo "Reproducibility test completed successfully"
echo "=========================================="

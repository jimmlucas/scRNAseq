#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(
  cd "$(dirname "${BASH_SOURCE[0]}")/../.." &&
  pwd
)"

cd "$PROJECT_ROOT"

INPUT_OBJECT="${INPUT_OBJECT:-results_test_repro/02_doublet_filter/seu_doublet_filtered.rds}"
APPLIED_POLICY="${APPLIED_POLICY:-results_test_repro/02_doublet_filter/applied_doublet_policy.yaml}"
OUTDIR="${OUTDIR:-results_test_repro}"
EXECUTION_ID="${EXECUTION_ID:-POST_DOUBLET_QC_REPRO_001}"
CONTAINER_IMAGE="${CONTAINER_IMAGE:-ghcr.io/jimmlucas/bioc-rnaseq-test01new2@sha256:b515735a14cfeb287ae628cfb3e9591492486afa521c9f0b2baba8fee2707b94}"

export DOCKER_DEFAULT_PLATFORM="${DOCKER_DEFAULT_PLATFORM:-linux/amd64}"

echo "=========================================="
echo "Post-doublet QC readiness test"
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
  echo "Doublet-filtered inputs not found; running doublet filter test first."
  tests/doublet_filter/run_test.sh
fi

if [[ ! -f "$INPUT_OBJECT" ]]; then
  echo "ERROR: input object not found: $INPUT_OBJECT" >&2
  exit 1
fi

if [[ ! -f "$APPLIED_POLICY" ]]; then
  echo "ERROR: applied policy not found: $APPLIED_POLICY" >&2
  exit 1
fi

echo "Pulling the pinned R container image..."

docker pull "$CONTAINER_IMAGE"

echo
echo "Running Nextflow post-doublet QC summary phase..."

nextflow run main.nf \
  --mode post_doublet_qc_summary_phase \
  --input_object "$INPUT_OBJECT" \
  --applied_doublet_policy "$APPLIED_POLICY" \
  --project_id M15214_TEST \
  --project_name CART_scRNAseq_TEST \
  --biological_run_id L13712_TEST \
  --execution_id "$EXECUTION_ID" \
  --outdir "$OUTDIR" \
  --r_image "$CONTAINER_IMAGE" \
  -profile docker,local_low_memory \
  -resume

OUTPUT_DIR="$OUTDIR/03_post_doublet_qc"

echo
echo "Validating generated readiness outputs..."

docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v "$PROJECT_ROOT:$PROJECT_ROOT" \
  -w "$PROJECT_ROOT" \
  "$CONTAINER_IMAGE" \
  Rscript \
  tests/post_doublet_qc/validate_post_doublet_qc.R \
  "$OUTPUT_DIR"

echo
echo "=========================================="
echo "Post-doublet QC readiness test completed successfully"
echo "=========================================="

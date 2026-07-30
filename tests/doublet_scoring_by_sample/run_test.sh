#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

cd "${ROOT_DIR}"

mkdir -p \
  "${TMP_DIR}/sample_objects/sample_a" \
  "${TMP_DIR}/sample_objects/sample_b" \
  "${TMP_DIR}/config"

touch \
  "${TMP_DIR}/sample_objects/sample_a/seu_qc_filtered.sample_a.rds" \
  "${TMP_DIR}/sample_objects/sample_b/seu_qc_filtered.sample_b.rds"

cat > "${TMP_DIR}/manifest.tsv" <<'TSV'
sample	sample_safe_id	sample_column	cells	features	object_path
sample_a	sample_a	sample	0	0	sample_objects/sample_a/seu_qc_filtered.sample_a.rds
sample_b	sample_b	sample	0	0	sample_objects/sample_b/seu_qc_filtered.sample_b.rds
TSV

cat > "${TMP_DIR}/config/doublet_config.yaml" <<'YAML'
project:
  id: M15214_TEST
doublet:
  expected_rate: 0.075
  pN: 0.25
  pcs: 10
  min_cells: 50
  nfeatures: 1000
  sweep_cores: 1
YAML

Rscript -e "invisible(parse(file='scripts/01g_merge_doublet_detection_scoring.R')); cat('Merge script parse: VALID\\n')"

nextflow run main.nf \
  --mode doublet_detection_scoring_by_sample_phase \
  --sample_manifest "${TMP_DIR}/manifest.tsv" \
  --sample_objects "${TMP_DIR}/sample_objects" \
  --doublet_config "${TMP_DIR}/config/doublet_config.yaml" \
  --project_id M15214_TEST \
  --project_name CART_scRNAseq_TEST \
  --biological_run_id L13712_TEST \
  --execution_id DOUBLET_BY_SAMPLE_STUB_001 \
  --outdir "${TMP_DIR}/results" \
  -profile docker,local_low_memory \
  -stub-run

test -f "${TMP_DIR}/results/01g_doublet_detection_scoring/seu_doublet_scored.rds"
test -f "${TMP_DIR}/results/01g_doublet_detection_scoring/01g_doublet_detection_summary.csv"
test -f "${TMP_DIR}/results/01g_doublet_detection_scoring/01g_doublet_detection_by_sample.csv"
test -f "${TMP_DIR}/results/01g_doublet_detection_scoring/01g_doublet_detection_cell_scores.csv"
test -f "${TMP_DIR}/results/01g_doublet_detection_scoring/01g_doublet_detection_pk_sweep.csv"
test -f "${TMP_DIR}/results/01g_doublet_detection_scoring/01g_doublet_detection_performance.csv"

test -f "${TMP_DIR}/results/01g_doublet_detection_scoring/per_sample/sample_a/seu_doublet_scored.sample_a.rds"
test -f "${TMP_DIR}/results/01g_doublet_detection_scoring/per_sample/sample_b/seu_doublet_scored.sample_b.rds"

echo "Doublet scoring by sample tests: VALID"

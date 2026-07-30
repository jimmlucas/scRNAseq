#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

R_BIN="${R_BIN:-Rscript}"

cd "${ROOT_DIR}"

"${R_BIN}" -e "invisible(parse(file='scripts/01f_split_seurat_by_sample.R')); cat('Split script parse: VALID\\n')"

nextflow run main.nf \
  --mode split_seurat_by_sample_phase \
  --input_object tests/doublet_scoring/fixtures/seu_qc_filtered.rds \
  --project_id M15214_TEST \
  --project_name CART_scRNAseq_TEST \
  --biological_run_id L13712_TEST \
  --execution_id SPLIT_SEURAT_BY_SAMPLE_STUB_001 \
  --outdir "${TMP_DIR}/results" \
  -profile docker,local_low_memory \
  -stub-run

test -d "${TMP_DIR}/results/01f_split_seurat_by_sample/sample_objects"
test -f "${TMP_DIR}/results/01f_split_seurat_by_sample/01f_split_by_sample_manifest.tsv"
test -f "${TMP_DIR}/results/01f_split_seurat_by_sample/01f_split_by_sample_summary.csv"
test -f "${TMP_DIR}/results/01f_split_seurat_by_sample/01f_sessionInfo.txt"

grep -q $'sample\tsample_safe_id\tsample_column\tcells\tfeatures\tobject_path' \
  "${TMP_DIR}/results/01f_split_seurat_by_sample/01f_split_by_sample_manifest.tsv"

echo "Split Seurat by sample tests: VALID"

#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(
  cd "$(dirname "${BASH_SOURCE[0]}")/../.." &&
  pwd
)"

cd "$PROJECT_ROOT"

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

command -v nextflow >/dev/null 2>&1 || {
  echo "ERROR: nextflow is not installed or not in PATH." >&2
  exit 1
}

if rg -n 'params\.(r_image|python_image)\.docker' modules; then
  echo "ERROR: container configuration treats an image string as a map." >&2
  exit 1
fi

BBUDY_SLURM_QUEUE=compute \
BBUDY_SLURM_ACCOUNT=scrna \
BBUDY_WORK_DIR="$TMP_DIR/work" \
BBUDY_QUEUE_SIZE=7 \
BBUDY_DOUBLET_CPUS=1 \
  nextflow config -profile slurm . >"$TMP_DIR/slurm.config"

nextflow config -profile singularity . >"$TMP_DIR/singularity.config"
nextflow config -profile apptainer . >"$TMP_DIR/apptainer.config"

grep -q "executor = 'slurm'" "$TMP_DIR/slurm.config"
grep -q "queue = 'compute'" "$TMP_DIR/slurm.config"
grep -q "clusterOptions = '--account=scrna'" "$TMP_DIR/slurm.config"
grep -q "queueSize = 7" "$TMP_DIR/slurm.config"
grep -q "workDir = '$TMP_DIR/work'" "$TMP_DIR/slurm.config"
grep -A3 "DOUBLET_DETECTION_SCORING" "$TMP_DIR/slurm.config" \
  | grep -q "cpus = 1"

if grep -q "maxForks" "$TMP_DIR/slurm.config"; then
  echo "ERROR: SLURM profile unexpectedly limits process parallelism." >&2
  exit 1
fi

grep -A4 "withLabel:very_high_memory" "$TMP_DIR/slurm.config" \
  | grep -q "72.h"
grep -q "singularity" "$TMP_DIR/singularity.config"
grep -q "apptainer" "$TMP_DIR/apptainer.config"

for config_file in \
  "$TMP_DIR/slurm.config" \
  "$TMP_DIR/singularity.config" \
  "$TMP_DIR/apptainer.config"
do
  grep -q \
    "container = 'docker://ghcr.io/jimmlucas/scrnaseq-control-python@sha256:" \
    "$config_file"
done

if NXF_OFFLINE=true nextflow run main.nf --mode compile_check \
  >"$TMP_DIR/compile.stdout" \
  2>"$TMP_DIR/compile.stderr"
then
  echo "ERROR: the compile-only invalid mode unexpectedly succeeded." >&2
  exit 1
fi

if grep -q "Script compilation failed" "$TMP_DIR/compile.stderr"; then
  cat "$TMP_DIR/compile.stderr" >&2
  exit 1
fi

grep -q "Invalid mode: compile_check" "$TMP_DIR/compile.stdout"

grep -q -- '--cores ${task.cpus}' \
  modules/doublet_detection_scoring.nf
grep -q -- '--cores ${task.cpus}' \
  modules/doublet_detection_score_sample.nf

echo "Server execution profiles: VALID"

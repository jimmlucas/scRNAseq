#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(
  cd "$(dirname "${BASH_SOURCE[0]}")/../.." &&
  pwd
)"

cd "$PROJECT_ROOT"

INPUT_OBJECT="${INPUT_OBJECT:-results_test_repro/02_doublet_filter/seu_doublet_filtered.rds}"
OUTDIR="${OUTDIR:-results_test_repro}"
EXECUTION_ID="${EXECUTION_ID:-NORMALIZATION_REPRO_001}"
CONTAINER_IMAGE="${CONTAINER_IMAGE:-ghcr.io/jimmlucas/bioc-rnaseq-test01new2@sha256:b515735a14cfeb287ae628cfb3e9591492486afa521c9f0b2baba8fee2707b94}"
TMP_DIR="$(mktemp -d)"

export DOCKER_DEFAULT_PLATFORM="${DOCKER_DEFAULT_PLATFORM:-linux/amd64}"

cleanup() {
  rm -rf "$TMP_DIR"
}

trap cleanup EXIT

echo "=========================================="
echo "Normalization reproducibility test"
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
  echo "Doublet-filtered object not found; running doublet filter fixture first."
  tests/doublet_filter/run_test.sh
fi

if [[ ! -f "$INPUT_OBJECT" ]]; then
  echo "ERROR: input object not found after doublet filter: $INPUT_OBJECT" >&2
  exit 1
fi

PROPOSED_POLICY="$OUTDIR/control/normalization_policy/normalization_policy.yaml"

if [[ ! -f "$PROPOSED_POLICY" ]]; then
  echo "Normalization policy not found; running policy resolver fixture first."
  tests/normalization_policy/run_test.sh
fi

if [[ ! -f "$PROPOSED_POLICY" ]]; then
  echo "ERROR: proposed normalization policy not found: $PROPOSED_POLICY" >&2
  exit 1
fi

REVIEWED_POLICY="$TMP_DIR/reviewed_normalization_policy.yaml"
APPROVED_POLICY="$TMP_DIR/approved_normalization_policy.yaml"

python3 - "$PROJECT_ROOT" "$PROPOSED_POLICY" "$REVIEWED_POLICY" <<'PY'
from __future__ import annotations

import copy
import importlib.util
import sys
from pathlib import Path

root = Path(sys.argv[1])
proposed_policy = Path(sys.argv[2])
reviewed_policy = Path(sys.argv[3])
validator_path = root / "scripts/control/validate_normalization_approval.py"

spec = importlib.util.spec_from_file_location(
    "validate_normalization_approval",
    validator_path,
)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

document = module.read_yaml(
    proposed_policy,
    "Proposed normalization policy",
)
document = copy.deepcopy(document)
approval = document["normalization_policy"]["approval"]
approval["status"] = "approved"
approval["reviewed_by"] = "jimmy"
approval["approved_at"] = "2026-07-27T12:00:00+02:00"

module.write_yaml(reviewed_policy, document)
PY

python3 scripts/control/validate_normalization_approval.py \
  --plan configs/analysis_plan_test.yaml \
  --policy "$REVIEWED_POLICY" \
  --validated-policy "$APPROVED_POLICY" \
  --summary "$TMP_DIR/normalization_policy_approval.summary.tsv" \
  --session-info "$TMP_DIR/normalization_policy_approval.session.txt"

echo
echo "Pulling the pinned R container image..."

docker pull "$CONTAINER_IMAGE"

echo
echo "Running Nextflow normalization phase..."

nextflow run main.nf \
  --mode normalization_phase \
  --input_object "$INPUT_OBJECT" \
  --approved_normalization_policy "$APPROVED_POLICY" \
  --project_id M15214_TEST \
  --project_name CART_scRNAseq_TEST \
  --biological_run_id L13712_TEST \
  --execution_id "$EXECUTION_ID" \
  --outdir "$OUTDIR" \
  --r_image "$CONTAINER_IMAGE" \
  -profile docker,local_low_memory \
  -resume

OUTPUT_DIR="$OUTDIR/04_normalization"

echo
echo "Validating generated outputs..."

docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v "$PROJECT_ROOT:$PROJECT_ROOT" \
  -w "$PROJECT_ROOT" \
  "$CONTAINER_IMAGE" \
  Rscript \
  tests/normalization/validate_normalization.R \
  "$INPUT_OBJECT" \
  "$OUTPUT_DIR"

python3 - "$PROJECT_ROOT" "$TMP_DIR" "$APPROVED_POLICY" "$INPUT_OBJECT" <<'PY'
from __future__ import annotations

import copy
import importlib.util
import subprocess
import sys
from pathlib import Path

root = Path(sys.argv[1])
tmp = Path(sys.argv[2])
approved_policy = Path(sys.argv[3])
input_object = Path(sys.argv[4])
validator_path = root / "scripts/control/validate_normalization_approval.py"
script = root / "scripts/04_normalization.R"

spec = importlib.util.spec_from_file_location(
    "validate_normalization_approval",
    validator_path,
)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

approved = module.read_yaml(
    approved_policy,
    "Approved normalization policy",
)


def write_case(name: str, document: dict) -> Path:
    path = tmp / f"{name}.yaml"
    module.write_yaml(path, document)
    return path


cases = {}

proposed = copy.deepcopy(approved)
proposed["normalization_policy"]["approval"]["status"] = "proposed"
cases["proposed"] = write_case("proposed", proposed)

invalid_validation = copy.deepcopy(approved)
invalid_validation["normalization_policy"]["validation"]["status"] = "invalid"
cases["invalid_validation"] = write_case(
    "invalid_validation",
    invalid_validation,
)

already_applied = copy.deepcopy(approved)
already_applied["normalization_policy"]["application"]["status"] = "applied"
already_applied["normalization_policy"]["application"][
    "normalization_applied"
] = True
cases["already_applied"] = write_case("already_applied", already_applied)

wrong_project = copy.deepcopy(approved)
wrong_project["normalization_policy"]["identity"]["project_id"] = "WRONG"
cases["wrong_project"] = write_case("wrong_project", wrong_project)

invalid_scope = copy.deepcopy(approved)
invalid_scope["normalization_policy"]["strategy"]["scope"] = "global"
cases["invalid_scope"] = write_case("invalid_scope", invalid_scope)

for case_name, policy_path in cases.items():
    out_dir = tmp / f"direct_{case_name}"
    out_dir.mkdir()
    command = [
        "docker",
        "run",
        "--rm",
        "-u",
        f"{subprocess.check_output(['id', '-u'], text=True).strip()}:"
        f"{subprocess.check_output(['id', '-g'], text=True).strip()}",
        "-v",
        f"{root}:{root}",
        "-v",
        f"{tmp}:{tmp}",
        "-w",
        str(out_dir),
        "ghcr.io/jimmlucas/bioc-rnaseq-test01new2@sha256:"
        "b515735a14cfeb287ae628cfb3e9591492486afa521c9f0b2baba8fee2707b94",
        "Rscript",
        str(script),
        "--input",
        str(root / input_object),
        "--policy",
        str(policy_path),
        "--output",
        "seu_normalized.rds",
        "--project-id",
        "M15214_TEST",
        "--biological-run-id",
        "L13712_TEST",
    ]
    result = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode == 0:
        print(f"Expected {case_name} to fail", file=sys.stderr)
        print(result.stdout, file=sys.stderr)
        raise SystemExit(1)

print("Normalization negative policy checks: VALID")
PY

echo
echo "Normalization reproducibility test: VALID"

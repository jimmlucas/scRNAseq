#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TMP_DIR="$(mktemp -d)"

cleanup() {
    rm -rf "${TMP_DIR}"
}

trap cleanup EXIT

VALIDATOR="${ROOT_DIR}/scripts/control/validate_normalization_approval.py"
PLAN="${ROOT_DIR}/configs/analysis_plan_test.yaml"
PROPOSED_POLICY="${ROOT_DIR}/results_test_repro/control/normalization_policy/normalization_policy.yaml"

python3 - "${ROOT_DIR}" "${TMP_DIR}" <<'PY'
from __future__ import annotations

import copy
import importlib.util
import sys
from pathlib import Path

root = Path(sys.argv[1])
tmp = Path(sys.argv[2])
validator_path = root / "scripts/control/validate_normalization_approval.py"

spec = importlib.util.spec_from_file_location(
    "validate_normalization_approval",
    validator_path,
)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

proposed = module.read_yaml(
    root
    / "results_test_repro/control/normalization_policy/normalization_policy.yaml",
    "Proposed normalization policy",
)


def reviewed(document: dict, status: str = "approved") -> dict:
    candidate = copy.deepcopy(document)
    approval = candidate["normalization_policy"]["approval"]
    approval["status"] = status
    approval["reviewed_by"] = "jimmy"
    approval["approved_at"] = "2026-07-27T12:00:00+02:00"
    return candidate


def write_case(name: str, document: dict) -> None:
    module.write_yaml(tmp / f"{name}.yaml", document)


approved = reviewed(proposed)
write_case("approved", approved)
write_case("proposed", proposed)

modified_and_approved = reviewed(proposed, "modified_and_approved")
modified_and_approved["normalization_policy"]["integration"][
    "required"
] = False
modified_and_approved["normalization_policy"]["integration"][
    "override_justification"
] = "pilot single-condition run; defer integration until cohort expansion"
write_case("modified_and_approved_with_override", modified_and_approved)

sctransform = reviewed(proposed)
sctransform["normalization_policy"]["strategy"]["method"] = "SCTransform"
sctransform["normalization_policy"]["strategy"]["regressors"] = "percent_mito"
write_case("sctransform", sctransform)

missing_reviewer = copy.deepcopy(approved)
missing_reviewer["normalization_policy"]["approval"]["reviewed_by"] = None
write_case("missing_reviewer", missing_reviewer)

missing_approved_at = copy.deepcopy(approved)
missing_approved_at["normalization_policy"]["approval"]["approved_at"] = None
write_case("missing_approved_at", missing_approved_at)

invalid_method = copy.deepcopy(approved)
invalid_method["normalization_policy"]["strategy"]["method"] = "MagicNormalize"
write_case("invalid_method", invalid_method)

invalid_scope = copy.deepcopy(approved)
invalid_scope["normalization_policy"]["strategy"]["scope"] = "global"
write_case("invalid_scope", invalid_scope)

invalid_nfeatures = copy.deepcopy(approved)
invalid_nfeatures["normalization_policy"]["strategy"]["nfeatures"] = 0
write_case("invalid_nfeatures", invalid_nfeatures)

invalid_log_regressor = copy.deepcopy(approved)
invalid_log_regressor["normalization_policy"]["strategy"][
    "regressors"
] = "percent_mito"
write_case("invalid_log_regressor", invalid_log_regressor)

integration_override_missing = copy.deepcopy(approved)
integration_override_missing["normalization_policy"]["integration"][
    "required"
] = False
write_case("integration_override_missing", integration_override_missing)

already_applied = copy.deepcopy(approved)
already_applied["normalization_policy"]["application"]["status"] = "applied"
already_applied["normalization_policy"]["application"][
    "normalization_applied"
] = True
write_case("already_applied", already_applied)

wrong_project = copy.deepcopy(approved)
wrong_project["normalization_policy"]["identity"]["project_id"] = "WRONG"
write_case("wrong_project", wrong_project)
PY

run_validator() {
    local case_name="$1"
    local out_dir="${TMP_DIR}/out_${case_name}"

    python3 "${VALIDATOR}" \
        --plan "${PLAN}" \
        --policy "${TMP_DIR}/${case_name}.yaml" \
        --validated-policy "${out_dir}/approved_normalization_policy.yaml" \
        --summary "${out_dir}/normalization_policy_approval.summary.tsv" \
        --session-info "${out_dir}/normalization_policy_approval.session.txt"
}

expect_pass() {
    local case_name="$1"

    if ! run_validator "${case_name}" >"${TMP_DIR}/${case_name}.stdout" 2>"${TMP_DIR}/${case_name}.stderr"; then
        echo "Expected ${case_name} to pass" >&2
        cat "${TMP_DIR}/${case_name}.stderr" >&2
        exit 1
    fi
}

expect_fail() {
    local case_name="$1"

    if run_validator "${case_name}" >"${TMP_DIR}/${case_name}.stdout" 2>"${TMP_DIR}/${case_name}.stderr"; then
        echo "Expected ${case_name} to fail" >&2
        cat "${TMP_DIR}/${case_name}.stdout" >&2
        exit 1
    fi
}

expect_pass "approved"
expect_pass "modified_and_approved_with_override"
expect_pass "sctransform"

for case_name in \
    proposed \
    missing_reviewer \
    missing_approved_at \
    invalid_method \
    invalid_scope \
    invalid_nfeatures \
    invalid_log_regressor \
    integration_override_missing \
    already_applied \
    wrong_project
do
    expect_fail "${case_name}"
done

python3 - "${ROOT_DIR}" "${TMP_DIR}" <<'PY'
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

root = Path(sys.argv[1])
tmp = Path(sys.argv[2])
validator_path = root / "scripts/control/validate_normalization_approval.py"

spec = importlib.util.spec_from_file_location(
    "validate_normalization_approval",
    validator_path,
)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

document = module.read_yaml(
    tmp / "out_approved/approved_normalization_policy.yaml",
    "Validated normalization policy",
)
policy = document["normalization_policy"]

assert policy["approval"]["status"] == "approved"
assert policy["validation"]["status"] == "valid"
assert policy["validation"]["validator"] == "validate_normalization_approval.py"
assert policy["validation"]["validated_for_normalization"] is True
assert policy["application"]["status"] == "not_applied"
assert policy["application"]["normalization_applied"] is False
assert policy["strategy"]["method"] == "LogNormalize"
assert policy["strategy"]["variable_features_method"] == "vst"
assert int(policy["strategy"]["nfeatures"]) == 2000
assert float(policy["strategy"]["scale_factor"]) == 10000
assert policy["strategy"]["scope"] == "per_sample"
assert policy["strategy"]["assay"] == "RNA"
assert policy["strategy"]["regressors"] == "none"
assert policy["integration"]["required"] is True

summary = (
    tmp
    / "out_approved/normalization_policy_approval.summary.tsv"
).read_text(encoding="utf-8").splitlines()
assert len(summary) == 2
assert "validation_status" in summary[0]
assert "\tvalid\tapproved\t" in summary[1]

print("Normalization approval validation tests: VALID")
PY

nextflow run "${ROOT_DIR}/main.nf" \
    -stub-run \
    --mode normalization_approval_phase \
    --analysis_plan "${PLAN}" \
    --approved_normalization_policy "${TMP_DIR}/approved.yaml" \
    --project_id M15214_TEST \
    --project_name CART_scRNAseq_TEST \
    --biological_run_id L13712_TEST \
    --outdir "${TMP_DIR}/nf_stub"

echo "Normalization approval phase stub run: VALID"

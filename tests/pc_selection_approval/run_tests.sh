#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TMP_DIR="$(mktemp -d)"

cleanup() {
    rm -rf "${TMP_DIR}"
}

trap cleanup EXIT

VALIDATOR="${ROOT_DIR}/scripts/control/validate_pc_selection_approval.py"
PLAN="${ROOT_DIR}/configs/analysis_plan_test.yaml"

python3 - "${ROOT_DIR}" "${TMP_DIR}" <<'PY'
from __future__ import annotations

import copy
import importlib.util
import sys
from pathlib import Path

root = Path(sys.argv[1])
tmp = Path(sys.argv[2])
validator_path = root / "scripts/control/validate_pc_selection_approval.py"

spec = importlib.util.spec_from_file_location(
    "validate_pc_selection_approval",
    validator_path,
)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

proposed = module.read_yaml(
    root / "results_test_repro/control/pc_selection_policy/pc_selection_policy.yaml",
    "Proposed PC-selection policy",
)


def reviewed(document: dict, status: str = "approved") -> dict:
    candidate = copy.deepcopy(document)
    approval = candidate["pc_selection_policy"]["approval"]
    approval["status"] = status
    approval["reviewed_by"] = "jimmy"
    approval["approved_at"] = "2026-07-27T12:00:00+02:00"
    return candidate


def write_case(name: str, document: dict) -> None:
    module.write_yaml(tmp / f"{name}.yaml", document)


approved = reviewed(proposed)
write_case("approved", approved)
write_case("proposed", proposed)

modified = reviewed(proposed, "modified_and_approved")
modified["pc_selection_policy"]["selected_pcs"] = 18
write_case("modified_and_approved", modified)

changed_without_modified_status = reviewed(proposed)
changed_without_modified_status["pc_selection_policy"]["selected_pcs"] = 18
write_case("changed_without_modified_status", changed_without_modified_status)

missing_reviewer = copy.deepcopy(approved)
missing_reviewer["pc_selection_policy"]["approval"]["reviewed_by"] = None
write_case("missing_reviewer", missing_reviewer)

missing_approved_at = copy.deepcopy(approved)
missing_approved_at["pc_selection_policy"]["approval"]["approved_at"] = None
write_case("missing_approved_at", missing_approved_at)

zero_selected = copy.deepcopy(modified)
zero_selected["pc_selection_policy"]["selected_pcs"] = 0
write_case("zero_selected", zero_selected)

too_many_selected = copy.deepcopy(modified)
too_many_selected["pc_selection_policy"]["selected_pcs"] = 31
write_case("too_many_selected", too_many_selected)

bad_recommendation = copy.deepcopy(approved)
bad_recommendation["pc_selection_policy"]["recommendation"][
    "recommended_pcs"
] = 31
write_case("bad_recommendation", bad_recommendation)

already_applied = copy.deepcopy(approved)
already_applied["pc_selection_policy"]["application"]["status"] = "applied"
already_applied["pc_selection_policy"]["application"][
    "pc_selection_applied"
] = True
write_case("already_applied", already_applied)

wrong_project = copy.deepcopy(approved)
wrong_project["pc_selection_policy"]["identity"]["project_id"] = "WRONG"
write_case("wrong_project", wrong_project)
PY

run_validator() {
    local case_name="$1"
    local out_dir="${TMP_DIR}/out_${case_name}"

    python3 "${VALIDATOR}" \
        --plan "${PLAN}" \
        --policy "${TMP_DIR}/${case_name}.yaml" \
        --validated-policy "${out_dir}/approved_pc_selection_policy.yaml" \
        --summary "${out_dir}/pc_selection_policy_approval.summary.tsv" \
        --session-info "${out_dir}/pc_selection_policy_approval.session.txt"
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
expect_pass "modified_and_approved"

for case_name in \
    proposed \
    changed_without_modified_status \
    missing_reviewer \
    missing_approved_at \
    zero_selected \
    too_many_selected \
    bad_recommendation \
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
validator_path = root / "scripts/control/validate_pc_selection_approval.py"

spec = importlib.util.spec_from_file_location(
    "validate_pc_selection_approval",
    validator_path,
)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

document = module.read_yaml(
    tmp / "out_modified_and_approved/approved_pc_selection_policy.yaml",
    "Validated PC-selection policy",
)
policy = document["pc_selection_policy"]

assert policy["proposal"]["recommended_pcs"] == 21
assert policy["selected_pcs"] == 18
assert policy["approval"]["status"] == "modified_and_approved"
assert policy["validation"]["status"] == "valid"
assert policy["validation"]["validator"] == "validate_pc_selection_approval.py"
assert policy["validation"]["validated_for_pc_selection"] is True
assert policy["validation"]["proposed_pcs"] == 21
assert policy["validation"]["selected_pcs"] == 18
assert policy["validation"]["max_computed_pcs"] == 30
assert policy["application"]["status"] == "not_applied"
assert policy["application"]["pc_selection_applied"] is False

summary = (
    tmp
    / "out_modified_and_approved/pc_selection_policy_approval.summary.tsv"
).read_text(encoding="utf-8").splitlines()
assert len(summary) == 2
assert "selected_pcs" in summary[0]
assert "\tvalid\tmodified_and_approved\t" in summary[1]

print("PC-selection approval validation tests: VALID")
PY

nextflow run "${ROOT_DIR}/main.nf" \
    -stub-run \
    --mode pc_selection_approval_phase \
    --analysis_plan "${PLAN}" \
    --approved_pc_selection_policy "${TMP_DIR}/approved.yaml" \
    --project_id M15214_TEST \
    --project_name CART_scRNAseq_TEST \
    --biological_run_id L13712_TEST \
    --outdir "${TMP_DIR}/nf_stub"

echo "PC-selection approval phase stub run: VALID"

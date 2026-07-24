#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TMP_DIR="$(mktemp -d)"

cleanup() {
    rm -rf "${TMP_DIR}"
}

trap cleanup EXIT

VALIDATOR="${ROOT_DIR}/scripts/control/validate_doublet_approval.py"
PLAN="${ROOT_DIR}/configs/analysis_plan_test.yaml"
APPROVED_FIXTURE="${ROOT_DIR}/tests/doublet_approval/fixtures/approved_policy.yaml"

python3 - "${ROOT_DIR}" "${TMP_DIR}" <<'PY'
from __future__ import annotations

import copy
import importlib.util
import sys
from pathlib import Path

root = Path(sys.argv[1])
tmp = Path(sys.argv[2])
validator_path = root / "scripts/control/validate_doublet_approval.py"

spec = importlib.util.spec_from_file_location(
    "validate_doublet_approval",
    validator_path,
)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

approved = module.read_yaml(
    root / "tests/doublet_approval/fixtures/approved_policy.yaml",
    "Approved policy fixture",
)


def write_case(name: str, document: dict) -> None:
    module.write_yaml(tmp / f"{name}.yaml", document)


write_case("approved", approved)

proposed = copy.deepcopy(approved)
proposed["doublet_policy"]["approval"]["status"] = "proposed"
proposed["doublet_policy"]["approval"]["reviewed_by"] = None
proposed["doublet_policy"]["approval"]["approved_at"] = None
for sample in proposed["doublet_policy"]["samples"].values():
    sample["decision_status"] = "unreviewed"
write_case("proposed", proposed)

missing_reviewer = copy.deepcopy(approved)
missing_reviewer["doublet_policy"]["approval"]["reviewed_by"] = None
write_case("missing_reviewer", missing_reviewer)

missing_approved_at = copy.deepcopy(approved)
missing_approved_at["doublet_policy"]["approval"]["approved_at"] = None
write_case("missing_approved_at", missing_approved_at)

missing_sample = copy.deepcopy(approved)
missing_sample["doublet_policy"]["samples"].pop("wta-5")
write_case("missing_sample", missing_sample)

new_sample = copy.deepcopy(approved)
new_sample["doublet_policy"]["samples"]["wta-extra"] = copy.deepcopy(
    new_sample["doublet_policy"]["samples"]["wta-1"]
)
write_case("new_sample", new_sample)

changed_counts = copy.deepcopy(approved)
changed_counts["doublet_policy"]["samples"]["wta-1"][
    "proposed_doublet_count"
] = 16
write_case("changed_counts", changed_counts)

invalid_action = copy.deepcopy(approved)
invalid_action["doublet_policy"]["samples"]["wta-1"][
    "proposed_action"
] = "delete_everything"
write_case("invalid_action", invalid_action)

wrong_project = copy.deepcopy(approved)
wrong_project["doublet_policy"]["identity"]["project_id"] = "WRONG_PROJECT"
write_case("wrong_project", wrong_project)

wrong_run = copy.deepcopy(approved)
wrong_run["doublet_policy"]["identity"][
    "biological_run_id"
] = "WRONG_RUN"
write_case("wrong_run", wrong_run)
PY

run_validator() {
    local case_name="$1"
    local out_dir="${TMP_DIR}/out_${case_name}"

    python3 "${VALIDATOR}" \
        --plan "${PLAN}" \
        --policy "${TMP_DIR}/${case_name}.yaml" \
        --validated-policy "${out_dir}/approved_doublet_policy.yaml" \
        --summary "${out_dir}/doublet_policy_approval.summary.tsv" \
        --session-info "${out_dir}/doublet_policy_approval.session.txt"
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

for case_name in \
    proposed \
    missing_reviewer \
    missing_approved_at \
    missing_sample \
    new_sample \
    changed_counts \
    invalid_action \
    wrong_project \
    wrong_run
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
validator_path = root / "scripts/control/validate_doublet_approval.py"

spec = importlib.util.spec_from_file_location(
    "validate_doublet_approval",
    validator_path,
)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

document = module.read_yaml(
    tmp / "out_approved/approved_doublet_policy.yaml",
    "Validated approved policy",
)
policy = document["doublet_policy"]

assert policy["approval"]["status"] == "approved"
assert policy["validation"]["status"] == "valid"
assert policy["validation"]["validator"] == "validate_doublet_approval.py"
assert policy["validation"]["validated_for_doublet_filter"] is True
assert policy["application"]["status"] == "not_applied"
assert policy["application"]["filtering_applied"] is False
assert int(policy["application"]["cells_removed"]) == 0

sample = policy["samples"]["wta-1"]
assert int(sample["cells_scored"]) == 199
assert int(sample["proposed_doublet_count"]) == 15
assert int(sample["proposed_singlet_count"]) == 184
assert float(sample["pK"]) == 0.3
assert int(sample["nExp_raw"]) == 15
assert int(sample["nExp_adjusted"]) == 15

print("Doublet approval validation tests: VALID")
PY

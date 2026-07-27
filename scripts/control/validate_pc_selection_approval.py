#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import platform
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

try:
    import yaml
except ModuleNotFoundError:  # pragma: no cover - exercised on minimal hosts.
    yaml = None


APPROVED_STATUSES = {
    "approved",
    "modified_and_approved",
}

VALIDATOR_NAME = "validate_pc_selection_approval.py"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate a human-reviewed PC-selection policy before "
            "downstream neighborhood, UMAP, clustering or integration phases "
            "can consume selected PCs."
        )
    )

    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--policy", required=True, type=Path)
    parser.add_argument("--validated-policy", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--session-info", required=True, type=Path)

    return parser.parse_args()


def parse_scalar(value: str) -> Any:
    value = value.strip()

    if value in {"", "null", "Null", "NULL", "~"}:
        return None

    if value in {"true", "True", "TRUE", "yes", "Yes", "YES"}:
        return True

    if value in {"false", "False", "FALSE", "no", "No", "NO"}:
        return False

    if (
        len(value) >= 2
        and value[0] == value[-1]
        and value[0] in {"'", '"'}
    ):
        return value[1:-1]

    try:
        parsed_int = int(value)
    except ValueError:
        parsed_int = None

    if parsed_int is not None:
        return parsed_int

    try:
        parsed_float = float(value)
    except ValueError:
        return value

    if not math.isfinite(parsed_float):
        return value

    return parsed_float


def read_simple_yaml(path: Path) -> dict[str, Any]:
    root: dict[str, Any] = {}
    stack: list[tuple[int, dict[str, Any]]] = [(-1, root)]

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip()

            if not line.strip() or line.lstrip().startswith("#"):
                continue

            stripped = line.lstrip()

            if stripped.startswith("- "):
                continue

            indent = len(line) - len(stripped)

            if ":" not in stripped:
                raise ValueError(
                    f"Unsupported YAML syntax at {path}:{line_number}"
                )

            key, value = stripped.split(":", 1)
            key = key.strip()

            if not key:
                raise ValueError(f"Empty YAML key at {path}:{line_number}")

            while stack and indent <= stack[-1][0]:
                stack.pop()

            if not stack:
                raise ValueError(
                    f"Invalid YAML indentation at {path}:{line_number}"
                )

            parent = stack[-1][1]

            if value.strip() == "":
                child: dict[str, Any] = {}
                parent[key] = child
                stack.append((indent, child))
            else:
                parent[key] = parse_scalar(value)

    return root


def read_yaml(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")

    if yaml is not None:
        with path.open("r", encoding="utf-8") as handle:
            document = yaml.safe_load(handle)
    else:
        document = read_simple_yaml(path)

    if not isinstance(document, dict):
        raise ValueError(f"{label} root must be a mapping: {path}")

    return document


def yaml_scalar(value: Any) -> str:
    if value is None:
        return "null"

    if isinstance(value, bool):
        return "true" if value else "false"

    if isinstance(value, (int, float)):
        return str(value)

    text = str(value)

    if text == "":
        return '""'

    escaped = text.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def write_yaml_value(handle, key: str, value: Any, indent: int = 0) -> None:
    prefix = " " * indent

    if isinstance(value, dict):
        handle.write(f"{prefix}{key}:\n")
        for child_key, child_value in value.items():
            write_yaml_value(handle, str(child_key), child_value, indent + 2)
        return

    handle.write(f"{prefix}{key}: {yaml_scalar(value)}\n")


def write_yaml(path: Path, document: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for key, value in document.items():
            write_yaml_value(handle, key, value)


def get_mapping(document: dict[str, Any], key: str, label: str) -> dict[str, Any]:
    value = document.get(key)

    if not isinstance(value, dict):
        raise ValueError(f"{label}.{key} must be a mapping.")

    return value


def require_non_empty_string(value: Any, label: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{label} must be a non-empty string.")

    value = value.strip()

    if not value:
        raise ValueError(f"{label} must be a non-empty string.")

    return value


def validate_timestamp(value: Any, label: str) -> str:
    timestamp = require_non_empty_string(value, label)

    try:
        datetime.fromisoformat(timestamp.replace("Z", "+00:00"))
    except ValueError as error:
        raise ValueError(f"{label} must use ISO-8601 format.") from error

    return timestamp


def parse_positive_int(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be a positive integer.")

    try:
        parsed_float = float(str(value))
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be a positive integer.") from error

    if not parsed_float.is_integer():
        raise ValueError(f"{label} must be an integer.")

    parsed = int(parsed_float)

    if parsed <= 0:
        raise ValueError(f"{label} must be positive.")

    return parsed


def validate_identity(
    plan: dict[str, Any],
    policy: dict[str, Any],
) -> tuple[str, str, str]:
    plan_section = get_mapping(plan, "plan", "analysis_plan")
    experiment = get_mapping(plan, "experiment", "analysis_plan")
    identity = get_mapping(policy, "identity", "pc_selection_policy")

    plan_id = require_non_empty_string(
        plan_section.get("plan_id"),
        "analysis_plan.plan.plan_id",
    )
    expected_project_id = require_non_empty_string(
        experiment.get("project_id"),
        "analysis_plan.experiment.project_id",
    )
    expected_run_id = require_non_empty_string(
        experiment.get("biological_run_id"),
        "analysis_plan.experiment.biological_run_id",
    )

    observed_project_id = str(identity.get("project_id", "")).strip()
    observed_run_id = str(identity.get("biological_run_id", "")).strip()

    if observed_project_id != expected_project_id:
        raise ValueError("Policy project_id does not match analysis plan.")

    if observed_run_id != expected_run_id:
        raise ValueError(
            "Policy biological_run_id does not match analysis plan."
        )

    return plan_id, expected_project_id, expected_run_id


def validate_approval(policy: dict[str, Any]) -> tuple[str, str, str]:
    approval = get_mapping(policy, "approval", "pc_selection_policy")

    if approval.get("required", True) is not True:
        raise ValueError("approval.required must be true for PC selection.")

    status = str(approval.get("status", "")).strip().lower()

    if status not in APPROVED_STATUSES:
        raise ValueError(
            "PC-selection policy is not approved. Accepted statuses: "
            + ", ".join(sorted(APPROVED_STATUSES))
            + f". Observed status: {status!r}"
        )

    reviewed_by = require_non_empty_string(
        approval.get("reviewed_by"),
        "approval.reviewed_by",
    )
    approved_at = validate_timestamp(
        approval.get("approved_at"),
        "approval.approved_at",
    )

    return status, reviewed_by, approved_at


def validate_input(policy: dict[str, Any]) -> dict[str, Any]:
    input_section = get_mapping(policy, "input", "pc_selection_policy")
    canonical = dict(input_section)

    canonical["cells"] = parse_positive_int(
        input_section.get("cells"),
        "input.cells",
    )
    canonical["features"] = parse_positive_int(
        input_section.get("features"),
        "input.features",
    )
    canonical["variable_features"] = parse_positive_int(
        input_section.get("variable_features"),
        "input.variable_features",
    )
    canonical["max_computed_pcs"] = parse_positive_int(
        input_section.get("max_computed_pcs"),
        "input.max_computed_pcs",
    )

    return canonical


def validate_recommendation(policy: dict[str, Any]) -> dict[str, Any]:
    recommendation = get_mapping(
        policy,
        "recommendation",
        "pc_selection_policy",
    )
    canonical = dict(recommendation)

    canonical["recommended_pcs"] = parse_positive_int(
        recommendation.get("recommended_pcs"),
        "recommendation.recommended_pcs",
    )
    canonical["max_computed_pcs"] = parse_positive_int(
        recommendation.get("max_computed_pcs"),
        "recommendation.max_computed_pcs",
    )

    if canonical["recommended_pcs"] > canonical["max_computed_pcs"]:
        raise ValueError(
            "recommendation.recommended_pcs cannot exceed max_computed_pcs."
        )

    return canonical


def validate_selected_pcs(
    policy: dict[str, Any],
    recommendation: dict[str, Any],
    approval_status: str,
) -> int:
    proposed_pcs = int(recommendation["recommended_pcs"])
    max_computed_pcs = int(recommendation["max_computed_pcs"])

    selected_pcs = policy.get("selected_pcs", proposed_pcs)
    selected_pcs = parse_positive_int(selected_pcs, "selected_pcs")

    if selected_pcs > max_computed_pcs:
        raise ValueError("selected_pcs cannot exceed max_computed_pcs.")

    if selected_pcs != proposed_pcs and approval_status != "modified_and_approved":
        raise ValueError(
            "selected_pcs differs from the proposed value and requires "
            "approval.status modified_and_approved."
        )

    return selected_pcs


def validate_application(policy: dict[str, Any]) -> dict[str, Any]:
    application = get_mapping(policy, "application", "pc_selection_policy")
    status = str(application.get("status", "")).strip().lower()

    if status != "not_applied":
        raise ValueError("application.status must be not_applied.")

    if application.get("pc_selection_applied", False) is not False:
        raise ValueError("application.pc_selection_applied must be false.")

    canonical = dict(application)
    canonical["status"] = status
    canonical["pc_selection_applied"] = False
    return canonical


def main() -> int:
    args = parse_args()

    try:
        plan_document = read_yaml(args.plan, "Validated analysis plan")
        policy_document = read_yaml(args.policy, "PC-selection policy")
        policy = get_mapping(
            policy_document,
            "pc_selection_policy",
            "policy_document",
        )

        plan_id, project_id, biological_run_id = validate_identity(
            plan_document,
            policy,
        )
        approval_status, reviewed_by, approved_at = validate_approval(policy)
        canonical_input = validate_input(policy)
        canonical_recommendation = validate_recommendation(policy)

        if (
            canonical_input["max_computed_pcs"]
            != canonical_recommendation["max_computed_pcs"]
        ):
            raise ValueError(
                "input.max_computed_pcs does not match "
                "recommendation.max_computed_pcs."
            )

        selected_pcs = validate_selected_pcs(
            policy,
            canonical_recommendation,
            approval_status,
        )
        canonical_application = validate_application(policy)

    except (FileNotFoundError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    validated_at = (
        datetime.now()
        .astimezone()
        .replace(microsecond=0)
        .isoformat()
    )

    canonical_policy = dict(policy)
    canonical_policy["input"] = canonical_input
    canonical_policy["recommendation"] = canonical_recommendation
    canonical_policy["proposal"] = {
        "recommended_pcs": canonical_recommendation["recommended_pcs"],
        "method": canonical_recommendation.get("method", ""),
        "max_computed_pcs": canonical_recommendation["max_computed_pcs"],
    }
    canonical_policy["selected_pcs"] = selected_pcs
    canonical_policy["approval"] = {
        **policy["approval"],
        "required": True,
        "status": approval_status,
        "reviewed_by": reviewed_by,
        "approved_at": approved_at,
    }
    canonical_policy["validation"] = {
        "status": "valid",
        "validated_at": validated_at,
        "validator": VALIDATOR_NAME,
        "validated_for_pc_selection": True,
        "proposed_pcs": canonical_recommendation["recommended_pcs"],
        "selected_pcs": selected_pcs,
        "max_computed_pcs": canonical_recommendation["max_computed_pcs"],
    }
    canonical_policy["application"] = canonical_application

    output_document = {
        "pc_selection_policy": canonical_policy
    }

    args.validated_policy.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.session_info.parent.mkdir(parents=True, exist_ok=True)

    write_yaml(args.validated_policy, output_document)

    with args.summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "plan_id",
                "project_id",
                "biological_run_id",
                "validation_status",
                "approval_status",
                "reviewed_by",
                "approved_at",
                "validated_at",
                "proposed_pcs",
                "selected_pcs",
                "max_computed_pcs",
            ]
        )
        writer.writerow(
            [
                plan_id,
                project_id,
                biological_run_id,
                "valid",
                approval_status,
                reviewed_by,
                approved_at,
                validated_at,
                canonical_recommendation["recommended_pcs"],
                selected_pcs,
                canonical_recommendation["max_computed_pcs"],
            ]
        )

    with args.session_info.open("w", encoding="utf-8") as handle:
        handle.write(f"python_version={platform.python_version()}\n")
        handle.write(f"platform={platform.platform()}\n")
        handle.write("validation_status=valid\n")
        handle.write(f"approval_status={approval_status}\n")
        handle.write(f"reviewed_by={reviewed_by}\n")
        handle.write(f"validated_at={validated_at}\n")
        handle.write(f"validator={VALIDATOR_NAME}\n")
        handle.write(
            f"proposed_pcs={canonical_recommendation['recommended_pcs']}\n"
        )
        handle.write(f"selected_pcs={selected_pcs}\n")

    print("PC-selection approval validated successfully.")
    print(f"Approval status: {approval_status}")
    print(f"Proposed PCs: {canonical_recommendation['recommended_pcs']}")
    print(f"Selected PCs: {selected_pcs}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

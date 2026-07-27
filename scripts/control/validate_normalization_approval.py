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

ALLOWED_METHODS = {
    "LogNormalize",
    "SCTransform",
}

ALLOWED_VARIABLE_FEATURE_METHODS = {
    "vst",
}

ALLOWED_REGRESSORS = {
    "percent_mito",
    "nCount_RNA",
    "cell_cycle",
}

VALIDATOR_NAME = "validate_normalization_approval.py"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate a human-reviewed normalization policy before "
            "allowing NORMALIZATION to transform the Seurat object."
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
                raise ValueError(
                    f"Empty YAML key at {path}:{line_number}"
                )

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


def parse_positive_number(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be positive.")

    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be positive.") from error

    if parsed <= 0 or not math.isfinite(parsed):
        raise ValueError(f"{label} must be a positive finite number.")

    return parsed


def normalise_regressors(value: Any) -> str:
    if value is None:
        return "none"

    if isinstance(value, list):
        values = [str(item).strip() for item in value if str(item).strip()]
    else:
        text = str(value).strip()
        if text.lower() in {"", "none", "null"}:
            return "none"
        values = [part.strip() for part in text.split(",") if part.strip()]

    unknown = sorted(set(values) - ALLOWED_REGRESSORS)
    if unknown:
        raise ValueError(
            "strategy.regressors contains unsupported values: "
            + ", ".join(unknown)
        )

    return ",".join(values) if values else "none"


def validate_identity(
    plan: dict[str, Any],
    policy: dict[str, Any],
) -> tuple[str, str, str]:
    plan_section = get_mapping(plan, "plan", "analysis_plan")
    experiment = get_mapping(plan, "experiment", "analysis_plan")
    identity = get_mapping(policy, "identity", "normalization_policy")

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
    approval = get_mapping(policy, "approval", "normalization_policy")

    if approval.get("required", True) is not True:
        raise ValueError(
            "approval.required must be true for normalization policy."
        )

    status = str(approval.get("status", "")).strip().lower()

    if status not in APPROVED_STATUSES:
        raise ValueError(
            "Normalization policy is not approved. Accepted statuses: "
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
    input_section = get_mapping(policy, "input", "normalization_policy")
    readiness = str(input_section.get("readiness", "")).strip().lower()

    if readiness != "ready":
        raise ValueError("input.readiness must be ready.")

    canonical_input = dict(input_section)
    canonical_input["cells"] = parse_positive_int(
        input_section.get("cells"),
        "input.cells",
    )
    canonical_input["features"] = parse_positive_int(
        input_section.get("features"),
        "input.features",
    )
    canonical_input["samples"] = parse_positive_int(
        input_section.get("samples"),
        "input.samples",
    )
    canonical_input["readiness"] = readiness

    return canonical_input


def validate_strategy(policy: dict[str, Any]) -> dict[str, Any]:
    strategy = get_mapping(policy, "strategy", "normalization_policy")
    method = require_non_empty_string(
        strategy.get("method"),
        "strategy.method",
    )

    if method not in ALLOWED_METHODS:
        raise ValueError(
            "strategy.method must be one of: "
            + ", ".join(sorted(ALLOWED_METHODS))
        )

    scope = require_non_empty_string(strategy.get("scope"), "strategy.scope")
    if scope != "per_sample":
        raise ValueError("strategy.scope must be per_sample.")

    assay = require_non_empty_string(strategy.get("assay"), "strategy.assay")
    if assay != "RNA":
        raise ValueError("strategy.assay must be RNA for this MVP.")

    variable_method = require_non_empty_string(
        strategy.get("variable_features_method"),
        "strategy.variable_features_method",
    )
    if variable_method not in ALLOWED_VARIABLE_FEATURE_METHODS:
        raise ValueError("strategy.variable_features_method must be vst.")

    nfeatures = parse_positive_int(
        strategy.get("nfeatures"),
        "strategy.nfeatures",
    )
    regressors = normalise_regressors(strategy.get("regressors"))

    canonical = dict(strategy)
    canonical["method"] = method
    canonical["variable_features_method"] = variable_method
    canonical["nfeatures"] = nfeatures
    canonical["scope"] = scope
    canonical["assay"] = assay
    canonical["regressors"] = regressors

    if method == "LogNormalize":
        canonical["scale_factor"] = parse_positive_number(
            strategy.get("scale_factor"),
            "strategy.scale_factor",
        )

        if regressors != "none":
            raise ValueError(
                "LogNormalize policy must keep strategy.regressors as none."
            )

    if method == "SCTransform":
        if strategy.get("scale_factor") is not None:
            canonical["scale_factor"] = parse_positive_number(
                strategy.get("scale_factor"),
                "strategy.scale_factor",
            )

    return canonical


def validate_integration(policy: dict[str, Any], sample_count: int) -> dict[str, Any]:
    integration = get_mapping(policy, "integration", "normalization_policy")
    required = integration.get("required")

    if not isinstance(required, bool):
        raise ValueError("integration.required must be boolean.")

    derived_required = sample_count > 1
    canonical = dict(integration)
    canonical["required"] = required

    if required != derived_required:
        justification = str(
            integration.get("override_justification", "")
        ).strip()

        if not justification:
            raise ValueError(
                "integration.required differs from the derived value "
                "and requires integration.override_justification."
            )

        canonical["override_justification"] = justification

    return canonical


def validate_application(policy: dict[str, Any]) -> dict[str, Any]:
    application = get_mapping(policy, "application", "normalization_policy")
    status = str(application.get("status", "")).strip().lower()

    if status != "not_applied":
        raise ValueError("application.status must be not_applied.")

    if application.get("normalization_applied", False) is not False:
        raise ValueError("application.normalization_applied must be false.")

    canonical = dict(application)
    canonical["status"] = status
    canonical["normalization_applied"] = False
    return canonical


def main() -> int:
    args = parse_args()

    try:
        plan_document = read_yaml(args.plan, "Validated analysis plan")
        policy_document = read_yaml(args.policy, "Normalization policy")
        policy = get_mapping(
            policy_document,
            "normalization_policy",
            "policy_document",
        )

        plan_id, project_id, biological_run_id = validate_identity(
            plan_document,
            policy,
        )
        approval_status, reviewed_by, approved_at = validate_approval(policy)
        canonical_input = validate_input(policy)
        canonical_strategy = validate_strategy(policy)
        canonical_integration = validate_integration(
            policy,
            int(canonical_input["samples"]),
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
    canonical_policy["strategy"] = canonical_strategy
    canonical_policy["integration"] = canonical_integration
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
        "validated_for_normalization": True,
        "method": canonical_strategy["method"],
        "scope": canonical_strategy["scope"],
    }
    canonical_policy["application"] = canonical_application

    output_document = {
        "normalization_policy": canonical_policy
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
                "method",
                "scope",
                "integration_required",
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
                canonical_strategy["method"],
                canonical_strategy["scope"],
                canonical_integration["required"],
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

    print("Normalization approval validated successfully.")
    print(f"Approval status: {approval_status}")
    print(f"Method: {canonical_strategy['method']}")
    print(f"Scope: {canonical_strategy['scope']}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

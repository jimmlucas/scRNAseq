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

APPROVED_DECISION_STATUSES = {
    "approved",
    "modified_and_approved",
}

ALLOWED_ACTIONS = {
    "exclude_algorithmic_doublets",
    "keep_all_cells",
    "no_filter_not_run",
    "no_filter_no_doublets_detected",
}

ALLOWED_APPLICATION_ACTIONS = {
    "remove_proposed_doublets",
    "keep_all",
    "exclude_sample_from_filtering",
}

PROPOSED_TO_APPLICATION_ACTION = {
    "exclude_algorithmic_doublets": "remove_proposed_doublets",
    "keep_all_cells": "keep_all",
    "no_filter_no_doublets_detected": "keep_all",
    "no_filter_not_run": "exclude_sample_from_filtering",
}

VALIDATOR_NAME = "validate_doublet_approval.py"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate a human-reviewed doublet policy before "
            "allowing DOUBLET_FILTER to remove cells."
        )
    )

    parser.add_argument(
        "--plan",
        required=True,
        type=Path,
        help="Validated BBUDY analysis plan YAML.",
    )
    parser.add_argument(
        "--policy",
        required=True,
        type=Path,
        help="Human-reviewed doublet policy YAML.",
    )
    parser.add_argument(
        "--validated-policy",
        required=True,
        type=Path,
        help="Canonical validated doublet policy YAML output.",
    )
    parser.add_argument(
        "--summary",
        required=True,
        type=Path,
        help="Validation summary TSV output.",
    )
    parser.add_argument(
        "--session-info",
        required=True,
        type=Path,
        help="Runtime session information output.",
    )

    return parser.parse_args()


def parse_scalar(value: str) -> Any:
    value = value.strip()

    if value in {"", "null", "Null", "NULL", "~"}:
        return None

    if value in {"true", "True", "TRUE"}:
        return True

    if value in {"false", "False", "FALSE"}:
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


def read_yaml(
    path: Path,
    label: str,
) -> dict[str, Any]:
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


def write_yaml_value(
    handle,
    key: str,
    value: Any,
    indent: int = 0,
) -> None:
    prefix = " " * indent

    if isinstance(value, dict):
        handle.write(f"{prefix}{key}:\n")

        for child_key, child_value in value.items():
            write_yaml_value(
                handle,
                str(child_key),
                child_value,
                indent + 2,
            )

        return

    handle.write(f"{prefix}{key}: {yaml_scalar(value)}\n")


def write_yaml(
    path: Path,
    document: dict[str, Any],
) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for key, value in document.items():
            write_yaml_value(handle, key, value)


def get_mapping(
    document: dict[str, Any],
    key: str,
    label: str,
) -> dict[str, Any]:
    value = document.get(key)

    if not isinstance(value, dict):
        raise ValueError(f"{label}.{key} must be a mapping.")

    return value


def require_non_empty_string(
    value: Any,
    label: str,
) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{label} must be a non-empty string.")

    value = value.strip()

    if not value:
        raise ValueError(f"{label} must be a non-empty string.")

    return value


def validate_timestamp(
    value: Any,
    label: str,
) -> str:
    timestamp = require_non_empty_string(value, label)
    normalised = timestamp.replace("Z", "+00:00")

    try:
        datetime.fromisoformat(normalised)
    except ValueError as error:
        raise ValueError(f"{label} must use ISO-8601 format.") from error

    return timestamp


def parse_int(
    value: Any,
    label: str,
) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be an integer.")

    try:
        parsed = int(float(str(value)))
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be an integer.") from error

    if parsed < 0:
        raise ValueError(f"{label} cannot be negative.")

    return parsed


def parse_optional_float(
    value: Any,
    label: str,
) -> float | None:
    if value is None:
        return None

    if isinstance(value, bool):
        raise ValueError(f"{label} must be numeric or null.")

    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be numeric or null.") from error

    if not math.isfinite(parsed):
        raise ValueError(f"{label} must be finite.")

    return parsed


def validate_identity(
    plan: dict[str, Any],
    policy: dict[str, Any],
) -> tuple[str, str, str]:
    plan_section = get_mapping(plan, "plan", "analysis_plan")
    experiment = get_mapping(plan, "experiment", "analysis_plan")
    identity = get_mapping(policy, "identity", "doublet_policy")

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

    comparisons = (
        ("project_id", expected_project_id, observed_project_id),
        ("biological_run_id", expected_run_id, observed_run_id),
    )

    for label, expected, observed in comparisons:
        if expected != observed:
            raise ValueError(
                f"Doublet policy {label} does not match the analysis "
                f"plan: expected={expected!r}, observed={observed!r}"
            )

    return plan_id, expected_project_id, expected_run_id


def validate_approval(
    policy: dict[str, Any],
) -> tuple[str, str, str]:
    approval = get_mapping(policy, "approval", "doublet_policy")

    if approval.get("required", True) is not True:
        raise ValueError(
            "approval.required must be true for scientific runs."
        )

    status = str(approval.get("status", "")).strip().lower()

    if status not in APPROVED_STATUSES:
        raise ValueError(
            "Doublet policy is not approved. Accepted statuses: "
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


def validate_cell_class_counts(
    sample: str,
    sample_policy: dict[str, Any],
) -> dict[str, int]:
    counts = get_mapping(
        sample_policy,
        "cell_class_counts",
        f"samples.{sample}",
    )

    canonical_counts: dict[str, int] = {}

    for class_name in ("Doublet", "Singlet", "NotRun", "Other"):
        canonical_counts[class_name] = parse_int(
            counts.get(class_name, 0),
            f"samples.{sample}.cell_class_counts.{class_name}",
        )

    return canonical_counts


def validate_samples(
    policy: dict[str, Any],
) -> tuple[dict[str, dict[str, Any]], int, int]:
    samples = get_mapping(policy, "samples", "doublet_policy")

    if not samples:
        raise ValueError("doublet_policy.samples cannot be empty.")

    canonical_samples: dict[str, dict[str, Any]] = {}
    total_doublets = 0
    total_scored = 0

    for sample, sample_policy in samples.items():
        if not isinstance(sample, str) or not sample:
            raise ValueError("Sample names must be non-empty strings.")

        if not isinstance(sample_policy, dict):
            raise ValueError(f"samples.{sample} must be a mapping.")

        detection_status = require_non_empty_string(
            sample_policy.get("detection_status"),
            f"samples.{sample}.detection_status",
        )
        decision_status = str(
            sample_policy.get("decision_status", "")
        ).strip().lower()

        if decision_status not in APPROVED_DECISION_STATUSES:
            raise ValueError(
                f"samples.{sample}.decision_status must be reviewed. "
                "Accepted statuses: "
                + ", ".join(sorted(APPROVED_DECISION_STATUSES))
                + f". Observed status: {decision_status!r}"
            )

        proposed_action = require_non_empty_string(
            sample_policy.get("proposed_action"),
            f"samples.{sample}.proposed_action",
        )

        if proposed_action not in ALLOWED_ACTIONS:
            raise ValueError(
                f"samples.{sample}.proposed_action is not supported: "
                f"{proposed_action!r}"
            )

        application_action = str(
            sample_policy.get(
                "application_action",
                PROPOSED_TO_APPLICATION_ACTION[proposed_action],
            )
        ).strip()

        if application_action not in ALLOWED_APPLICATION_ACTIONS:
            raise ValueError(
                f"samples.{sample}.application_action is not supported: "
                f"{application_action!r}"
            )

        cells_scored = parse_int(
            sample_policy.get("cells_scored"),
            f"samples.{sample}.cells_scored",
        )
        proposed_doublets = parse_int(
            sample_policy.get("proposed_doublet_count"),
            f"samples.{sample}.proposed_doublet_count",
        )
        proposed_singlets = parse_int(
            sample_policy.get("proposed_singlet_count"),
            f"samples.{sample}.proposed_singlet_count",
        )
        n_exp_adjusted = parse_int(
            sample_policy.get("nExp_adjusted"),
            f"samples.{sample}.nExp_adjusted",
        )

        if proposed_doublets + proposed_singlets > cells_scored:
            raise ValueError(
                f"samples.{sample} proposes more classified cells than "
                "were scored."
            )

        if (
            application_action == "remove_proposed_doublets"
            and detection_status != "completed"
        ):
            raise ValueError(
                f"samples.{sample} cannot exclude doublets because "
                f"detection_status={detection_status!r}."
            )

        if (
            application_action == "remove_proposed_doublets"
            and proposed_doublets <= 0
        ):
            raise ValueError(
                f"samples.{sample} cannot exclude zero proposed doublets."
            )

        cell_class_counts = validate_cell_class_counts(
            sample,
            sample_policy,
        )

        if cell_class_counts["Doublet"] != proposed_doublets:
            raise ValueError(
                f"samples.{sample}.cell_class_counts.Doublet does not "
                "match proposed_doublet_count."
            )

        canonical_sample = dict(sample_policy)
        canonical_sample["detection_status"] = detection_status
        canonical_sample["decision_status"] = decision_status
        canonical_sample["proposed_action"] = proposed_action
        canonical_sample["application_action"] = application_action
        canonical_sample["cells_scored"] = cells_scored
        canonical_sample["proposed_doublet_count"] = proposed_doublets
        canonical_sample["proposed_singlet_count"] = proposed_singlets
        canonical_sample["expected_doublet_rate"] = parse_optional_float(
            sample_policy.get("expected_doublet_rate"),
            f"samples.{sample}.expected_doublet_rate",
        )
        canonical_sample["pN"] = parse_optional_float(
            sample_policy.get("pN"),
            f"samples.{sample}.pN",
        )
        canonical_sample["pK"] = parse_optional_float(
            sample_policy.get("pK"),
            f"samples.{sample}.pK",
        )
        canonical_sample["nExp_raw"] = parse_int(
            sample_policy.get("nExp_raw"),
            f"samples.{sample}.nExp_raw",
        )
        canonical_sample["homotypic_prop"] = parse_optional_float(
            sample_policy.get("homotypic_prop"),
            f"samples.{sample}.homotypic_prop",
        )
        canonical_sample["nExp_adjusted"] = n_exp_adjusted
        canonical_sample["cell_class_counts"] = cell_class_counts

        canonical_samples[sample] = canonical_sample
        total_doublets += proposed_doublets
        total_scored += cells_scored

    return canonical_samples, total_scored, total_doublets


def validate_summary(
    policy: dict[str, Any],
    sample_count: int,
    total_scored: int,
    total_doublets: int,
) -> dict[str, Any]:
    summary = get_mapping(policy, "summary", "doublet_policy")
    cells_removed = parse_int(
        summary.get("cells_removed_by_scoring", 0),
        "summary.cells_removed_by_scoring",
    )

    if cells_removed != 0:
        raise ValueError(
            "Doublet scoring must not remove cells before approval."
        )

    observed_samples = parse_int(
        summary.get("samples"),
        "summary.samples",
    )
    observed_cells = parse_int(
        summary.get("cells_scored"),
        "summary.cells_scored",
    )
    observed_doublets = parse_int(
        summary.get("proposed_doublets"),
        "summary.proposed_doublets",
    )

    if observed_samples != sample_count:
        raise ValueError(
            "summary.samples does not match doublet_policy.samples."
        )

    if observed_cells != total_scored:
        raise ValueError(
            "summary.cells_scored does not match sample totals."
        )

    if observed_doublets != total_doublets:
        raise ValueError(
            "summary.proposed_doublets does not match sample totals."
        )

    canonical_summary = dict(summary)
    canonical_summary["cells_scored"] = total_scored
    canonical_summary["samples"] = sample_count
    canonical_summary["proposed_doublets"] = total_doublets
    canonical_summary["cells_removed_by_scoring"] = 0

    return canonical_summary


def main() -> int:
    args = parse_args()

    try:
        plan_document = read_yaml(args.plan, "Validated analysis plan")
        policy_document = read_yaml(args.policy, "Doublet policy")
        policy = get_mapping(
            policy_document,
            "doublet_policy",
            "policy_document",
        )

        plan_id, project_id, biological_run_id = validate_identity(
            plan_document,
            policy,
        )
        approval_status, reviewed_by, approved_at = validate_approval(
            policy
        )
        canonical_samples, total_scored, total_doublets = validate_samples(
            policy
        )
        canonical_summary = validate_summary(
            policy,
            len(canonical_samples),
            total_scored,
            total_doublets,
        )

    except (FileNotFoundError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    except Exception as error:
        if yaml is not None and error.__class__.__name__.endswith("YAMLError"):
            print(f"ERROR: {error}", file=sys.stderr)
            return 1
        raise

    validated_at = (
        datetime.now()
        .astimezone()
        .replace(microsecond=0)
        .isoformat()
    )

    canonical_policy = dict(policy)
    canonical_policy["summary"] = canonical_summary
    canonical_policy["samples"] = canonical_samples
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
        "validated_for_doublet_filter": True,
        "sample_count": len(canonical_samples),
        "cells_scored": total_scored,
        "approved_doublet_count": total_doublets,
    }
    canonical_policy["application"] = {
        "status": "not_applied",
        "filtering_applied": False,
        "cells_removed": 0,
    }

    output_document = {
        "doublet_policy": canonical_policy
    }

    args.validated_policy.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.session_info.parent.mkdir(parents=True, exist_ok=True)

    write_yaml(args.validated_policy, output_document)

    with args.summary.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
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
                "sample_count",
                "cells_scored",
                "approved_doublet_count",
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
                len(canonical_samples),
                total_scored,
                total_doublets,
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

    print("Doublet policy approval validated successfully.")
    print(f"Approval status: {approval_status}")
    print(f"Reviewed by: {reviewed_by}")
    print(f"Approved doublets: {total_doublets}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

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

import yaml


APPROVED_STATUSES = {
    "approved",
    "modified_and_approved",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate a human-reviewed QC threshold policy before "
            "allowing QC_FILTER to remove cells."
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
        help="Human-reviewed QC threshold policy YAML.",
    )
    parser.add_argument(
        "--validated-policy",
        required=True,
        type=Path,
        help="Canonical validated policy YAML output.",
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


def read_yaml(
    path: Path,
    label: str,
) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(
            f"{label} not found: {path}"
        )

    with path.open(
        "r",
        encoding="utf-8",
    ) as handle:
        document = yaml.safe_load(handle)

    if not isinstance(document, dict):
        raise ValueError(
            f"{label} root must be a mapping: {path}"
        )

    return document


def get_mapping(
    document: dict[str, Any],
    key: str,
    label: str,
) -> dict[str, Any]:
    value = document.get(key)

    if not isinstance(value, dict):
        raise ValueError(
            f"{label}.{key} must be a mapping."
        )

    return value


def parse_optional_float(
    value: Any,
    label: str,
) -> float | None:
    if value is None:
        return None

    if isinstance(value, bool):
        raise ValueError(
            f"{label} must be numeric or null."
        )

    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"{label} must be numeric or null."
        ) from error

    if not math.isfinite(parsed):
        raise ValueError(
            f"{label} must be finite."
        )

    return parsed


def require_non_empty_string(
    value: Any,
    label: str,
) -> str:
    if not isinstance(value, str):
        raise ValueError(
            f"{label} must be a non-empty string."
        )

    value = value.strip()

    if not value:
        raise ValueError(
            f"{label} must be a non-empty string."
        )

    return value


def validate_timestamp(
    value: Any,
    label: str,
) -> str:
    timestamp = require_non_empty_string(
        value,
        label,
    )

    normalised = timestamp.replace(
        "Z",
        "+00:00",
    )

    try:
        datetime.fromisoformat(normalised)
    except ValueError as error:
        raise ValueError(
            f"{label} must use ISO-8601 format."
        ) from error

    return timestamp


def validate_identity(
    plan: dict[str, Any],
    policy: dict[str, Any],
) -> tuple[str, str, str]:
    plan_section = get_mapping(
        plan,
        "plan",
        "analysis_plan",
    )
    experiment = get_mapping(
        plan,
        "experiment",
        "analysis_plan",
    )

    expected_plan_id = str(
        plan_section.get(
            "plan_id",
            "",
        )
    )
    expected_project_id = str(
        experiment.get(
            "project_id",
            "",
        )
    )
    expected_run_id = str(
        experiment.get(
            "biological_run_id",
            "",
        )
    )

    policy_plan_id = str(
        policy.get(
            "plan_id",
            "",
        )
    )
    policy_project_id = str(
        policy.get(
            "project_id",
            "",
        )
    )
    policy_run_id = str(
        policy.get(
            "biological_run_id",
            "",
        )
    )

    comparisons = (
        (
            "plan_id",
            expected_plan_id,
            policy_plan_id,
        ),
        (
            "project_id",
            expected_project_id,
            policy_project_id,
        ),
        (
            "biological_run_id",
            expected_run_id,
            policy_run_id,
        ),
    )

    for label, expected, observed in comparisons:
        if expected != observed:
            raise ValueError(
                f"Policy {label} does not match the analysis plan: "
                f"expected={expected!r}, observed={observed!r}"
            )

    return (
        expected_plan_id,
        expected_project_id,
        expected_run_id,
    )


def validate_approval(
    policy: dict[str, Any],
) -> tuple[str, str, str]:
    approval = get_mapping(
        policy,
        "approval",
        "qc_threshold_policy",
    )

    required = approval.get(
        "required",
        True,
    )

    if required is not True:
        raise ValueError(
            "approval.required must be true for scientific runs."
        )

    status = str(
        approval.get(
            "status",
            "",
        )
    ).strip().lower()

    if status not in APPROVED_STATUSES:
        raise ValueError(
            "QC threshold policy is not approved. "
            "Accepted statuses: "
            + ", ".join(
                sorted(APPROVED_STATUSES)
            )
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

    return (
        status,
        reviewed_by,
        approved_at,
    )


def validate_thresholds(
    policy: dict[str, Any],
) -> tuple[dict[str, Any], int]:
    thresholds = get_mapping(
        policy,
        "thresholds",
        "qc_threshold_policy",
    )

    if not thresholds:
        raise ValueError(
            "qc_threshold_policy.thresholds cannot be empty."
        )

    enabled_count = 0
    canonical_thresholds: dict[str, Any] = {}

    for metric_name, threshold in thresholds.items():
        if not isinstance(metric_name, str) or not metric_name:
            raise ValueError(
                "Threshold metric names must be non-empty strings."
            )

        if not isinstance(threshold, dict):
            raise ValueError(
                f"Threshold definition for {metric_name} "
                "must be a mapping."
            )

        enabled = threshold.get(
            "enabled",
            False,
        )

        if not isinstance(enabled, bool):
            raise ValueError(
                f"thresholds.{metric_name}.enabled "
                "must be boolean."
            )

        lower = parse_optional_float(
            threshold.get("min"),
            f"thresholds.{metric_name}.min",
        )

        upper = parse_optional_float(
            threshold.get("max"),
            f"thresholds.{metric_name}.max",
        )

        if (
            lower is not None
            and upper is not None
            and lower >= upper
        ):
            raise ValueError(
                f"Invalid bounds for {metric_name}: "
                f"min={lower}, max={upper}"
            )

        if enabled and lower is None and upper is None:
            raise ValueError(
                f"Enabled threshold {metric_name} "
                "must define min and/or max."
            )

        if enabled:
            enabled_count += 1

        canonical_threshold = dict(threshold)
        canonical_threshold["enabled"] = enabled
        canonical_threshold["min"] = lower
        canonical_threshold["max"] = upper

        canonical_thresholds[
            metric_name
        ] = canonical_threshold

    if enabled_count == 0:
        raise ValueError(
            "The approved policy contains no enabled thresholds."
        )

    return canonical_thresholds, enabled_count


def main() -> int:
    args = parse_args()

    try:
        plan_document = read_yaml(
            args.plan,
            "Validated analysis plan",
        )
        policy_document = read_yaml(
            args.policy,
            "QC threshold policy",
        )

        policy = get_mapping(
            policy_document,
            "qc_threshold_policy",
            "policy_document",
        )

        (
            plan_id,
            project_id,
            biological_run_id,
        ) = validate_identity(
            plan_document,
            policy,
        )

        (
            approval_status,
            reviewed_by,
            approved_at,
        ) = validate_approval(policy)

        (
            canonical_thresholds,
            enabled_threshold_count,
        ) = validate_thresholds(policy)

    except (
        FileNotFoundError,
        ValueError,
        yaml.YAMLError,
    ) as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )
        return 1

    canonical_policy = dict(policy)
    canonical_policy["thresholds"] = canonical_thresholds
    canonical_policy["approval"] = {
        **policy["approval"],
        "required": True,
        "status": approval_status,
        "reviewed_by": reviewed_by,
        "approved_at": approved_at,
    }
    canonical_policy["validation"] = {
        "status": "valid",
        "validated_for_qc_filter": True,
        "enabled_threshold_count": (
            enabled_threshold_count
        ),
    }

    output_document = {
        "qc_threshold_policy": canonical_policy
    }

    args.validated_policy.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    args.summary.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    args.session_info.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with args.validated_policy.open(
        "w",
        encoding="utf-8",
    ) as handle:
        yaml.safe_dump(
            output_document,
            handle,
            sort_keys=False,
            allow_unicode=True,
        )

    with args.summary.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
        )

        writer.writerow(
            [
                "plan_id",
                "project_id",
                "biological_run_id",
                "validation_status",
                "approval_status",
                "reviewed_by",
                "approved_at",
                "enabled_threshold_count",
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
                enabled_threshold_count,
            ]
        )

    with args.session_info.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            f"python_version={platform.python_version()}\n"
        )
        handle.write(
            f"platform={platform.platform()}\n"
        )
        handle.write(
            "validation_status=valid\n"
        )
        handle.write(
            f"approval_status={approval_status}\n"
        )
        handle.write(
            f"reviewed_by={reviewed_by}\n"
        )

    print(
        "QC threshold approval validated successfully."
    )
    print(
        f"Approval status: {approval_status}"
    )
    print(
        f"Reviewed by: {reviewed_by}"
    )
    print(
        "Enabled thresholds: "
        f"{enabled_threshold_count}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

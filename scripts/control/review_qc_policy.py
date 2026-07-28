#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import copy
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml


DECISION_TO_STATUS = {
    "approve": "approved",
    "modify-and-approve": "modified_and_approved",
    "reject": "rejected",
    "needs-revision": "needs_revision",
}

APPROVING_DECISIONS = {
    "approve",
    "modify-and-approve",
}

ALLOWED_THRESHOLD_FIELDS = {
    "enabled",
    "min",
    "max",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare a human-reviewed QC threshold policy from a "
            "proposed policy and its impact diagnostics. This command "
            "does not validate or apply the policy."
        )
    )

    parser.add_argument(
        "--results-dir",
        type=Path,
        help=(
            "Pipeline results directory. Used to locate "
            "control/qc_thresholds/qc_threshold_policy.yaml and "
            "control/qc_threshold_impact artifacts."
        ),
    )
    parser.add_argument(
        "--policy-yaml",
        type=Path,
        help="Explicit proposed QC threshold policy YAML.",
    )
    parser.add_argument(
        "--impact-dir",
        type=Path,
        help="Explicit QC threshold impact artifact directory.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help=(
            "Directory for reviewed policy outputs. Defaults to "
            "configs/real_validation/<project_id>/review when that "
            "project directory exists, otherwise "
            "<results-dir>/control/qc_threshold_review."
        ),
    )
    parser.add_argument(
        "--decision",
        required=True,
        choices=sorted(DECISION_TO_STATUS),
        help="Human review decision.",
    )
    parser.add_argument(
        "--reviewer",
        required=True,
        help="Reviewer name or identifier.",
    )
    parser.add_argument(
        "--approved-at",
        help=(
            "ISO-8601 timestamp. Defaults to current UTC time. Used "
            "for approving decisions; recorded as reviewed_at for all "
            "decisions."
        ),
    )
    parser.add_argument(
        "--reason",
        default="",
        help="Free-text reason for the review decision.",
    )
    parser.add_argument(
        "--set",
        action="append",
        default=[],
        metavar="METRIC.FIELD=VALUE",
        help=(
            "Modify an allowed threshold field. Accepted fields are "
            "enabled, min and max. Examples: percent.mito.max=15, "
            "thresholds.nFeature_RNA.min=200."
        ),
    )

    return parser.parse_args()


def read_yaml(path: Path, label: str) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")

    with path.open("r", encoding="utf-8") as handle:
        document = yaml.safe_load(handle)

    if not isinstance(document, dict):
        raise ValueError(f"{label} root must be a mapping: {path}")

    return document


def require_mapping(
    document: dict[str, Any],
    key: str,
    label: str,
) -> dict[str, Any]:
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


def iso_timestamp(value: str | None) -> str:
    if value:
        normalised = value.replace("Z", "+00:00")
        try:
            datetime.fromisoformat(normalised)
        except ValueError as error:
            raise ValueError(
                "--approved-at must use ISO-8601 format."
            ) from error
        return value

    return datetime.now(timezone.utc).isoformat()


def resolve_policy_path(args: argparse.Namespace) -> Path:
    if args.policy_yaml:
        return args.policy_yaml

    if not args.results_dir:
        raise ValueError(
            "Either --results-dir or --policy-yaml must be provided."
        )

    return (
        args.results_dir
        / "control"
        / "qc_thresholds"
        / "qc_threshold_policy.yaml"
    )


def resolve_impact_dir(
    args: argparse.Namespace,
    policy_path: Path,
) -> Path | None:
    if args.impact_dir:
        return args.impact_dir

    if args.results_dir:
        return (
            args.results_dir
            / "control"
            / "qc_threshold_impact"
        )

    candidate = (
        policy_path.parent.parent
        / "qc_threshold_impact"
    )

    if candidate.is_dir():
        return candidate

    return None


def resolve_output_dir(
    args: argparse.Namespace,
    policy: dict[str, Any],
) -> Path:
    if args.output_dir:
        return args.output_dir

    project_id = require_non_empty_string(
        policy.get("project_id"),
        "qc_threshold_policy.project_id",
    )

    config_project_dir = (
        Path("configs")
        / "real_validation"
        / project_id
    )

    if config_project_dir.is_dir():
        return config_project_dir / "review"

    if args.results_dir:
        return (
            args.results_dir
            / "control"
            / "qc_threshold_review"
        )

    return Path("review") / project_id


def parse_scalar(value: str) -> Any:
    trimmed = value.strip()
    lowered = trimmed.lower()

    if lowered in {"true", "false"}:
        return lowered == "true"

    if lowered in {"null", "none"}:
        return None

    try:
        return int(trimmed)
    except ValueError:
        pass

    try:
        return float(trimmed)
    except ValueError:
        return trimmed


def split_assignment(assignment: str) -> tuple[str, str, Any]:
    if "=" not in assignment:
        raise ValueError(
            f"Invalid --set assignment {assignment!r}; expected key=value."
        )

    raw_key, raw_value = assignment.split("=", 1)
    key = raw_key.strip()

    if not key:
        raise ValueError("Invalid --set assignment with empty key.")

    parts = key.split(".")

    if parts[0] == "thresholds":
        parts = parts[1:]

    if len(parts) < 2:
        raise ValueError(
            f"Invalid --set key {key!r}; expected metric.field."
        )

    field = parts[-1]
    metric = ".".join(parts[:-1])

    if field not in ALLOWED_THRESHOLD_FIELDS:
        raise ValueError(
            f"Cannot modify threshold field {field!r}. Allowed fields: "
            + ", ".join(sorted(ALLOWED_THRESHOLD_FIELDS))
        )

    return metric, field, parse_scalar(raw_value)


def validate_threshold_value(
    field: str,
    value: Any,
    assignment: str,
) -> Any:
    if field == "enabled":
        if not isinstance(value, bool):
            raise ValueError(
                f"{assignment}: enabled must be true or false."
            )
        return value

    if value is None:
        return None

    if isinstance(value, bool):
        raise ValueError(
            f"{assignment}: {field} must be numeric or null."
        )

    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"{assignment}: {field} must be numeric or null."
        ) from error

    return parsed


def apply_changes(
    policy: dict[str, Any],
    assignments: list[str],
) -> list[dict[str, Any]]:
    thresholds = require_mapping(
        policy,
        "thresholds",
        "qc_threshold_policy",
    )
    changes = []

    for assignment in assignments:
        metric, field, raw_value = split_assignment(assignment)

        if metric not in thresholds:
            raise ValueError(
                f"Cannot modify unknown threshold metric: {metric}"
            )

        threshold = thresholds[metric]
        if not isinstance(threshold, dict):
            raise ValueError(
                f"thresholds.{metric} must be a mapping."
            )

        value = validate_threshold_value(
            field,
            raw_value,
            assignment,
        )
        previous = threshold.get(field)
        threshold[field] = value

        changes.append(
            {
                "field": f"thresholds.{metric}.{field}",
                "proposed": previous,
                "reviewed": value,
            }
        )

    for metric, threshold in thresholds.items():
        if not isinstance(threshold, dict):
            continue

        minimum = threshold.get("min")
        maximum = threshold.get("max")

        if (
            minimum is not None
            and maximum is not None
            and float(minimum) >= float(maximum)
        ):
            raise ValueError(
                f"Invalid reviewed bounds for {metric}: "
                f"min={minimum}, max={maximum}"
            )

    return changes


def read_summary_tsv(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}

    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if not rows:
        return {}

    return dict(rows[0])


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []

    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def collect_impact_summary(
    impact_dir: Path | None,
) -> dict[str, Any]:
    if not impact_dir:
        return {
            "status": "not_available",
        }

    summary = read_summary_tsv(
        impact_dir / "qc_threshold_impact.summary.tsv"
    )
    by_sample = read_csv_rows(
        impact_dir / "qc_threshold_impact_by_sample.csv"
    )
    by_metric = read_csv_rows(
        impact_dir / "qc_threshold_impact_by_metric.csv"
    )
    lower_tail = read_csv_rows(
        impact_dir / "qc_lower_tail_summary.csv"
    )

    if not summary:
        return {
            "status": "not_available",
            "impact_dir": str(impact_dir),
        }

    top_metrics = [
        row
        for row in by_metric
        if row.get("sample") == "ALL"
        and row.get("Group") == "ALL"
    ]
    top_metrics.sort(
        key=lambda row: int(row.get("flagged_cells", "0") or 0),
        reverse=True,
    )

    return {
        "status": "available",
        "impact_dir": str(impact_dir),
        "cell_count": int(summary.get("cell_count", "0") or 0),
        "flagged_any": int(summary.get("flagged_any", "0") or 0),
        "retained": int(summary.get("retained", "0") or 0),
        "flagged_fraction": float(
            summary.get("flagged_fraction", "0") or 0
        ),
        "by_sample": by_sample,
        "top_metrics": top_metrics,
        "lower_tail": lower_tail,
    }


def write_review_summary(
    path: Path,
    policy: dict[str, Any],
    status: str,
    reviewer: str,
    timestamp: str,
    change_count: int,
    impact: dict[str, Any],
) -> None:
    fieldnames = [
        "plan_id",
        "project_id",
        "biological_run_id",
        "review_status",
        "reviewed_by",
        "reviewed_at",
        "change_count",
        "impact_status",
        "cell_count",
        "flagged_any",
        "retained",
        "flagged_fraction",
    ]

    row = {
        "plan_id": policy.get("plan_id", ""),
        "project_id": policy.get("project_id", ""),
        "biological_run_id": policy.get("biological_run_id", ""),
        "review_status": status,
        "reviewed_by": reviewer,
        "reviewed_at": timestamp,
        "change_count": change_count,
        "impact_status": impact.get("status", "not_available"),
        "cell_count": impact.get("cell_count", ""),
        "flagged_any": impact.get("flagged_any", ""),
        "retained": impact.get("retained", ""),
        "flagged_fraction": impact.get("flagged_fraction", ""),
    }

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(row)


def main() -> int:
    args = parse_args()

    try:
        reviewer = require_non_empty_string(
            args.reviewer,
            "--reviewer",
        )

        if args.decision == "modify-and-approve" and not args.set:
            raise ValueError(
                "--decision modify-and-approve requires at least one --set."
            )

        if args.decision != "modify-and-approve" and args.set:
            raise ValueError(
                "--set is only accepted with "
                "--decision modify-and-approve."
            )

        policy_path = resolve_policy_path(args)
        policy_document = read_yaml(
            policy_path,
            "QC threshold policy",
        )
        proposed_policy = require_mapping(
            policy_document,
            "qc_threshold_policy",
            "policy_document",
        )
        reviewed_policy = copy.deepcopy(proposed_policy)

        impact_dir = resolve_impact_dir(args, policy_path)
        impact = collect_impact_summary(impact_dir)

        reviewed_at = iso_timestamp(args.approved_at)
        review_status = DECISION_TO_STATUS[args.decision]
        changes = apply_changes(
            reviewed_policy,
            args.set,
        )

        approval = dict(
            reviewed_policy.get("approval", {})
        )
        approval.update(
            {
                "required": True,
                "status": review_status,
                "reviewed_by": reviewer,
                "reviewed_at": reviewed_at,
            }
        )

        if args.decision in APPROVING_DECISIONS:
            approval["approved_at"] = reviewed_at
            approval["approved_by"] = reviewer
        else:
            approval["approved_at"] = None
            approval["approved_by"] = None

        reviewed_policy["approval"] = approval
        reviewed_policy["review"] = {
            "decision": args.decision,
            "status": review_status,
            "reviewed_by": reviewer,
            "reviewed_at": reviewed_at,
            "reason": args.reason.strip(),
            "changes": changes,
            "impact_snapshot": impact,
            "produced_by": "review_qc_policy.py",
        }
        reviewed_policy["application"] = {
            **dict(reviewed_policy.get("application", {})),
            "status": "not_applied",
        }

        if args.decision == "approve" and changes:
            raise ValueError(
                "approve cannot include threshold changes."
            )

        output_dir = resolve_output_dir(
            args,
            reviewed_policy,
        )
        output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        reviewed_policy_path = (
            output_dir
            / "qc_threshold_policy.reviewed.yaml"
        )
        summary_path = (
            output_dir
            / "qc_threshold_policy.review.summary.tsv"
        )
        session_path = (
            output_dir
            / "qc_threshold_policy.review.session.txt"
        )

        with reviewed_policy_path.open(
            "w",
            encoding="utf-8",
        ) as handle:
            yaml.safe_dump(
                {
                    "qc_threshold_policy": reviewed_policy,
                },
                handle,
                sort_keys=False,
                allow_unicode=True,
            )

        write_review_summary(
            summary_path,
            reviewed_policy,
            review_status,
            reviewer,
            reviewed_at,
            len(changes),
            impact,
        )

        with session_path.open(
            "w",
            encoding="utf-8",
        ) as handle:
            handle.write(
                f"python_version={platform.python_version()}\n"
            )
            handle.write(f"platform={platform.platform()}\n")
            handle.write("process=review_qc_policy.py\n")
            handle.write(f"review_status={review_status}\n")
            handle.write(f"reviewed_by={reviewer}\n")
            handle.write(f"change_count={len(changes)}\n")

    except (
        FileNotFoundError,
        ValueError,
        yaml.YAMLError,
    ) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    print("QC threshold policy reviewed.")
    print(f"Decision: {args.decision}")
    print(f"Review status: {review_status}")
    print(f"Reviewer: {reviewer}")
    print(f"Reviewed policy: {reviewed_policy_path}")

    if impact.get("status") == "available":
        print(
            "Impact: "
            f"{impact['flagged_any']} / {impact['cell_count']} "
            f"flagged ({impact['flagged_fraction']:.4%}); "
            f"{impact['retained']} retained"
        )
    else:
        print("Impact: not available")

    if args.decision in APPROVING_DECISIONS:
        print(
            "Next step: run qc_threshold_approval_phase with "
            f"--approved_qc_threshold_policy {reviewed_policy_path}"
        )
    else:
        print(
            "Next step: revise thresholds or rerun QC threshold resolution; "
            "this policy is not valid for filtering."
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

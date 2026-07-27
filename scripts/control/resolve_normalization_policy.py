#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import platform
import sys
from pathlib import Path
from typing import Any

try:
    import yaml
except ModuleNotFoundError:  # pragma: no cover - used on minimal hosts.
    yaml = None


DEFAULT_METHOD = "LogNormalize"
DEFAULT_VARIABLE_FEATURES_METHOD = "vst"
DEFAULT_NFEATURES = 2000
DEFAULT_SCALE_FACTOR = 10000
DEFAULT_SCOPE = "per_sample"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Resolve a proposed normalization policy from the "
            "post-doublet QC readiness contract."
        )
    )

    parser.add_argument(
        "--readiness",
        required=True,
        type=Path,
        help="Post-doublet analysis readiness YAML.",
    )
    parser.add_argument(
        "--summary",
        required=True,
        type=Path,
        help="Post-doublet QC global summary CSV.",
    )
    parser.add_argument(
        "--by-sample",
        required=True,
        type=Path,
        help="Post-doublet QC per-sample summary CSV.",
    )
    parser.add_argument(
        "--warnings",
        required=True,
        type=Path,
        help="Post-doublet QC warnings TSV.",
    )
    parser.add_argument(
        "--policy-yaml",
        required=True,
        type=Path,
        help="Output proposed normalization policy YAML.",
    )
    parser.add_argument(
        "--policy-tsv",
        required=True,
        type=Path,
        help="Output proposed normalization policy TSV.",
    )
    parser.add_argument(
        "--resolution-summary",
        required=True,
        type=Path,
        help="Output resolution summary TSV.",
    )
    parser.add_argument(
        "--session-info",
        required=True,
        type=Path,
        help="Output runtime session information.",
    )
    parser.add_argument(
        "--project-id",
        default="",
        help="Stable project identifier.",
    )
    parser.add_argument(
        "--project-name",
        default="",
        help="Project display name.",
    )
    parser.add_argument(
        "--biological-run-id",
        default="",
        help="Stable biological run identifier.",
    )

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
        return float(value)
    except ValueError:
        return value


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
        raise ValueError(f"{label} root must be a mapping.")

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

    if isinstance(value, list):
        handle.write(f"{prefix}{key}:\n")
        for item in value:
            if isinstance(item, dict):
                handle.write(f"{prefix}  -\n")
                for child_key, child_value in item.items():
                    write_yaml_value(
                        handle,
                        str(child_key),
                        child_value,
                        indent + 4,
                    )
            else:
                handle.write(f"{prefix}  - {yaml_scalar(item)}\n")
        return

    handle.write(f"{prefix}{key}: {yaml_scalar(value)}\n")


def write_yaml(path: Path, document: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for key, value in document.items():
            write_yaml_value(handle, key, value)


def read_delimited_rows(
    path: Path,
    label: str,
    delimiter: str = ",",
) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)

        if reader.fieldnames is None:
            raise ValueError(f"{label} has no header: {path}")

        return list(reader)


def metric_values(rows: list[dict[str, str]]) -> dict[str, str]:
    values: dict[str, str] = {}

    for row in rows:
        metric = str(row.get("metric", "")).strip()
        if metric:
            values[metric] = str(row.get("value", "")).strip()

    return values


def parse_int(value: Any, label: str) -> int:
    try:
        parsed = int(float(str(value)))
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be an integer.") from error

    if parsed < 0:
        raise ValueError(f"{label} cannot be negative.")

    return parsed


def get_mapping(document: dict[str, Any], key: str, label: str) -> dict[str, Any]:
    value = document.get(key)

    if not isinstance(value, dict):
        raise ValueError(f"{label}.{key} must be a mapping.")

    return value


def normalize_modalities(value: Any) -> list[str]:
    if value is None:
        return []

    if isinstance(value, list):
        return [str(item) for item in value]

    return [str(value)]


def serialize_modalities(value: Any) -> str:
    return ",".join(normalize_modalities(value))


def resolve_policy(
    readiness_document: dict[str, Any],
    summary: dict[str, str],
    by_sample_rows: list[dict[str, str]],
    warning_rows: list[dict[str, str]],
    args: argparse.Namespace,
) -> dict[str, Any]:
    readiness = get_mapping(
        readiness_document,
        "analysis_readiness",
        "readiness_document",
    )
    dataset = get_mapping(
        readiness,
        "dataset",
        "analysis_readiness",
    )

    readiness_status = str(readiness.get("status", "")).strip().lower()

    if readiness_status != "ready":
        raise ValueError(
            "Normalization policy can only be resolved from a ready "
            f"dataset. Observed status: {readiness_status!r}"
        )

    cells = parse_int(dataset.get("cells"), "dataset.cells")
    features = parse_int(dataset.get("features"), "dataset.features")
    samples = parse_int(dataset.get("samples"), "dataset.samples")
    modalities = serialize_modalities(dataset.get("modalities"))
    warning_count = len(warning_rows)

    if cells != parse_int(summary.get("cells"), "summary.cells"):
        raise ValueError("Readiness cells do not match summary cells.")

    if features != parse_int(summary.get("features"), "summary.features"):
        raise ValueError(
            "Readiness features do not match summary features."
        )

    if samples != len(by_sample_rows):
        raise ValueError(
            "Readiness sample count does not match by-sample rows."
        )

    if warning_count != parse_int(summary.get("warning_count", 0), "warning_count"):
        raise ValueError(
            "Warning count does not match post-doublet QC summary."
        )

    if cells <= 0 or features <= 0 or samples <= 0:
        raise ValueError(
            "Normalization requires positive cells, features and samples."
        )

    integration_required = samples > 1
    integration_reason = (
        "multiple_samples_detected"
        if integration_required
        else "single_sample_dataset"
    )

    return {
        "normalization_policy": {
            "schema_version": "1.0",
            "status": "proposed",
            "identity": {
                "project_id": args.project_id,
                "project_name": args.project_name,
                "biological_run_id": args.biological_run_id,
            },
            "input": {
                "cells": cells,
                "features": features,
                "samples": samples,
                "modalities": modalities,
                "readiness": readiness_status,
                "warning_count": warning_count,
            },
            "strategy": {
                "method": DEFAULT_METHOD,
                "variable_features_method": (
                    DEFAULT_VARIABLE_FEATURES_METHOD
                ),
                "nfeatures": DEFAULT_NFEATURES,
                "scale_factor": DEFAULT_SCALE_FACTOR,
                "scope": DEFAULT_SCOPE,
                "assay": "RNA",
                "regressors": "none",
            },
            "integration": {
                "required": integration_required,
                "reason": integration_reason,
            },
            "rationale": {
                "selected_method_reason": (
                    "stable_default_for_mvp"
                ),
                "scope_reason": (
                    "per_sample_normalization_limits_memory_and_supports_resume"
                ),
                "integration_reason": integration_reason,
            },
            "alternatives": {
                "SCTransform": {
                    "method": "SCTransform",
                    "status": "available",
                    "tradeoff": (
                        "higher_memory_and_modeling_complexity"
                    ),
                }
            },
            "review": {
                "required": True,
                "reason": (
                    "normalization_choices_affect_downstream_structure"
                ),
            },
            "approval": {
                "required": True,
                "status": "proposed",
                "reviewed_by": None,
                "approved_at": None,
            },
            "application": {
                "status": "not_applied",
                "normalization_applied": False,
            },
        }
    }


def write_policy_tsv(path: Path, policy: dict[str, Any]) -> None:
    root = policy["normalization_policy"]
    strategy = root["strategy"]
    integration = root["integration"]

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "project_id",
                "biological_run_id",
                "readiness",
                "cells",
                "features",
                "samples",
                "method",
                "variable_features_method",
                "nfeatures",
                "scale_factor",
                "scope",
                "assay",
                "integration_required",
                "approval_status",
                "application_status",
            ]
        )
        writer.writerow(
            [
                root["identity"]["project_id"],
                root["identity"]["biological_run_id"],
                root["input"]["readiness"],
                root["input"]["cells"],
                root["input"]["features"],
                root["input"]["samples"],
                strategy["method"],
                strategy["variable_features_method"],
                strategy["nfeatures"],
                strategy["scale_factor"],
                strategy["scope"],
                strategy["assay"],
                integration["required"],
                root["approval"]["status"],
                root["application"]["status"],
            ]
        )


def write_resolution_summary(path: Path, policy: dict[str, Any]) -> None:
    root = policy["normalization_policy"]

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "project_id",
                "biological_run_id",
                "readiness",
                "cells",
                "features",
                "samples",
                "method",
                "scope",
                "integration_required",
                "approval_status",
            ]
        )
        writer.writerow(
            [
                root["identity"]["project_id"],
                root["identity"]["biological_run_id"],
                root["input"]["readiness"],
                root["input"]["cells"],
                root["input"]["features"],
                root["input"]["samples"],
                root["strategy"]["method"],
                root["strategy"]["scope"],
                root["integration"]["required"],
                root["approval"]["status"],
            ]
        )


def main() -> int:
    args = parse_args()

    try:
        readiness_document = read_yaml(
            args.readiness,
            "Analysis readiness",
        )
        summary_rows = read_delimited_rows(
            args.summary,
            "Post-doublet QC summary",
        )
        by_sample_rows = read_delimited_rows(
            args.by_sample,
            "Post-doublet QC by-sample summary",
        )
        warning_rows = read_delimited_rows(
            args.warnings,
            "Post-doublet QC warnings",
            delimiter="\t",
        )
        policy = resolve_policy(
            readiness_document,
            metric_values(summary_rows),
            by_sample_rows,
            warning_rows,
            args,
        )

    except (FileNotFoundError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    args.policy_yaml.parent.mkdir(parents=True, exist_ok=True)
    args.policy_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.resolution_summary.parent.mkdir(parents=True, exist_ok=True)
    args.session_info.parent.mkdir(parents=True, exist_ok=True)

    write_yaml(args.policy_yaml, policy)
    write_policy_tsv(args.policy_tsv, policy)
    write_resolution_summary(args.resolution_summary, policy)

    root = policy["normalization_policy"]

    with args.session_info.open("w", encoding="utf-8") as handle:
        handle.write(f"python_version={platform.python_version()}\n")
        handle.write(f"platform={platform.platform()}\n")
        handle.write("resolver=NORMALIZATION_POLICY_RESOLUTION\n")
        handle.write("approval_status=proposed\n")
        handle.write(f"method={root['strategy']['method']}\n")
        handle.write(f"scope={root['strategy']['scope']}\n")

    print("Normalization policy resolved: proposed")
    print(f"Method: {root['strategy']['method']}")
    print(f"Scope: {root['strategy']['scope']}")
    print(f"Samples: {root['input']['samples']}")
    print(f"Integration required: {root['integration']['required']}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

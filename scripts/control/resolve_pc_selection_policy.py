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


DEFAULT_CUMULATIVE_VARIANCE_TARGET = 80.0
DEFAULT_MIN_PCS = 10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Resolve a proposed PC-selection policy from dimensionality "
            "reduction metrics."
        )
    )
    parser.add_argument(
        "--dimensionality-summary",
        required=True,
        type=Path,
        help="DIMENSIONALITY summary CSV.",
    )
    parser.add_argument(
        "--variable-features",
        required=True,
        type=Path,
        help="Variable feature summary CSV.",
    )
    parser.add_argument(
        "--elbow-metrics",
        required=True,
        type=Path,
        help="PCA elbow metrics CSV.",
    )
    parser.add_argument(
        "--policy-yaml",
        required=True,
        type=Path,
        help="Output proposed PC-selection policy YAML.",
    )
    parser.add_argument(
        "--policy-tsv",
        required=True,
        type=Path,
        help="Output proposed PC-selection policy TSV.",
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
    parser.add_argument(
        "--cumulative-variance-target",
        type=float,
        default=DEFAULT_CUMULATIVE_VARIANCE_TARGET,
        help="Cumulative explained-variance target for the proposed PCs.",
    )
    parser.add_argument(
        "--min-pcs",
        type=int,
        default=DEFAULT_MIN_PCS,
        help="Minimum number of PCs to propose when enough PCs exist.",
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


def parse_float(value: Any, label: str) -> float:
    try:
        parsed = float(str(value))
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be numeric.") from error

    if not parsed >= 0:
        raise ValueError(f"{label} cannot be negative.")

    return parsed


def pc_index(value: Any, label: str) -> int:
    text = str(value).strip()

    if text.startswith("PC_"):
        text = text[3:]

    parsed = parse_int(text, label)

    if parsed <= 0:
        raise ValueError(f"{label} must identify a positive PC.")

    return parsed


def resolve_recommended_pcs(
    elbow_rows: list[dict[str, str]],
    target: float,
    min_pcs: int,
) -> tuple[int, str, float]:
    if target <= 0 or target > 100:
        raise ValueError("cumulative variance target must be in (0, 100].")

    if min_pcs <= 0:
        raise ValueError("min_pcs must be positive.")

    parsed_rows = []

    for row in elbow_rows:
        parsed_rows.append(
            {
                "pc": pc_index(row.get("pc"), "elbow.pc"),
                "cumulative": parse_float(
                    row.get("cumulative_variance_percent"),
                    "elbow.cumulative_variance_percent",
                ),
                "variance_percent": parse_float(
                    row.get("variance_percent"),
                    "elbow.variance_percent",
                ),
            }
        )

    if len(parsed_rows) < 2:
        raise ValueError("At least two PCA metric rows are required.")

    max_computed = max(row["pc"] for row in parsed_rows)

    if len({row["pc"] for row in parsed_rows}) != len(parsed_rows):
        raise ValueError("PCA metric rows contain duplicated PCs.")

    sorted_rows = sorted(parsed_rows, key=lambda row: row["pc"])

    recommended = max_computed
    achieved = sorted_rows[-1]["cumulative"]

    for row in sorted_rows:
        if row["cumulative"] >= target:
            recommended = row["pc"]
            achieved = row["cumulative"]
            break

    recommended = max(min(recommended, max_computed), min(min_pcs, max_computed))
    selected = next(row for row in sorted_rows if row["pc"] == recommended)
    achieved = selected["cumulative"]

    return (
        recommended,
        "cumulative_variance_target_with_minimum_pc_floor",
        achieved,
    )


def resolve_policy(
    summary: dict[str, str],
    variable_feature_rows: list[dict[str, str]],
    elbow_rows: list[dict[str, str]],
    args: argparse.Namespace,
) -> dict[str, Any]:
    cells_input = parse_int(summary.get("cells_input"), "summary.cells_input")
    cells_output = parse_int(summary.get("cells_output"), "summary.cells_output")
    features_input = parse_int(
        summary.get("features_input"),
        "summary.features_input",
    )
    features_output = parse_int(
        summary.get("features_output"),
        "summary.features_output",
    )
    npcs_computed = parse_int(
        summary.get("npcs_computed"),
        "summary.npcs_computed",
    )
    variable_features = parse_int(
        summary.get("variable_features"),
        "summary.variable_features",
    )

    if cells_input != cells_output:
        raise ValueError("Dimensionality changed cell count.")

    if features_input != features_output:
        raise ValueError("Dimensionality changed feature count.")

    if str(summary.get("scaling_applied", "")).lower() != "true":
        raise ValueError("Dimensionality summary does not mark scaling applied.")

    if str(summary.get("pca_applied", "")).lower() != "true":
        raise ValueError("Dimensionality summary does not mark PCA applied.")

    if str(summary.get("normalization_application_status", "")).lower() != "applied":
        raise ValueError("Normalization must be applied before PC selection.")

    if npcs_computed != len(elbow_rows):
        raise ValueError("npcs_computed does not match elbow metric rows.")

    if variable_features != len(variable_feature_rows):
        raise ValueError(
            "variable_features does not match variable feature rows."
        )

    if cells_input <= 0 or features_input <= 0:
        raise ValueError("PC selection requires positive cells and features.")

    recommended_pcs, method, achieved_cumulative = resolve_recommended_pcs(
        elbow_rows,
        args.cumulative_variance_target,
        args.min_pcs,
    )

    if recommended_pcs > npcs_computed:
        raise ValueError("Recommended PCs cannot exceed computed PCs.")

    return {
        "pc_selection_policy": {
            "schema_version": "1.0",
            "status": "proposed",
            "identity": {
                "project_id": args.project_id,
                "project_name": args.project_name,
                "biological_run_id": args.biological_run_id,
            },
            "input": {
                "cells": cells_input,
                "features": features_input,
                "variable_features": variable_features,
                "max_computed_pcs": npcs_computed,
                "assay": summary.get("assay", ""),
                "normalization_method": summary.get(
                    "normalization_method",
                    "",
                ),
            },
            "recommendation": {
                "recommended_pcs": recommended_pcs,
                "method": method,
                "cumulative_variance_target": args.cumulative_variance_target,
                "achieved_cumulative_variance": achieved_cumulative,
                "min_pcs": args.min_pcs,
                "max_computed_pcs": npcs_computed,
            },
            "rationale": {
                "selected_method_reason": (
                    "uses_published_pca_elbow_metrics_without_running_"
                    "neighbors_umap_clustering_or_integration"
                ),
                "review_reason": (
                    "pc_count_affects_neighbors_umap_clustering_and_"
                    "integration"
                ),
            },
            "allowed_downstream_uses": {
                "neighbors": True,
                "umap": True,
                "clustering": True,
                "integration": True,
            },
            "review": {
                "required": True,
                "reason": (
                    "selected_pcs_affect_downstream_manifold_and_cluster_"
                    "structure"
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
                "pc_selection_applied": False,
            },
        }
    }


def write_policy_tsv(path: Path, policy: dict[str, Any]) -> None:
    root = policy["pc_selection_policy"]
    recommendation = root["recommendation"]

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "project_id",
                "biological_run_id",
                "cells",
                "features",
                "variable_features",
                "recommended_pcs",
                "max_computed_pcs",
                "method",
                "achieved_cumulative_variance",
                "approval_status",
                "application_status",
            ]
        )
        writer.writerow(
            [
                root["identity"]["project_id"],
                root["identity"]["biological_run_id"],
                root["input"]["cells"],
                root["input"]["features"],
                root["input"]["variable_features"],
                recommendation["recommended_pcs"],
                recommendation["max_computed_pcs"],
                recommendation["method"],
                recommendation["achieved_cumulative_variance"],
                root["approval"]["status"],
                root["application"]["status"],
            ]
        )


def write_resolution_summary(path: Path, policy: dict[str, Any]) -> None:
    root = policy["pc_selection_policy"]
    recommendation = root["recommendation"]

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "project_id",
                "biological_run_id",
                "recommended_pcs",
                "max_computed_pcs",
                "method",
                "approval_status",
            ]
        )
        writer.writerow(
            [
                root["identity"]["project_id"],
                root["identity"]["biological_run_id"],
                recommendation["recommended_pcs"],
                recommendation["max_computed_pcs"],
                recommendation["method"],
                root["approval"]["status"],
            ]
        )


def main() -> int:
    args = parse_args()

    try:
        summary_rows = read_delimited_rows(
            args.dimensionality_summary,
            "Dimensionality summary",
        )
        variable_feature_rows = read_delimited_rows(
            args.variable_features,
            "Variable feature summary",
        )
        elbow_rows = read_delimited_rows(
            args.elbow_metrics,
            "PCA elbow metrics",
        )
        policy = resolve_policy(
            metric_values(summary_rows),
            variable_feature_rows,
            elbow_rows,
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

    root = policy["pc_selection_policy"]
    recommendation = root["recommendation"]

    with args.session_info.open("w", encoding="utf-8") as handle:
        handle.write(f"python_version={platform.python_version()}\n")
        handle.write(f"platform={platform.platform()}\n")
        handle.write("resolver=PC_SELECTION_POLICY_RESOLUTION\n")
        handle.write("approval_status=proposed\n")
        handle.write(f"recommended_pcs={recommendation['recommended_pcs']}\n")
        handle.write(f"max_computed_pcs={recommendation['max_computed_pcs']}\n")

    print("PC-selection policy resolved: proposed")
    print(f"Recommended PCs: {recommendation['recommended_pcs']}")
    print(f"Max computed PCs: {recommendation['max_computed_pcs']}")
    print(f"Method: {recommendation['method']}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

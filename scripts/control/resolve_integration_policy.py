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


ALLOWED_METHODS = {
    "none",
    "Harmony",
}

PLANNED_METHODS = {
    "SeuratIntegration",
    "WNN",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Resolve a proposed integration policy from approved PC "
            "selection, readiness and normalization contracts."
        )
    )
    parser.add_argument(
        "--pca-object",
        required=True,
        type=Path,
        help="PCA Seurat object path used as provenance input.",
    )
    parser.add_argument(
        "--dimensionality-summary",
        required=True,
        type=Path,
        help="DIMENSIONALITY summary CSV.",
    )
    parser.add_argument(
        "--approved-pc-selection-policy",
        required=True,
        type=Path,
        help="Validated approved PC-selection policy YAML.",
    )
    parser.add_argument(
        "--analysis-readiness",
        required=True,
        type=Path,
        help="Post-doublet analysis readiness YAML.",
    )
    parser.add_argument(
        "--applied-normalization-policy",
        required=True,
        type=Path,
        help="Applied normalization policy YAML.",
    )
    parser.add_argument(
        "--policy-yaml",
        required=True,
        type=Path,
        help="Output proposed integration policy YAML.",
    )
    parser.add_argument(
        "--policy-tsv",
        required=True,
        type=Path,
        help="Output proposed integration policy TSV.",
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
    parser.add_argument("--project-id", default="")
    parser.add_argument("--project-name", default="")
    parser.add_argument("--biological-run-id", default="")
    parser.add_argument(
        "--default-batch-variable",
        default="sample",
        help="Default batch variable proposed when integration is required.",
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
    stack: list[tuple[int, Any, dict[str, Any] | None, str | None]] = [
        (-1, root, None, None)
    ]

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip()

            if not line.strip() or line.lstrip().startswith("#"):
                continue

            stripped = line.lstrip()

            if stripped in {"---", "..."}:
                continue

            indent = len(line) - len(stripped)

            if stripped.startswith("- "):
                while stack and indent <= stack[-1][0]:
                    stack.pop()

                if not stack:
                    raise ValueError(
                        f"Invalid YAML list indentation at {path}:{line_number}"
                    )

                current_indent, current, parent, parent_key = stack[-1]

                if isinstance(current, dict):
                    if current or parent is None or parent_key is None:
                        raise ValueError(
                            f"Unsupported YAML list syntax at {path}:{line_number}"
                        )

                    replacement: list[Any] = []
                    parent[parent_key] = replacement
                    stack[-1] = (
                        current_indent,
                        replacement,
                        parent,
                        parent_key,
                    )
                    current = replacement

                if not isinstance(current, list):
                    raise ValueError(
                        f"Unsupported YAML list syntax at {path}:{line_number}"
                    )

                item = stripped[2:].strip()
                current.append(parse_scalar(item))
                continue


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

            parent_container = stack[-1][1]

            if not isinstance(parent_container, dict):
                raise ValueError(
                    f"Invalid YAML mapping inside list at {path}:{line_number}"
                )

            if value.strip() == "":
                child: dict[str, Any] = {}
                parent_container[key] = child
                stack.append((indent, child, parent_container, key))
            else:
                parent_container[key] = parse_scalar(value)

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


def read_delimited_rows(path: Path, label: str) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)

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


def get_mapping(document: dict[str, Any], key: str, label: str) -> dict[str, Any]:
    value = document.get(key)

    if not isinstance(value, dict):
        raise ValueError(f"{label}.{key} must be a mapping.")

    return value


def parse_int(value: Any, label: str) -> int:
    try:
        parsed = int(float(str(value)))
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be an integer.") from error

    if parsed < 0:
        raise ValueError(f"{label} cannot be negative.")

    return parsed


def normalize_modalities(value: Any) -> list[str]:
    if value is None:
        return []

    if isinstance(value, list):
        return [str(item) for item in value]

    return [part.strip() for part in str(value).split(",") if part.strip()]


def validate_identity(root: dict[str, Any], args: argparse.Namespace) -> None:
    identity = get_mapping(root, "identity", "policy")
    project_id = str(identity.get("project_id", "")).strip()
    run_id = str(identity.get("biological_run_id", "")).strip()

    if args.project_id and project_id != args.project_id:
        raise ValueError("Policy project_id does not match execution metadata.")

    if args.biological_run_id and run_id != args.biological_run_id:
        raise ValueError(
            "Policy biological_run_id does not match execution metadata."
        )


def validate_pc_policy(document: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    root = get_mapping(document, "pc_selection_policy", "pc_selection_document")
    validate_identity(root, args)

    approval = get_mapping(root, "approval", "pc_selection_policy")
    validation = get_mapping(root, "validation", "pc_selection_policy")
    application = get_mapping(root, "application", "pc_selection_policy")

    if str(approval.get("status", "")).strip().lower() not in {
        "approved",
        "modified_and_approved",
    }:
        raise ValueError("PC-selection policy must be approved.")

    if str(validation.get("status", "")).strip().lower() != "valid":
        raise ValueError("PC-selection policy validation.status must be valid.")

    if validation.get("validated_for_pc_selection") is not True:
        raise ValueError(
            "PC-selection policy has not been validated for PC selection."
        )

    if str(application.get("status", "")).strip().lower() != "not_applied":
        raise ValueError("PC-selection policy application.status must be not_applied.")

    selected_pcs = parse_int(root.get("selected_pcs"), "selected_pcs")
    max_computed_pcs = parse_int(
        validation.get("max_computed_pcs"),
        "validation.max_computed_pcs",
    )

    if selected_pcs <= 0:
        raise ValueError("selected_pcs must be positive.")

    if selected_pcs > max_computed_pcs:
        raise ValueError("selected_pcs cannot exceed max_computed_pcs.")

    return root


def validate_normalization_policy(
    document: dict[str, Any],
    args: argparse.Namespace,
) -> dict[str, Any]:
    root = get_mapping(
        document,
        "normalization_policy",
        "normalization_document",
    )
    validate_identity(root, args)

    validation = get_mapping(root, "validation", "normalization_policy")
    application = get_mapping(root, "application", "normalization_policy")

    if str(validation.get("status", "")).strip().lower() != "valid":
        raise ValueError("Normalization policy validation.status must be valid.")

    if str(application.get("status", "")).strip().lower() != "applied":
        raise ValueError("Normalization policy application.status must be applied.")

    return root


def validate_readiness(document: dict[str, Any]) -> dict[str, Any]:
    root = get_mapping(document, "analysis_readiness", "readiness_document")
    status = str(root.get("status", "")).strip().lower()

    if status != "ready":
        raise ValueError("Analysis readiness must be ready.")

    return root


def choose_strategy(
    samples: int,
    modalities: list[str],
    normalization_root: dict[str, Any],
    args: argparse.Namespace,
) -> tuple[bool, str, str, str, str]:
    normalization_integration = get_mapping(
        normalization_root,
        "integration",
        "normalization_policy",
    )

    derived_required = samples > 1
    upstream_required = normalization_integration.get("required")
    required = derived_required or upstream_required is True

    if not required:
        return (
            False,
            "none",
            "pca",
            "pca",
            "single_sample_or_no_batch_structure_detected",
        )

    if "RNA" not in modalities:
        return (
            False,
            "none",
            "pca",
            "pca",
            "rna_modality_not_available_for_harmony_mvp",
        )

    batch_variable = args.default_batch_variable.strip() or "sample"

    return (
        True,
        "Harmony",
        "pca",
        "harmony",
        f"multiple_samples_detected_use_{batch_variable}_as_initial_batch_variable",
    )


def resolve_policy(
    pca_object: Path,
    dimensionality_summary: dict[str, str],
    pc_policy: dict[str, Any],
    readiness: dict[str, Any],
    normalization_root: dict[str, Any],
    args: argparse.Namespace,
) -> dict[str, Any]:
    if not pca_object.exists():
        raise FileNotFoundError(f"PCA object not found: {pca_object}")

    dataset = get_mapping(readiness, "dataset", "analysis_readiness")
    samples = parse_int(dataset.get("samples"), "dataset.samples")
    cells = parse_int(dataset.get("cells"), "dataset.cells")
    features = parse_int(dataset.get("features"), "dataset.features")
    modalities = normalize_modalities(dataset.get("modalities"))

    summary_cells = parse_int(
        dimensionality_summary.get("cells_output"),
        "dimensionality_summary.cells_output",
    )
    summary_features = parse_int(
        dimensionality_summary.get("features_output"),
        "dimensionality_summary.features_output",
    )
    npcs_computed = parse_int(
        dimensionality_summary.get("npcs_computed"),
        "dimensionality_summary.npcs_computed",
    )

    if cells != summary_cells:
        raise ValueError("Readiness cells do not match dimensionality cells.")

    if features != summary_features:
        raise ValueError(
            "Readiness features do not match dimensionality features."
        )

    selected_pcs = parse_int(pc_policy.get("selected_pcs"), "selected_pcs")

    if selected_pcs > npcs_computed:
        raise ValueError("selected_pcs cannot exceed computed PCA dimensions.")

    required, method, input_reduction, output_reduction, reason = choose_strategy(
        samples,
        modalities,
        normalization_root,
        args,
    )

    batch_variable = (
        args.default_batch_variable.strip()
        if required and method == "Harmony"
        else None
    )

    return {
        "integration_policy": {
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
                "modalities": ",".join(modalities),
                "pca_object": str(pca_object),
                "selected_pcs": selected_pcs,
                "max_computed_pcs": npcs_computed,
                "pc_selection_validation_status": "valid",
                "normalization_application_status": "applied",
            },
            "strategy": {
                "integration_required": required,
                "method": method,
                "input_reduction": input_reduction,
                "output_reduction": output_reduction,
                "dimensions_to_use": selected_pcs,
                "batch_variable": batch_variable,
                "grouping_variable": None,
                "multimodal_strategy": "none",
            },
            "rationale": {
                "reason": reason,
                "selected_pcs_source": "approved_pc_selection_policy",
                "method_scope": "mvp_allowed_methods_none_or_harmony",
            },
            "evidence": {
                "samples_detected": samples,
                "modalities_detected": ",".join(modalities),
                "selected_pcs": selected_pcs,
                "max_computed_pcs": npcs_computed,
                "normalization_method": normalization_root.get(
                    "strategy",
                    {},
                ).get("method", ""),
                "integration_required_reason": reason,
                "batch_variable_source": (
                    "default_mvp_sample_metadata"
                    if batch_variable
                    else "not_required"
                ),
                "detected_causistics": "not_available_in_current_inputs",
            },
            "allowed_methods": {
                "available": sorted(ALLOWED_METHODS),
                "planned": sorted(PLANNED_METHODS),
            },
            "review": {
                "required": True,
                "reason": (
                    "integration_strategy_changes_the_downstream_"
                    "representation_used_for_neighbors_umap_and_clustering"
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
                "integration_applied": False,
            },
        }
    }


def write_policy_tsv(path: Path, policy: dict[str, Any]) -> None:
    root = policy["integration_policy"]
    strategy = root["strategy"]

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "project_id",
                "biological_run_id",
                "samples",
                "selected_pcs",
                "integration_required",
                "method",
                "input_reduction",
                "output_reduction",
                "batch_variable",
                "approval_status",
                "application_status",
            ]
        )
        writer.writerow(
            [
                root["identity"]["project_id"],
                root["identity"]["biological_run_id"],
                root["input"]["samples"],
                root["input"]["selected_pcs"],
                strategy["integration_required"],
                strategy["method"],
                strategy["input_reduction"],
                strategy["output_reduction"],
                strategy["batch_variable"],
                root["approval"]["status"],
                root["application"]["status"],
            ]
        )


def write_resolution_summary(path: Path, policy: dict[str, Any]) -> None:
    root = policy["integration_policy"]
    strategy = root["strategy"]

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "project_id",
                "biological_run_id",
                "integration_required",
                "method",
                "dimensions_to_use",
                "approval_status",
            ]
        )
        writer.writerow(
            [
                root["identity"]["project_id"],
                root["identity"]["biological_run_id"],
                strategy["integration_required"],
                strategy["method"],
                strategy["dimensions_to_use"],
                root["approval"]["status"],
            ]
        )


def main() -> int:
    args = parse_args()

    try:
        dimensionality_summary = metric_values(
            read_delimited_rows(
                args.dimensionality_summary,
                "Dimensionality summary",
            )
        )
        pc_policy_document = read_yaml(
            args.approved_pc_selection_policy,
            "Approved PC-selection policy",
        )
        readiness_document = read_yaml(
            args.analysis_readiness,
            "Analysis readiness",
        )
        normalization_document = read_yaml(
            args.applied_normalization_policy,
            "Applied normalization policy",
        )

        pc_policy = validate_pc_policy(pc_policy_document, args)
        readiness = validate_readiness(readiness_document)
        normalization_root = validate_normalization_policy(
            normalization_document,
            args,
        )
        policy = resolve_policy(
            args.pca_object,
            dimensionality_summary,
            pc_policy,
            readiness,
            normalization_root,
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

    root = policy["integration_policy"]
    strategy = root["strategy"]

    with args.session_info.open("w", encoding="utf-8") as handle:
        handle.write(f"python_version={platform.python_version()}\n")
        handle.write(f"platform={platform.platform()}\n")
        handle.write("resolver=INTEGRATION_POLICY_RESOLUTION\n")
        handle.write("approval_status=proposed\n")
        handle.write(f"method={strategy['method']}\n")
        handle.write(
            f"integration_required={strategy['integration_required']}\n"
        )
        handle.write(f"dimensions_to_use={strategy['dimensions_to_use']}\n")

    print("Integration policy resolved: proposed")
    print(f"Method: {strategy['method']}")
    print(f"Integration required: {strategy['integration_required']}")
    print(f"Dimensions to use: {strategy['dimensions_to_use']}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

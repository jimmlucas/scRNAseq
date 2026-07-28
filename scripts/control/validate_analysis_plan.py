#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Any

import yaml


ALLOWED_EXPERIMENT_TYPES = {
    "scRNAseq",
    "snRNAseq",
    "CITEseq",
    "multiome",
    "spatial",
    "unknown",
}

ALLOWED_MATERIALS = {
    "cells",
    "nuclei",
    "mixed",
    "unknown",
}

ALLOWED_ORGANISMS = {
    "human",
    "mouse",
    "other",
    "unknown",
}

ALLOWED_MATRIX_TYPES = {
    "raw",
    "filtered",
    "seurat_object",
    "single_cell_experiment",
    "unknown",
}

ALLOWED_THRESHOLD_MODES = {
    "fixed",
    "data_driven",
    "adaptive",
    "hybrid",
}

DEPRECATED_THRESHOLD_MODE_ALIASES = {
    "adaptive": "data_driven",
}

ROUTING_MODULES = {
    "input_discovery",
    "load_data",
    "qc_metrics",
    "qc_visualization",
    "qc_threshold_resolution",
    "qc_filter",
    "doublet_detection",
    "sample_resolution",
    "merge_samples",
    "qc_scores",
    "normalization",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate a single-cell analysis plan."
    )
    parser.add_argument(
        "--plan",
        required=True,
        type=Path,
        help="Input analysis plan in YAML format.",
    )
    parser.add_argument(
        "--validated-plan",
        required=True,
        type=Path,
        help="Normalized validated plan output.",
    )
    parser.add_argument(
        "--summary",
        required=True,
        type=Path,
        help="Validation summary TSV.",
    )
    return parser.parse_args()


def require_mapping(
    value: Any,
    name: str,
    errors: list[str],
) -> dict[str, Any]:
    if not isinstance(value, dict):
        errors.append(f"Section '{name}' must be a mapping.")
        return {}

    return value


def require_value(
    mapping: dict[str, Any],
    key: str,
    section: str,
    errors: list[str],
) -> Any:
    value = mapping.get(key)

    if value is None or value == "":
        errors.append(f"Missing required value: {section}.{key}")

    return value


def validate_allowed(
    value: Any,
    allowed: set[str],
    field: str,
    errors: list[str],
) -> None:
    if value is not None and value not in allowed:
        errors.append(
            f"Invalid value for {field}: '{value}'. "
            f"Allowed values: {', '.join(sorted(allowed))}"
        )


def main() -> int:
    args = parse_args()

    if not args.plan.exists():
        print(
            f"Analysis plan not found: {args.plan}",
            file=sys.stderr,
        )
        return 2

    with args.plan.open("r", encoding="utf-8") as handle:
        document = yaml.safe_load(handle)

    errors: list[str] = []
    warnings: list[str] = []

    if not isinstance(document, dict):
        errors.append("The analysis plan root must be a mapping.")
        document = {}

    plan = require_mapping(
        document.get("plan"),
        "plan",
        errors,
    )
    experiment = require_mapping(
        document.get("experiment"),
        "experiment",
        errors,
    )
    inputs = require_mapping(
        document.get("inputs"),
        "inputs",
        errors,
    )
    routing = require_mapping(
        document.get("routing"),
        "routing",
        errors,
    )
    qc = require_mapping(
        document.get("qc"),
        "qc",
        errors,
    )

    schema_version = require_value(
        plan,
        "schema_version",
        "plan",
        errors,
    )
    plan_id = require_value(
        plan,
        "plan_id",
        "plan",
        errors,
    )

    experiment_type = require_value(
        experiment,
        "experiment_type",
        "experiment",
        errors,
    )
    material = require_value(
        experiment,
        "material",
        "experiment",
        errors,
    )
    organism = require_value(
        experiment,
        "organism",
        "experiment",
        errors,
    )
    matrix_type = require_value(
        inputs,
        "matrix_type",
        "inputs",
        errors,
    )
    threshold_mode = require_value(
        qc,
        "threshold_mode",
        "qc",
        errors,
    )

    validate_allowed(
        experiment_type,
        ALLOWED_EXPERIMENT_TYPES,
        "experiment.experiment_type",
        errors,
    )
    validate_allowed(
        material,
        ALLOWED_MATERIALS,
        "experiment.material",
        errors,
    )
    validate_allowed(
        organism,
        ALLOWED_ORGANISMS,
        "experiment.organism",
        errors,
    )
    validate_allowed(
        matrix_type,
        ALLOWED_MATRIX_TYPES,
        "inputs.matrix_type",
        errors,
    )
    validate_allowed(
        threshold_mode,
        ALLOWED_THRESHOLD_MODES,
        "qc.threshold_mode",
        errors,
    )

    if threshold_mode in DEPRECATED_THRESHOLD_MODE_ALIASES:
        warnings.append(
            "qc.threshold_mode='adaptive' is deprecated; "
            "use 'data_driven'."
        )

    unknown_routes = set(routing) - ROUTING_MODULES

    if unknown_routes:
        errors.append(
            "Unknown routing modules: "
            + ", ".join(sorted(unknown_routes))
        )

    for module_name in ROUTING_MODULES:
        module_value = routing.get(module_name, False)

        if not isinstance(module_value, bool):
            errors.append(
                f"routing.{module_name} must be true or false."
            )

    if routing.get("qc_filter") and not routing.get("qc_metrics"):
        errors.append(
            "routing.qc_filter requires routing.qc_metrics=true."
        )

    if (
        routing.get("qc_threshold_resolution")
        and not routing.get("qc_metrics")
    ):
        errors.append(
            "routing.qc_threshold_resolution requires "
            "routing.qc_metrics=true."
        )

    if (
        routing.get("sample_resolution")
        and not routing.get("load_data")
    ):
        errors.append(
            "routing.sample_resolution requires routing.load_data=true."
        )

    if routing.get("merge_samples") and not routing.get("load_data"):
        errors.append(
            "routing.merge_samples requires routing.load_data=true."
        )

    if routing.get("normalization") and not routing.get("qc_filter"):
        warnings.append(
            "Normalization is enabled while QC filtering is disabled."
        )

    if schema_version != "1.0":
        warnings.append(
            f"Plan schema version '{schema_version}' "
            "differs from supported version '1.0'."
        )

    validation_status = "valid" if not errors else "invalid"

    normalized_document = dict(document)
    normalized_document.setdefault("validation", {})
    normalized_document["validation"] = {
        "status": validation_status,
        "errors": errors,
        "warnings": warnings,
    }

    args.validated_plan.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    args.summary.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with args.validated_plan.open(
        "w",
        encoding="utf-8",
    ) as handle:
        yaml.safe_dump(
            normalized_document,
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
                "status",
                "error_count",
                "warning_count",
            ]
        )
        writer.writerow(
            [
                plan_id or "",
                validation_status,
                len(errors),
                len(warnings),
            ]
        )

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)

        return 1

    for warning in warnings:
        print(f"WARNING: {warning}", file=sys.stderr)

    print(
        f"Analysis plan '{plan_id}' validated successfully."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

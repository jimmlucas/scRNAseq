#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import platform
import sys
from pathlib import Path
from typing import Any

import yaml


SUPPORTED_THRESHOLD_MODES = {
    "fixed",
    "data_driven",
    "hybrid",
}


METRIC_COLUMNS = {
    "nFeature_RNA": "nFeature_RNA",
    "nCount_RNA": "nCount_RNA",
    "percent.mito": "percent.mito",
    "percent.ribo": "percent.ribo",
    "percent.globin": "percent.globin",
    "log10GenesPerUMI": "log10GenesPerUMI",
}


# Metrics that may generate filtering bounds automatically.
# Ribosomal and globin percentages remain diagnostic unless the
# analysis plan provides explicit thresholds for them.
AUTO_FILTER_METRICS = {
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mito",
    "log10GenesPerUMI",
}


DEFAULT_BOUNDS = {
    "cells": {
        "nFeature_RNA": {
            "min": 200.0,
            "max": None,
        },
        "nCount_RNA": {
            "min": 500.0,
            "max": None,
        },
        "percent.mito": {
            "min": None,
            "max": 20.0,
        },
        "percent.ribo": {
            "min": None,
            "max": None,
        },
        "percent.globin": {
            "min": None,
            "max": None,
        },
        "log10GenesPerUMI": {
            "min": 0.80,
            "max": None,
        },
    },
    "nuclei": {
        "nFeature_RNA": {
            "min": 100.0,
            "max": None,
        },
        "nCount_RNA": {
            "min": 200.0,
            "max": None,
        },
        "percent.mito": {
            "min": None,
            "max": 10.0,
        },
        "percent.ribo": {
            "min": None,
            "max": None,
        },
        "percent.globin": {
            "min": None,
            "max": None,
        },
        "log10GenesPerUMI": {
            "min": 0.75,
            "max": None,
        },
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Resolve a proposed QC threshold policy from a validated "
            "analysis plan and observed QC metrics."
        )
    )

    parser.add_argument(
        "--plan",
        required=True,
        type=Path,
        help="Validated analysis plan YAML.",
    )
    parser.add_argument(
        "--global-metrics",
        required=True,
        type=Path,
        help="Global QC metrics CSV.",
    )
    parser.add_argument(
        "--sample-metrics",
        required=True,
        type=Path,
        help="Per-sample QC metrics CSV.",
    )
    parser.add_argument(
        "--cell-metadata",
        required=True,
        type=Path,
        help="Cell-level QC metadata CSV.",
    )
    parser.add_argument(
        "--policy-yaml",
        required=True,
        type=Path,
        help="Output proposed QC policy YAML.",
    )
    parser.add_argument(
        "--policy-tsv",
        required=True,
        type=Path,
        help="Output proposed QC policy TSV.",
    )
    parser.add_argument(
        "--summary",
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

    return parser.parse_args()


def read_yaml(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")

    with path.open("r", encoding="utf-8") as handle:
        document = yaml.safe_load(handle)

    if not isinstance(document, dict):
        raise ValueError(f"{label} root must be a mapping: {path}")

    return document


def read_csv_rows(
    path: Path,
    label: str,
) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")

    with path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            raise ValueError(f"{label} has no header: {path}")

        return list(reader)


def parse_float(
    value: Any,
    default: float | None = None,
) -> float | None:
    if value is None:
        return default

    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return default

    if not math.isfinite(parsed):
        return default

    return parsed


def percentile(
    values: list[float],
    probability: float,
) -> float | None:
    if not values:
        return None

    ordered = sorted(values)

    if len(ordered) == 1:
        return ordered[0]

    position = probability * (len(ordered) - 1)
    lower_index = math.floor(position)
    upper_index = math.ceil(position)

    if lower_index == upper_index:
        return ordered[lower_index]

    fraction = position - lower_index

    return (
        ordered[lower_index]
        + (
            ordered[upper_index]
            - ordered[lower_index]
        )
        * fraction
    )


def iqr_bounds(
    values: list[float],
    multiplier: float = 3.0,
) -> tuple[float | None, float | None]:
    q1 = percentile(values, 0.25)
    q3 = percentile(values, 0.75)

    if q1 is None or q3 is None:
        return None, None

    spread = q3 - q1

    return (
        q1 - multiplier * spread,
        q3 + multiplier * spread,
    )


def collect_column(
    rows: list[dict[str, str]],
    column: str,
) -> list[float]:
    values: list[float] = []

    for row in rows:
        value = parse_float(row.get(column))

        if value is not None:
            values.append(value)

    return values


def get_nested(
    document: dict[str, Any],
    path: list[str],
    default: Any = None,
) -> Any:
    current: Any = document

    for key in path:
        if not isinstance(current, dict):
            return default

        if key not in current:
            return default

        current = current[key]

    return current


def resolve_material(
    plan: dict[str, Any],
) -> str:
    material = str(
        get_nested(
            plan,
            ["experiment", "material"],
            "cells",
        )
    ).strip().lower()

    if material in {
        "nucleus",
        "nuclei",
        "snrna",
        "snrna-seq",
        "snrnaseq",
    }:
        return "nuclei"

    return "cells"


def resolve_threshold_mode(
    plan: dict[str, Any],
) -> str:
    mode = str(
        get_nested(
            plan,
            ["qc", "threshold_mode"],
            "hybrid",
        )
    ).strip().lower()

    if mode not in SUPPORTED_THRESHOLD_MODES:
        raise ValueError(
            "Unsupported qc.threshold_mode: "
            f"{mode}. Supported values: "
            + ", ".join(sorted(SUPPORTED_THRESHOLD_MODES))
        )

    return mode


def resolve_explicit_thresholds(
    plan: dict[str, Any],
) -> dict[str, dict[str, float | None]]:
    explicit = get_nested(
        plan,
        ["qc", "thresholds"],
        {},
    )

    if not isinstance(explicit, dict):
        return {}

    resolved: dict[str, dict[str, float | None]] = {}

    for metric_name, bounds in explicit.items():
        if metric_name not in METRIC_COLUMNS:
            continue

        if not isinstance(bounds, dict):
            continue

        resolved[metric_name] = {
            "min": parse_float(bounds.get("min")),
            "max": parse_float(bounds.get("max")),
        }

    return resolved


def propose_data_driven_bounds(
    metric_name: str,
    values: list[float],
) -> dict[str, float | None]:
    lower_iqr, upper_iqr = iqr_bounds(
        values,
        multiplier=3.0,
    )

    if metric_name in {
        "percent.mito",
        "percent.ribo",
        "percent.globin",
    }:
        return {
            "min": None,
            "max": upper_iqr,
        }

    if metric_name == "log10GenesPerUMI":
        return {
            "min": lower_iqr,
            "max": None,
        }

    return {
        "min": lower_iqr,
        "max": upper_iqr,
    }


def combine_hybrid_bounds(
    default_bounds: dict[str, float | None],
    data_bounds: dict[str, float | None],
) -> dict[str, float | None]:
    default_min = default_bounds.get("min")
    default_max = default_bounds.get("max")
    data_min = data_bounds.get("min")
    data_max = data_bounds.get("max")

    resolved_min: float | None
    resolved_max: float | None

    if default_min is None:
        resolved_min = data_min
    elif data_min is None:
        resolved_min = default_min
    else:
        resolved_min = max(
            default_min,
            data_min,
        )

    if default_max is None:
        resolved_max = data_max
    elif data_max is None:
        resolved_max = default_max
    else:
        resolved_max = min(
            default_max,
            data_max,
        )

    return {
        "min": resolved_min,
        "max": resolved_max,
    }


def normalise_bounds(
    metric_name: str,
    bounds: dict[str, float | None],
) -> dict[str, float | None]:
    lower = bounds.get("min")
    upper = bounds.get("max")

    if metric_name in {
        "nFeature_RNA",
        "nCount_RNA",
        "percent.mito",
        "percent.ribo",
        "percent.globin",
    }:
        if lower is not None:
            lower = max(0.0, lower)

        if upper is not None:
            upper = max(0.0, upper)

    if (
        lower is not None
        and upper is not None
        and lower >= upper
    ):
        raise ValueError(
            f"Invalid resolved bounds for {metric_name}: "
            f"min={lower}, max={upper}"
        )

    return {
        "min": lower,
        "max": upper,
    }


def count_flagged_cells(
    values: list[float],
    lower: float | None,
    upper: float | None,
) -> int:
    flagged = 0

    for value in values:
        if lower is not None and value < lower:
            flagged += 1
            continue

        if upper is not None and value > upper:
            flagged += 1

    return flagged


def main() -> int:
    args = parse_args()

    try:
        plan = read_yaml(
            args.plan,
            "Validated analysis plan",
        )
        global_metrics = read_csv_rows(
            args.global_metrics,
            "Global QC metrics",
        )
        sample_metrics = read_csv_rows(
            args.sample_metrics,
            "Per-sample QC metrics",
        )
        cell_rows = read_csv_rows(
            args.cell_metadata,
            "Cell-level QC metadata",
        )

        threshold_mode = resolve_threshold_mode(plan)
        material = resolve_material(plan)
        explicit_thresholds = resolve_explicit_thresholds(plan)

    except (
        FileNotFoundError,
        ValueError,
        yaml.YAMLError,
    ) as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )
        return 2

    plan_section = plan.get("plan", {})
    experiment_section = plan.get("experiment", {})

    if not isinstance(plan_section, dict):
        plan_section = {}

    if not isinstance(experiment_section, dict):
        experiment_section = {}

    plan_id = str(
        plan_section.get(
            "plan_id",
            "",
        )
    )

    project_id = str(
        experiment_section.get(
            "project_id",
            "",
        )
    )

    biological_run_id = str(
        experiment_section.get(
            "biological_run_id",
            "",
        )
    )

    default_policy = DEFAULT_BOUNDS[material]
    resolved_thresholds: dict[str, dict[str, Any]] = {}
    total_cells = len(cell_rows)

    for metric_name, column_name in METRIC_COLUMNS.items():
        values = collect_column(
            cell_rows,
            column_name,
        )

        if not values:
            continue

        data_bounds = propose_data_driven_bounds(
            metric_name,
            values,
        )
        default_bounds = default_policy.get(
            metric_name,
            {
                "min": None,
                "max": None,
            },
        )

        if metric_name in explicit_thresholds:
            final_bounds = explicit_thresholds[metric_name]
            source = "analysis_plan_explicit"
            rationale = (
                "Explicit bounds supplied in qc.thresholds."
            )

        elif metric_name not in AUTO_FILTER_METRICS:
            final_bounds = {
                "min": None,
                "max": None,
            }
            source = "observational_only"
            rationale = (
                "Metric retained for diagnostics but not used for "
                "automatic filtering unless explicitly configured."
            )

        elif threshold_mode == "fixed":
            final_bounds = default_bounds
            source = "material_default"
            rationale = (
                f"Fixed default policy for material={material}."
            )

        elif threshold_mode == "data_driven":
            final_bounds = data_bounds
            source = "observed_distribution"
            rationale = (
                "Bounds proposed from the observed distribution "
                "using 3×IQR."
            )

        else:
            final_bounds = combine_hybrid_bounds(
                default_bounds,
                data_bounds,
            )
            source = "hybrid_default_and_distribution"
            rationale = (
                "Conservative combination of material defaults "
                "and observed 3×IQR bounds."
            )

        final_bounds = normalise_bounds(
            metric_name,
            final_bounds,
        )

        lower = final_bounds.get("min")
        upper = final_bounds.get("max")

        flagged_cells = count_flagged_cells(
            values,
            lower,
            upper,
        )

        flagged_fraction = (
            flagged_cells / total_cells
            if total_cells > 0
            else 0.0
        )

        resolved_thresholds[metric_name] = {
            "enabled": (
                lower is not None
                or upper is not None
            ),
            "min": lower,
            "max": upper,
            "source": source,
            "rationale": rationale,
            "observed": {
                "minimum": min(values),
                "q25": percentile(values, 0.25),
                "median": percentile(values, 0.50),
                "q75": percentile(values, 0.75),
                "maximum": max(values),
            },
            "estimated_flagged_cells": flagged_cells,
            "estimated_flagged_fraction": flagged_fraction,
        }

    policy_document = {
        "qc_threshold_policy": {
            "schema_version": "1.0",
            "plan_id": plan_id,
            "project_id": project_id,
            "biological_run_id": biological_run_id,
            "threshold_mode": threshold_mode,
            "material": material,
            "approval": {
                "required": True,
                "status": "proposed",
                "approved_by": None,
                "approved_at": None,
            },
            "application_semantics": {
                "minimum_operator": ">=",
                "maximum_operator": "<=",
                "combine_metrics_with": "AND",
                "missing_metric_action": "error",
            },
            "thresholds": resolved_thresholds,
            "provenance": {
                "global_metrics_file": args.global_metrics.name,
                "sample_metrics_file": args.sample_metrics.name,
                "cell_metadata_file": args.cell_metadata.name,
                "global_metric_rows": len(global_metrics),
                "sample_metric_rows": len(sample_metrics),
                "cell_rows": total_cells,
            },
        }
    }

    args.policy_yaml.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    args.policy_tsv.parent.mkdir(
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

    with args.policy_yaml.open(
        "w",
        encoding="utf-8",
    ) as handle:
        yaml.safe_dump(
            policy_document,
            handle,
            sort_keys=False,
            allow_unicode=True,
        )

    with args.policy_tsv.open(
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
                "metric",
                "enabled",
                "min",
                "max",
                "source",
                "estimated_flagged_cells",
                "estimated_flagged_fraction",
                "approval_status",
            ]
        )

        for metric_name, threshold in resolved_thresholds.items():
            writer.writerow(
                [
                    metric_name,
                    str(
                        threshold["enabled"]
                    ).lower(),
                    threshold["min"],
                    threshold["max"],
                    threshold["source"],
                    threshold[
                        "estimated_flagged_cells"
                    ],
                    threshold[
                        "estimated_flagged_fraction"
                    ],
                    "proposed",
                ]
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
                "threshold_mode",
                "material",
                "approval_status",
                "metric_count",
                "cell_count",
            ]
        )

        writer.writerow(
            [
                plan_id,
                project_id,
                biological_run_id,
                threshold_mode,
                material,
                "proposed",
                len(resolved_thresholds),
                total_cells,
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
            f"threshold_mode={threshold_mode}\n"
        )
        handle.write(
            f"material={material}\n"
        )
        handle.write(
            "approval_status=proposed\n"
        )

    print(
        "QC threshold policy resolved: proposed"
    )
    print(
        f"Resolved metrics: {len(resolved_thresholds)}"
    )
    print(
        f"Observed cells: {total_cells}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

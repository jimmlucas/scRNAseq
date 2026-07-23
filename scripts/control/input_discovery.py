#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Inspect single-cell project inputs without loading "
            "expression matrices."
        )
    )

    parser.add_argument(
        "--plan",
        required=True,
        type=Path,
        help="Validated analysis plan YAML.",
    )
    parser.add_argument(
        "--routing",
        required=True,
        type=Path,
        help="Resolved analysis routing YAML.",
    )
    parser.add_argument(
        "--project-dir",
        required=True,
        type=Path,
        help="Pipeline project directory.",
    )
    parser.add_argument(
        "--inventory",
        required=True,
        type=Path,
        help="Output input inventory YAML.",
    )
    parser.add_argument(
        "--sample-summary",
        required=True,
        type=Path,
        help="Output sample-level TSV.",
    )
    parser.add_argument(
        "--summary",
        required=True,
        type=Path,
        help="Output discovery summary TSV.",
    )

    return parser.parse_args()


def read_yaml(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")

    with path.open("r", encoding="utf-8") as handle:
        document = yaml.safe_load(handle)

    if not isinstance(document, dict):
        raise ValueError(f"{label} root must be a mapping.")

    return document


def resolve_project_path(
    value: str | None,
    project_dir: Path,
) -> Path | None:
    if value is None or str(value).strip() == "":
        return None

    path = Path(str(value)).expanduser()

    if not path.is_absolute():
        path = project_dir / path

    return path.resolve()


def detect_delimiter(path: Path) -> str:
    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        sample = handle.read(8192)

    try:
        dialect = csv.Sniffer().sniff(
            sample,
            delimiters=",;\t",
        )
        return dialect.delimiter
    except csv.Error:
        return ","


def read_targets(
    targets_path: Path,
) -> tuple[list[str], list[dict[str, str]]]:
    delimiter = detect_delimiter(targets_path)

    with targets_path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter=delimiter,
        )

        columns = reader.fieldnames or []
        rows = [
            {
                str(key): (
                    value.strip()
                    if isinstance(value, str)
                    else ""
                )
                for key, value in row.items()
                if key is not None
            }
            for row in reader
        ]

    return columns, rows


def non_empty_unique(
    rows: list[dict[str, str]],
    column: str | None,
) -> list[str]:
    if not column:
        return []

    return sorted(
        {
            row.get(column, "").strip()
            for row in rows
            if row.get(column, "").strip()
        }
    )


def duplicate_values(
    rows: list[dict[str, str]],
    column: str | None,
) -> list[str]:
    if not column:
        return []

    counts = Counter(
        row.get(column, "").strip()
        for row in rows
        if row.get(column, "").strip()
    )

    return sorted(
        value
        for value, count in counts.items()
        if count > 1
    )


def count_directory_entries(path: Path | None) -> dict[str, Any]:
    result: dict[str, Any] = {
        "exists": False,
        "is_directory": False,
        "entry_count": 0,
        "file_count": 0,
        "directory_count": 0,
    }

    if path is None:
        return result

    result["exists"] = path.exists()
    result["is_directory"] = path.is_dir()

    if not path.is_dir():
        return result

    entries = list(path.iterdir())

    result["entry_count"] = len(entries)
    result["file_count"] = sum(
        entry.is_file()
        for entry in entries
    )
    result["directory_count"] = sum(
        entry.is_dir()
        for entry in entries
    )

    return result


def normalize_modalities(value: Any) -> list[str]:
    if value is None:
        return []

    if isinstance(value, str):
        return [value]

    if isinstance(value, list):
        return [
            str(item)
            for item in value
            if str(item).strip()
        ]

    return []


def main() -> int:
    args = parse_args()
    project_dir = args.project_dir.resolve()

    try:
        plan_document = read_yaml(
            args.plan,
            "Validated analysis plan",
        )
        routing_document = read_yaml(
            args.routing,
            "Resolved analysis routing",
        )
    except (FileNotFoundError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

    validation = plan_document.get("validation", {})
    if not isinstance(validation, dict):
        validation = {}

    if validation.get("status") != "valid":
        print(
            "ERROR: Input discovery requires a valid analysis plan.",
            file=sys.stderr,
        )
        return 1

    routes = routing_document.get("routes", {})
    if not isinstance(routes, dict):
        routes = {}

    discovery_route = routes.get(
        "input_discovery",
        {},
    )
    if not isinstance(discovery_route, dict):
        discovery_route = {}

    requested = bool(
        discovery_route.get("requested", False)
    )

    plan = plan_document.get("plan", {})
    experiment = plan_document.get("experiment", {})
    inputs = plan_document.get("inputs", {})

    if not isinstance(plan, dict):
        plan = {}

    if not isinstance(experiment, dict):
        experiment = {}

    if not isinstance(inputs, dict):
        inputs = {}

    plan_id = str(plan.get("plan_id", ""))
    project_id = str(experiment.get("project_id", ""))
    biological_run_id = str(
        experiment.get("biological_run_id", "")
    )

    inventory: dict[str, Any] = {
        "input_discovery": {
            "schema_version": "1.0",
            "plan_id": plan_id,
            "project_id": project_id,
            "biological_run_id": biological_run_id,
            "requested": requested,
        }
    }

    if not requested:
        inventory["input_discovery"].update(
            {
                "status": "skipped",
                "reason": "not_requested",
            }
        )

        args.inventory.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        with args.inventory.open(
            "w",
            encoding="utf-8",
        ) as handle:
            yaml.safe_dump(
                inventory,
                handle,
                sort_keys=False,
                allow_unicode=True,
            )

        with args.sample_summary.open(
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
                    "row_number",
                    "sample_id",
                    "biological_sample_id",
                    "group",
                    "run_id",
                    "lane_id",
                    "technical_replicate_id",
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
                    "status",
                    "sample_count",
                    "biological_sample_count",
                    "missing_column_count",
                    "warning_count",
                    "error_count",
                ]
            )
            writer.writerow(
                [
                    plan_id,
                    "skipped",
                    0,
                    0,
                    0,
                    0,
                    0,
                ]
            )

        print("Input discovery skipped: not requested.")
        return 0

    load_config_value = inputs.get("load_config")
    load_config_path = resolve_project_path(
        str(load_config_value)
        if load_config_value is not None
        else None,
        project_dir,
    )

    errors: list[str] = []
    warnings: list[str] = []

    runtime_config: dict[str, Any] = {}

    if load_config_path is None:
        errors.append("inputs.load_config is not defined.")
    elif not load_config_path.exists():
        errors.append(
            f"Load configuration not found: {load_config_path}"
        )
    else:
        try:
            runtime_config = read_yaml(
                load_config_path,
                "Load configuration",
            )
        except (FileNotFoundError, ValueError) as error:
            errors.append(str(error))

    config_project = runtime_config.get("project", {})
    config_input = runtime_config.get("input", {})

    if not isinstance(config_project, dict):
        config_project = {}

    if not isinstance(config_input, dict):
        config_input = {}

    targets_path = resolve_project_path(
        config_input.get("targets"),
        project_dir,
    )
    count_dir_path = resolve_project_path(
        config_input.get("count_dir"),
        project_dir,
    )

    sample_column = (
        inputs.get("sample_column")
        or config_project.get("sample_column")
    )
    biological_sample_column = (
        inputs.get("biological_sample_column")
        or sample_column
    )
    group_column = (
        inputs.get("group_column")
        or config_project.get("group_column")
    )
    run_column = inputs.get("run_column")
    lane_column = inputs.get("lane_column")
    replicate_column = inputs.get("replicate_column")

    declared_columns = {
        "sample_column": sample_column,
        "biological_sample_column": biological_sample_column,
        "group_column": group_column,
        "run_column": run_column,
        "lane_column": lane_column,
        "replicate_column": replicate_column,
    }

    targets_columns: list[str] = []
    targets_rows: list[dict[str, str]] = []

    if targets_path is None:
        errors.append(
            "input.targets is not defined in the load configuration."
        )
    elif not targets_path.exists():
        errors.append(
            f"Targets file not found: {targets_path}"
        )
    elif not targets_path.is_file():
        errors.append(
            f"Targets path is not a file: {targets_path}"
        )
    else:
        try:
            targets_columns, targets_rows = read_targets(
                targets_path
            )
        except (OSError, csv.Error) as error:
            errors.append(
                f"Unable to read targets file: {error}"
            )

    required_columns = {
        key: value
        for key, value in declared_columns.items()
        if value is not None
    }

    missing_columns = sorted(
        {
            str(column)
            for column in required_columns.values()
            if str(column) not in targets_columns
        }
    )

    if missing_columns:
        errors.append(
            "Missing declared target columns: "
            + ", ".join(missing_columns)
        )

    samples = non_empty_unique(
        targets_rows,
        str(sample_column)
        if sample_column
        else None,
    )
    biological_samples = non_empty_unique(
        targets_rows,
        str(biological_sample_column)
        if biological_sample_column
        else None,
    )
    groups = non_empty_unique(
        targets_rows,
        str(group_column)
        if group_column
        else None,
    )
    runs = non_empty_unique(
        targets_rows,
        str(run_column)
        if run_column
        else None,
    )
    lanes = non_empty_unique(
        targets_rows,
        str(lane_column)
        if lane_column
        else None,
    )
    replicates = non_empty_unique(
        targets_rows,
        str(replicate_column)
        if replicate_column
        else None,
    )

    duplicated_sample_ids = duplicate_values(
        targets_rows,
        str(sample_column)
        if sample_column
        else None,
    )

    if duplicated_sample_ids:
        warnings.append(
            "Duplicated sample identifiers detected: "
            + ", ".join(duplicated_sample_ids)
        )

    modalities = normalize_modalities(
        experiment.get("modalities")
    )

    multimodal = len(modalities) > 1
    technical_replicates_detected = bool(
        replicates
        or (
            biological_samples
            and len(samples) > len(biological_samples)
        )
        or duplicated_sample_ids
    )

    merge_required = len(samples) > 1
    sample_resolution_required = (
        technical_replicates_detected
    )

    count_dir_inventory = count_directory_entries(
        count_dir_path
    )

    if count_dir_path is None:
        errors.append(
            "input.count_dir is not defined in the load configuration."
        )
    elif not count_dir_path.exists():
        errors.append(
            f"Count directory not found: {count_dir_path}"
        )
    elif not count_dir_path.is_dir():
        errors.append(
            f"Count directory path is not a directory: "
            f"{count_dir_path}"
        )

    status = "ready" if not errors else "blocked"

    inventory["input_discovery"].update(
        {
            "status": status,
            "experiment": {
                "experiment_type": experiment.get(
                    "experiment_type"
                ),
                "organism": experiment.get("organism"),
                "material": experiment.get("material"),
                "tissue": experiment.get("tissue"),
                "platform": experiment.get("platform"),
                "matrix_type": inputs.get("matrix_type"),
            },
            "modalities": {
                "declared": modalities,
                "multimodal": multimodal,
            },
            "configuration": {
                "load_config": {
                    "path": (
                        str(load_config_path)
                        if load_config_path
                        else None
                    ),
                    "exists": bool(
                        load_config_path
                        and load_config_path.exists()
                    ),
                },
            },
            "inputs": {
                "targets_file": {
                    "path": (
                        str(targets_path)
                        if targets_path
                        else None
                    ),
                    "exists": bool(
                        targets_path
                        and targets_path.exists()
                    ),
                    "row_count": len(targets_rows),
                    "columns": targets_columns,
                },
                "count_dir": {
                    "path": (
                        str(count_dir_path)
                        if count_dir_path
                        else None
                    ),
                    **count_dir_inventory,
                },
            },
            "metadata": {
                "declared_columns": declared_columns,
                "missing_columns": missing_columns,
            },
            "samples": {
                "row_count": len(targets_rows),
                "sample_count": len(samples),
                "biological_sample_count": len(
                    biological_samples
                ),
                "group_count": len(groups),
                "run_count": len(runs),
                "lane_count": len(lanes),
                "technical_replicate_count": len(replicates),
                "sample_ids": samples,
                "biological_sample_ids": biological_samples,
                "groups": groups,
                "runs": runs,
                "lanes": lanes,
                "technical_replicates": replicates,
                "duplicated_sample_ids": duplicated_sample_ids,
            },
            "routing_observations": {
                "sample_resolution_required": (
                    sample_resolution_required
                ),
                "merge_required": merge_required,
                "multimodal_qc_required": multimodal,
            },
            "warnings": warnings,
            "errors": errors,
        }
    )

    args.inventory.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    args.sample_summary.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    args.summary.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with args.inventory.open(
        "w",
        encoding="utf-8",
    ) as handle:
        yaml.safe_dump(
            inventory,
            handle,
            sort_keys=False,
            allow_unicode=True,
        )

    with args.sample_summary.open(
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
                "row_number",
                "sample_id",
                "biological_sample_id",
                "group",
                "run_id",
                "lane_id",
                "technical_replicate_id",
            ]
        )

        for index, row in enumerate(
            targets_rows,
            start=1,
        ):
            writer.writerow(
                [
                    index,
                    row.get(str(sample_column), "")
                    if sample_column
                    else "",
                    row.get(
                        str(biological_sample_column),
                        "",
                    )
                    if biological_sample_column
                    else "",
                    row.get(str(group_column), "")
                    if group_column
                    else "",
                    row.get(str(run_column), "")
                    if run_column
                    else "",
                    row.get(str(lane_column), "")
                    if lane_column
                    else "",
                    row.get(str(replicate_column), "")
                    if replicate_column
                    else "",
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
                "status",
                "sample_count",
                "biological_sample_count",
                "missing_column_count",
                "warning_count",
                "error_count",
            ]
        )

        writer.writerow(
            [
                plan_id,
                status,
                len(samples),
                len(biological_samples),
                len(missing_columns),
                len(warnings),
                len(errors),
            ]
        )

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)

        return 1

    for warning in warnings:
        print(f"WARNING: {warning}", file=sys.stderr)

    print(
        "Input discovery completed: "
        f"{len(samples)} samples, "
        f"{len(biological_samples)} biological samples."
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

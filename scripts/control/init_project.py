#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
import csv
import sys
import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Initialize a scRNAseq project and generate runtime contracts."
    )

    parser.add_argument("--project-id", required=True)
    parser.add_argument("--biological-run-id", required=True)
    parser.add_argument("--experiment-type", required=True)
    parser.add_argument("--organism", required=True)
    parser.add_argument("--material", required=True)

    parser.add_argument(
        "--targets",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--count-dir",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--sample-column",
        default="SampleName",
    )

    parser.add_argument(
        "--shortname-column",
        default="ShortName",
    )

    parser.add_argument(
        "--group-column",
        default="Group",
    )

    parser.add_argument(
        "--threshold-mode",
        default="data_driven",
        choices=[
            "fixed",
            "data_driven",
            "adaptive",
            "hybrid",
        ],
    )

    return parser.parse_args()


def validate_targets(
    targets_path: Path,
    required_columns: list[str],
) -> None:

    if not targets_path.is_file():
        raise FileNotFoundError(
            f"Targets file not found: {targets_path}"
        )

    with targets_path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        reader = csv.DictReader(handle)

        columns = reader.fieldnames or []

        missing = [
            column
            for column in required_columns
            if column not in columns
        ]

        if missing:
            raise ValueError(
                "Missing required columns in targets: "
                + ", ".join(missing)
            )

        rows = list(reader)

    if not rows:
        raise ValueError("Targets file contains no samples.")

    sample_ids = [
        row[required_columns[1]].strip()
        for row in rows
    ]

    if any(not sample_id for sample_id in sample_ids):
        raise ValueError(
            f"Empty values detected in {required_columns[1]}."
        )

    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError(
            f"Duplicate values detected in {required_columns[1]}."
        )


def validate_sample_directories(
    count_dir: Path,
    targets_path: Path,
    shortname_column: str,
) -> None:

    if not count_dir.is_dir():
        raise FileNotFoundError(
            f"Count directory not found: {count_dir}"
        )

    with targets_path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        reader = csv.DictReader(handle)

        for row in reader:

            sample_id = row[shortname_column].strip()

            sample_dir = count_dir / sample_id

            if not sample_dir.is_dir():
                raise FileNotFoundError(
                    f"Sample directory not found: {sample_dir}"
                )

            rds_files = list(
                sample_dir.rglob("*.rds")
            )

            if len(rds_files) == 0:
                raise FileNotFoundError(
                    f"No RDS found for sample {sample_id}: "
                    f"{sample_dir}"
                )

            if len(rds_files) > 1:
                raise ValueError(
                    f"Multiple RDS files found for sample "
                    f"{sample_id}: {sample_dir}"
                )


def relative_or_absolute(
    path: Path,
    project_dir: Path,
) -> str:

    path = path.resolve()

    try:
        return str(
            path.relative_to(project_dir.resolve())
        )
    except ValueError:
        return str(path)


def main() -> int:

    args = parse_args()

    project_dir = Path.cwd().resolve()

    targets = args.targets.resolve()
    count_dir = args.count_dir.resolve()
    output_dir = args.output_dir.resolve()

    required_columns = [
        args.sample_column,
        args.shortname_column,
        args.group_column,
    ]

    try:

        validate_targets(
            targets,
            required_columns,
        )

        validate_sample_directories(
            count_dir,
            targets,
            args.shortname_column,
        )

    except (FileNotFoundError, ValueError) as error:

        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )

        return 1

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    targets_value = relative_or_absolute(
        targets,
        project_dir,
    )

    count_dir_value = relative_or_absolute(
        count_dir,
        project_dir,
    )

    load_config_path = (
        output_dir
        / "load_config.yaml"
    )

    analysis_plan_path = (
        output_dir
        / "analysis_plan.yaml"
    )

    load_config_value = relative_or_absolute(
        load_config_path,
        project_dir,
    )

    load_config = {
        "project": {
            "name": args.project_id,
            "organism": args.organism,
            "sample_column": args.sample_column,
            "shortname_column": args.shortname_column,
            "group_column": args.group_column,
        },
        "input": {
            "targets": targets_value,
            "count_dir": count_dir_value,
        },
        "output": {
            "data_dir": count_dir_value,
            "results_dir": (
                f"results/{args.project_id}/00_load_data"
            ),
            "merged_object": "seu.rds",
        },
        "qc": {
            "min_features": 0,
            "max_features": 100000,
            "min_counts": 0,
            "max_counts": 1000000,
            "max_percent_mito": 100,
            "max_percent_ribo": 100,
            "max_percent_globin": 100,
        },
        "doublet": {
            "expected_rate": 0.075,
            "pN": 0.25,
            "pcs": 15,
            "min_cells": 500,
            "nfeatures": 1500,
        },
        "scoring": {
            "scale_factor": 10000,
        },
    }

    analysis_plan = {
        "plan": {
            "schema_version": "1.0",
            "plan_id": (
                f"{args.project_id}_"
                f"{args.biological_run_id}_analysis"
            ),
        },
        "experiment": {
            "project_id": args.project_id,
            "biological_run_id": args.biological_run_id,
            "experiment_type": args.experiment_type,
            "material": args.material,
            "organism": args.organism,
            "modalities": [
                "RNA",
            ],
        },
        "inputs": {
            "matrix_type": "seurat_object",
            "load_config": load_config_value,
            "sample_column": args.sample_column,
            "biological_sample_column": args.sample_column,
            "group_column": args.group_column,
        },
        "routing": {
            "input_discovery": True,
            "load_data": True,
            "qc_metrics": True,
            "qc_visualization": False,
            "qc_threshold_resolution": True,
            "qc_filter": False,
            "doublet_detection": False,
            "sample_resolution": False,
            "merge_samples": False,
            "qc_scores": False,
            "normalization": False,
        },
        "qc": {
            "threshold_mode": args.threshold_mode,
        },
    }

    with load_config_path.open(
        "w",
        encoding="utf-8",
    ) as handle:

        yaml.safe_dump(
            load_config,
            handle,
            sort_keys=False,
            allow_unicode=True,
        )

    with analysis_plan_path.open(
        "w",
        encoding="utf-8",
    ) as handle:

        yaml.safe_dump(
            analysis_plan,
            handle,
            sort_keys=False,
            allow_unicode=True,
        )

    print("Project initialized successfully.")
    print(f"Project ID: {args.project_id}")
    print(
        "Biological run ID: "
        f"{args.biological_run_id}"
    )
    print(f"Targets: {targets_value}")
    print(f"Count directory: {count_dir_value}")
    print(f"Load config: {load_config_path}")
    print(f"Analysis plan: {analysis_plan_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

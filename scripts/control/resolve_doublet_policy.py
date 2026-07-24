#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import platform
import sys
from pathlib import Path
from typing import Any


COMPLETED_STATUS = "completed"
PROPOSED_ACTION = "exclude_algorithmic_doublets"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Resolve a proposed doublet policy from technical "
            "doublet scoring artifacts."
        )
    )

    parser.add_argument(
        "--summary",
        required=True,
        type=Path,
        help="Doublet scoring summary CSV.",
    )
    parser.add_argument(
        "--by-sample",
        required=True,
        type=Path,
        help="Doublet scoring per-sample CSV.",
    )
    parser.add_argument(
        "--cell-scores",
        required=True,
        type=Path,
        help="Doublet scoring cell-level CSV.",
    )
    parser.add_argument(
        "--policy-yaml",
        required=True,
        type=Path,
        help="Output proposed doublet policy YAML.",
    )
    parser.add_argument(
        "--policy-tsv",
        required=True,
        type=Path,
        help="Output proposed doublet policy TSV.",
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


def read_summary_values(path: Path) -> dict[str, str]:
    rows = read_csv_rows(path, "Doublet scoring summary")
    values: dict[str, str] = {}

    for row in rows:
        metric = str(row.get("metric", "")).strip()

        if metric:
            values[metric] = str(row.get("value", "")).strip()

    return values


def require_columns(
    rows: list[dict[str, str]],
    columns: list[str],
    label: str,
) -> None:
    if not rows:
        return

    missing = [
        column
        for column in columns
        if column not in rows[0]
    ]

    if missing:
        raise ValueError(
            f"{label} is missing required columns: "
            + ", ".join(missing)
        )


def parse_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def parse_float(value: Any, default: float | None = None) -> float | None:
    try:
        parsed = float(str(value))
    except (TypeError, ValueError):
        return default

    return parsed


def count_cell_classes(
    cell_rows: list[dict[str, str]],
) -> dict[str, dict[str, int]]:
    counts: dict[str, dict[str, int]] = {}

    for row in cell_rows:
        sample = str(row.get("sample", "")).strip()
        doublet_class = str(row.get("doublet_class", "")).strip()

        if not sample:
            sample = "unknown"

        if sample not in counts:
            counts[sample] = {
                "Doublet": 0,
                "Singlet": 0,
                "NotRun": 0,
                "Other": 0,
            }

        if doublet_class in counts[sample]:
            counts[sample][doublet_class] += 1
        else:
            counts[sample]["Other"] += 1

    return counts


def proposed_action_for_status(
    detection_status: str,
    proposed_doublets: int,
) -> str:
    if detection_status != COMPLETED_STATUS:
        return "no_filter_not_run"

    if proposed_doublets <= 0:
        return "no_filter_no_doublets_detected"

    return PROPOSED_ACTION


def build_sample_policy(
    by_sample_rows: list[dict[str, str]],
    class_counts: dict[str, dict[str, int]],
) -> dict[str, dict[str, Any]]:
    samples: dict[str, dict[str, Any]] = {}

    for row in by_sample_rows:
        sample = str(row.get("sample", "")).strip()

        if not sample:
            raise ValueError(
                "Doublet by-sample table contains a row without sample."
            )

        detection_status = str(
            row.get("doublet_detection_status", "")
        ).strip()

        decision_status = str(
            row.get("doublet_decision_status", "")
        ).strip()

        approval_status = str(
            row.get("approval_status", "")
        ).strip()

        proposed_doublets = parse_int(
            row.get("proposed_doublets")
        )

        cells_scored = parse_int(
            row.get("cells_scored")
        )

        if decision_status != "unreviewed":
            raise ValueError(
                f"Sample {sample} is not unreviewed: "
                f"{decision_status}"
            )

        if approval_status != "proposed":
            raise ValueError(
                f"Sample {sample} approval status is not proposed: "
                f"{approval_status}"
            )

        samples[sample] = {
            "detection_status": detection_status,
            "decision_status": decision_status,
            "proposed_action": proposed_action_for_status(
                detection_status,
                proposed_doublets,
            ),
            "cells_scored": cells_scored,
            "proposed_doublet_count": proposed_doublets,
            "proposed_singlet_count": parse_int(
                row.get("proposed_singlets")
            ),
            "expected_doublet_rate": parse_float(
                row.get("expected_doublet_rate")
            ),
            "pN": parse_float(row.get("pN")),
            "pK": parse_float(row.get("pK")),
            "nExp_raw": parse_int(row.get("nExp_raw")),
            "homotypic_prop": parse_float(
                row.get("homotypic_prop")
            ),
            "nExp_adjusted": parse_int(
                row.get("nExp_adjusted")
            ),
            "cell_class_counts": class_counts.get(
                sample,
                {},
            ),
        }

    return samples


def write_policy_tsv(
    path: Path,
    samples: dict[str, dict[str, Any]],
) -> None:
    with path.open(
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
                "sample",
                "detection_status",
                "decision_status",
                "proposed_action",
                "cells_scored",
                "proposed_doublet_count",
                "proposed_singlet_count",
                "pK",
                "nExp_adjusted",
                "approval_status",
            ]
        )

        for sample, policy in samples.items():
            writer.writerow(
                [
                    sample,
                    policy["detection_status"],
                    policy["decision_status"],
                    policy["proposed_action"],
                    policy["cells_scored"],
                    policy["proposed_doublet_count"],
                    policy["proposed_singlet_count"],
                    policy["pK"],
                    policy["nExp_adjusted"],
                    "proposed",
                ]
            )


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


def write_policy_yaml(
    path: Path,
    policy_document: dict[str, Any],
) -> None:
    with path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        for key, value in policy_document.items():
            write_yaml_value(
                handle,
                key,
                value,
                indent=0,
            )


def main() -> int:
    args = parse_args()

    try:
        summary_values = read_summary_values(args.summary)
        by_sample_rows = read_csv_rows(
            args.by_sample,
            "Doublet by-sample scores",
        )
        cell_rows = read_csv_rows(
            args.cell_scores,
            "Doublet cell scores",
        )

        require_columns(
            by_sample_rows,
            [
                "sample",
                "cells_scored",
                "pK",
                "nExp_adjusted",
                "proposed_singlets",
                "proposed_doublets",
                "doublet_detection_status",
                "doublet_decision_status",
                "approval_status",
            ],
            "Doublet by-sample scores",
        )

        require_columns(
            cell_rows,
            [
                "cell_id",
                "sample",
                "doublet_class",
                "is_doublet",
                "doublet_detection_status",
                "doublet_decision_status",
            ],
            "Doublet cell scores",
        )

        class_counts = count_cell_classes(cell_rows)
        sample_policy = build_sample_policy(
            by_sample_rows,
            class_counts,
        )

        total_cells = parse_int(
            summary_values.get("cells_output")
        )
        proposed_doublets = parse_int(
            summary_values.get("proposed_doublets")
        )
        cells_removed = parse_int(
            summary_values.get("cells_removed")
        )

        if cells_removed != 0:
            raise ValueError(
                "Doublet scoring artifacts indicate cells were removed. "
                f"cells_removed={cells_removed}"
            )

    except (
        FileNotFoundError,
        ValueError,
    ) as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )
        return 2

    policy_document = {
        "doublet_policy": {
            "schema_version": "1.0",
            "identity": {
                "project_id": args.project_id,
                "project_name": args.project_name,
                "biological_run_id": args.biological_run_id,
            },
            "method": {
                "name": "DoubletFinder",
                "source_artifact": "seu_doublet_scored.rds",
                "source_summary": args.summary.name,
                "source_by_sample": args.by_sample.name,
                "source_cell_scores": args.cell_scores.name,
            },
            "application_semantics": {
                "scope": "per_sample",
                "proposed_doublet_class": "Doublet",
                "not_run_action": "do_not_filter",
                "missing_cell_action": "error",
                "approval_required_before_filtering": True,
            },
            "summary": {
                "cells_scored": total_cells,
                "samples": len(sample_policy),
                "proposed_doublets": proposed_doublets,
                "cells_removed_by_scoring": cells_removed,
            },
            "samples": sample_policy,
            "approval": {
                "required": True,
                "status": "proposed",
                "reviewed_by": None,
                "approved_at": None,
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
    args.resolution_summary.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    args.session_info.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    write_policy_yaml(
        args.policy_yaml,
        policy_document,
    )

    write_policy_tsv(
        args.policy_tsv,
        sample_policy,
    )

    with args.resolution_summary.open(
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
                "project_id",
                "biological_run_id",
                "samples",
                "cells_scored",
                "proposed_doublets",
                "approval_status",
            ]
        )
        writer.writerow(
            [
                args.project_id,
                args.biological_run_id,
                len(sample_policy),
                total_cells,
                proposed_doublets,
                "proposed",
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
            "approval_status=proposed\n"
        )
        handle.write(
            "resolver=DOUBLET_POLICY_RESOLUTION\n"
        )

    print("Doublet policy resolved: proposed")
    print(f"Samples: {len(sample_policy)}")
    print(f"Cells scored: {total_cells}")
    print(f"Proposed doublets: {proposed_doublets}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

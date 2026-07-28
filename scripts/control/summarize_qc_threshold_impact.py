#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import itertools
import platform
from pathlib import Path
from typing import Any

import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize the simulated impact of a proposed QC threshold "
            "policy without filtering cells."
        )
    )
    parser.add_argument(
        "--policy-yaml",
        required=True,
        type=Path,
        help="Proposed QC threshold policy YAML.",
    )
    parser.add_argument(
        "--cell-metadata",
        required=True,
        type=Path,
        help="Cell-level QC metadata CSV.",
    )
    parser.add_argument(
        "--by-sample",
        required=True,
        type=Path,
        help="Per-sample QC metrics CSV.",
    )
    parser.add_argument(
        "--impact-yaml",
        required=True,
        type=Path,
        help="Output impact summary YAML.",
    )
    parser.add_argument(
        "--by-sample-output",
        required=True,
        type=Path,
        help="Output impact-by-sample CSV.",
    )
    parser.add_argument(
        "--by-metric-output",
        required=True,
        type=Path,
        help="Output impact-by-metric CSV.",
    )
    parser.add_argument(
        "--overlap-output",
        required=True,
        type=Path,
        help="Output metric overlap CSV.",
    )
    parser.add_argument(
        "--lower-tail-by-sample-output",
        required=True,
        type=Path,
        help="Output lower-tail diagnostic by sample CSV.",
    )
    parser.add_argument(
        "--lower-tail-summary-output",
        required=True,
        type=Path,
        help="Output global lower-tail diagnostic CSV.",
    )
    parser.add_argument(
        "--summary",
        required=True,
        type=Path,
        help="Output one-row summary TSV.",
    )
    parser.add_argument(
        "--session-info",
        required=True,
        type=Path,
        help="Output runtime session information.",
    )
    parser.add_argument(
        "--project-id",
        default=None,
        help="Optional expected project_id for identity validation.",
    )
    parser.add_argument(
        "--biological-run-id",
        default=None,
        help="Optional expected biological_run_id for identity validation.",
    )
    parser.add_argument(
        "--plan-id",
        default=None,
        help="Optional expected plan_id for identity validation.",
    )
    return parser.parse_args()


def read_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Policy YAML not found: {path}")

    document = yaml.safe_load(path.read_text(encoding="utf-8"))

    if not isinstance(document, dict):
        raise ValueError(f"Policy YAML root must be a mapping: {path}")

    policy = document.get("qc_threshold_policy")

    if not isinstance(policy, dict):
        raise ValueError("Missing qc_threshold_policy mapping.")

    return policy


def read_csv_rows(path: Path, label: str) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            raise ValueError(f"{label} has no header: {path}")

        return list(reader)


def parse_float(value: Any) -> float | None:
    if value is None or value == "":
        return None

    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def resolve_sample_column(rows: list[dict[str, str]]) -> str:
    if not rows:
        raise ValueError("Cell metadata contains no rows.")

    columns = rows[0].keys()

    for candidate in (
        "sample",
        "Sample",
        "SampleName",
        "orig.ident",
        "ShortName",
    ):
        if candidate in columns:
            return candidate

    raise ValueError(
        "Unable to infer sample column from cell metadata. "
        "Expected one of: sample, Sample, SampleName, orig.ident, ShortName."
    )


def resolve_group_column(rows: list[dict[str, str]]) -> str | None:
    if not rows:
        return None

    columns = rows[0].keys()

    for candidate in (
        "Group",
        "group",
        "condition",
        "Condition",
    ):
        if candidate in columns:
            return candidate

    return None


def require_identity(
    policy: dict[str, Any],
    expected_key: str,
    expected_value: str | None,
) -> None:
    if expected_value is None or expected_value == "":
        return

    observed_value = str(policy.get(expected_key, ""))

    if observed_value != expected_value:
        raise ValueError(
            f"Policy {expected_key} mismatch: "
            f"expected '{expected_value}', observed '{observed_value}'."
        )


def resolve_application_semantics(
    policy: dict[str, Any],
) -> dict[str, str]:
    semantics = policy.get("application_semantics")

    if not isinstance(semantics, dict):
        raise ValueError(
            "qc_threshold_policy.application_semantics must be a mapping."
        )

    minimum_operator = str(
        semantics.get("minimum_operator", "")
    )
    maximum_operator = str(
        semantics.get("maximum_operator", "")
    )
    combine_metrics_with = str(
        semantics.get("combine_metrics_with", "")
    ).upper()
    missing_metric_action = str(
        semantics.get("missing_metric_action", "")
    ).lower()

    if minimum_operator != ">=":
        raise ValueError(
            "Unsupported minimum operator. Expected >=."
        )

    if maximum_operator != "<=":
        raise ValueError(
            "Unsupported maximum operator. Expected <=."
        )

    if combine_metrics_with != "AND":
        raise ValueError(
            "Unsupported metric combination. Expected AND."
        )

    if missing_metric_action != "error":
        raise ValueError(
            "Unsupported missing_metric_action. Expected error."
        )

    return {
        "minimum_operator": minimum_operator,
        "maximum_operator": maximum_operator,
        "combine_metrics_with": combine_metrics_with,
        "missing_metric_action": missing_metric_action,
    }


def resolve_thresholds(
    policy: dict[str, Any],
) -> dict[str, dict[str, float | None]]:
    thresholds = policy.get("thresholds")

    if not isinstance(thresholds, dict):
        raise ValueError("qc_threshold_policy.thresholds must be a mapping.")

    resolved: dict[str, dict[str, float | None]] = {}

    for metric, rule in thresholds.items():
        if not isinstance(rule, dict):
            continue

        if not bool(rule.get("enabled")):
            continue

        resolved[metric] = {
            "min": parse_float(rule.get("min")),
            "max": parse_float(rule.get("max")),
        }

        if (
            resolved[metric]["min"] is None
            and resolved[metric]["max"] is None
        ):
            raise ValueError(
                f"Enabled threshold has no bounds: {metric}"
            )

        if (
            resolved[metric]["min"] is not None
            and resolved[metric]["max"] is not None
            and resolved[metric]["min"] >= resolved[metric]["max"]
        ):
            raise ValueError(
                f"Invalid bounds for {metric}: "
                f"min={resolved[metric]['min']}, "
                f"max={resolved[metric]['max']}"
            )

    if not resolved:
        raise ValueError("No enabled thresholds found in QC policy.")

    return resolved


def is_flagged(
    value: float | None,
    bounds: dict[str, float | None],
) -> bool:
    if value is None:
        return False

    minimum = bounds.get("min")
    maximum = bounds.get("max")

    if minimum is not None and value < minimum:
        return True

    if maximum is not None and value > maximum:
        return True

    return False


def ratio(
    numerator: int,
    denominator: int,
) -> float:
    if denominator == 0:
        return 0.0

    return numerator / denominator


def write_csv(
    path: Path,
    rows: list[dict[str, Any]],
    fieldnames: list[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_summary_tsv(
    path: Path,
    row: dict[str, Any],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "plan_id",
        "project_id",
        "biological_run_id",
        "threshold_mode",
        "approval_status",
        "cell_count",
        "flagged_any",
        "retained",
        "flagged_fraction",
        "enabled_metric_count",
    ]

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow({key: row.get(key, "") for key in fieldnames})


def lower_tail_thresholds(
    columns: set[str],
) -> dict[str, list[float]]:
    probes = {
        "nFeature_RNA": [
            100.0,
            200.0,
            500.0,
        ],
        "nCount_RNA": [
            500.0,
            1000.0,
        ],
        "log10GenesPerUMI": [
            0.80,
            0.85,
        ],
    }

    return {
        metric: thresholds
        for metric, thresholds in probes.items()
        if metric in columns
    }


def main() -> int:
    args = parse_args()

    policy = read_yaml(args.policy_yaml)
    require_identity(policy, "project_id", args.project_id)
    require_identity(
        policy,
        "biological_run_id",
        args.biological_run_id,
    )
    require_identity(policy, "plan_id", args.plan_id)
    application_semantics = resolve_application_semantics(policy)
    cell_rows = read_csv_rows(args.cell_metadata, "Cell metadata")
    sample_rows = read_csv_rows(args.by_sample, "Per-sample metrics")
    thresholds = resolve_thresholds(policy)

    sample_col = resolve_sample_column(cell_rows)
    group_col = resolve_group_column(cell_rows)

    known_groups: dict[str, str] = {}

    for row in sample_rows:
        sample = row.get("sample") or row.get("Sample") or row.get(sample_col)
        group = row.get("Group") or row.get("group") or ""

        if sample:
            known_groups[sample] = group

    by_sample: dict[str, dict[str, Any]] = {}
    by_metric: dict[tuple[str, str], dict[str, Any]] = {}
    metric_sets: dict[str, set[str]] = {
        metric: set()
        for metric in thresholds
    }
    sample_metric_sets: dict[tuple[str, str], set[str]] = {}
    flagged_any_cells: set[str] = set()

    cell_columns = set(cell_rows[0])

    for metric in thresholds:
        if metric not in cell_columns:
            raise ValueError(
                f"Enabled QC metric is missing from metadata: {metric}"
            )

    lower_tail_probes = lower_tail_thresholds(cell_columns)
    lower_tail_by_sample: dict[tuple[str, str, float], dict[str, Any]] = {}
    lower_tail_global: dict[tuple[str, float], dict[str, Any]] = {}

    for index, row in enumerate(cell_rows):
        cell_id = (
            row.get("cell")
            or row.get("barcode")
            or row.get("Cell")
            or f"cell_{index + 1}"
        )
        sample = row.get(sample_col, "")
        group = (
            row.get(group_col, "")
            if group_col is not None
            else known_groups.get(sample, "")
        )

        sample_record = by_sample.setdefault(
            sample,
            {
                "sample": sample,
                "Group": group,
                "cells_before": 0,
                "flagged_any": 0,
                "retained": 0,
            },
        )
        sample_record["cells_before"] += 1

        metric_flags: dict[str, bool] = {}

        for metric, probe_thresholds in lower_tail_probes.items():
            raw_value = row.get(metric)
            value = parse_float(raw_value)

            if value is None:
                continue

            for probe_threshold in probe_thresholds:
                global_key = (metric, probe_threshold)
                global_record = lower_tail_global.setdefault(
                    global_key,
                    {
                        "metric": metric,
                        "threshold": probe_threshold,
                        "cells_evaluated": 0,
                        "cells_below_threshold": 0,
                    },
                )
                global_record["cells_evaluated"] += 1

                sample_key = (sample, metric, probe_threshold)
                sample_tail_record = lower_tail_by_sample.setdefault(
                    sample_key,
                    {
                        "sample": sample,
                        "Group": group,
                        "metric": metric,
                        "threshold": probe_threshold,
                        "cells_evaluated": 0,
                        "cells_below_threshold": 0,
                    },
                )
                sample_tail_record["cells_evaluated"] += 1

                if value < probe_threshold:
                    global_record["cells_below_threshold"] += 1
                    sample_tail_record["cells_below_threshold"] += 1

        for metric, bounds in thresholds.items():
            raw_value = row.get(metric)
            value = parse_float(raw_value)

            if raw_value is None or raw_value == "" or value is None:
                if application_semantics["missing_metric_action"] == "error":
                    raise ValueError(
                        "Enabled QC metric contains missing or non-numeric "
                        f"values: {metric}"
                    )

            flagged = is_flagged(value, bounds)
            metric_flags[metric] = flagged

            metric_key = (sample, metric)
            metric_record = by_metric.setdefault(
                metric_key,
                {
                    "sample": sample,
                    "Group": group,
                    "metric": metric,
                    "cells_before": 0,
                    "flagged_cells": 0,
                },
            )
            metric_record["cells_before"] += 1

            if flagged:
                metric_record["flagged_cells"] += 1
                metric_sets[metric].add(cell_id)
                sample_metric_sets.setdefault(
                    metric_key,
                    set(),
                ).add(cell_id)

        if any(metric_flags.values()):
            sample_record["flagged_any"] += 1
            flagged_any_cells.add(cell_id)

    by_sample_rows = []

    for sample in sorted(by_sample):
        row = by_sample[sample]
        cells_before = int(row["cells_before"])
        flagged_any = int(row["flagged_any"])
        retained = cells_before - flagged_any
        row["retained"] = retained
        row["flagged_fraction"] = ratio(flagged_any, cells_before)
        by_sample_rows.append(row)

    by_metric_rows = []

    for sample, metric in sorted(by_metric):
        row = by_metric[(sample, metric)]
        cells_before = int(row["cells_before"])
        flagged_cells = int(row["flagged_cells"])
        metric_cell_set = sample_metric_sets.get(
            (sample, metric),
            set(),
        )
        other_metric_cells = set().union(
            *[
                sample_metric_sets.get((sample, other_metric), set())
                for other_metric in thresholds
                if other_metric != metric
            ]
        )
        unique_flagged_cells = len(
            metric_cell_set - other_metric_cells
        )
        overlapped_flagged_cells = (
            flagged_cells - unique_flagged_cells
        )
        row["flagged_fraction"] = ratio(flagged_cells, cells_before)
        row["unique_flagged_cells"] = unique_flagged_cells
        row["unique_flagged_fraction"] = ratio(
            unique_flagged_cells,
            flagged_cells,
        )
        row["overlapped_flagged_cells"] = overlapped_flagged_cells
        by_metric_rows.append(row)

    global_metric_rows = []

    for metric in sorted(thresholds):
        metric_cell_set = metric_sets[metric]
        other_metric_cells = set().union(
            *[
                metric_sets[other_metric]
                for other_metric in thresholds
                if other_metric != metric
            ]
        )
        flagged_cells = len(metric_cell_set)
        unique_flagged_cells = len(
            metric_cell_set - other_metric_cells
        )
        global_metric_rows.append(
            {
                "sample": "ALL",
                "Group": "ALL",
                "metric": metric,
                "cells_before": len(cell_rows),
                "flagged_cells": flagged_cells,
                "flagged_fraction": ratio(
                    flagged_cells,
                    len(cell_rows),
                ),
                "unique_flagged_cells": unique_flagged_cells,
                "unique_flagged_fraction": ratio(
                    unique_flagged_cells,
                    flagged_cells,
                ),
                "overlapped_flagged_cells": (
                    flagged_cells - unique_flagged_cells
                ),
            }
        )

    by_metric_rows = global_metric_rows + by_metric_rows

    lower_tail_summary_rows = []

    for metric, probe_threshold in sorted(lower_tail_global):
        row = lower_tail_global[(metric, probe_threshold)]
        cells_evaluated = int(row["cells_evaluated"])
        cells_below = int(row["cells_below_threshold"])
        row["fraction_below_threshold"] = ratio(
            cells_below,
            cells_evaluated,
        )
        lower_tail_summary_rows.append(row)

    lower_tail_by_sample_rows = []

    for sample, metric, probe_threshold in sorted(lower_tail_by_sample):
        row = lower_tail_by_sample[
            (sample, metric, probe_threshold)
        ]
        cells_evaluated = int(row["cells_evaluated"])
        cells_below = int(row["cells_below_threshold"])
        row["fraction_below_threshold"] = ratio(
            cells_below,
            cells_evaluated,
        )
        lower_tail_by_sample_rows.append(row)

    overlap_rows = []

    for metric_a, metric_b in itertools.combinations(
        sorted(thresholds),
        2,
    ):
        set_a = metric_sets[metric_a]
        set_b = metric_sets[metric_b]
        flagged_both = len(set_a & set_b)
        flagged_union = len(set_a | set_b)
        overlap_rows.append(
            {
                "metric_a": metric_a,
                "metric_b": metric_b,
                "metric_a_flagged": len(set_a),
                "metric_b_flagged": len(set_b),
                "flagged_both": flagged_both,
                "flagged_union": flagged_union,
                "overlap_fraction_of_metric_a": ratio(
                    flagged_both,
                    len(set_a),
                ),
                "overlap_fraction_of_metric_b": ratio(
                    flagged_both,
                    len(set_b),
                ),
            }
        )

    cell_count = len(cell_rows)
    flagged_any = len(flagged_any_cells)
    retained = cell_count - flagged_any

    summary_row = {
        "plan_id": policy.get("plan_id", ""),
        "project_id": policy.get("project_id", ""),
        "biological_run_id": policy.get("biological_run_id", ""),
        "threshold_mode": policy.get("threshold_mode", ""),
        "approval_status": (
            policy.get("approval", {}).get("status")
            if isinstance(policy.get("approval"), dict)
            else ""
        ),
        "cell_count": cell_count,
        "flagged_any": flagged_any,
        "retained": retained,
        "flagged_fraction": ratio(flagged_any, cell_count),
        "enabled_metric_count": len(thresholds),
    }

    impact_document = {
        "qc_threshold_impact": {
            "schema_version": "1.0",
            "plan_id": summary_row["plan_id"],
            "project_id": summary_row["project_id"],
            "biological_run_id": summary_row["biological_run_id"],
            "threshold_mode": summary_row["threshold_mode"],
            "approval_status": summary_row["approval_status"],
            "application_status": "not_applied",
            "cell_count": cell_count,
            "flagged_any": flagged_any,
            "retained": retained,
            "flagged_fraction": summary_row["flagged_fraction"],
            "enabled_metrics": sorted(thresholds),
            "lower_tail_diagnostics": {
                "status": "available"
                if lower_tail_probes
                else "not_available",
                "metrics": sorted(lower_tail_probes),
            },
            "inputs": {
                "policy_yaml": str(args.policy_yaml),
                "cell_metadata": str(args.cell_metadata),
                "by_sample": str(args.by_sample),
            },
        }
    }

    args.impact_yaml.parent.mkdir(parents=True, exist_ok=True)
    args.impact_yaml.write_text(
        yaml.safe_dump(
            impact_document,
            sort_keys=False,
            allow_unicode=True,
        ),
        encoding="utf-8",
    )

    write_csv(
        args.by_sample_output,
        by_sample_rows,
        [
            "sample",
            "Group",
            "cells_before",
            "flagged_any",
            "retained",
            "flagged_fraction",
        ],
    )
    write_csv(
        args.by_metric_output,
        by_metric_rows,
        [
            "sample",
            "Group",
            "metric",
            "cells_before",
            "flagged_cells",
            "flagged_fraction",
            "unique_flagged_cells",
            "unique_flagged_fraction",
            "overlapped_flagged_cells",
        ],
    )
    write_csv(
        args.overlap_output,
        overlap_rows,
        [
            "metric_a",
            "metric_b",
            "metric_a_flagged",
            "metric_b_flagged",
            "flagged_both",
            "flagged_union",
            "overlap_fraction_of_metric_a",
            "overlap_fraction_of_metric_b",
        ],
    )
    write_csv(
        args.lower_tail_summary_output,
        lower_tail_summary_rows,
        [
            "metric",
            "threshold",
            "cells_evaluated",
            "cells_below_threshold",
            "fraction_below_threshold",
        ],
    )
    write_csv(
        args.lower_tail_by_sample_output,
        lower_tail_by_sample_rows,
        [
            "sample",
            "Group",
            "metric",
            "threshold",
            "cells_evaluated",
            "cells_below_threshold",
            "fraction_below_threshold",
        ],
    )
    write_summary_tsv(args.summary, summary_row)

    args.session_info.write_text(
        "\n".join(
            [
                f"python={platform.python_version()}",
                f"platform={platform.platform()}",
                "process=QC_THRESHOLD_IMPACT_SUMMARY",
                f"threshold_mode={summary_row['threshold_mode']}",
                f"cell_count={cell_count}",
                f"flagged_any={flagged_any}",
                f"enabled_metric_count={len(thresholds)}",
                "lower_tail_metrics="
                + ",".join(sorted(lower_tail_probes)),
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print("QC threshold impact summarized.")
    print(f"Cells: {cell_count}")
    print(f"Flagged by any enabled metric: {flagged_any}")
    print(f"Retained if applied: {retained}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

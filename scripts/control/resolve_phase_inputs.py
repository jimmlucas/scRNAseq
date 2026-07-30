#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import shlex
import sys
from pathlib import Path
from typing import Any


ARTIFACTS = {
    "loaded_object": "00_load_data/seu.rds",
    "qc_metrics_object": "01a_qc_metrics/seu_qc_metrics.rds",
    "qc_threshold_policy": (
        "control/qc_thresholds/qc_threshold_policy.yaml"
    ),
    "qc_threshold_impact": (
        "control/qc_threshold_impact/qc_threshold_impact.yaml"
    ),
    "qc_threshold_impact_summary": (
        "control/qc_threshold_impact/qc_threshold_impact.summary.tsv"
    ),
    "approved_qc_threshold_policy": (
        "control/qc_threshold_approval/approved_qc_threshold_policy.yaml"
    ),
    "qc_filtered_object": (
        "01e_qc_filter_policy/seu_qc_filtered.rds"
    ),
    "qc_filter_summary": (
        "01e_qc_filter_policy/01e_qc_filter_summary.csv"
    ),
    "sample_objects": (
        "01f_split_seurat_by_sample/sample_objects"
    ),
    "sample_manifest": (
        "01f_split_seurat_by_sample/01f_split_by_sample_manifest.tsv"
    ),
    "split_by_sample_summary": (
        "01f_split_seurat_by_sample/01f_split_by_sample_summary.csv"
    ),
    "doublet_scored_object": (
        "01g_doublet_detection_scoring/seu_doublet_scored.rds"
    ),
    "doublet_detection_summary": (
        "01g_doublet_detection_scoring/01g_doublet_detection_summary.csv"
    ),
    "doublet_detection_by_sample": (
        "01g_doublet_detection_scoring/01g_doublet_detection_by_sample.csv"
    ),
    "doublet_detection_cell_scores": (
        "01g_doublet_detection_scoring/01g_doublet_detection_cell_scores.csv"
    ),
    "doublet_detection_performance": (
        "01g_doublet_detection_scoring/01g_doublet_detection_performance.csv"
    ),
    "doublet_policy": "control/doublet_policy/doublet_policy.yaml",
    "approved_doublet_policy": (
        "control/doublet_policy_approval/approved_doublet_policy.yaml"
    ),
    "doublet_filtered_object": (
        "02_doublet_filter/seu_doublet_filtered.rds"
    ),
    "applied_doublet_policy": (
        "02_doublet_filter/applied_doublet_policy.yaml"
    ),
    "analysis_readiness": (
        "03_post_doublet_qc/03_analysis_readiness.yaml"
    ),
    "normalization_policy": (
        "control/normalization_policy/normalization_policy.yaml"
    ),
    "approved_normalization_policy": (
        "control/normalization_policy_approval/"
        "approved_normalization_policy.yaml"
    ),
    "normalized_object": "04_normalization/seu_normalized.rds",
    "applied_normalization_policy": (
        "04_normalization/applied_normalization_policy.yaml"
    ),
    "pca_object": "05_dimensionality/seu_pca.rds",
    "pca_elbow_metrics": (
        "05_dimensionality/05_pca_elbow_metrics.csv"
    ),
    "dimensionality_summary": (
        "05_dimensionality/05_dimensionality_summary.csv"
    ),
    "variable_features_summary": (
        "04_normalization/04_variable_features.csv"
    ),
    "pc_selection_policy": (
        "control/pc_selection_policy/pc_selection_policy.yaml"
    ),
    "approved_pc_selection_policy": (
        "control/pc_selection_policy_approval/"
        "approved_pc_selection_policy.yaml"
    ),
    "integration_policy": (
        "control/integration_policy/integration_policy.yaml"
    ),
    "doublet_config": None,
}


PHASE_INPUTS = {
    "qc_filter_policy_phase": {
        "input_object": "qc_metrics_object",
        "approved_qc_threshold_policy": "approved_qc_threshold_policy",
    },
    "split_seurat_by_sample_phase": {
        "input_object": "qc_filtered_object",
    },
    "doublet_detection_scoring_phase": {
        "input_object": "qc_filtered_object",
        "doublet_config": "doublet_config",
    },
    "doublet_detection_scoring_by_sample_phase": {
        "sample_manifest": "sample_manifest",
        "sample_objects": "sample_objects",
        "doublet_config": "doublet_config",
    },
    "doublet_policy_resolution_phase": {
        "summary": "doublet_detection_summary",
        "by_sample": "doublet_detection_by_sample",
        "cell_scores": "doublet_detection_cell_scores",
    },
    "doublet_filter_phase": {
        "input_object": "doublet_scored_object",
        "approved_doublet_policy": "approved_doublet_policy",
    },
    "post_doublet_qc_summary_phase": {
        "input_object": "doublet_filtered_object",
        "applied_doublet_policy": "applied_doublet_policy",
    },
    "normalization_policy_resolution_phase": {
        "input_object": "doublet_filtered_object",
        "readiness": "analysis_readiness",
        "applied_doublet_policy": "applied_doublet_policy",
    },
    "normalization_phase": {
        "input_object": "doublet_filtered_object",
        "approved_normalization_policy": "approved_normalization_policy",
    },
    "dimensionality_phase": {
        "input_object": "normalized_object",
        "applied_normalization_policy": "applied_normalization_policy",
    },
    "pc_selection_policy_resolution_phase": {
        "pca_elbow_metrics": "pca_elbow_metrics",
        "dimensionality_summary": "dimensionality_summary",
        "variable_features_summary": "variable_features_summary",
    },
    "integration_policy_resolution_phase": {
        "input_object": "pca_object",
        "approved_pc_selection_policy": "approved_pc_selection_policy",
        "analysis_readiness": "analysis_readiness",
        "applied_normalization_policy": "applied_normalization_policy",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Resolve canonical pipeline artifact paths from an outdir. "
            "This is a control-plane helper; it does not execute Nextflow."
        )
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Pipeline output directory.",
    )
    parser.add_argument(
        "--phase",
        choices=sorted(PHASE_INPUTS),
        help="Return only the input paths required by one phase.",
    )
    parser.add_argument(
        "--artifact",
        action="append",
        choices=sorted(ARTIFACTS),
        help="Return one or more specific artifact paths.",
    )
    parser.add_argument(
        "--project-id",
        help=(
            "Project id used to resolve conventional config paths such "
            "as configs/real_validation/<project_id>/doublet_config.yaml."
        ),
    )
    parser.add_argument(
        "--format",
        choices=[
            "json",
            "shell",
            "tsv",
        ],
        default="shell",
        help="Output format.",
    )
    parser.add_argument(
        "--absolute",
        action="store_true",
        help="Print absolute paths.",
    )
    parser.add_argument(
        "--require-existing",
        action="store_true",
        help="Fail if any resolved path does not exist.",
    )

    return parser.parse_args()


def path_for_artifact(
    artifact: str,
    outdir: Path,
    project_id: str | None,
) -> Path:
    if artifact == "doublet_config":
        if project_id:
            return (
                Path("configs")
                / "real_validation"
                / project_id
                / "doublet_config.yaml"
            )

        outdir_name = outdir.name
        return (
            Path("configs")
            / "real_validation"
            / outdir_name
            / "doublet_config.yaml"
        )

    if artifact not in ARTIFACTS:
        raise ValueError(f"Unknown artifact: {artifact}")

    return outdir / ARTIFACTS[artifact]


def display_path(
    path: Path,
    absolute: bool,
) -> str:
    if absolute:
        return str(path.resolve())

    return str(path)


def resolve_requested(
    args: argparse.Namespace,
) -> dict[str, dict[str, Any]]:
    outdir = args.outdir
    requested: dict[str, str]

    if args.phase:
        requested = PHASE_INPUTS[args.phase]
    elif args.artifact:
        requested = {
            artifact: artifact
            for artifact in args.artifact
        }
    else:
        requested = {
            artifact: artifact
            for artifact in ARTIFACTS
        }

    resolved = {}

    for parameter, artifact in requested.items():
        path = path_for_artifact(
            artifact,
            outdir,
            args.project_id,
        )
        exists = path.exists()

        resolved[parameter] = {
            "artifact": artifact,
            "path": display_path(path, args.absolute),
            "exists": exists,
        }

    return resolved


def print_shell(
    resolved: dict[str, dict[str, Any]],
) -> None:
    for parameter, record in resolved.items():
        variable = parameter.upper()
        print(
            f"{variable}={shlex.quote(str(record['path']))}"
        )


def print_tsv(
    resolved: dict[str, dict[str, Any]],
) -> None:
    print("parameter\tartifact\tpath\texists")

    for parameter, record in resolved.items():
        print(
            "\t".join(
                [
                    parameter,
                    str(record["artifact"]),
                    str(record["path"]),
                    str(record["exists"]).lower(),
                ]
            )
        )


def main() -> int:
    args = parse_args()

    try:
        resolved = resolve_requested(args)

        missing = [
            (parameter, record)
            for parameter, record in resolved.items()
            if not record["exists"]
        ]

        if args.format == "json":
            print(
                json.dumps(
                    {
                        "outdir": display_path(
                            args.outdir,
                            args.absolute,
                        ),
                        "phase": args.phase,
                        "paths": resolved,
                    },
                    indent=2,
                    sort_keys=True,
                )
            )
        elif args.format == "tsv":
            print_tsv(resolved)
        else:
            print_shell(resolved)

        if args.require_existing and missing:
            print(
                "ERROR: Missing required resolved path(s):",
                file=sys.stderr,
            )
            for parameter, record in missing:
                print(
                    f"- {parameter} ({record['artifact']}): "
                    f"{record['path']}",
                    file=sys.stderr,
                )
            return 1

    except ValueError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

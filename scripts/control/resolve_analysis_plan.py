#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Any

import yaml


ROUTING_TO_PIPELINE_MODULE = {
    "input_discovery": "INPUT_DISCOVERY",
    "load_data": "LOAD_DATA",
    "qc_metrics": "QC_METRICS",
    "qc_visualization": "QC_VISUALIZATION",
    "qc_threshold_resolution": "QC_THRESHOLD_RESOLUTION",
    "qc_filter": "QC_FILTER",
    "doublet_detection": "QC_DOUBLET_DETECTION",
    "sample_resolution": "SAMPLE_RESOLUTION",
    "merge_samples": "MERGE_SAMPLES",
    "qc_scores": "QC_SCORES",
    "normalization": "NORMALIZATION",
}


# This records the state of the plan-driven pipeline architecture.
#
# available:
#   A usable Nextflow implementation currently exists.
#
# legacy_available:
#   A module exists, but still belongs to the old architecture and
#   must be refactored before being considered native to the new flow.
#
# planned:
#   The module is part of the target design but is not implemented.
MODULE_IMPLEMENTATION_STATE = {
    "INPUT_DISCOVERY": "available",
    "LOAD_DATA": "available",
    "QC_METRICS": "available",
    "QC_VISUALIZATION": "legacy_available",
    "QC_THRESHOLD_RESOLUTION": "available",
    "QC_FILTER": "legacy_available",
    "QC_DOUBLET_DETECTION": "legacy_available",
    "SAMPLE_RESOLUTION": "planned",
    "MERGE_SAMPLES": "planned",
    "QC_SCORES": "legacy_available",
    "NORMALIZATION": "planned",
}


APPROVAL_ROUTE_MAP = {
    "qc_threshold_resolution": "thresholds",
    "sample_resolution": "sample_resolution",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Resolve requested analysis routes against "
            "the current Nextflow implementation."
        )
    )

    parser.add_argument(
        "--plan",
        required=True,
        type=Path,
        help="Validated analysis plan YAML.",
    )
    parser.add_argument(
        "--module-catalog",
        required=True,
        type=Path,
        help="Pipeline module catalog YAML.",
    )
    parser.add_argument(
        "--routing-output",
        required=True,
        type=Path,
        help="Resolved routing YAML.",
    )
    parser.add_argument(
        "--summary",
        required=True,
        type=Path,
        help="Resolved routing TSV.",
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


def resolve_approval(
    route_name: str,
    qc_section: dict[str, Any],
    plan_status: str,
) -> tuple[str, bool]:
    approval_key = APPROVAL_ROUTE_MAP.get(route_name)

    if approval_key is None:
        return "not_required", True

    approval_policy = qc_section.get("approval_policy", {})
    if not isinstance(approval_policy, dict):
        approval_policy = {}

    requested_policy = approval_policy.get(
        approval_key,
        "automatic",
    )

    if requested_policy == "automatic":
        return "automatic", True

    if requested_policy == "review":
        return "review", plan_status == "approved"

    if requested_policy == "required":
        return "required", plan_status == "approved"

    return str(requested_policy), False


def main() -> int:
    args = parse_args()

    try:
        plan_document = read_yaml(
            args.plan,
            "Validated analysis plan",
        )
        catalog_document = read_yaml(
            args.module_catalog,
            "Module catalog",
        )
    except (FileNotFoundError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

    validation = plan_document.get("validation", {})
    if not isinstance(validation, dict):
        validation = {}

    if validation.get("status") != "valid":
        print(
            "ERROR: Analysis plan has not passed validation.",
            file=sys.stderr,
        )
        return 1

    plan_section = plan_document.get("plan", {})
    routing_section = plan_document.get("routing", {})
    qc_section = plan_document.get("qc", {})

    if not isinstance(plan_section, dict):
        plan_section = {}

    if not isinstance(routing_section, dict):
        print(
            "ERROR: routing section must be a mapping.",
            file=sys.stderr,
        )
        return 1

    if not isinstance(qc_section, dict):
        qc_section = {}

    plan_id = str(plan_section.get("plan_id", ""))
    plan_status = str(
        plan_section.get("approval_status", "draft")
    )

    catalog_modules = catalog_document.get("modules", {})
    catalog_future = catalog_document.get(
        "future_modules",
        {},
    )

    if not isinstance(catalog_modules, dict):
        catalog_modules = {}

    if not isinstance(catalog_future, dict):
        catalog_future = {}

    resolved_routes: dict[str, dict[str, Any]] = {}
    blocking_routes: list[str] = []
    executable_routes: list[str] = []

    for route_name, module_name in ROUTING_TO_PIPELINE_MODULE.items():
        requested = routing_section.get(route_name, False)

        if not isinstance(requested, bool):
            print(
                f"ERROR: routing.{route_name} must be boolean.",
                file=sys.stderr,
            )
            return 1

        implementation_state = MODULE_IMPLEMENTATION_STATE.get(
            module_name,
            "unknown",
        )

        catalog_state = (
            "implemented"
            if module_name in catalog_modules
            else (
                "future"
                if module_name in catalog_future
                else "not_catalogued"
            )
        )

        approval_policy, approval_satisfied = resolve_approval(
            route_name=route_name,
            qc_section=qc_section,
            plan_status=plan_status,
        )

        native_available = implementation_state == "available"
        legacy_available = (
            implementation_state == "legacy_available"
        )

        executable = (
            requested
            and native_available
            and approval_satisfied
        )

        if not requested:
            reason = "not_requested"

        elif not approval_satisfied:
            reason = "approval_pending"

        elif implementation_state == "planned":
            reason = "module_not_implemented"

        elif legacy_available:
            reason = "legacy_module_requires_refactor"

        elif implementation_state == "available":
            reason = "ready"

        else:
            reason = "unknown_module_state"

        blocking = requested and not executable

        if executable:
            executable_routes.append(route_name)

        if blocking:
            blocking_routes.append(route_name)

        resolved_routes[route_name] = {
            "module": module_name,
            "requested": requested,
            "implementation_state": implementation_state,
            "catalog_state": catalog_state,
            "approval_policy": approval_policy,
            "approval_satisfied": approval_satisfied,
            "executable": executable,
            "blocking": blocking,
            "reason": reason,
        }

    routing_status = (
        "ready"
        if not blocking_routes
        else "blocked"
    )

    output_document = {
        "routing_resolution": {
            "schema_version": "1.0",
            "plan_id": plan_id,
            "plan_approval_status": plan_status,
            "status": routing_status,
            "executable_routes": executable_routes,
            "blocking_routes": blocking_routes,
        },
        "routes": resolved_routes,
    }

    args.routing_output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    args.summary.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with args.routing_output.open(
        "w",
        encoding="utf-8",
    ) as handle:
        yaml.safe_dump(
            output_document,
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
                "route",
                "module",
                "requested",
                "implementation_state",
                "catalog_state",
                "approval_policy",
                "approval_satisfied",
                "executable",
                "blocking",
                "reason",
            ]
        )

        for route_name, route in resolved_routes.items():
            writer.writerow(
                [
                    route_name,
                    route["module"],
                    str(route["requested"]).lower(),
                    route["implementation_state"],
                    route["catalog_state"],
                    route["approval_policy"],
                    str(
                        route["approval_satisfied"]
                    ).lower(),
                    str(route["executable"]).lower(),
                    str(route["blocking"]).lower(),
                    route["reason"],
                ]
            )

    print(
        f"Analysis plan routing resolved: {routing_status}"
    )

    if blocking_routes:
        print(
            "Requested routes currently blocked: "
            + ", ".join(blocking_routes)
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

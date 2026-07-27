#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path
from typing import Any

import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Resolve and validate the LOAD_DATA configuration "
            "from an analysis plan."
        )
    )

    parser.add_argument(
        "--plan",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--routing",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--inventory",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--project-dir",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--project-id",
        required=True,
    )
    parser.add_argument(
        "--biological-run-id",
        required=True,
    )
    parser.add_argument(
        "--resolved-config",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--summary",
        required=True,
        type=Path,
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


def resolve_path(
    value: Any,
    project_dir: Path,
) -> Path | None:
    if value is None or str(value).strip() == "":
        return None

    path = Path(str(value)).expanduser()

    if not path.is_absolute():
        path = project_dir / path

    return path.resolve()


def write_summary(
    path: Path,
    *,
    plan_id: str,
    status: str,
    requested: bool,
    executable: bool,
    reason: str,
    resolved_config: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")

        writer.writerow(
            [
                "plan_id",
                "status",
                "requested",
                "executable",
                "reason",
                "resolved_config",
            ]
        )

        writer.writerow(
            [
                plan_id,
                status,
                str(requested).lower(),
                str(executable).lower(),
                reason,
                resolved_config,
            ]
        )


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
            "Resolved routing",
        )
        inventory_document = read_yaml(
            args.inventory,
            "Input inventory",
        )
    except (FileNotFoundError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

    validation = plan_document.get("validation", {})
    if not isinstance(validation, dict):
        validation = {}

    if validation.get("status") != "valid":
        print(
            "ERROR: LOAD_DATA requires a valid analysis plan.",
            file=sys.stderr,
        )
        return 1

    plan = plan_document.get("plan", {})
    experiment = plan_document.get("experiment", {})
    inputs = plan_document.get("inputs", {})
    routes = routing_document.get("routes", {})
    discovery = inventory_document.get(
        "input_discovery",
        {},
    )

    if not isinstance(plan, dict):
        plan = {}

    if not isinstance(experiment, dict):
        experiment = {}

    if not isinstance(inputs, dict):
        inputs = {}

    if not isinstance(routes, dict):
        routes = {}

    if not isinstance(discovery, dict):
        discovery = {}

    plan_id = str(plan.get("plan_id", ""))

    plan_project_id = str(
        experiment.get("project_id", "")
    )
    plan_biological_run_id = str(
        experiment.get("biological_run_id", "")
    )

    if plan_project_id != args.project_id:
        print(
            "ERROR: Project identity mismatch: "
            f"plan={plan_project_id}, "
            f"execution={args.project_id}",
            file=sys.stderr,
        )
        return 1

    if plan_biological_run_id != args.biological_run_id:
        print(
            "ERROR: Biological run identity mismatch: "
            f"plan={plan_biological_run_id}, "
            f"execution={args.biological_run_id}",
            file=sys.stderr,
        )
        return 1

    load_route = routes.get("load_data", {})
    if not isinstance(load_route, dict):
        load_route = {}

    requested = bool(
        load_route.get("requested", False)
    )
    executable = bool(
        load_route.get("executable", False)
    )

    if not requested:
        write_summary(
            args.summary,
            plan_id=plan_id,
            status="skipped",
            requested=False,
            executable=False,
            reason="not_requested",
            resolved_config="",
        )

        print("LOAD_DATA skipped: not requested.")
        return 0

    if not executable:
        reason = str(
            load_route.get(
                "reason",
                "route_not_executable",
            )
        )

        write_summary(
            args.summary,
            plan_id=plan_id,
            status="blocked",
            requested=True,
            executable=False,
            reason=reason,
            resolved_config="",
        )

        print(
            f"ERROR: LOAD_DATA route is blocked: {reason}",
            file=sys.stderr,
        )
        return 1

    discovery_status = str(
        discovery.get("status", "")
    )

    if discovery_status != "ready":
        write_summary(
            args.summary,
            plan_id=plan_id,
            status="blocked",
            requested=True,
            executable=True,
            reason="input_discovery_not_ready",
            resolved_config="",
        )

        print(
            "ERROR: LOAD_DATA requires input discovery "
            "status 'ready'.",
            file=sys.stderr,
        )
        return 1

    load_config_path = resolve_path(
        inputs.get("load_config"),
        project_dir,
    )

    if load_config_path is None:
        print(
            "ERROR: inputs.load_config is not defined.",
            file=sys.stderr,
        )
        return 1

    if not load_config_path.is_file():
        print(
            f"ERROR: Load configuration not found: "
            f"{load_config_path}",
            file=sys.stderr,
        )
        return 1

    try:
        load_config = read_yaml(
            load_config_path,
            "Load configuration",
        )
    except (FileNotFoundError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    config_project = load_config.get("project", {})
    if not isinstance(config_project, dict):
        config_project = {}

    config_organism = str(
        config_project.get("organism", "")
    )
    plan_organism = str(
        experiment.get("organism", "")
    )

    if (
        config_organism
        and plan_organism
        and config_organism != plan_organism
    ):
        print(
            "ERROR: Organism mismatch between analysis plan "
            f"and load configuration: plan={plan_organism}, "
            f"config={config_organism}",
            file=sys.stderr,
        )
        return 1

    args.resolved_config.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    shutil.copyfile(
        load_config_path,
        args.resolved_config,
    )

    write_summary(
        args.summary,
        plan_id=plan_id,
        status="ready",
        requested=True,
        executable=True,
        reason="ready",
        resolved_config=str(load_config_path),
    )

    print(
        f"LOAD_DATA configuration resolved: "
        f"{load_config_path}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

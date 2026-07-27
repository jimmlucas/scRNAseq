#!/usr/bin/env python3

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path


def load_resolver(root: Path):
    resolver_path = root / "scripts/control/resolve_pc_selection_policy.py"
    spec = importlib.util.spec_from_file_location(
        "resolve_pc_selection_policy",
        resolver_path,
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def main() -> int:
    root = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    resolver = load_resolver(root)

    policy = resolver.read_yaml(
        output_dir / "pc_selection_policy.yaml",
        "PC-selection policy",
    )["pc_selection_policy"]

    assert policy["status"] == "proposed"
    assert int(policy["input"]["cells"]) == 363
    assert int(policy["input"]["features"]) == 27522
    assert int(policy["input"]["variable_features"]) == 3036
    assert int(policy["input"]["max_computed_pcs"]) == 30

    recommendation = policy["recommendation"]
    assert int(recommendation["recommended_pcs"]) == 21
    assert int(recommendation["max_computed_pcs"]) == 30
    assert recommendation["method"] == (
        "cumulative_variance_target_with_minimum_pc_floor"
    )
    assert float(recommendation["cumulative_variance_target"]) == 80.0
    assert float(recommendation["achieved_cumulative_variance"]) >= 80.0

    assert policy["review"]["required"] is True
    assert policy["approval"]["status"] == "proposed"
    assert policy["application"]["status"] == "not_applied"
    assert policy["application"]["pc_selection_applied"] is False

    with (output_dir / "pc_selection_policy.tsv").open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 1
    row = rows[0]
    assert int(row["recommended_pcs"]) == 21
    assert int(row["max_computed_pcs"]) == 30
    assert row["approval_status"] == "proposed"
    assert row["application_status"] == "not_applied"

    print("PC-selection policy contract: VALID")
    print("Recommended PCs:", recommendation["recommended_pcs"])
    print("Max computed PCs:", recommendation["max_computed_pcs"])
    print("Method:", recommendation["method"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

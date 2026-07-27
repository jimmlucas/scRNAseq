#!/usr/bin/env python3

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path


def load_resolver(root: Path):
    resolver_path = root / "scripts/control/resolve_integration_policy.py"
    spec = importlib.util.spec_from_file_location(
        "resolve_integration_policy",
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
        output_dir / "integration_policy.yaml",
        "Integration policy",
    )["integration_policy"]

    assert policy["status"] == "proposed"
    assert int(policy["input"]["cells"]) == 363
    assert int(policy["input"]["features"]) == 27522
    assert int(policy["input"]["samples"]) == 2
    assert int(policy["input"]["selected_pcs"]) == 21
    assert int(policy["input"]["max_computed_pcs"]) == 30

    strategy = policy["strategy"]
    assert strategy["integration_required"] is True
    assert strategy["method"] == "Harmony"
    assert strategy["input_reduction"] == "pca"
    assert strategy["output_reduction"] == "harmony"
    assert int(strategy["dimensions_to_use"]) == 21
    assert int(strategy["dimensions_to_use"]) > 0
    assert int(strategy["dimensions_to_use"]) <= int(policy["input"]["max_computed_pcs"])
    assert strategy["output_reduction"] != strategy["input_reduction"]
    assert strategy["batch_variable"] == "sample"
    assert strategy["batch_variable"]
    assert strategy["multimodal_strategy"] == "none"

    evidence = policy["evidence"]
    assert int(evidence["samples_detected"]) == 2
    assert int(evidence["selected_pcs"]) == 21
    assert evidence["normalization_method"] == "LogNormalize"
    assert evidence["batch_variable_source"] == "default_mvp_sample_metadata"

    assert policy["approval"]["status"] == "proposed"
    assert policy["application"]["status"] == "not_applied"
    assert policy["application"]["integration_applied"] is False
    assert "Harmony" in policy["allowed_methods"]["available"]
    assert "none" in policy["allowed_methods"]["available"]
    assert "WNN" in policy["allowed_methods"]["planned"]

    with (output_dir / "integration_policy.tsv").open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 1
    row = rows[0]
    assert row["method"] == "Harmony"
    assert row["integration_required"] == "True"
    assert int(row["selected_pcs"]) == 21
    assert row["approval_status"] == "proposed"
    assert row["application_status"] == "not_applied"

    print("Integration policy contract: VALID")
    print("Method:", strategy["method"])
    print("Integration required:", strategy["integration_required"])
    print("Dimensions:", strategy["dimensions_to_use"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

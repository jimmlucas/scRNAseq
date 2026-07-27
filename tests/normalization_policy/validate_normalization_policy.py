#!/usr/bin/env python3

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path


def load_resolver(root: Path):
    resolver_path = root / "scripts/control/resolve_normalization_policy.py"
    spec = importlib.util.spec_from_file_location(
        "resolve_normalization_policy",
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
        output_dir / "normalization_policy.yaml",
        "Normalization policy",
    )["normalization_policy"]

    assert policy["status"] == "proposed"
    assert int(policy["input"]["cells"]) == 363
    assert int(policy["input"]["features"]) == 27522
    assert int(policy["input"]["samples"]) == 2
    assert policy["input"]["readiness"] == "ready"

    strategy = policy["strategy"]
    assert strategy["method"] == "LogNormalize"
    assert strategy["variable_features_method"] == "vst"
    assert int(strategy["nfeatures"]) == 2000
    assert int(strategy["scale_factor"]) == 10000
    assert strategy["scope"] == "per_sample"
    assert strategy["assay"] == "RNA"

    assert policy["integration"]["required"] is True
    assert policy["integration"]["reason"] == "multiple_samples_detected"
    assert policy["review"]["required"] is True
    assert policy["approval"]["status"] == "proposed"
    assert policy["application"]["status"] == "not_applied"
    assert policy["application"]["normalization_applied"] is False

    with (output_dir / "normalization_policy.tsv").open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 1
    row = rows[0]
    assert row["method"] == "LogNormalize"
    assert row["scope"] == "per_sample"
    assert row["integration_required"] == "True"
    assert row["approval_status"] == "proposed"
    assert row["application_status"] == "not_applied"

    print("Normalization policy contract: VALID")
    print("Method: LogNormalize")
    print("Variable features: vst")
    print("Scope: per_sample")
    print("Integration required: true")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

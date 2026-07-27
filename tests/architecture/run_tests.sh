#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

python3 - "${ROOT_DIR}" <<'PY'
from __future__ import annotations

import re
import sys
from pathlib import Path

root = Path(sys.argv[1])

def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


engine_paths = [
    root / "main.nf",
    root / "nextflow.config",
    root / "Dockerfile.control-python",
    root / "control",
    root / "modules",
    root / "scripts/control",
    root / "subworkflow",
]

for path in engine_paths:
    paths = [path]
    if path.is_dir():
        paths = [
            child
            for child in path.rglob("*")
            if child.is_file()
            and child.suffix in {".nf", ".py", ".yaml", ".yml"}
        ]

    for child in paths:
        text = read_text(child)
        if re.search(r"\bBBUDY\b|\bbbudy\b|\bbudy\b", text):
            raise AssertionError(
                f"Active pipeline engine contains supervisor-specific "
                f"naming: {child.relative_to(root)}"
            )

main_text = read_text(root / "main.nf")

assert "analysis_plan_qc_process" in main_text
assert "bbudy_qc_process" not in main_text
assert "subworkflow/control/analysis_plan_qc_process" in main_text

assert (root / "subworkflow/control/analysis_plan_qc_process.nf").exists()
assert not (root / "subworkflow/bbudy").exists()

module_catalog_text = read_text(root / "control/module_catalog.yaml")
execution_modes_text = read_text(root / "control/execution_modes.yaml")
causistics_catalog_text = read_text(root / "control/causistics_catalog.yaml")
legacy_patterns_text = read_text(
    root / "docs/analysis_patterns_from_legacy_rmds.md"
)

assert 'id: scrna_seq' in module_catalog_text
assert 'orchestrator_awareness: none' in execution_modes_text
assert 'causistics_contract: structured_evidence_flags' in execution_modes_text
assert 'analysis_plan_qc_process' in execution_modes_text
assert 'normalization_approval_phase' in execution_modes_text
assert 'scRNA-seq:' in execution_modes_text
assert 'snRNA-seq:' in execution_modes_text
assert 'CITE-seq:' in execution_modes_text
assert execution_modes_text.count('requires_human_review: true') >= 3

for causistic_id in [
    "high_globin_content",
    "repeated_sample_across_runs",
    "sampletag_multiplet",
    "ab_marker_removed",
    "marker_based_doublet_signal",
    "integration_required",
    "multimodal_wnn_applicable",
    "pseudobulk_design_required",
    "bcr_metadata_available",
]:
    assert causistic_id in causistics_catalog_text

for expected_phrase in [
    "Mouse scRNA Patterns",
    "Human PBMC CITE-seq Patterns",
    "Analysis Type Profiles",
]:
    assert expected_phrase in legacy_patterns_text

print("Architecture contract: VALID")
PY

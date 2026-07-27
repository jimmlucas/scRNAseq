# Doublet detection scoring regression test

This test reproduces the validated `DOUBLET_DETECTION_SCORING`
execution using a fixed Seurat input object and a pinned Docker image.

## Validated environment

- Pipeline branch: current QC orchestration feature branch
- Docker image:
  `ghcr.io/jimmlucas/bioc-rnaseq-test01new2:v0.3.0`
- Docker digest:
  `sha256:b515735a14cfeb287ae628cfb3e9591492486afa521c9f0b2baba8fee2707b94`
- DoubletFinder: `2.0.6`
- Expected input cells: `392`
- Expected input features: `27,522`
- Expected samples: `2`

## Run

From the repository root:

```bash
tests/doublet_scoring/run_test.sh
```

## Expected result

- 392 cells preserved
- 27,522 features preserved
- 0 cells removed
- 29 algorithmic doublet proposals
- 363 algorithmic singlet proposals
- all cells marked as `completed`
- all decisions marked as `unreviewed`

The scoring phase annotates cells but does not remove them.

## Optional variables

```bash
OUTDIR=my_results \
EXECUTION_ID=EXTERNAL_MACHINE_TEST_001 \
tests/doublet_scoring/run_test.sh
```

The container is pinned by digest in `run_test.sh`.

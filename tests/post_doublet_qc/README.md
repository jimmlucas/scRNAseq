# Post-Doublet QC Readiness Test

This test exercises `POST_DOUBLET_QC_SUMMARY` as the non-destructive
checkpoint between QC/doublet handling and downstream analysis.

It verifies:

- final cells, features and samples;
- per-sample final cell counts;
- doublets removed by policy;
- post-filter QC metrics;
- empty warning table for the current fixture;
- `analysis_readiness.status: ready`;
- recommended next phase: `normalization_policy_resolution`.

Run:

```bash
tests/post_doublet_qc/run_test.sh
```

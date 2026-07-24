# Doublet Filter Regression Test

This test exercises `DOUBLET_FILTER` as a policy-application step.

It checks that filtering:

- consumes a doublet-scored Seurat object;
- consumes an approved, valid and not-yet-applied doublet policy;
- does not recalculate doublets;
- removes only cells authorized by per-sample `application_action`;
- preserves features;
- records removed cells in an audit table;
- writes a separate `applied_doublet_policy.yaml`.

Run:

```bash
tests/doublet_filter/run_test.sh
```

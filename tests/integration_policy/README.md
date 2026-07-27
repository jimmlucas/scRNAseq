# Integration Policy Resolution Test

This test validates `INTEGRATION_POLICY_RESOLUTION`.

It consumes:

- `seu_pca.rds`
- `05_dimensionality_summary.csv`
- `approved_pc_selection_policy.yaml`
- `03_analysis_readiness.yaml`
- `applied_normalization_policy.yaml`

and verifies that the resolver proposes an integration policy without running
Harmony, Seurat integration, neighbors, UMAP or clustering.

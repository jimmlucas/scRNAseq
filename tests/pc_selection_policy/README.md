# PC Selection Policy Resolution Test

This test validates `PC_SELECTION_POLICY_RESOLUTION`.

It consumes:

- `05_dimensionality_summary.csv`
- `05_variable_features_summary.csv`
- `05_pca_elbow_metrics.csv`

and verifies that the resolver proposes a reviewable PC-selection policy
without running neighbors, UMAP, clustering or integration.

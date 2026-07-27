# Dimensionality Phase Test

This test validates the `DIMENSIONALITY` technical application phase.

It consumes:

- `seu_normalized.rds`
- `applied_normalization_policy.yaml`

and verifies that scaling and PCA are applied without changing cell or feature
counts. This phase does not approve or decide integration; it publishes PCA
evidence for later readiness and integration policy steps.

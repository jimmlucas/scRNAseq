# Normalization Phase Test

This test validates the `NORMALIZATION` scientific application phase.

It consumes:

- `seu_doublet_filtered.rds`
- `approved_normalization_policy.yaml`

and verifies that normalization applies the approved policy without changing
cell or feature counts.


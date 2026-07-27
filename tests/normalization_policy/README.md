# Normalization Policy Resolution Test

This test exercises `NORMALIZATION_POLICY_RESOLUTION` as a proposal-only
control-plane step.

It verifies that the resolver:

- consumes post-doublet readiness outputs;
- proposes `LogNormalize` + `vst`;
- uses per-sample scope;
- marks integration as required for multiple samples;
- leaves approval as `proposed`;
- leaves application as `not_applied`;
- rejects a non-ready readiness contract.

Run:

```bash
tests/normalization_policy/run_test.sh
```

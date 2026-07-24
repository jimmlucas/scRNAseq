# Doublet Approval Validation Test

This test exercises `VALIDATE_DOUBLET_APPROVAL` without running the full
Nextflow graph.

It validates that a human-reviewed doublet policy can be normalized for a
future filtering step, while unsafe or inconsistent policies fail before any
cell-removal module can consume them.

Checked cases:

- approved policy passes;
- proposed policy fails;
- missing reviewer fails;
- missing approval timestamp fails;
- missing or unexpected sample fails through summary/sample consistency;
- changed doublet counts fail;
- unsupported action fails;
- mismatched project or biological run identity fails.

Run:

```bash
tests/doublet_approval/run_tests.sh
```

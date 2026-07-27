# Architecture Contract Test

This test protects the separation between the scientific Nextflow engine and
any external supervisor.

It verifies that active pipeline code uses neutral control-plane naming, that
the plan-driven mode is `analysis_plan_qc_process`, and that the public
contract files are parseable.


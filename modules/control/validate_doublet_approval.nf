nextflow.enable.dsl=2


process VALIDATE_DOUBLET_APPROVAL {

    tag "${plan_meta.project_id ?: 'project'}:${plan_meta.biological_run_id ?: 'run'}"

    label 'lightweight'

    publishDir(
        "${params.outdir}/control/doublet_policy_approval",
        mode: 'copy',
        overwrite: true
    )

    input:

    tuple val(plan_meta), path(validated_plan)
    tuple val(policy_meta), path(reviewed_policy)
    path validator_script

    output:

    tuple val(output_meta), path("approved_doublet_policy.yaml"),
        emit: approved_policy

    tuple val(summary_meta), path("doublet_policy_approval.summary.tsv"),
        emit: summary

    path "doublet_policy_approval.session.txt",
        emit: session_info

    script:

    if (
        plan_meta.project_id != policy_meta.project_id
        || plan_meta.biological_run_id != policy_meta.biological_run_id
    ) {
        error(
            "VALIDATE_DOUBLET_APPROVAL received incompatible " +
            "metadata between the analysis plan and reviewed policy."
        )
    }

    output_meta = plan_meta + [
        stage           : 'doublet_policy_approval',
        artifact_type   : 'approved_doublet_policy',
        approval_status : 'approved'
    ]

    summary_meta = plan_meta + [
        stage         : 'doublet_policy_approval',
        artifact_type : 'doublet_policy_approval_summary'
    ]

    """
    python3 ${validator_script} \
        --plan ${validated_plan} \
        --policy ${reviewed_policy} \
        --validated-policy approved_doublet_policy.yaml \
        --summary doublet_policy_approval.summary.tsv \
        --session-info doublet_policy_approval.session.txt
    """

    stub:

    output_meta = plan_meta + [
        stage           : 'doublet_policy_approval',
        artifact_type   : 'approved_doublet_policy',
        approval_status : 'approved'
    ]

    summary_meta = plan_meta + [
        stage         : 'doublet_policy_approval',
        artifact_type : 'doublet_policy_approval_summary'
    ]

    """
    cat > approved_doublet_policy.yaml <<'YAML'
    doublet_policy:
      schema_version: "1.0"
      identity:
        project_id: "stub_project"
        project_name: "stub_project"
        biological_run_id: "stub_run"
      summary:
        cells_scored: 1
        samples: 1
        proposed_doublets: 1
        cells_removed_by_scoring: 0
      samples:
        stub_sample:
          detection_status: "completed"
          decision_status: "approved"
          proposed_action: "exclude_algorithmic_doublets"
          cells_scored: 1
          proposed_doublet_count: 1
          proposed_singlet_count: 0
          expected_doublet_rate: 0.075
          pN: 0.25
          pK: 0.3
          nExp_raw: 1
          homotypic_prop: 0.0
          nExp_adjusted: 1
          cell_class_counts:
            Doublet: 1
            Singlet: 0
            NotRun: 0
            Other: 0
      approval:
        required: true
        status: "approved"
        reviewed_by: "stub_reviewer"
        approved_at: "2026-01-01T00:00:00+00:00"
      validation:
        status: "valid"
        validated_at: "2026-01-01T00:00:00+00:00"
        validator: "validate_doublet_approval.py"
        validated_for_doublet_filter: true
        sample_count: 1
        cells_scored: 1
        approved_doublet_count: 1
      application:
        status: "not_applied"
        filtering_applied: false
        cells_removed: 0
    YAML

    cat > doublet_policy_approval.summary.tsv <<'TSV'
    plan_id	project_id	biological_run_id	validation_status	approval_status	reviewed_by	approved_at	validated_at	sample_count	cells_scored	approved_doublet_count
    stub_plan	stub_project	stub_run	valid	approved	stub_reviewer	2026-01-01T00:00:00+00:00	2026-01-01T00:00:00+00:00	1	1	1
    TSV

    cat > doublet_policy_approval.session.txt <<'TXT'
    stub=true
    process=VALIDATE_DOUBLET_APPROVAL
    validation_status=valid
    approval_status=approved
    TXT
    """
}

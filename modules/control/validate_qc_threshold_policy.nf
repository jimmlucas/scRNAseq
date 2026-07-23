nextflow.enable.dsl=2


process VALIDATE_QC_THRESHOLD_APPROVAL {

    tag "${plan_meta.project_id ?: 'project'}:${plan_meta.biological_run_id ?: 'run'}"

    label 'lightweight'

    publishDir(
        "${params.outdir}/control/qc_threshold_approval",
        mode: 'copy',
        overwrite: true
    )

    input:

    tuple val(plan_meta), path(validated_plan)
    tuple val(policy_meta), path(reviewed_policy)
    path validator_script

    output:

    tuple val(output_meta), path("approved_qc_threshold_policy.yaml"),
        emit: approved_policy

    tuple val(summary_meta), path("qc_threshold_approval.summary.tsv"),
        emit: summary

    path "qc_threshold_approval.session.txt",
        emit: session_info

    script:

    if (
        plan_meta.project_id != policy_meta.project_id
        || plan_meta.biological_run_id != policy_meta.biological_run_id
    ) {
        error(
            "VALIDATE_QC_THRESHOLD_APPROVAL received incompatible " +
            "metadata between the analysis plan and reviewed policy."
        )
    }

    output_meta = plan_meta + [
        stage           : 'qc_threshold_approval',
        artifact_type   : 'approved_qc_threshold_policy',
        approval_status : 'approved'
    ]

    summary_meta = plan_meta + [
        stage         : 'qc_threshold_approval',
        artifact_type : 'qc_threshold_approval_summary'
    ]

    """
    python3 ${validator_script} \
        --plan ${validated_plan} \
        --policy ${reviewed_policy} \
        --validated-policy approved_qc_threshold_policy.yaml \
        --summary qc_threshold_approval.summary.tsv \
        --session-info qc_threshold_approval.session.txt
    """

    stub:

    output_meta = plan_meta + [
        stage           : 'qc_threshold_approval',
        artifact_type   : 'approved_qc_threshold_policy',
        approval_status : 'approved'
    ]

    summary_meta = plan_meta + [
        stage         : 'qc_threshold_approval',
        artifact_type : 'qc_threshold_approval_summary'
    ]

    """
    cat > approved_qc_threshold_policy.yaml <<'YAML'
    qc_threshold_policy:
      schema_version: "1.0"
      plan_id: "stub_plan"
      project_id: "stub_project"
      biological_run_id: "stub_run"
      approval:
        required: true
        status: "approved"
        reviewed_by: "stub_reviewer"
        approved_at: "2026-01-01T00:00:00+00:00"
      thresholds:
        nFeature_RNA:
          enabled: true
          min: 200
          max: null
      validation:
        status: "valid"
        validated_for_qc_filter: true
        enabled_threshold_count: 1
    YAML

    cat > qc_threshold_approval.summary.tsv <<'TSV'
    plan_id	project_id	biological_run_id	validation_status	approval_status	reviewed_by	approved_at	enabled_threshold_count
    stub_plan	stub_project	stub_run	valid	approved	stub_reviewer	2026-01-01T00:00:00+00:00	1
    TSV

    cat > qc_threshold_approval.session.txt <<'TXT'
    stub=true
    process=VALIDATE_QC_THRESHOLD_APPROVAL
    validation_status=valid
    approval_status=approved
    TXT
    """
}

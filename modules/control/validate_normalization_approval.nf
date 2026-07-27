nextflow.enable.dsl=2


process VALIDATE_NORMALIZATION_APPROVAL {

    tag "${plan_meta.project_id ?: 'project'}:${plan_meta.biological_run_id ?: 'run'}"

    label 'lightweight'

    publishDir(
        "${params.outdir}/control/normalization_policy_approval",
        mode: 'copy',
        overwrite: true
    )

    input:

    tuple val(plan_meta), path(validated_plan)
    tuple val(policy_meta), path(reviewed_policy)
    path validator_script

    output:

    tuple val(output_meta), path("approved_normalization_policy.yaml"),
        emit: approved_policy

    tuple val(summary_meta), path("normalization_policy_approval.summary.tsv"),
        emit: summary

    path "normalization_policy_approval.session.txt",
        emit: session_info

    script:

    if (
        plan_meta.project_id != policy_meta.project_id
        || plan_meta.biological_run_id != policy_meta.biological_run_id
    ) {
        error(
            "VALIDATE_NORMALIZATION_APPROVAL received incompatible " +
            "metadata between the analysis plan and reviewed policy."
        )
    }

    output_meta = plan_meta + [
        stage           : 'normalization_policy_approval',
        artifact_type   : 'approved_normalization_policy',
        approval_status : 'approved'
    ]

    summary_meta = plan_meta + [
        stage         : 'normalization_policy_approval',
        artifact_type : 'normalization_policy_approval_summary'
    ]

    """
    python3 ${validator_script} \
        --plan ${validated_plan} \
        --policy ${reviewed_policy} \
        --validated-policy approved_normalization_policy.yaml \
        --summary normalization_policy_approval.summary.tsv \
        --session-info normalization_policy_approval.session.txt
    """

    stub:

    output_meta = plan_meta + [
        stage           : 'normalization_policy_approval',
        artifact_type   : 'approved_normalization_policy',
        approval_status : 'approved'
    ]

    summary_meta = plan_meta + [
        stage         : 'normalization_policy_approval',
        artifact_type : 'normalization_policy_approval_summary'
    ]

    """
    cat > approved_normalization_policy.yaml <<'YAML'
    normalization_policy:
      schema_version: "1.0"
      status: "proposed"
      identity:
        project_id: "stub_project"
        project_name: "stub_project"
        biological_run_id: "stub_run"
      input:
        cells: 1
        features: 1
        samples: 1
        readiness: "ready"
      strategy:
        method: "LogNormalize"
        variable_features_method: "vst"
        nfeatures: 2000
        scale_factor: 10000
        scope: "per_sample"
        assay: "RNA"
        regressors: "none"
      integration:
        required: false
      approval:
        required: true
        status: "approved"
        reviewed_by: "stub_reviewer"
        approved_at: "2026-01-01T00:00:00+00:00"
      validation:
        status: "valid"
        validated_at: "2026-01-01T00:00:00+00:00"
        validator: "validate_normalization_approval.py"
        validated_for_normalization: true
      application:
        status: "not_applied"
        normalization_applied: false
    YAML

    cat > normalization_policy_approval.summary.tsv <<'TSV'
    plan_id	project_id	biological_run_id	validation_status	approval_status	reviewed_by	approved_at	validated_at	method	scope	integration_required
    stub_plan	stub_project	stub_run	valid	approved	stub_reviewer	2026-01-01T00:00:00+00:00	2026-01-01T00:00:00+00:00	LogNormalize	per_sample	false
    TSV

    cat > normalization_policy_approval.session.txt <<'TXT'
    stub=true
    process=VALIDATE_NORMALIZATION_APPROVAL
    validation_status=valid
    approval_status=approved
    TXT
    """
}

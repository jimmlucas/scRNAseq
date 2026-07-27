nextflow.enable.dsl=2


process VALIDATE_PC_SELECTION_APPROVAL {

    tag "${plan_meta.project_id ?: 'project'}:${plan_meta.biological_run_id ?: 'run'}"

    label 'lightweight'

    publishDir(
        "${params.outdir}/control/pc_selection_policy_approval",
        mode: 'copy',
        overwrite: true
    )

    input:

    tuple val(plan_meta), path(validated_plan)
    tuple val(policy_meta), path(reviewed_policy)
    path validator_script

    output:

    tuple val(output_meta), path("approved_pc_selection_policy.yaml"),
        emit: approved_policy

    tuple val(summary_meta), path("pc_selection_policy_approval.summary.tsv"),
        emit: summary

    path "pc_selection_policy_approval.session.txt",
        emit: session_info

    script:

    if (
        plan_meta.project_id != policy_meta.project_id
        || plan_meta.biological_run_id != policy_meta.biological_run_id
    ) {
        error(
            "VALIDATE_PC_SELECTION_APPROVAL received incompatible " +
            "metadata between the analysis plan and reviewed policy."
        )
    }

    output_meta = plan_meta + [
        stage           : 'pc_selection_policy_approval',
        artifact_type   : 'approved_pc_selection_policy',
        approval_status : 'approved'
    ]

    summary_meta = plan_meta + [
        stage         : 'pc_selection_policy_approval',
        artifact_type : 'pc_selection_policy_approval_summary'
    ]

    """
    python3 ${validator_script} \
        --plan ${validated_plan} \
        --policy ${reviewed_policy} \
        --validated-policy approved_pc_selection_policy.yaml \
        --summary pc_selection_policy_approval.summary.tsv \
        --session-info pc_selection_policy_approval.session.txt
    """

    stub:

    output_meta = plan_meta + [
        stage           : 'pc_selection_policy_approval',
        artifact_type   : 'approved_pc_selection_policy',
        approval_status : 'approved'
    ]

    summary_meta = plan_meta + [
        stage         : 'pc_selection_policy_approval',
        artifact_type : 'pc_selection_policy_approval_summary'
    ]

    """
    cat > approved_pc_selection_policy.yaml <<'YAML'
    pc_selection_policy:
      schema_version: "1.0"
      status: "proposed"
      identity:
        project_id: "stub_project"
        project_name: "stub_project"
        biological_run_id: "stub_run"
      input:
        cells: 1
        features: 1
        variable_features: 1
        max_computed_pcs: 30
      recommendation:
        recommended_pcs: 10
        max_computed_pcs: 30
        method: "cumulative_variance_target_with_minimum_pc_floor"
      proposal:
        recommended_pcs: 10
        max_computed_pcs: 30
        method: "cumulative_variance_target_with_minimum_pc_floor"
      selected_pcs: 10
      approval:
        required: true
        status: "approved"
        reviewed_by: "stub_reviewer"
        approved_at: "2026-01-01T00:00:00+00:00"
      validation:
        status: "valid"
        validated_at: "2026-01-01T00:00:00+00:00"
        validator: "validate_pc_selection_approval.py"
        validated_for_pc_selection: true
        proposed_pcs: 10
        selected_pcs: 10
        max_computed_pcs: 30
      application:
        status: "not_applied"
        pc_selection_applied: false
    YAML

    cat > pc_selection_policy_approval.summary.tsv <<'TSV'
    plan_id	project_id	biological_run_id	validation_status	approval_status	reviewed_by	approved_at	validated_at	proposed_pcs	selected_pcs	max_computed_pcs
    stub_plan	stub_project	stub_run	valid	approved	stub_reviewer	2026-01-01T00:00:00+00:00	2026-01-01T00:00:00+00:00	10	10	30
    TSV

    cat > pc_selection_policy_approval.session.txt <<'TXT'
    stub=true
    process=VALIDATE_PC_SELECTION_APPROVAL
    validation_status=valid
    approval_status=approved
    TXT
    """
}

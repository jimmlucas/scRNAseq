nextflow.enable.dsl=2


process PC_SELECTION_POLICY_RESOLUTION {

    tag "${meta.project_id ?: 'project'}:${meta.biological_run_id ?: 'run'}"

    label 'lightweight'

    publishDir(
        "${params.outdir}/control/pc_selection_policy",
        mode: 'copy',
        overwrite: true
    )

    input:

    val meta
    path dimensionality_summary
    path variable_features
    path elbow_metrics
    path resolver_script

    output:

    tuple val(policy_meta), path("pc_selection_policy.yaml"),
        emit: policy_yaml

    tuple val(policy_meta), path("pc_selection_policy.tsv"),
        emit: policy_tsv

    tuple val(summary_meta), path("pc_selection_policy_resolution.summary.tsv"),
        emit: summary

    path "pc_selection_policy_resolution.session.txt",
        emit: session_info

    script:

    policy_meta = meta + [
        stage          : 'pc_selection_policy_resolution',
        artifact_type  : 'pc_selection_policy',
        approval_status: 'proposed'
    ]

    summary_meta = meta + [
        stage         : 'pc_selection_policy_resolution',
        artifact_type : 'pc_selection_policy_resolution'
    ]

    """
    python3 ${resolver_script} \
        --dimensionality-summary ${dimensionality_summary} \
        --variable-features ${variable_features} \
        --elbow-metrics ${elbow_metrics} \
        --policy-yaml pc_selection_policy.yaml \
        --policy-tsv pc_selection_policy.tsv \
        --resolution-summary pc_selection_policy_resolution.summary.tsv \
        --session-info pc_selection_policy_resolution.session.txt \
        --project-id '${meta.project_id ?: ""}' \
        --project-name '${meta.project_name ?: ""}' \
        --biological-run-id '${meta.biological_run_id ?: ""}'
    """

    stub:

    policy_meta = meta + [
        stage          : 'pc_selection_policy_resolution',
        artifact_type  : 'pc_selection_policy',
        approval_status: 'proposed'
    ]

    summary_meta = meta + [
        stage         : 'pc_selection_policy_resolution',
        artifact_type : 'pc_selection_policy_resolution'
    ]

    """
    cat > pc_selection_policy.yaml <<'YAML'
    pc_selection_policy:
      schema_version: "1.0"
      status: "proposed"
      recommendation:
        recommended_pcs: 10
        max_computed_pcs: 30
        method: "cumulative_variance_target_with_minimum_pc_floor"
      approval:
        required: true
        status: "proposed"
      application:
        status: "not_applied"
        pc_selection_applied: false
    YAML

    cat > pc_selection_policy.tsv <<'TSV'
    project_id	biological_run_id	cells	features	variable_features	recommended_pcs	max_computed_pcs	method	achieved_cumulative_variance	approval_status	application_status
    stub_project	stub_run	0	0	0	10	30	cumulative_variance_target_with_minimum_pc_floor	0	proposed	not_applied
    TSV

    cat > pc_selection_policy_resolution.summary.tsv <<'TSV'
    project_id	biological_run_id	recommended_pcs	max_computed_pcs	method	approval_status
    stub_project	stub_run	10	30	cumulative_variance_target_with_minimum_pc_floor	proposed
    TSV

    cat > pc_selection_policy_resolution.session.txt <<'TXT'
    stub=true
    process=PC_SELECTION_POLICY_RESOLUTION
    approval_status=proposed
    TXT
    """
}

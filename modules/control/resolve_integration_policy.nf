nextflow.enable.dsl=2


process INTEGRATION_POLICY_RESOLUTION {

    tag "${meta.project_id ?: 'project'}:${meta.biological_run_id ?: 'run'}"

    label 'lightweight'

    publishDir(
        "${params.outdir}/control/integration_policy",
        mode: 'copy',
        overwrite: true
    )

    input:

    val meta
    path pca_object
    path dimensionality_summary
    path approved_pc_selection_policy
    path analysis_readiness
    path applied_normalization_policy
    path resolver_script

    output:

    tuple val(policy_meta), path("integration_policy.yaml"),
        emit: policy_yaml

    tuple val(policy_meta), path("integration_policy.tsv"),
        emit: policy_tsv

    tuple val(summary_meta), path("integration_policy_resolution.summary.tsv"),
        emit: summary

    path "integration_policy_resolution.session.txt",
        emit: session_info

    script:

    policy_meta = meta + [
        stage          : 'integration_policy_resolution',
        artifact_type  : 'integration_policy',
        approval_status: 'proposed'
    ]

    summary_meta = meta + [
        stage         : 'integration_policy_resolution',
        artifact_type : 'integration_policy_resolution'
    ]

    """
    python3 ${resolver_script} \
        --pca-object ${pca_object} \
        --dimensionality-summary ${dimensionality_summary} \
        --approved-pc-selection-policy ${approved_pc_selection_policy} \
        --analysis-readiness ${analysis_readiness} \
        --applied-normalization-policy ${applied_normalization_policy} \
        --policy-yaml integration_policy.yaml \
        --policy-tsv integration_policy.tsv \
        --resolution-summary integration_policy_resolution.summary.tsv \
        --session-info integration_policy_resolution.session.txt \
        --project-id '${meta.project_id ?: ""}' \
        --project-name '${meta.project_name ?: ""}' \
        --biological-run-id '${meta.biological_run_id ?: ""}'
    """

    stub:

    policy_meta = meta + [
        stage          : 'integration_policy_resolution',
        artifact_type  : 'integration_policy',
        approval_status: 'proposed'
    ]

    summary_meta = meta + [
        stage         : 'integration_policy_resolution',
        artifact_type : 'integration_policy_resolution'
    ]

    """
    cat > integration_policy.yaml <<'YAML'
    integration_policy:
      schema_version: "1.0"
      status: "proposed"
      strategy:
        integration_required: true
        method: "Harmony"
        input_reduction: "pca"
        output_reduction: "harmony"
        dimensions_to_use: 10
        batch_variable: "sample"
      approval:
        required: true
        status: "proposed"
      application:
        status: "not_applied"
        integration_applied: false
    YAML

    cat > integration_policy.tsv <<'TSV'
    project_id	biological_run_id	samples	selected_pcs	integration_required	method	input_reduction	output_reduction	batch_variable	approval_status	application_status
    stub_project	stub_run	2	10	true	Harmony	pca	harmony	sample	proposed	not_applied
    TSV

    cat > integration_policy_resolution.summary.tsv <<'TSV'
    project_id	biological_run_id	integration_required	method	dimensions_to_use	approval_status
    stub_project	stub_run	true	Harmony	10	proposed
    TSV

    cat > integration_policy_resolution.session.txt <<'TXT'
    stub=true
    process=INTEGRATION_POLICY_RESOLUTION
    approval_status=proposed
    TXT
    """
}

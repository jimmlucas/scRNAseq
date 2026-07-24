nextflow.enable.dsl=2


process DOUBLET_POLICY_RESOLUTION {

    tag "${meta.project_id ?: 'project'}:${meta.biological_run_id ?: 'run'}"

    label 'lightweight'

    publishDir(
        "${params.outdir}/control/doublet_policy",
        mode: 'copy',
        overwrite: true
    )

    input:

    val meta
    path doublet_summary
    path doublet_by_sample
    path doublet_cell_scores
    path resolver_script

    output:

    tuple val(policy_meta), path("doublet_policy.yaml"),
        emit: policy_yaml

    tuple val(policy_meta), path("doublet_policy.tsv"),
        emit: policy_tsv

    tuple val(summary_meta), path("doublet_policy_resolution.summary.tsv"),
        emit: summary

    path "doublet_policy_resolution.session.txt",
        emit: session_info

    script:

    policy_meta = meta + [
        stage          : 'doublet_policy_resolution',
        artifact_type  : 'doublet_detection_policy',
        approval_status: 'proposed'
    ]

    summary_meta = meta + [
        stage         : 'doublet_policy_resolution',
        artifact_type : 'doublet_policy_resolution'
    ]

    """
    python3 ${resolver_script} \
        --summary ${doublet_summary} \
        --by-sample ${doublet_by_sample} \
        --cell-scores ${doublet_cell_scores} \
        --policy-yaml doublet_policy.yaml \
        --policy-tsv doublet_policy.tsv \
        --resolution-summary doublet_policy_resolution.summary.tsv \
        --session-info doublet_policy_resolution.session.txt \
        --project-id '${meta.project_id ?: ""}' \
        --project-name '${meta.project_name ?: ""}' \
        --biological-run-id '${meta.biological_run_id ?: ""}'
    """

    stub:

    policy_meta = meta + [
        stage          : 'doublet_policy_resolution',
        artifact_type  : 'doublet_detection_policy',
        approval_status: 'proposed'
    ]

    summary_meta = meta + [
        stage         : 'doublet_policy_resolution',
        artifact_type : 'doublet_policy_resolution'
    ]

    """
    cat > doublet_policy.yaml <<'YAML'
    doublet_policy:
      schema_version: "1.0"
      approval:
        required: true
        status: "proposed"
        reviewed_by: null
        approved_at: null
      samples: {}
    YAML

    cat > doublet_policy.tsv <<'TSV'
    sample	detection_status	decision_status	proposed_action	cells_scored	proposed_doublet_count	proposed_singlet_count	pK	nExp_adjusted	approval_status
    TSV

    cat > doublet_policy_resolution.summary.tsv <<'TSV'
    project_id	biological_run_id	samples	cells_scored	proposed_doublets	approval_status
    stub_project	stub_run	0	0	0	proposed
    TSV

    cat > doublet_policy_resolution.session.txt <<'TXT'
    stub=true
    process=DOUBLET_POLICY_RESOLUTION
    approval_status=proposed
    TXT
    """
}

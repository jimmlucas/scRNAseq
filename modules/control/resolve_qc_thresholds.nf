nextflow.enable.dsl=2


process QC_THRESHOLD_RESOLUTION {

    tag "${plan_meta.project_id ?: 'project'}:${plan_meta.biological_run_id ?: 'run'}"

    label 'lightweight'

    publishDir(
        "${params.outdir}/control/qc_thresholds",
        mode: 'copy',
        overwrite: true
    )

    input:

    tuple val(plan_meta), path(validated_plan)
    path global_metrics
    path sample_metrics
    path cell_metadata
    path resolver_script

    output:

    tuple val(policy_meta), path("qc_threshold_policy.yaml"),
        emit: policy_yaml

    tuple val(policy_meta), path("qc_threshold_policy.tsv"),
        emit: policy_tsv

    tuple val(summary_meta), path("qc_threshold_resolution.summary.tsv"),
        emit: summary

    path "qc_threshold_resolution.session.txt",
        emit: session_info

    script:

    policy_meta = plan_meta + [
        stage         : 'qc_threshold_resolution',
        artifact_type : 'qc_threshold_policy',
        approval_status: 'proposed'
    ]

    summary_meta = plan_meta + [
        stage         : 'qc_threshold_resolution',
        artifact_type : 'qc_threshold_resolution'
    ]

    """
    python3 ${resolver_script} \
        --plan ${validated_plan} \
        --global-metrics ${global_metrics} \
        --sample-metrics ${sample_metrics} \
        --cell-metadata ${cell_metadata} \
        --policy-yaml qc_threshold_policy.yaml \
        --policy-tsv qc_threshold_policy.tsv \
        --summary qc_threshold_resolution.summary.tsv \
        --session-info qc_threshold_resolution.session.txt
    """

    stub:

    policy_meta = plan_meta + [
        stage         : 'qc_threshold_resolution',
        artifact_type : 'qc_threshold_policy',
        approval_status: 'proposed'
    ]

    summary_meta = plan_meta + [
        stage         : 'qc_threshold_resolution',
        artifact_type : 'qc_threshold_resolution'
    ]

    """
    cat > qc_threshold_policy.yaml <<'YAML'
    qc_threshold_policy:
      schema_version: "1.0"
      plan_id: "stub_plan"
      threshold_mode: "hybrid"
      material: "cells"
      approval:
        required: true
        status: "proposed"
        approved_by: null
        approved_at: null
      thresholds: {}
    YAML

    cat > qc_threshold_policy.tsv <<'TSV'
    metric	enabled	min	max	source	estimated_flagged_cells	estimated_flagged_fraction	approval_status
    TSV

    cat > qc_threshold_resolution.summary.tsv <<'TSV'
    plan_id	project_id	biological_run_id	threshold_mode	material	approval_status	metric_count	cell_count
    stub_plan	stub_project	stub_run	hybrid	cells	proposed	0	0
    TSV

    cat > qc_threshold_resolution.session.txt <<'TXT'
    stub=true
    process=QC_THRESHOLD_RESOLUTION
    approval_status=proposed
    TXT
    """
}

nextflow.enable.dsl=2


process QC_VISUALIZATION_POLICY {

    tag "${before_meta.project_id ?: 'project'}:${before_meta.biological_run_id ?: 'run'}"

    label 'medium_memory'

    publishDir(
        "${params.outdir}/01f_qc_visualization_policy",
        mode: 'copy',
        overwrite: true
    )

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        "docker://${params.r_image}" :
        params.r_image}"

    input:

    tuple val(before_meta), path(before_object)
    tuple val(after_meta), path(after_object)
    path cell_decisions
    tuple val(policy_meta), path(approved_policy)
    path qc_visualization_script

    output:

    tuple val(summary_meta), path("01f_qc_visualization_summary.csv"),
        emit: summary

    tuple val(distribution_meta), path("01f_qc_metric_distributions.csv"),
        emit: metric_distributions

    tuple val(by_sample_meta), path("01f_qc_by_sample_before_after.csv"),
        emit: by_sample

    tuple val(retention_meta), path("01f_qc_sample_retention.csv"),
        emit: sample_retention

    tuple val(exclusion_meta), path("01f_qc_exclusion_reasons.csv"),
        emit: exclusion_reasons

    tuple val(cell_metadata_meta), path("01f_qc_cell_metadata_before_after.csv"),
        emit: cell_metadata

    tuple val(plots_meta), path("01f_qc_*.png"),
        emit: plots

    tuple val(session_meta), path("01f_sessionInfo.txt"),
        emit: session_info

    script:

    if (
        before_meta.project_id != after_meta.project_id
        || before_meta.project_id != policy_meta.project_id
        || before_meta.biological_run_id != after_meta.biological_run_id
        || before_meta.biological_run_id != policy_meta.biological_run_id
    ) {
        error(
            "QC_VISUALIZATION_POLICY received incompatible metadata."
        )
    }

    output_base_meta = before_meta + [
        stage                : 'qc_visualization',
        source_stage         : after_meta.stage ?: '',
        source_artifact_type : after_meta.artifact_type ?: ''
    ]

    summary_meta = output_base_meta + [
        artifact_type: 'qc_visualization_summary'
    ]

    distribution_meta = output_base_meta + [
        artifact_type: 'qc_metric_distributions'
    ]

    by_sample_meta = output_base_meta + [
        artifact_type: 'per_sample_metrics'
    ]

    retention_meta = output_base_meta + [
        artifact_type: 'qc_sample_retention'
    ]

    exclusion_meta = output_base_meta + [
        artifact_type: 'qc_exclusion_reasons'
    ]

    cell_metadata_meta = output_base_meta + [
        artifact_type: 'cell_metadata'
    ]

    plots_meta = output_base_meta + [
        artifact_type: 'qc_visualization_report'
    ]

    session_meta = output_base_meta + [
        artifact_type: 'session_info'
    ]

    """
    Rscript ${qc_visualization_script} \
        --before-object ${before_object} \
        --after-object ${after_object} \
        --cell-decisions ${cell_decisions} \
        --policy ${approved_policy}
    """

    stub:

    output_base_meta = before_meta + [
        stage                : 'qc_visualization',
        source_stage         : after_meta.stage ?: '',
        source_artifact_type : after_meta.artifact_type ?: ''
    ]

    summary_meta = output_base_meta + [
        artifact_type: 'qc_visualization_summary'
    ]

    distribution_meta = output_base_meta + [
        artifact_type: 'qc_metric_distributions'
    ]

    by_sample_meta = output_base_meta + [
        artifact_type: 'per_sample_metrics'
    ]

    retention_meta = output_base_meta + [
        artifact_type: 'qc_sample_retention'
    ]

    exclusion_meta = output_base_meta + [
        artifact_type: 'qc_exclusion_reasons'
    ]

    cell_metadata_meta = output_base_meta + [
        artifact_type: 'cell_metadata'
    ]

    plots_meta = output_base_meta + [
        artifact_type: 'qc_visualization_report'
    ]

    session_meta = output_base_meta + [
        artifact_type: 'session_info'
    ]

    """
    cat > 01f_qc_visualization_summary.csv <<'CSV'
    metric,value
    cells_before,0
    cells_after,0
    cells_removed,0
    fraction_removed,0
    features_before,0
    features_after,0
    approved_threshold_count,0
    CSV

    cat > 01f_qc_metric_distributions.csv <<'CSV'
    qc_stage,metric,cells,minimum,q25,median,mean,q75,maximum
    CSV

    cat > 01f_qc_by_sample_before_after.csv <<'CSV'
    qc_stage,sample,Group,cells
    CSV

    cat > 01f_qc_sample_retention.csv <<'CSV'
    sample,Group,cells_before,cells_after,cells_removed,fraction_retained
    CSV

    cat > 01f_qc_exclusion_reasons.csv <<'CSV'
    qc_exclusion_reason,cells
    CSV

    cat > 01f_qc_cell_metadata_before_after.csv <<'CSV'
    cell_id,qc_stage
    CSV

    touch 01f_qc_scatter_after_filtering.png

    cat > 01f_sessionInfo.txt <<'TXT'
    stub=true
    process=QC_VISUALIZATION_POLICY
    TXT
    """
}

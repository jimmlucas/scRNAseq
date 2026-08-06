nextflow.enable.dsl=2


process QC_FILTER_POLICY {

    tag "${object_meta.project_id ?: 'project'}:${object_meta.biological_run_id ?: 'run'}"

    label 'high_memory'

    publishDir(
        "${params.outdir}/01e_qc_filter_policy",
        mode: 'copy',
        overwrite: true
    )

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        "docker://${params.r_image}" :
        params.r_image}"

    input:

    tuple val(object_meta), path(qc_metrics_object)
    tuple val(policy_meta), path(approved_policy)
    path qc_filter_script

    output:

    tuple val(output_meta), path("seu_qc_filtered.rds"),
        emit: filtered_object

    tuple val(summary_meta), path("01e_qc_filter_summary.csv"),
        emit: summary

    path "01e_qc_filter_by_threshold.csv",
        emit: by_threshold

    path "01e_qc_filter_by_sample.csv",
        emit: by_sample

    path "01e_qc_filter_cell_decisions.csv",
        emit: cell_decisions

    path "01e_sessionInfo.txt",
        emit: session_info

    script:

    if (
        object_meta.project_id != policy_meta.project_id
        || object_meta.biological_run_id != policy_meta.biological_run_id
    ) {
        error(
            "QC_FILTER_POLICY received incompatible metadata between " +
            "QC metrics object and approved threshold policy."
        )
    }

    if (
        policy_meta.approval_status != 'approved'
        && policy_meta.approval_status != 'modified_and_approved'
    ) {
        error(
            "QC_FILTER_POLICY requires an approved threshold policy."
        )
    }

    output_meta = object_meta + [
        stage         : 'qc_filter',
        artifact_type : 'qc_filtered_seurat',
        policy_status : policy_meta.approval_status
    ]

    summary_meta = object_meta + [
        stage         : 'qc_filter',
        artifact_type : 'metrics'
    ]

    """
    Rscript ${qc_filter_script} \
        --input ${qc_metrics_object} \
        --policy ${approved_policy} \
        --output seu_qc_filtered.rds
    """

    stub:

    output_meta = object_meta + [
        stage         : 'qc_filter',
        artifact_type : 'qc_filtered_seurat',
        policy_status : 'approved'
    ]

    summary_meta = object_meta + [
        stage         : 'qc_filter',
        artifact_type : 'metrics'
    ]

    """
    touch seu_qc_filtered.rds

    cat > 01e_qc_filter_summary.csv <<'CSV'
    metric,value
    cells_before_filtering,0
    cells_after_filtering,0
    cells_removed,0
    fraction_removed,0
    features_before_filtering,0
    features_after_filtering,0
    enabled_threshold_count,0
    CSV

    cat > 01e_qc_filter_by_threshold.csv <<'CSV'
    metric,min,max,cells_below_min,cells_above_max,cells_flagged_for_metric
    CSV

    cat > 01e_qc_filter_by_sample.csv <<'CSV'
    sample,cells_before,cells_after,cells_removed,fraction_removed
    CSV

    cat > 01e_qc_filter_cell_decisions.csv <<'CSV'
    cell_id,qc_pass,qc_fail,qc_exclusion_reason
    CSV

    cat > 01e_sessionInfo.txt <<'TXT'
    stub=true
    process=QC_FILTER_POLICY
    TXT
    """
}

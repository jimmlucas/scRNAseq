nextflow.enable.dsl=2


process QC_METRICS {

    tag "${meta.project_id ?: 'project'}:${meta.biological_run_id ?: 'run'}"

    label 'high_memory'

    publishDir(
        "${params.outdir}/01a_qc_metrics",
        mode: 'copy',
        overwrite: true
    )

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        "docker://${params.r_image}" :
        params.r_image}"

    input:

    tuple val(meta), path(seurat_object)
    tuple val(config_meta), path(config_file)
    path qc_metrics_script

    output:

    tuple val(output_meta), path("seu_qc_metrics.rds"),
        emit: metrics_object

    path "01a_qc_metrics_summary.csv",
        emit: summary

    path "01a_qc_metrics_by_sample.csv",
        emit: by_sample

    path "01a_qc_metrics_cell_metadata.csv",
        emit: cell_metadata

    path "01a_sessionInfo.txt",
        emit: session_info

    script:

    if (
        meta.project_id != config_meta.project_id
        || meta.biological_run_id != config_meta.biological_run_id
    ) {
        error(
            "QC_METRICS received incompatible metadata between " +
            "Seurat object and configuration."
        )
    }

    output_meta = meta + [
        stage         : 'qc_metrics',
        artifact_type : 'qc_metrics_seurat'
    ]

    """
    Rscript ${qc_metrics_script} \
        --input ${seurat_object} \
        --config ${config_file} \
        --output seu_qc_metrics.rds
    """

    stub:

    output_meta = meta + [
        stage         : 'qc_metrics',
        artifact_type : 'qc_metrics_seurat'
    ]

    """
    touch seu_qc_metrics.rds

    cat > 01a_qc_metrics_summary.csv <<'CSV'
    metric,value
    cells,0
    rna_features,0
    CSV

    cat > 01a_qc_metrics_by_sample.csv <<'CSV'
    sample,Group,cells
    stub,stub,0
    CSV

    cat > 01a_qc_metrics_cell_metadata.csv <<'CSV'
    cell_id,sample,Group
    CSV

    cat > 01a_sessionInfo.txt <<'TXT'
    stub=true
    process=QC_METRICS
    TXT
    """
}

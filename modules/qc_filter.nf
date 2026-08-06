process QC_FILTER {

    tag "qc_filter"
    label "high_memory"

    publishDir "${params.outdir}/01_qc_filter", mode: 'copy'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        "docker://${params.r_image}" :
        params.r_image}"

    input:
    tuple val(meta), path(seurat_object)
    path config_file
    path qc_filter_script

    output:
    tuple val(output_meta), path("seu_qc_filtered.rds"), emit: filtered_object
    path "01_qc_summary.csv", emit: summary
    path "01_sessionInfo.txt", emit: session_info

    script:
    output_meta = meta + [
        stage         : 'qc_filter',
        artifact_type : 'filtered_seurat'
    ]

    """
    Rscript ${qc_filter_script} \
        --input ${seurat_object} \
        --config ${config_file} \
        --output seu_qc_filtered.rds
    """
}

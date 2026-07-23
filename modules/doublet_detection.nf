process QC_DOUBLET_DETECTION {

    tag "doublet_detection"
    label "very_high_memory"

    publishDir "${params.outdir}/01c_doublet_detection",
        mode: 'copy',
        overwrite: true

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.r_image.docker}" :
        params.r_image}"

    input:
    tuple val(meta), path(seurat_object)
    path config_file
    path doublet_detection_script

    output:
    tuple val(singlet_meta), path("seu_singlets.rds"), emit: singlet_object
    tuple val(summary_meta), path("01c_doublet_summary.csv"), emit: summary
    tuple val(by_sample_meta), path("01c_doublet_by_sample.csv"), emit: by_sample
    tuple val(metadata_meta), path("01c_doublet_metadata.csv"), emit: metadata
    tuple val(session_info_meta), path("01c_sessionInfo.txt"), emit: session_info

    script:

    /*
     * Preserve the scientific identity of the filtered Seurat object
     * and record the provenance of every generated artifact.
     */
    def output_base_meta = meta + [
        stage                : "doublet_detection",
        source_stage         : meta.stage ?: "",
        source_artifact_type : meta.artifact_type ?: ""
    ]

    singlet_meta = output_base_meta + [
        artifact_type: "singlet_seurat"
    ]

    summary_meta = output_base_meta + [
        artifact_type: "doublet_summary_report"
    ]

    by_sample_meta = output_base_meta + [
        artifact_type: "doublet_by_sample_report"
    ]

    metadata_meta = output_base_meta + [
        artifact_type: "doublet_cell_metadata"
    ]

    session_info_meta = output_base_meta + [
        artifact_type: "session_info"
    ]

    """
    Rscript ${doublet_detection_script} \
        --input ${seurat_object} \
        --config ${config_file} \
        --output seu_singlets.rds
    """
}
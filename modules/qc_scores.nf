process QC_SCORES {

    tag "qc_scores"
    label "high_memory"

    publishDir "${params.outdir}/01d_qc_scores",
        mode: 'copy',
        overwrite: true

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        "docker://${params.r_image}" :
        params.r_image}"

    input:
    tuple val(meta), path(seurat_object)
    path config_file
    path qc_scores_script

    output:
    tuple val(scored_meta), path("seu_qc_scored.rds"), emit: scored_object
    tuple val(summary_meta), path("01d_qc_scores_summary.csv"), emit: summary
    tuple val(phase_summary_meta), path("01d_cell_cycle_phase_summary.csv"), emit: phase_summary
    tuple val(metadata_meta), path("01d_qc_scores_metadata.csv"), emit: metadata
    tuple val(session_info_meta), path("01d_sessionInfo.txt"), emit: session_info

    script:

    /*
     * Preserve the scientific identity and provenance of the singlet
     * Seurat object while assigning a contract to each output artifact.
     */
    def output_base_meta = meta + [
        stage                : "qc_scores",
        source_stage         : meta.stage ?: "",
        source_artifact_type : meta.artifact_type ?: ""
    ]

    scored_meta = output_base_meta + [
        artifact_type: "scored_seurat"
    ]

    summary_meta = output_base_meta + [
        artifact_type: "qc_scores_summary_report"
    ]

    phase_summary_meta = output_base_meta + [
        artifact_type: "cell_cycle_phase_summary_report"
    ]

    metadata_meta = output_base_meta + [
        artifact_type: "qc_scores_cell_metadata"
    ]

    session_info_meta = output_base_meta + [
        artifact_type: "session_info"
    ]

    """
    Rscript ${qc_scores_script} \
        --input ${seurat_object} \
        --config ${config_file} \
        --output seu_qc_scored.rds
    """
}

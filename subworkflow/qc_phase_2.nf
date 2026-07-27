nextflow.enable.dsl = 2

include { QC_VISUALIZATION }     from '../modules/qc_visualization'
include { QC_DOUBLET_DETECTION } from '../modules/doublet_detection'
include { QC_SCORES }            from '../modules/qc_scores'


workflow QC_PHASE_2 {

    take:
    scientific_meta_ch
    filtered_object_ch
    config_ch

    main:

    qc_visualization_script_ch = channel.value(
        file(
            "${projectDir}/scripts/01b_qc_visualization.R",
            checkIfExists: true
        )
    )

    doublet_detection_script_ch = channel.value(
        file(
            "${projectDir}/scripts/01c_doublet_detection.R",
            checkIfExists: true
        )
    )

    qc_scores_script_ch = channel.value(
        file(
            "${projectDir}/scripts/01d_qc_scores.R",
            checkIfExists: true
        )
    )

    /*
     * Standalone checkpoint reconstruction.
     *
     * qc_phase_2 starts from a previously published QC_FILTER artifact.
     */
    filtered_artifact_ch =
        scientific_meta_ch
            .combine(filtered_object_ch)
            .map { meta, filtered_object ->

                def checkpoint_meta = meta + [
                    stage        : "qc_filter",
                    artifact_type: "filtered_seurat"
                ]

                tuple(checkpoint_meta, filtered_object)
            }

    /*
     * Both branches consume the same scientific artifact:
     *
     * tuple(meta, filtered_seurat.rds)
     */
    QC_VISUALIZATION(
        filtered_artifact_ch,
        config_ch,
        qc_visualization_script_ch
    )

    QC_DOUBLET_DETECTION(
        filtered_artifact_ch,
        config_ch,
        doublet_detection_script_ch
    )

    /*
     * Direct contract-preserving connection:
     *
     * tuple(meta, singlet_seurat.rds)
     */
    QC_SCORES(
        QC_DOUBLET_DETECTION.out.singlet_object,
        config_ch,
        qc_scores_script_ch
    )

    emit:

    visual_summary      = QC_VISUALIZATION.out.summary
    visual_by_sample    = QC_VISUALIZATION.out.by_sample
    visual_by_group     = QC_VISUALIZATION.out.by_group
    visual_metadata     = QC_VISUALIZATION.out.cell_metadata
    visual_plots        = QC_VISUALIZATION.out.plots
    visual_session_info = QC_VISUALIZATION.out.session_info

    singlet_object       = QC_DOUBLET_DETECTION.out.singlet_object
    doublet_summary      = QC_DOUBLET_DETECTION.out.summary
    doublet_by_sample    = QC_DOUBLET_DETECTION.out.by_sample
    doublet_metadata     = QC_DOUBLET_DETECTION.out.metadata
    doublet_session_info = QC_DOUBLET_DETECTION.out.session_info

    scored_object       = QC_SCORES.out.scored_object
    score_summary       = QC_SCORES.out.summary
    score_phase_summary = QC_SCORES.out.phase_summary
    score_metadata      = QC_SCORES.out.metadata
    score_session_info  = QC_SCORES.out.session_info
}
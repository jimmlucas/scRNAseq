nextflow.enable.dsl = 2

include { LOAD_DATA }            from '../modules/load'
include { QC_FILTER }            from '../modules/qc_filter'
include { QC_VISUALIZATION }     from '../modules/qc_visualization'
include { QC_DOUBLET_DETECTION } from '../modules/doublet_detection'
include { QC_SCORES }            from '../modules/qc_scores'


workflow qc_process {

    take:
    scientific_meta_ch
    config_ch

    main:

    /*
     * Scripts are explicit process inputs.
     *
     * Any scientific script change therefore invalidates only the
     * relevant process cache.
     */
    load_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/00-Loading_data.R",
            checkIfExists: true
        )
    )

    qc_filter_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/01_qc_filter.R",
            checkIfExists: true
        )
    )

    qc_visualization_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/01b_qc_visualization.R",
            checkIfExists: true
        )
    )

    doublet_detection_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/01c_doublet_detection.R",
            checkIfExists: true
        )
    )

    qc_scores_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/01d_qc_scores.R",
            checkIfExists: true
        )
    )

    /*
     * Complete contract-preserving QC workflow.
     *
     * Every scientific artifact travels as:
     *
     * tuple(meta, path)
     */
    LOAD_DATA(
        scientific_meta_ch,
        config_ch,
        load_script_ch
    )

    QC_FILTER(
        LOAD_DATA.out.seurat_object,
        config_ch,
        qc_filter_script_ch
    )

    /*
     * Both branches consume the same filtered scientific artifact:
     *
     * tuple(meta, filtered_seurat.rds)
     */
    QC_VISUALIZATION(
        QC_FILTER.out.filtered_object,
        config_ch,
        qc_visualization_script_ch
    )

    QC_DOUBLET_DETECTION(
        QC_FILTER.out.filtered_object,
        config_ch,
        doublet_detection_script_ch
    )

    /*
     * Direct contract-preserving transition:
     *
     * tuple(meta, singlet_seurat.rds)
     *     ->
     * tuple(meta, scored_seurat.rds)
     */
    QC_SCORES(
        QC_DOUBLET_DETECTION.out.singlet_object,
        config_ch,
        qc_scores_script_ch
    )

    emit:

    loaded_object     = LOAD_DATA.out.seurat_object
    load_summary      = LOAD_DATA.out.summary
    load_session_info = LOAD_DATA.out.session_info

    filtered_object = QC_FILTER.out.filtered_object
    qc_summary      = QC_FILTER.out.summary
    qc_session_info = QC_FILTER.out.session_info

    qc_visual_summary      = QC_VISUALIZATION.out.summary
    qc_visual_by_sample    = QC_VISUALIZATION.out.by_sample
    qc_visual_by_group     = QC_VISUALIZATION.out.by_group
    qc_visual_metadata     = QC_VISUALIZATION.out.cell_metadata
    qc_visual_plots        = QC_VISUALIZATION.out.plots
    qc_visual_session_info = QC_VISUALIZATION.out.session_info

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

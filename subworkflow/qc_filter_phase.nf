nextflow.enable.dsl=2

include { QC_FILTER } from '../modules/qc_filter'


workflow QC_FILTER_PHASE {

    take:
    scientific_meta_ch
    loaded_object_ch
    qc_filter_config_ch

    main:

    qc_filter_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/01_qc_filter.R",
            checkIfExists: true
        )
    )

    /*
     * Reconstruct the scientific artifact contract when starting
     * directly from a previously published loaded Seurat object.
     */
    loaded_artifact_ch = scientific_meta_ch
        .combine(loaded_object_ch)
        .map { meta, seurat_object ->
            def checkpoint_meta = meta + [
                stage         : 'load_data',
                artifact_type : 'loaded_seurat'
            ]

            tuple(checkpoint_meta, seurat_object)
        }

    QC_FILTER(
        loaded_artifact_ch,
        qc_filter_config_ch,
        qc_filter_script_ch
    )

    emit:

    filtered_object = QC_FILTER.out.filtered_object
    summary = QC_FILTER.out.summary
    session_info = QC_FILTER.out.session_info
}
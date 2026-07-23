nextflow.enable.dsl=2

include { LOAD_DATA } from '../modules/load'
include { QC_FILTER } from '../modules/qc_filter'


workflow QC_PHASE_1 {

    take:
    scientific_meta_ch
    load_config_ch
    qc_filter_config_ch

    main:

    load_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/00-Loading_data.R",
            checkIfExists: true
        )
    )

    load_config_artifact_ch = scientific_meta_ch
        .combine(load_config_ch)
        .map { meta, config_file ->
            tuple(
                meta + [
                    stage       : 'load_config',
                    artifact    : 'resolved_load_config'
                ],
                config_file
            )
        }

    LOAD_DATA(
        load_config_artifact_ch,
        load_script_ch
    )

    qc_filter_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/01_qc_filter.R",
            checkIfExists: true
        )
    )

    QC_FILTER(
        LOAD_DATA.out.seurat_object,
        qc_filter_config_ch,
        qc_filter_script_ch
    )

    emit:

    loaded_object = LOAD_DATA.out.seurat_object
    load_summary = LOAD_DATA.out.summary
    load_session_info = LOAD_DATA.out.session_info

    filtered_object = QC_FILTER.out.filtered_object
    qc_summary = QC_FILTER.out.summary
    qc_session_info = QC_FILTER.out.session_info
}
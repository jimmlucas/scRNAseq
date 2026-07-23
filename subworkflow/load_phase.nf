nextflow.enable.dsl=2

include { LOAD_DATA } from '../modules/load'


workflow LOAD_PHASE {

    take:
    scientific_meta_ch
    load_config_ch

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

    emit:

    seurat_object = LOAD_DATA.out.seurat_object
    summary = LOAD_DATA.out.summary
    session_info = LOAD_DATA.out.session_info
}
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

    LOAD_DATA(
        scientific_meta_ch,
        load_config_ch,
        load_script_ch
    )

    emit:

    seurat_object = LOAD_DATA.out.seurat_object
    summary = LOAD_DATA.out.summary
    session_info = LOAD_DATA.out.session_info
}
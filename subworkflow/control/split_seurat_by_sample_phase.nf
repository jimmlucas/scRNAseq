nextflow.enable.dsl=2

include { SPLIT_SEURAT_BY_SAMPLE } from '../../modules/split_seurat_by_sample'


workflow SPLIT_SEURAT_BY_SAMPLE_PHASE {

    take:

    scientific_meta_ch
    input_object_ch

    main:

    split_script_ch = channel.value(
        file(
            "${projectDir}/scripts/01f_split_seurat_by_sample.R",
            checkIfExists: true
        )
    )

    input_artifact_ch = scientific_meta_ch
        .combine(input_object_ch)
        .map { meta, input_object ->
            tuple(
                meta + [
                    stage         : 'qc_filter',
                    artifact_type : 'qc_filtered_seurat'
                ],
                input_object
            )
        }

    SPLIT_SEURAT_BY_SAMPLE(
        input_artifact_ch,
        split_script_ch
    )

    emit:

    sample_objects =
        SPLIT_SEURAT_BY_SAMPLE.out.sample_objects

    manifest =
        SPLIT_SEURAT_BY_SAMPLE.out.manifest

    summary =
        SPLIT_SEURAT_BY_SAMPLE.out.summary

    session_info =
        SPLIT_SEURAT_BY_SAMPLE.out.session_info
}

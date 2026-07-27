nextflow.enable.dsl=2

include { DOUBLET_DETECTION_SCORING } from '../../modules/doublet_detection_scoring'


workflow DOUBLET_DETECTION_SCORING_PHASE {

    take:

    scientific_meta_ch
    filtered_object_ch
    config_ch

    main:

    doublet_detection_script_ch = channel.value(
        file(
            "${projectDir}/scripts/01g_doublet_detection_scoring.R",
            checkIfExists: true
        )
    )

    filtered_artifact_ch = scientific_meta_ch
        .combine(filtered_object_ch)
        .map { meta, filtered_object ->
            tuple(
                meta + [
                    stage         : 'qc_filter',
                    artifact_type : 'qc_filtered_seurat'
                ],
                filtered_object
            )
        }

    config_artifact_ch = scientific_meta_ch
        .combine(config_ch)
        .map { meta, config_file ->
            tuple(
                meta + [
                    stage         : 'doublet_detection_config',
                    artifact_type : 'runtime_config'
                ],
                config_file
            )
        }

    DOUBLET_DETECTION_SCORING(
        filtered_artifact_ch,
        config_artifact_ch,
        doublet_detection_script_ch
    )

    emit:

    scored_object =
        DOUBLET_DETECTION_SCORING.out.scored_object

    summary =
        DOUBLET_DETECTION_SCORING.out.summary

    by_sample =
        DOUBLET_DETECTION_SCORING.out.by_sample

    cell_scores =
        DOUBLET_DETECTION_SCORING.out.cell_scores

    pk_sweep =
        DOUBLET_DETECTION_SCORING.out.pk_sweep

    session_info =
        DOUBLET_DETECTION_SCORING.out.session_info
}

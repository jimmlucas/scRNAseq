nextflow.enable.dsl=2

include { DOUBLET_DETECTION_SCORE_SAMPLE } from '../../modules/doublet_detection_score_sample'
include { MERGE_DOUBLET_DETECTION_SCORING } from '../../modules/merge_doublet_detection_scoring'


workflow DOUBLET_DETECTION_SCORING_BY_SAMPLE_PHASE {

    take:

    scientific_meta_ch
    sample_manifest_ch
    sample_objects_ch
    config_ch

    main:

    doublet_detection_script_ch = channel.value(
        file(
            "${projectDir}/scripts/01g_doublet_detection_scoring.R",
            checkIfExists: true
        )
    )

    merge_script_ch = channel.value(
        file(
            "${projectDir}/scripts/01g_merge_doublet_detection_scoring.R",
            checkIfExists: true
        )
    )

    sample_input_ch = sample_manifest_ch
        .splitCsv(
            header: true,
            sep: '\t'
        )
        .combine(sample_objects_ch)
        .combine(scientific_meta_ch)
        .map { row, sample_objects_root, meta ->
            relative_object = row.object_path.toString()
            relative_object = relative_object.replaceFirst(
                /^sample_objects\//,
                ''
            )

            sample_safe_id = row.sample_safe_id.toString()

            tuple(
                meta + [
                    stage          : 'split_by_sample',
                    artifact_type  : 'sample_seurat_object',
                    sample         : row.sample.toString(),
                    sample_safe_id : sample_safe_id
                ],
                file(
                    "${sample_objects_root}/${relative_object}",
                    checkIfExists: true
                )
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

    sample_with_config_ch = sample_input_ch
        .combine(config_artifact_ch)
        .map { sample_meta, sample_object, config_meta, config_file ->
            tuple(
                sample_meta,
                sample_object,
                config_meta,
                config_file
            )
        }

    DOUBLET_DETECTION_SCORE_SAMPLE(
        sample_with_config_ch,
        doublet_detection_script_ch
    )

    MERGE_DOUBLET_DETECTION_SCORING(
        scientific_meta_ch,
        DOUBLET_DETECTION_SCORE_SAMPLE.out.scored_object.collect(),
        DOUBLET_DETECTION_SCORE_SAMPLE.out.by_sample.collect(),
        DOUBLET_DETECTION_SCORE_SAMPLE.out.cell_scores.collect(),
        DOUBLET_DETECTION_SCORE_SAMPLE.out.pk_sweep.collect(),
        DOUBLET_DETECTION_SCORE_SAMPLE.out.performance.collect(),
        merge_script_ch
    )

    emit:

    scored_object =
        MERGE_DOUBLET_DETECTION_SCORING.out.scored_object

    summary =
        MERGE_DOUBLET_DETECTION_SCORING.out.summary

    by_sample =
        MERGE_DOUBLET_DETECTION_SCORING.out.by_sample

    cell_scores =
        MERGE_DOUBLET_DETECTION_SCORING.out.cell_scores

    pk_sweep =
        MERGE_DOUBLET_DETECTION_SCORING.out.pk_sweep

    performance =
        MERGE_DOUBLET_DETECTION_SCORING.out.performance

    session_info =
        MERGE_DOUBLET_DETECTION_SCORING.out.session_info
}

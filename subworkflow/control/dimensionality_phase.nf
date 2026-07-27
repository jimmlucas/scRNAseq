nextflow.enable.dsl=2

include { DIMENSIONALITY } from '../../modules/dimensionality'


workflow DIMENSIONALITY_PHASE {

    take:

    scientific_meta_ch
    normalized_object_ch
    applied_policy_ch

    main:

    dimensionality_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/05_dimensionality.R",
            checkIfExists: true
        )
    )

    normalized_object_artifact_ch = scientific_meta_ch
        .combine(normalized_object_ch)
        .map { meta, normalized_object ->
            tuple(
                meta + [
                    stage         : 'normalization',
                    artifact_type : 'normalized_seurat'
                ],
                normalized_object
            )
        }

    applied_policy_artifact_ch = scientific_meta_ch
        .combine(applied_policy_ch)
        .map { meta, applied_policy ->
            tuple(
                meta + [
                    stage              : 'normalization',
                    artifact_type      : 'applied_normalization_policy',
                    application_status : 'applied'
                ],
                applied_policy
            )
        }

    DIMENSIONALITY(
        normalized_object_artifact_ch,
        applied_policy_artifact_ch,
        dimensionality_script_ch
    )

    emit:

    pca_object =
        DIMENSIONALITY.out.pca_object

    summary =
        DIMENSIONALITY.out.summary

    variable_features_summary =
        DIMENSIONALITY.out.variable_features_summary

    embeddings_preview =
        DIMENSIONALITY.out.embeddings_preview

    loadings_preview =
        DIMENSIONALITY.out.loadings_preview

    elbow_metrics =
        DIMENSIONALITY.out.elbow_metrics

    session_info =
        DIMENSIONALITY.out.session_info
}

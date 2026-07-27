nextflow.enable.dsl=2

include { NORMALIZATION } from '../../modules/normalization'


workflow NORMALIZATION_PHASE {

    take:

    scientific_meta_ch
    filtered_object_ch
    approved_policy_ch

    main:

    normalization_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/04_normalization.R",
            checkIfExists: true
        )
    )

    filtered_object_artifact_ch = scientific_meta_ch
        .combine(filtered_object_ch)
        .map { meta, filtered_object ->
            tuple(
                meta + [
                    stage         : 'doublet_filter',
                    artifact_type : 'singlet_seurat'
                ],
                filtered_object
            )
        }

    approved_policy_artifact_ch = scientific_meta_ch
        .combine(approved_policy_ch)
        .map { meta, approved_policy ->
            tuple(
                meta + [
                    stage           : 'normalization_policy_approval',
                    artifact_type   : 'approved_normalization_policy',
                    approval_status : 'approved'
                ],
                approved_policy
            )
        }

    NORMALIZATION(
        filtered_object_artifact_ch,
        approved_policy_artifact_ch,
        normalization_script_ch
    )

    emit:

    normalized_object =
        NORMALIZATION.out.normalized_object

    sample_objects =
        NORMALIZATION.out.sample_objects

    summary =
        NORMALIZATION.out.summary

    by_sample =
        NORMALIZATION.out.by_sample

    variable_features =
        NORMALIZATION.out.variable_features

    applied_policy =
        NORMALIZATION.out.applied_policy

    session_info =
        NORMALIZATION.out.session_info
}

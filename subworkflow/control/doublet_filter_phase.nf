nextflow.enable.dsl=2

include { DOUBLET_FILTER } from '../../modules/doublet_filter'


workflow DOUBLET_FILTER_PHASE {

    take:

    scientific_meta_ch
    scored_object_ch
    approved_policy_ch

    main:

    doublet_filter_script_ch = channel.value(
        file(
            "${projectDir}/scripts/02_doublet_filter.R",
            checkIfExists: true
        )
    )

    scored_object_artifact_ch = scientific_meta_ch
        .combine(scored_object_ch)
        .map { meta, scored_object ->
            tuple(
                meta + [
                    stage         : 'doublet_detection',
                    artifact_type : 'doublet_scored_seurat'
                ],
                scored_object
            )
        }

    approved_policy_artifact_ch = scientific_meta_ch
        .combine(approved_policy_ch)
        .map { meta, approved_policy ->
            tuple(
                meta + [
                    stage           : 'doublet_policy_approval',
                    artifact_type   : 'approved_doublet_policy',
                    approval_status : 'approved'
                ],
                approved_policy
            )
        }

    DOUBLET_FILTER(
        scored_object_artifact_ch,
        approved_policy_artifact_ch,
        doublet_filter_script_ch
    )

    emit:

    filtered_object =
        DOUBLET_FILTER.out.filtered_object

    summary =
        DOUBLET_FILTER.out.summary

    by_sample =
        DOUBLET_FILTER.out.by_sample

    removed_cells =
        DOUBLET_FILTER.out.removed_cells

    applied_policy =
        DOUBLET_FILTER.out.applied_policy

    session_info =
        DOUBLET_FILTER.out.session_info
}

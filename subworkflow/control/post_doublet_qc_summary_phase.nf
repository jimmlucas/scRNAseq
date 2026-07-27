nextflow.enable.dsl=2

include { POST_DOUBLET_QC_SUMMARY } from '../../modules/post_doublet_qc_summary'


workflow POST_DOUBLET_QC_SUMMARY_PHASE {

    take:

    scientific_meta_ch
    filtered_object_ch
    applied_policy_ch

    main:

    summary_script_ch = channel.value(
        file(
            "${projectDir}/scripts/03_post_doublet_qc_summary.R",
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

    applied_policy_artifact_ch = scientific_meta_ch
        .combine(applied_policy_ch)
        .map { meta, applied_policy ->
            tuple(
                meta + [
                    stage         : 'doublet_filter',
                    artifact_type : 'applied_doublet_policy',
                    policy_status : 'applied'
                ],
                applied_policy
            )
        }

    POST_DOUBLET_QC_SUMMARY(
        filtered_object_artifact_ch,
        applied_policy_artifact_ch,
        summary_script_ch
    )

    emit:

    summary =
        POST_DOUBLET_QC_SUMMARY.out.summary

    by_sample =
        POST_DOUBLET_QC_SUMMARY.out.by_sample

    warnings =
        POST_DOUBLET_QC_SUMMARY.out.warnings

    readiness =
        POST_DOUBLET_QC_SUMMARY.out.readiness

    session_info =
        POST_DOUBLET_QC_SUMMARY.out.session_info
}

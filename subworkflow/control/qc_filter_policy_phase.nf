nextflow.enable.dsl=2

include { QC_FILTER_POLICY } from '../../modules/qc_filter_policy'


workflow QC_FILTER_POLICY_PHASE {

    take:

    scientific_meta_ch
    qc_metrics_object_ch
    approved_policy_ch

    main:

    qc_filter_script_ch = channel.value(
        file(
            "${projectDir}/scripts/01e_qc_filter_policy.R",
            checkIfExists: true
        )
    )

    /*
     * Reconstruct artifact contracts when starting from
     * previously published checkpoint files.
     */

    qc_metrics_artifact_ch = scientific_meta_ch
        .combine(qc_metrics_object_ch)
        .map { meta, qc_metrics_object ->
            tuple(
                meta + [
                    stage         : 'qc_metrics',
                    artifact_type : 'qc_metrics_seurat'
                ],
                qc_metrics_object
            )
        }

    approved_policy_artifact_ch = scientific_meta_ch
        .combine(approved_policy_ch)
        .map { meta, approved_policy ->
            tuple(
                meta + [
                    stage           : 'qc_threshold_approval',
                    artifact_type   : 'approved_qc_threshold_policy',
                    approval_status : 'approved'
                ],
                approved_policy
            )
        }

    QC_FILTER_POLICY(
        qc_metrics_artifact_ch,
        approved_policy_artifact_ch,
        qc_filter_script_ch
    )

    emit:

    filtered_object =
        QC_FILTER_POLICY.out.filtered_object

    summary =
        QC_FILTER_POLICY.out.summary

    by_threshold =
        QC_FILTER_POLICY.out.by_threshold

    by_sample =
        QC_FILTER_POLICY.out.by_sample

    cell_decisions =
        QC_FILTER_POLICY.out.cell_decisions

    session_info =
        QC_FILTER_POLICY.out.session_info
}

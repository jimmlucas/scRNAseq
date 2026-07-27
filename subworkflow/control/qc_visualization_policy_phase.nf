nextflow.enable.dsl=2

include { QC_VISUALIZATION_POLICY } from '../../modules/qc_visualization_policy'


workflow QC_VISUALIZATION_POLICY_PHASE {

    take:

    scientific_meta_ch
    qc_metrics_object_ch
    filtered_object_ch
    cell_decisions_ch
    approved_policy_ch

    main:

    visualization_script_ch = channel.value(
        file(
            "${projectDir}/scripts/01f_qc_visualization_policy.R",
            checkIfExists: true
        )
    )

    before_artifact_ch = scientific_meta_ch
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

    after_artifact_ch = scientific_meta_ch
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

    QC_VISUALIZATION_POLICY(
        before_artifact_ch,
        after_artifact_ch,
        cell_decisions_ch,
        approved_policy_artifact_ch,
        visualization_script_ch
    )

    emit:

    summary =
        QC_VISUALIZATION_POLICY.out.summary

    metric_distributions =
        QC_VISUALIZATION_POLICY.out.metric_distributions

    by_sample =
        QC_VISUALIZATION_POLICY.out.by_sample

    sample_retention =
        QC_VISUALIZATION_POLICY.out.sample_retention

    exclusion_reasons =
        QC_VISUALIZATION_POLICY.out.exclusion_reasons

    cell_metadata =
        QC_VISUALIZATION_POLICY.out.cell_metadata

    plots =
        QC_VISUALIZATION_POLICY.out.plots

    session_info =
        QC_VISUALIZATION_POLICY.out.session_info
}

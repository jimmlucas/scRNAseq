nextflow.enable.dsl = 2

include { QC_VISUALIZATION } from '../modules/qc_visualization'


workflow QC_VISUALIZATION_PHASE {

    take:
    scientific_meta_ch
    filtered_object_ch
    config_ch

    main:

    qc_visualization_script_ch = channel.value(
        file(
            "${projectDir}/scripts/01b_qc_visualization.R",
            checkIfExists: true
        )
    )

    /*
     * Standalone checkpoint reconstruction.
     *
     * In qc_visualization_phase, the input RDS comes from a previously
     * published QC_FILTER result. Reconstruct its scientific artifact
     * contract before passing it to the module.
     */
    filtered_artifact_ch =
        scientific_meta_ch
            .combine(filtered_object_ch)
            .map { meta, filtered_object ->

                def checkpoint_meta = meta + [
                        stage        : "qc_filter",
                        artifact_type: "filtered_seurat"
                ]

                tuple(checkpoint_meta, filtered_object)
            }

    QC_VISUALIZATION(
        filtered_artifact_ch,
        config_ch,
        qc_visualization_script_ch
    )

    emit:

    summary =
        QC_VISUALIZATION.out.summary

    by_sample =
        QC_VISUALIZATION.out.by_sample

    by_group =
        QC_VISUALIZATION.out.by_group

    cell_metadata =
        QC_VISUALIZATION.out.cell_metadata

    plots =
        QC_VISUALIZATION.out.plots

    session_info =
        QC_VISUALIZATION.out.session_info
}
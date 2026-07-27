nextflow.enable.dsl = 2

include { QC_DOUBLET_DETECTION } from '../modules/doublet_detection'


workflow DOUBLET_PHASE {

    take:
    scientific_meta_ch
    filtered_object_ch
    config_ch

    main:

    doublet_detection_script_ch = channel.value(
        file(
            "${projectDir}/scripts/01c_doublet_detection.R",
            checkIfExists: true
        )
    )

    /*
     * Standalone checkpoint reconstruction.
     *
     * The RDS supplied to doublet_phase is expected to be a previously
     * published QC_FILTER artifact.
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

    QC_DOUBLET_DETECTION(
        filtered_artifact_ch,
        config_ch,
        doublet_detection_script_ch
    )

    emit:

    singlet_object = QC_DOUBLET_DETECTION.out.singlet_object
    summary        = QC_DOUBLET_DETECTION.out.summary
    by_sample      = QC_DOUBLET_DETECTION.out.by_sample
    metadata       = QC_DOUBLET_DETECTION.out.metadata
    session_info   = QC_DOUBLET_DETECTION.out.session_info
}
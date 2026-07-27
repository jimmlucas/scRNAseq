nextflow.enable.dsl = 2

include { QC_SCORES } from '../modules/qc_scores'


workflow QC_SCORES_PHASE {

    take:
    scientific_meta_ch
    singlet_object_ch
    config_ch

    main:

    qc_scores_script_ch = channel.value(
        file(
            "${projectDir}/scripts/01d_qc_scores.R",
            checkIfExists: true
        )
    )

    /*
     * Standalone checkpoint reconstruction.
     *
     * The supplied RDS is expected to be a previously published
     * QC_DOUBLET_DETECTION artifact.
     */
    singlet_artifact_ch =
        scientific_meta_ch
            .combine(singlet_object_ch)
            .map { meta, singlet_object ->

                def checkpoint_meta = meta + [
                    stage        : "doublet_detection",
                    artifact_type: "singlet_seurat"
                ]

                tuple(checkpoint_meta, singlet_object)
            }

    QC_SCORES(
        singlet_artifact_ch,
        config_ch,
        qc_scores_script_ch
    )

    emit:

    scored_object = QC_SCORES.out.scored_object
    summary        = QC_SCORES.out.summary
    phase_summary  = QC_SCORES.out.phase_summary
    metadata       = QC_SCORES.out.metadata
    session_info   = QC_SCORES.out.session_info
}
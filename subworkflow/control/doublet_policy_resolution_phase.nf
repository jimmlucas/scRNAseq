nextflow.enable.dsl=2

include { DOUBLET_POLICY_RESOLUTION } from '../../modules/control/resolve_doublet_policy'


workflow DOUBLET_POLICY_RESOLUTION_PHASE {

    take:

    scientific_meta_ch
    doublet_summary_ch
    doublet_by_sample_ch
    doublet_cell_scores_ch

    main:

    resolver_script_ch = channel.value(
        file(
            "${projectDir}/scripts/control/resolve_doublet_policy.py",
            checkIfExists: true
        )
    )

    DOUBLET_POLICY_RESOLUTION(
        scientific_meta_ch,
        doublet_summary_ch,
        doublet_by_sample_ch,
        doublet_cell_scores_ch,
        resolver_script_ch
    )

    emit:

    policy_yaml =
        DOUBLET_POLICY_RESOLUTION.out.policy_yaml

    policy_tsv =
        DOUBLET_POLICY_RESOLUTION.out.policy_tsv

    summary =
        DOUBLET_POLICY_RESOLUTION.out.summary

    session_info =
        DOUBLET_POLICY_RESOLUTION.out.session_info
}

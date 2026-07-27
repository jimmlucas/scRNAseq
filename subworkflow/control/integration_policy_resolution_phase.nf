nextflow.enable.dsl=2

include { INTEGRATION_POLICY_RESOLUTION } from '../../modules/control/resolve_integration_policy'


workflow INTEGRATION_POLICY_RESOLUTION_PHASE {

    take:

    scientific_meta_ch
    pca_object_ch
    dimensionality_summary_ch
    approved_pc_selection_policy_ch
    analysis_readiness_ch
    applied_normalization_policy_ch

    main:

    resolver_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/control/resolve_integration_policy.py",
            checkIfExists: true
        )
    )

    INTEGRATION_POLICY_RESOLUTION(
        scientific_meta_ch,
        pca_object_ch,
        dimensionality_summary_ch,
        approved_pc_selection_policy_ch,
        analysis_readiness_ch,
        applied_normalization_policy_ch,
        resolver_script_ch
    )

    emit:

    policy_yaml =
        INTEGRATION_POLICY_RESOLUTION.out.policy_yaml

    policy_tsv =
        INTEGRATION_POLICY_RESOLUTION.out.policy_tsv

    summary =
        INTEGRATION_POLICY_RESOLUTION.out.summary

    session_info =
        INTEGRATION_POLICY_RESOLUTION.out.session_info
}

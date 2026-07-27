nextflow.enable.dsl=2

include { NORMALIZATION_POLICY_RESOLUTION } from '../../modules/control/resolve_normalization_policy'


workflow NORMALIZATION_POLICY_RESOLUTION_PHASE {

    take:

    scientific_meta_ch
    readiness_ch
    post_doublet_summary_ch
    post_doublet_by_sample_ch
    post_doublet_warnings_ch

    main:

    resolver_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/control/resolve_normalization_policy.py",
            checkIfExists: true
        )
    )

    NORMALIZATION_POLICY_RESOLUTION(
        scientific_meta_ch,
        readiness_ch,
        post_doublet_summary_ch,
        post_doublet_by_sample_ch,
        post_doublet_warnings_ch,
        resolver_script_ch
    )

    emit:

    policy_yaml =
        NORMALIZATION_POLICY_RESOLUTION.out.policy_yaml

    policy_tsv =
        NORMALIZATION_POLICY_RESOLUTION.out.policy_tsv

    summary =
        NORMALIZATION_POLICY_RESOLUTION.out.summary

    session_info =
        NORMALIZATION_POLICY_RESOLUTION.out.session_info
}

nextflow.enable.dsl=2

include { PC_SELECTION_POLICY_RESOLUTION } from '../../modules/control/resolve_pc_selection_policy'


workflow PC_SELECTION_POLICY_RESOLUTION_PHASE {

    take:

    scientific_meta_ch
    dimensionality_summary_ch
    variable_features_ch
    elbow_metrics_ch

    main:

    resolver_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/control/resolve_pc_selection_policy.py",
            checkIfExists: true
        )
    )

    PC_SELECTION_POLICY_RESOLUTION(
        scientific_meta_ch,
        dimensionality_summary_ch,
        variable_features_ch,
        elbow_metrics_ch,
        resolver_script_ch
    )

    emit:

    policy_yaml =
        PC_SELECTION_POLICY_RESOLUTION.out.policy_yaml

    policy_tsv =
        PC_SELECTION_POLICY_RESOLUTION.out.policy_tsv

    summary =
        PC_SELECTION_POLICY_RESOLUTION.out.summary

    session_info =
        PC_SELECTION_POLICY_RESOLUTION.out.session_info
}

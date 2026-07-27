nextflow.enable.dsl=2

include { VALIDATE_ANALYSIS_PLAN } from '../../modules/control/validate_analysis_plan'
include { VALIDATE_NORMALIZATION_APPROVAL } from '../../modules/control/validate_normalization_approval'


workflow NORMALIZATION_APPROVAL_PHASE {

    take:

    plan_artifact_ch
    reviewed_policy_artifact_ch

    main:

    validation_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/control/validate_analysis_plan.py",
            checkIfExists: true
        )
    )

    approval_validator_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/control/validate_normalization_approval.py",
            checkIfExists: true
        )
    )

    VALIDATE_ANALYSIS_PLAN(
        plan_artifact_ch,
        validation_script_ch
    )

    VALIDATE_NORMALIZATION_APPROVAL(
        VALIDATE_ANALYSIS_PLAN.out.validated_plan,
        reviewed_policy_artifact_ch,
        approval_validator_script_ch
    )

    emit:

    validated_plan =
        VALIDATE_ANALYSIS_PLAN.out.validated_plan

    plan_validation_summary =
        VALIDATE_ANALYSIS_PLAN.out.validation_summary

    plan_validation_session_info =
        VALIDATE_ANALYSIS_PLAN.out.session_info

    approved_policy =
        VALIDATE_NORMALIZATION_APPROVAL.out.approved_policy

    approval_summary =
        VALIDATE_NORMALIZATION_APPROVAL.out.summary

    approval_session_info =
        VALIDATE_NORMALIZATION_APPROVAL.out.session_info
}

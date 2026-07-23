nextflow.enable.dsl=2

include {
    VALIDATE_ANALYSIS_PLAN
} from '../../modules/control/validate_analysis_plan'

include {
    RESOLVE_ANALYSIS_PLAN
} from '../../modules/control/resolve_analysis_plan'


workflow BBUDY_QC_PROCESS {

    take:

    scientific_meta_ch
    analysis_plan_ch

    main:

    validation_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/control/validate_analysis_plan.py",
            checkIfExists: true
        )
    )

    routing_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/control/resolve_analysis_plan.py",
            checkIfExists: true
        )
    )

    module_catalog_ch = Channel.value(
        file(
            "${projectDir}/control/module_catalog.yaml",
            checkIfExists: true
        )
    )

    plan_artifact_ch = scientific_meta_ch
        .combine(analysis_plan_ch)
        .map { meta, plan ->
            tuple(
                meta + [
                    stage         : 'analysis_plan',
                    artifact_type : 'analysis_plan_yaml'
                ],
                plan
            )
        }

    VALIDATE_ANALYSIS_PLAN(
        plan_artifact_ch,
        validation_script_ch
    )

    RESOLVE_ANALYSIS_PLAN(
        VALIDATE_ANALYSIS_PLAN.out.validated_plan,
        module_catalog_ch,
        routing_script_ch
    )

    emit:

    validated_plan =
        VALIDATE_ANALYSIS_PLAN.out.validated_plan

    validation_summary =
        VALIDATE_ANALYSIS_PLAN.out.validation_summary

    resolved_routing =
        RESOLVE_ANALYSIS_PLAN.out.resolved_routing

    routing_summary =
        RESOLVE_ANALYSIS_PLAN.out.routing_summary

    validation_session_info =
        VALIDATE_ANALYSIS_PLAN.out.session_info

    routing_session_info =
        RESOLVE_ANALYSIS_PLAN.out.session_info
}

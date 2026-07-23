nextflow.enable.dsl=2

include { VALIDATE_ANALYSIS_PLAN }   from '../../modules/control/validate_analysis_plan'
include { RESOLVE_ANALYSIS_PLAN }    from '../../modules/control/resolve_analysis_plan'
include { INPUT_DISCOVERY }          from '../../modules/control/input_discovery'
include { RESOLVE_LOAD_CONFIG }      from '../../modules/control/resolve_load_config'
include { LOAD_DATA }                from '../../modules/load'
include { QC_METRICS }               from '../../modules/qc_metrics'
include { QC_THRESHOLD_RESOLUTION }  from '../../modules/control/resolve_qc_thresholds'

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

    input_discovery_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/control/input_discovery.py",
            checkIfExists: true
        )
    )

    load_config_resolve_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/control/resolve_load_config.py",
            checkIfExists: true
        )
    )

    load_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/00-Loading_data.R",
            checkIfExists: true
        )
    )

    qc_metrics_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/01a_qc_metrics.R",
            checkIfExists: true
        )
    )

    qc_threshold_resolution_script_ch = Channel.value(
        file(
            "${projectDir}/scripts/control/resolve_qc_thresholds.py",
            checkIfExists: true
        )
    )

    VALIDATE_ANALYSIS_PLAN(
        plan_artifact_ch,
        validation_script_ch
    )

    RESOLVE_ANALYSIS_PLAN(
        VALIDATE_ANALYSIS_PLAN.out.validated_plan,
        module_catalog_ch,
        routing_script_ch
    )

    INPUT_DISCOVERY(
        VALIDATE_ANALYSIS_PLAN.out.validated_plan,
        RESOLVE_ANALYSIS_PLAN.out.resolved_routing,
        input_discovery_script_ch
    )

    RESOLVE_LOAD_CONFIG(
        VALIDATE_ANALYSIS_PLAN.out.validated_plan,
        RESOLVE_ANALYSIS_PLAN.out.resolved_routing,
        INPUT_DISCOVERY.out.inventory,
        load_config_resolve_script_ch
    )

    LOAD_DATA(
        RESOLVE_LOAD_CONFIG.out.load_config,
        load_script_ch
    )

    QC_METRICS(
        LOAD_DATA.out.seurat_object,
        RESOLVE_LOAD_CONFIG.out.load_config,
        qc_metrics_script_ch
    )

    QC_THRESHOLD_RESOLUTION(
        VALIDATE_ANALYSIS_PLAN.out.validated_plan,
        QC_METRICS.out.summary,
        QC_METRICS.out.by_sample,
        QC_METRICS.out.cell_metadata,
        qc_threshold_resolution_script_ch
    )

    emit:

    validated_plan = VALIDATE_ANALYSIS_PLAN.out.validated_plan
    validation_summary = VALIDATE_ANALYSIS_PLAN.out.validation_summary
    validation_session_info = VALIDATE_ANALYSIS_PLAN.out.session_info

    resolved_routing = RESOLVE_ANALYSIS_PLAN.out.resolved_routing
    resolved_load_config = RESOLVE_LOAD_CONFIG.out.load_config
    routing_summary = RESOLVE_ANALYSIS_PLAN.out.routing_summary
    routing_session_info = RESOLVE_ANALYSIS_PLAN.out.session_info

    load_config_summary = RESOLVE_LOAD_CONFIG.out.summary
    load_config_session_info = RESOLVE_LOAD_CONFIG.out.session_info
    loaded_seurat = LOAD_DATA.out.seurat_object
    load_summary = LOAD_DATA.out.summary
    load_session_info = LOAD_DATA.out.session_info

    input_inventory = INPUT_DISCOVERY.out.inventory
    input_sample_summary = INPUT_DISCOVERY.out.sample_summary
    input_discovery_summary = INPUT_DISCOVERY.out.summary
    input_discovery_session_info = INPUT_DISCOVERY.out.session_info

    qc_metrics_object = QC_METRICS.out.metrics_object
    qc_metrics_summary = QC_METRICS.out.summary
    qc_metrics_by_sample = QC_METRICS.out.by_sample
    qc_metrics_cell_metadata = QC_METRICS.out.cell_metadata
    qc_metrics_session_info = QC_METRICS.out.session_info

    qc_threshold_policy_yaml =
        QC_THRESHOLD_RESOLUTION.out.policy_yaml

    qc_threshold_policy_tsv =
        QC_THRESHOLD_RESOLUTION.out.policy_tsv

    qc_threshold_resolution_summary =
        QC_THRESHOLD_RESOLUTION.out.summary

    qc_threshold_resolution_session_info =
        QC_THRESHOLD_RESOLUTION.out.session_info

}
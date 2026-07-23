nextflow.enable.dsl=2


process RESOLVE_ANALYSIS_PLAN {

    tag "${meta.project_id ?: 'project'}:${meta.biological_run_id ?: 'run'}"

    label 'lightweight'

    publishDir(
        "${params.outdir}/control/analysis_plan",
        mode: 'copy',
        overwrite: true
    )

    input:

    tuple val(meta), path(validated_plan)
    path module_catalog
    path routing_script

    output:

    tuple val(output_meta), path("analysis_plan.routing.yaml"),
        emit: resolved_routing

    tuple val(output_meta), path("analysis_plan.routing.tsv"),
        emit: routing_summary

    path "analysis_plan_routing.session.txt",
        emit: session_info

    script:

    output_meta = meta + [
        stage         : 'plan_routing',
        artifact_type : 'resolved_analysis_routing'
    ]

    """
    python3 ${routing_script} \
        --plan ${validated_plan} \
        --module-catalog ${module_catalog} \
        --routing-output analysis_plan.routing.yaml \
        --summary analysis_plan.routing.tsv

    {
        python3 --version
        echo "resolver=${routing_script}"
        echo "validated_plan=${validated_plan}"
        echo "module_catalog=${module_catalog}"
    } > analysis_plan_routing.session.txt 2>&1
    """

    stub:

    output_meta = meta + [
        stage         : 'plan_routing',
        artifact_type : 'resolved_analysis_routing'
    ]

    """
    cat > analysis_plan.routing.yaml <<'YAML'
    routing_resolution:
      schema_version: "1.0"
      plan_id: "stub_plan"
      plan_approval_status: "approved"
      status: "ready"
      executable_routes:
        - load_data
      blocking_routes: []
    routes:
      load_data:
        module: "LOAD_DATA"
        requested: true
        implementation_state: "available"
        catalog_state: "implemented"
        approval_policy: "not_required"
        approval_satisfied: true
        executable: true
        blocking: false
        reason: "ready"
    YAML

    cat > analysis_plan.routing.tsv <<'TSV'
    route	module	requested	implementation_state	catalog_state	approval_policy	approval_satisfied	executable	blocking	reason
    load_data	LOAD_DATA	true	available	implemented	not_required	true	true	false	ready
    TSV

    echo "stub=true" > analysis_plan_routing.session.txt
    """
}

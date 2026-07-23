nextflow.enable.dsl=2


process RESOLVE_LOAD_CONFIG {

    tag "${plan_meta.project_id ?: 'project'}:${plan_meta.biological_run_id ?: 'run'}"

    label 'lightweight'

    publishDir(
        "${params.outdir}/control/load_config",
        mode: 'copy',
        overwrite: true
    )

    input:

    tuple val(plan_meta), path(validated_plan)
    tuple val(routing_meta), path(resolved_routing)
    tuple val(inventory_meta), path(input_inventory)
    path resolver_script

    output:

    tuple val(output_meta), path("load_config.resolved.yaml"),
        optional: true,
        emit: load_config

    tuple val(output_meta), path("load_config.resolution.tsv"),
        emit: summary

    path "load_config_resolution.session.txt",
        emit: session_info

    script:

    if (
        plan_meta.project_id != routing_meta.project_id
        || plan_meta.project_id != inventory_meta.project_id
        || plan_meta.biological_run_id != routing_meta.biological_run_id
        || plan_meta.biological_run_id != inventory_meta.biological_run_id
    ) {
        error(
            "RESOLVE_LOAD_CONFIG received incompatible metadata."
        )
    }

    output_meta = plan_meta + [
        stage         : 'load_config_resolution',
        artifact_type : 'resolved_load_config'
    ]

    """
    python3 ${resolver_script} \
        --plan ${validated_plan} \
        --routing ${resolved_routing} \
        --inventory ${input_inventory} \
        --project-dir ${projectDir} \
        --project-id '${plan_meta.project_id}' \
        --biological-run-id '${plan_meta.biological_run_id}' \
        --resolved-config load_config.resolved.yaml \
        --summary load_config.resolution.tsv

    {
        python3 --version
        echo "resolver=${resolver_script}"
        echo "validated_plan=${validated_plan}"
        echo "resolved_routing=${resolved_routing}"
        echo "input_inventory=${input_inventory}"
    } > load_config_resolution.session.txt 2>&1
    """

    stub:

    output_meta = plan_meta + [
        stage         : 'load_config_resolution',
        artifact_type : 'resolved_load_config'
    ]

    """
    cat > load_config.resolved.yaml <<'YAML'
    project:
      name: "stub"
      organism: "mouse"
    input:
      targets: "stub_targets.csv"
      count_dir: "stub_count"
    output:
      merged_object: "seu.rds"
    YAML

    cat > load_config.resolution.tsv <<'TSV'
    plan_id	status	requested	executable	reason	resolved_config
    stub_plan	ready	true	true	ready	stub_config.yaml
    TSV

    echo "stub=true" > load_config_resolution.session.txt
    """
}

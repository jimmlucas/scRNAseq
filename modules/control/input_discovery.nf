nextflow.enable.dsl=2


process INPUT_DISCOVERY {

    tag "${meta.project_id ?: 'project'}:${meta.biological_run_id ?: 'run'}"

    label 'lightweight'

    publishDir(
        "${params.outdir}/control/input_discovery",
        mode: 'copy',
        overwrite: true
    )

    input:

    tuple val(meta), path(validated_plan)
    tuple val(routing_meta), path(resolved_routing)
    path discovery_script

    output:

    tuple val(output_meta), path("input_discovery.inventory.yaml"),
        emit: inventory

    tuple val(output_meta), path("input_discovery.samples.tsv"),
        emit: sample_summary

    tuple val(output_meta), path("input_discovery.summary.tsv"),
        emit: summary

    path "input_discovery.session.txt",
        emit: session_info

    script:

    if (
        meta.project_id != routing_meta.project_id
        || meta.biological_run_id != routing_meta.biological_run_id
    ) {
        error(
            "INPUT_DISCOVERY received incompatible metadata: " +
            "plan=${meta}, routing=${routing_meta}"
        )
    }

    output_meta = meta + [
        stage         : 'input_discovery',
        artifact_type : 'input_inventory'
    ]

    """
    python3 ${discovery_script} \
        --plan ${validated_plan} \
        --routing ${resolved_routing} \
        --project-dir ${projectDir} \
        --inventory input_discovery.inventory.yaml \
        --sample-summary input_discovery.samples.tsv \
        --summary input_discovery.summary.tsv

    {
        python3 --version
        echo "discovery_script=${discovery_script}"
        echo "validated_plan=${validated_plan}"
        echo "resolved_routing=${resolved_routing}"
        echo "project_dir=${projectDir}"
    } > input_discovery.session.txt 2>&1
    """

    stub:

    output_meta = meta + [
        stage         : 'input_discovery',
        artifact_type : 'input_inventory'
    ]

    """
    cat > input_discovery.inventory.yaml <<'YAML'
    input_discovery:
      schema_version: "1.0"
      plan_id: "stub_plan"
      requested: true
      status: "ready"
      samples:
        sample_count: 0
        biological_sample_count: 0
      warnings: []
      errors: []
    YAML

    cat > input_discovery.samples.tsv <<'TSV'
    row_number	sample_id	biological_sample_id	group	run_id	lane_id	technical_replicate_id
    TSV

    cat > input_discovery.summary.tsv <<'TSV'
    plan_id	status	sample_count	biological_sample_count	missing_column_count	warning_count	error_count
    stub_plan	ready	0	0	0	0	0
    TSV

    echo "stub=true" > input_discovery.session.txt
    """
}

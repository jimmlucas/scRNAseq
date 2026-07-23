nextflow.enable.dsl=2


process VALIDATE_ANALYSIS_PLAN {

    tag "${meta.project_id ?: 'project'}:${meta.biological_run_id ?: 'run'}"

    label 'lightweight'

    publishDir(
        "${params.outdir}/control/analysis_plan",
        mode: 'copy',
        overwrite: true
    )

    input:

    tuple val(meta), path(analysis_plan)
    path validation_script

    output:

    tuple val(meta), path("analysis_plan.validated.yaml"),
        emit: validated_plan

    tuple val(meta), path("analysis_plan.validation.tsv"),
        emit: validation_summary

    path "analysis_plan_validation.session.txt",
        emit: session_info

    script:

    """
    python3 ${validation_script} \
        --plan ${analysis_plan} \
        --validated-plan analysis_plan.validated.yaml \
        --summary analysis_plan.validation.tsv

    {
        python3 --version
        echo "validator=${validation_script}"
        echo "input_plan=${analysis_plan}"
    } > analysis_plan_validation.session.txt 2>&1
    """

    stub:

    """
    cat > analysis_plan.validated.yaml <<'YAML'
    plan:
      schema_version: "1.0"
      plan_id: "stub_plan"
    validation:
      status: "valid"
      errors: []
      warnings: []
    YAML

    cat > analysis_plan.validation.tsv <<'TSV'
    plan_id	status	error_count	warning_count
    stub_plan	valid	0	0
    TSV

    echo "stub=true" > analysis_plan_validation.session.txt
    """
}

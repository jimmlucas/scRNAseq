nextflow.enable.dsl=2


process DOUBLET_FILTER {

    tag "${object_meta.project_id ?: 'project'}:${object_meta.biological_run_id ?: 'run'}"

    label 'high_memory'

    publishDir(
        "${params.outdir}/02_doublet_filter",
        mode: 'copy',
        overwrite: true
    )

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        "docker://${params.r_image}" :
        params.r_image}"

    input:

    tuple val(object_meta), path(scored_object)
    tuple val(policy_meta), path(approved_policy)
    path doublet_filter_script

    output:

    tuple val(output_meta), path("seu_doublet_filtered.rds"),
        emit: filtered_object

    tuple val(summary_meta), path("02_doublet_filter_summary.csv"),
        emit: summary

    path "02_doublet_filter_by_sample.csv",
        emit: by_sample

    path "02_doublet_filter_removed_cells.csv",
        emit: removed_cells

    tuple val(applied_policy_meta), path("applied_doublet_policy.yaml"),
        emit: applied_policy

    path "02_doublet_filter_sessionInfo.txt",
        emit: session_info

    script:

    if (
        object_meta.project_id != policy_meta.project_id
        || object_meta.biological_run_id != policy_meta.biological_run_id
    ) {
        error(
            "DOUBLET_FILTER received incompatible metadata between " +
            "doublet-scored object and approved doublet policy."
        )
    }

    if (policy_meta.approval_status != 'approved') {
        error(
            "DOUBLET_FILTER requires an approved doublet policy."
        )
    }

    output_meta = object_meta + [
        stage         : 'doublet_filter',
        artifact_type : 'singlet_seurat',
        policy_status : policy_meta.approval_status
    ]

    summary_meta = object_meta + [
        stage         : 'doublet_filter',
        artifact_type : 'metrics'
    ]

    applied_policy_meta = object_meta + [
        stage         : 'doublet_filter',
        artifact_type : 'applied_doublet_policy',
        policy_status : 'applied'
    ]

    """
    Rscript ${doublet_filter_script} \
        --input ${scored_object} \
        --policy ${approved_policy} \
        --output seu_doublet_filtered.rds \
        --project-id '${object_meta.project_id ?: ""}' \
        --biological-run-id '${object_meta.biological_run_id ?: ""}'
    """

    stub:

    output_meta = object_meta + [
        stage         : 'doublet_filter',
        artifact_type : 'singlet_seurat',
        policy_status : 'approved'
    ]

    summary_meta = object_meta + [
        stage         : 'doublet_filter',
        artifact_type : 'metrics'
    ]

    applied_policy_meta = object_meta + [
        stage         : 'doublet_filter',
        artifact_type : 'applied_doublet_policy',
        policy_status : 'applied'
    ]

    """
    touch seu_doublet_filtered.rds

    cat > 02_doublet_filter_summary.csv <<'CSV'
    metric,value
    cells_input,0
    cells_output,0
    cells_removed,0
    features_input,0
    features_output,0
    samples,0
    approval_status,approved
    validation_status,valid
    application_status,applied
    CSV

    cat > 02_doublet_filter_by_sample.csv <<'CSV'
    sample,action,cells_before,approved_doublets,cells_removed,cells_after,filtering_applied
    CSV

    cat > 02_doublet_filter_removed_cells.csv <<'CSV'
    cell_id,sample,action,doublet_class,is_doublet,doublet_detection_status,doublet_decision_status
    CSV

    cat > applied_doublet_policy.yaml <<'YAML'
    doublet_policy:
      approval:
        status: "approved"
      validation:
        status: "valid"
      application:
        status: "applied"
        filtering_applied: true
        cells_removed: 0
        applied_at: "2026-01-01T00:00:00Z"
    YAML

    cat > 02_doublet_filter_sessionInfo.txt <<'TXT'
    stub=true
    process=DOUBLET_FILTER
    TXT
    """
}

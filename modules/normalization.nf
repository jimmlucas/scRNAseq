nextflow.enable.dsl=2


process NORMALIZATION {

    tag "${object_meta.project_id ?: 'project'}:${object_meta.biological_run_id ?: 'run'}"

    label 'high_memory'

    publishDir(
        "${params.outdir}/04_normalization",
        mode: 'copy',
        overwrite: true
    )

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        "docker://${params.r_image}" :
        params.r_image}"

    input:

    tuple val(object_meta), path(filtered_object)
    tuple val(policy_meta), path(approved_policy)
    path normalization_script

    output:

    tuple val(output_meta), path("seu_normalized.rds"),
        emit: normalized_object

    tuple val(sample_objects_meta), path("normalized_samples"),
        emit: sample_objects

    tuple val(summary_meta), path("04_normalization_summary.csv"),
        emit: summary

    path "04_normalization_by_sample.csv",
        emit: by_sample

    path "04_variable_features.csv",
        emit: variable_features

    tuple val(applied_policy_meta), path("applied_normalization_policy.yaml"),
        emit: applied_policy

    path "04_normalization_sessionInfo.txt",
        emit: session_info

    script:

    if (
        object_meta.project_id != policy_meta.project_id
        || object_meta.biological_run_id != policy_meta.biological_run_id
    ) {
        error(
            "NORMALIZATION received incompatible metadata between " +
            "filtered object and approved normalization policy."
        )
    }

    if (policy_meta.approval_status != 'approved') {
        error(
            "NORMALIZATION requires an approved normalization policy."
        )
    }

    output_meta = object_meta + [
        stage         : 'normalization',
        artifact_type : 'normalized_seurat',
        policy_status : policy_meta.approval_status
    ]

    sample_objects_meta = object_meta + [
        stage         : 'normalization',
        artifact_type : 'normalized_sample_objects',
        policy_status : policy_meta.approval_status
    ]

    summary_meta = object_meta + [
        stage         : 'normalization',
        artifact_type : 'normalization_summary'
    ]

    applied_policy_meta = object_meta + [
        stage         : 'normalization',
        artifact_type : 'applied_normalization_policy',
        policy_status : 'applied'
    ]

    """
    Rscript ${normalization_script} \
        --input ${filtered_object} \
        --policy ${approved_policy} \
        --output seu_normalized.rds \
        --sample-output-dir normalized_samples \
        --project-id '${object_meta.project_id ?: ""}' \
        --biological-run-id '${object_meta.biological_run_id ?: ""}'
    """

    stub:

    output_meta = object_meta + [
        stage         : 'normalization',
        artifact_type : 'normalized_seurat',
        policy_status : 'approved'
    ]

    sample_objects_meta = object_meta + [
        stage         : 'normalization',
        artifact_type : 'normalized_sample_objects',
        policy_status : 'approved'
    ]

    summary_meta = object_meta + [
        stage         : 'normalization',
        artifact_type : 'normalization_summary'
    ]

    applied_policy_meta = object_meta + [
        stage         : 'normalization',
        artifact_type : 'applied_normalization_policy',
        policy_status : 'applied'
    ]

    """
    touch seu_normalized.rds
    mkdir -p normalized_samples
    touch normalized_samples/stub.normalized.rds

    cat > 04_normalization_summary.csv <<'CSV'
    metric,value
    cells_input,0
    cells_output,0
    features_input,0
    features_output,0
    samples,0
    method,LogNormalize
    scope,per_sample
    assay,RNA
    variable_features_method,vst
    nfeatures_requested,2000
    scale_factor,10000
    regressors,none
    approval_status,approved
    validation_status,valid
    application_status,applied
    normalized_object,seu_normalized.rds
    CSV

    cat > 04_normalization_by_sample.csv <<'CSV'
    sample,cells,features,method,assay,variable_features_method,variable_features,nfeatures_requested,scale_factor,regressors,output_file
    CSV

    cat > 04_variable_features.csv <<'CSV'
    sample,rank,feature
    CSV

    cat > applied_normalization_policy.yaml <<'YAML'
    normalization_policy:
      approval:
        status: "approved"
      validation:
        status: "valid"
        validated_for_normalization: true
      application:
        status: "applied"
        normalization_applied: true
        cells_normalized: 0
        features_preserved: 0
        samples_normalized: 0
        applied_at: "2026-01-01T00:00:00Z"
    YAML

    cat > 04_normalization_sessionInfo.txt <<'TXT'
    stub=true
    process=NORMALIZATION
    TXT
    """
}

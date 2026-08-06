nextflow.enable.dsl=2


process DIMENSIONALITY {

    tag "${object_meta.project_id ?: 'project'}:${object_meta.biological_run_id ?: 'run'}"

    label 'high_memory'

    publishDir(
        "${params.outdir}/05_dimensionality",
        mode: 'copy',
        overwrite: true
    )

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        "docker://${params.r_image}" :
        params.r_image}"

    input:

    tuple val(object_meta), path(normalized_object)
    tuple val(policy_meta), path(applied_policy)
    path dimensionality_script

    output:

    tuple val(output_meta), path("seu_pca.rds"),
        emit: pca_object

    tuple val(summary_meta), path("05_dimensionality_summary.csv"),
        emit: summary

    path "05_variable_features_summary.csv",
        emit: variable_features_summary

    path "05_pca_embeddings_preview.csv",
        emit: embeddings_preview

    path "05_pca_loadings_preview.csv",
        emit: loadings_preview

    path "05_pca_elbow_metrics.csv",
        emit: elbow_metrics

    path "05_dimensionality_sessionInfo.txt",
        emit: session_info

    script:

    if (
        object_meta.project_id != policy_meta.project_id
        || object_meta.biological_run_id != policy_meta.biological_run_id
    ) {
        error(
            "DIMENSIONALITY received incompatible metadata between " +
            "normalized object and applied normalization policy."
        )
    }

    if (policy_meta.application_status != 'applied') {
        error(
            "DIMENSIONALITY requires an applied normalization policy."
        )
    }

    output_meta = object_meta + [
        stage         : 'dimensionality',
        artifact_type : 'pca_seurat',
        reduction     : 'pca'
    ]

    summary_meta = object_meta + [
        stage         : 'dimensionality',
        artifact_type : 'dimensionality_summary',
        reduction     : 'pca'
    ]

    """
    Rscript ${dimensionality_script} \
        --input ${normalized_object} \
        --policy ${applied_policy} \
        --output seu_pca.rds \
        --project-id '${object_meta.project_id ?: ""}' \
        --biological-run-id '${object_meta.biological_run_id ?: ""}' \
        --npcs ${params.dimensionality_npcs}
    """

    stub:

    output_meta = object_meta + [
        stage         : 'dimensionality',
        artifact_type : 'pca_seurat',
        reduction     : 'pca'
    ]

    summary_meta = object_meta + [
        stage         : 'dimensionality',
        artifact_type : 'dimensionality_summary',
        reduction     : 'pca'
    ]

    """
    touch seu_pca.rds

    cat > 05_dimensionality_summary.csv <<'CSV'
    metric,value
    cells_input,0
    cells_output,0
    features_input,0
    features_output,0
    assay,RNA
    normalization_method,LogNormalize
    variable_features,0
    nfeatures_requested,2000
    npcs_requested,30
    npcs_computed,0
    scaling_applied,true
    pca_applied,true
    approval_status,approved
    validation_status,valid
    normalization_application_status,applied
    pca_object,seu_pca.rds
    CSV

    cat > 05_variable_features_summary.csv <<'CSV'
    rank,feature
    CSV

    cat > 05_pca_embeddings_preview.csv <<'CSV'
    cell,PC_1,PC_2
    CSV

    cat > 05_pca_loadings_preview.csv <<'CSV'
    feature,PC_1,PC_2
    CSV

    cat > 05_pca_elbow_metrics.csv <<'CSV'
    pc,stdev,variance,variance_percent,cumulative_variance_percent
    CSV

    cat > 05_dimensionality_sessionInfo.txt <<'TXT'
    stub=true
    process=DIMENSIONALITY
    TXT
    """
}

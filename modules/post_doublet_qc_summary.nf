nextflow.enable.dsl=2


process POST_DOUBLET_QC_SUMMARY {

    tag "${object_meta.project_id ?: 'project'}:${object_meta.biological_run_id ?: 'run'}"

    label 'high_memory'

    publishDir(
        "${params.outdir}/03_post_doublet_qc",
        mode: 'copy',
        overwrite: true
    )

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.r_image.docker}" :
        params.r_image}"

    input:

    tuple val(object_meta), path(filtered_object)
    tuple val(policy_meta), path(applied_policy)
    path summary_script

    output:

    tuple val(summary_meta), path("03_post_doublet_qc_summary.csv"),
        emit: summary

    tuple val(by_sample_meta), path("03_post_doublet_qc_by_sample.csv"),
        emit: by_sample

    tuple val(warnings_meta), path("03_post_doublet_qc_warnings.tsv"),
        emit: warnings

    tuple val(readiness_meta), path("03_analysis_readiness.yaml"),
        emit: readiness

    path "03_sessionInfo.txt",
        emit: session_info

    script:

    if (
        object_meta.project_id != policy_meta.project_id
        || object_meta.biological_run_id != policy_meta.biological_run_id
    ) {
        error(
            "POST_DOUBLET_QC_SUMMARY received incompatible metadata " +
            "between filtered object and applied doublet policy."
        )
    }

    if (policy_meta.policy_status != 'applied') {
        error(
            "POST_DOUBLET_QC_SUMMARY requires an applied doublet policy."
        )
    }

    summary_meta = object_meta + [
        stage         : 'post_doublet_qc',
        artifact_type : 'metrics'
    ]

    by_sample_meta = object_meta + [
        stage         : 'post_doublet_qc',
        artifact_type : 'per_sample_metrics'
    ]

    warnings_meta = object_meta + [
        stage         : 'post_doublet_qc',
        artifact_type : 'analysis_warnings'
    ]

    readiness_meta = object_meta + [
        stage         : 'post_doublet_qc',
        artifact_type : 'analysis_readiness'
    ]

    """
    Rscript ${summary_script} \
        --input ${filtered_object} \
        --applied-policy ${applied_policy} \
        --project-id '${object_meta.project_id ?: ""}' \
        --biological-run-id '${object_meta.biological_run_id ?: ""}'
    """

    stub:

    summary_meta = object_meta + [
        stage         : 'post_doublet_qc',
        artifact_type : 'metrics'
    ]

    by_sample_meta = object_meta + [
        stage         : 'post_doublet_qc',
        artifact_type : 'per_sample_metrics'
    ]

    warnings_meta = object_meta + [
        stage         : 'post_doublet_qc',
        artifact_type : 'analysis_warnings'
    ]

    readiness_meta = object_meta + [
        stage         : 'post_doublet_qc',
        artifact_type : 'analysis_readiness'
    ]

    """
    cat > 03_post_doublet_qc_summary.csv <<'CSV'
    metric,value
    cells,0
    features,0
    samples,0
    doublets_removed,0
    warning_count,0
    CSV

    cat > 03_post_doublet_qc_by_sample.csv <<'CSV'
    sample,action,cells_original,doublets_removed,cells_final,median_nFeature_RNA,median_nCount_RNA,median_percent_mito
    CSV

    cat > 03_post_doublet_qc_warnings.tsv <<'TSV'
    severity	code	message	sample	metric	value	threshold
    TSV

    cat > 03_analysis_readiness.yaml <<'YAML'
    analysis_readiness:
      status: "ready"
      qc:
        completed: true
      doublets:
        scoring_completed: true
        human_approval_validated: true
        filtering_applied: true
      dataset:
        cells: 0
        features: 0
        samples: 0
      next_phase:
        recommended: "normalization_policy_resolution"
    YAML

    cat > 03_sessionInfo.txt <<'TXT'
    stub=true
    process=POST_DOUBLET_QC_SUMMARY
    TXT
    """
}

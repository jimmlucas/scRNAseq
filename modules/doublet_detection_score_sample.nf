nextflow.enable.dsl=2


process DOUBLET_DETECTION_SCORE_SAMPLE {

    tag "${meta.project_id ?: 'project'}:${meta.biological_run_id ?: 'run'}:${meta.sample ?: 'sample'}"

    label 'very_high_memory'

    publishDir(
        { "${params.outdir}/01g_doublet_detection_scoring/per_sample/${meta.sample_safe_id}" },
        mode: 'copy',
        overwrite: true
    )

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        "docker://${params.r_image}" :
        params.r_image}"

    input:

    tuple val(meta), path(sample_object), val(config_meta), path(config_file)
    path doublet_detection_script

    output:

    path "seu_doublet_scored.${meta.sample_safe_id}.rds",
        emit: scored_object

    path "01g_doublet_detection_summary.${meta.sample_safe_id}.csv",
        emit: summary

    path "01g_doublet_detection_by_sample.${meta.sample_safe_id}.csv",
        emit: by_sample

    path "01g_doublet_detection_cell_scores.${meta.sample_safe_id}.csv",
        emit: cell_scores

    path "01g_doublet_detection_pk_sweep.${meta.sample_safe_id}.csv",
        emit: pk_sweep

    path "01g_doublet_detection_performance.${meta.sample_safe_id}.csv",
        emit: performance

    path "01g_sessionInfo.${meta.sample_safe_id}.txt",
        emit: session_info

    script:

    if (
        meta.project_id != config_meta.project_id
        || meta.biological_run_id != config_meta.biological_run_id
    ) {
        error(
            "DOUBLET_DETECTION_SCORE_SAMPLE received incompatible " +
            "metadata between sample object and configuration."
        )
    }

    """
    Rscript ${doublet_detection_script} \
        --input ${sample_object} \
        --config ${config_file} \
        --cores ${task.cpus} \
        --output seu_doublet_scored.${meta.sample_safe_id}.rds

    mv 01g_doublet_detection_summary.csv \
        01g_doublet_detection_summary.${meta.sample_safe_id}.csv
    mv 01g_doublet_detection_by_sample.csv \
        01g_doublet_detection_by_sample.${meta.sample_safe_id}.csv
    mv 01g_doublet_detection_cell_scores.csv \
        01g_doublet_detection_cell_scores.${meta.sample_safe_id}.csv
    mv 01g_doublet_detection_pk_sweep.csv \
        01g_doublet_detection_pk_sweep.${meta.sample_safe_id}.csv
    mv 01g_doublet_detection_performance.csv \
        01g_doublet_detection_performance.${meta.sample_safe_id}.csv
    mv 01g_sessionInfo.txt \
        01g_sessionInfo.${meta.sample_safe_id}.txt
    """

    stub:

    """
    touch seu_doublet_scored.${meta.sample_safe_id}.rds

    cat > 01g_doublet_detection_summary.${meta.sample_safe_id}.csv <<'CSV'
    metric,value
    cells_input,0
    cells_output,0
    features_input,0
    features_output,0
    cells_removed,0
    samples,1
    proposed_doublets,0
    proposed_doublet_fraction,0
    approval_status,proposed
    CSV

    cat > 01g_doublet_detection_by_sample.${meta.sample_safe_id}.csv <<CSV
    sample,cells_scored,expected_doublet_rate,pN,pK,nExp_raw,homotypic_prop,nExp_adjusted,proposed_singlets,proposed_doublets,doublet_detection_status,doublet_decision_status,approval_status
    ${meta.sample},0,0.075,0.25,NA,0,0,0,0,0,stub,unreviewed,proposed
    CSV

    cat > 01g_doublet_detection_cell_scores.${meta.sample_safe_id}.csv <<'CSV'
    cell_id,sample,doublet_score,doublet_class,is_doublet,doublet_pK,doublet_nExp_raw,doublet_nExp_adjusted,doublet_detection_status,doublet_decision_status,doublet_reviewed_class
    CSV

    cat > 01g_doublet_detection_pk_sweep.${meta.sample_safe_id}.csv <<'CSV'
    sample,pK,BCmetric
    CSV

    cat > 01g_doublet_detection_performance.${meta.sample_safe_id}.csv <<CSV
    sample,cells,nfeatures,pcs,sweep_cores,normalize_seconds,find_variable_features_seconds,scale_data_seconds,run_pca_seconds,param_sweep_seconds,doublet_finder_seconds,total_seconds
    ${meta.sample},0,0,0,1,0,0,0,0,0,0,0
    CSV

    cat > 01g_sessionInfo.${meta.sample_safe_id}.txt <<'TXT'
    stub=true
    process=DOUBLET_DETECTION_SCORE_SAMPLE
    TXT
    """
}

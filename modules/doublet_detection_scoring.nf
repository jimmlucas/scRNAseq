nextflow.enable.dsl=2


process DOUBLET_DETECTION_SCORING {

    tag "${meta.project_id ?: 'project'}:${meta.biological_run_id ?: 'run'}"

    label 'very_high_memory'

    publishDir(
        "${params.outdir}/01g_doublet_detection_scoring",
        mode: 'copy',
        overwrite: true
    )

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        "docker://${params.r_image}" :
        params.r_image}"

    input:

    tuple val(meta), path(filtered_object)
    tuple val(config_meta), path(config_file)
    path doublet_detection_script

    output:

    tuple val(output_meta), path("seu_doublet_scored.rds"),
        emit: scored_object

    tuple val(summary_meta), path("01g_doublet_detection_summary.csv"),
        emit: summary

    path "01g_doublet_detection_by_sample.csv",
        emit: by_sample

    path "01g_doublet_detection_cell_scores.csv",
        emit: cell_scores

    path "01g_doublet_detection_pk_sweep.csv",
        emit: pk_sweep

    path "01g_doublet_detection_performance.csv",
        emit: performance

    path "01g_sessionInfo.txt",
        emit: session_info

    script:

    if (
        meta.project_id != config_meta.project_id
        || meta.biological_run_id != config_meta.biological_run_id
    ) {
        error(
            "DOUBLET_DETECTION_SCORING received incompatible metadata " +
            "between QC-filtered object and configuration."
        )
    }

    output_meta = meta + [
        stage         : 'doublet_detection',
        artifact_type : 'doublet_scored_seurat',
        decision_status: 'unreviewed'
    ]

    summary_meta = meta + [
        stage         : 'doublet_detection',
        artifact_type : 'metrics',
        decision_status: 'proposed'
    ]

    """
    Rscript ${doublet_detection_script} \
        --input ${filtered_object} \
        --config ${config_file} \
        --cores ${task.cpus} \
        --output seu_doublet_scored.rds
    """

    stub:

    output_meta = meta + [
        stage         : 'doublet_detection',
        artifact_type : 'doublet_scored_seurat',
        decision_status: 'unreviewed'
    ]

    summary_meta = meta + [
        stage         : 'doublet_detection',
        artifact_type : 'metrics',
        decision_status: 'proposed'
    ]

    """
    touch seu_doublet_scored.rds

    cat > 01g_doublet_detection_summary.csv <<'CSV'
    metric,value
    cells_input,0
    cells_output,0
    features_input,0
    features_output,0
    cells_removed,0
    samples,0
    proposed_doublets,0
    proposed_doublet_fraction,0
    approval_status,proposed
    CSV

    cat > 01g_doublet_detection_by_sample.csv <<'CSV'
    sample,cells_scored,expected_doublet_rate,pN,pK,pk_mode,pk_source,pk_stable,pk_delta,selected_bcmetric,sweep_cells,sweep_fraction,sampling_method,sweep_attempts,nExp_raw,homotypic_prop,nExp_adjusted,proposed_singlets,proposed_doublets,doublet_detection_status,doublet_decision_status,approval_status
    CSV

    cat > 01g_doublet_detection_cell_scores.csv <<'CSV'
    cell_id,sample,doublet_score,doublet_class,is_doublet,doublet_pK,doublet_nExp_raw,doublet_nExp_adjusted,doublet_detection_status,doublet_decision_status,doublet_reviewed_class
    CSV

    cat > 01g_doublet_detection_pk_sweep.csv <<'CSV'
    sample,pk_mode,pK_numeric,BCmetric_numeric,sweep_iteration,sweep_seed,sweep_cells,sampling_method,selected,selected_final
    CSV

    cat > 01g_doublet_detection_performance.csv <<'CSV'
    sample,cells,nfeatures,pcs,sweep_cores,pk_mode,pk_source,pk_stable,pk_delta,selected_bcmetric,sweep_cells,sweep_fraction,sampling_method,sweep_attempts,normalize_seconds,find_variable_features_seconds,scale_data_seconds,run_pca_seconds,param_sweep_seconds,doublet_finder_seconds,total_seconds
    CSV

    cat > 01g_sessionInfo.txt <<'TXT'
    stub=true
    process=DOUBLET_DETECTION_SCORING
    TXT
    """
}

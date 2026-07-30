nextflow.enable.dsl=2


process SPLIT_SEURAT_BY_SAMPLE {

    tag "${meta.project_id ?: 'project'}:${meta.biological_run_id ?: 'run'}"

    label 'high_memory'

    publishDir(
        "${params.outdir}/01f_split_seurat_by_sample",
        mode: 'copy',
        overwrite: true
    )

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.r_image.docker}" :
        params.r_image}"

    input:

    tuple val(meta), path(input_object)
    path split_script

    output:

    tuple val(sample_objects_meta), path("sample_objects"),
        emit: sample_objects

    tuple val(manifest_meta), path("01f_split_by_sample_manifest.tsv"),
        emit: manifest

    tuple val(summary_meta), path("01f_split_by_sample_summary.csv"),
        emit: summary

    path "01f_sessionInfo.txt",
        emit: session_info

    script:

    sample_objects_meta = meta + [
        stage         : 'split_by_sample',
        artifact_type : 'sample_seurat_objects'
    ]

    manifest_meta = meta + [
        stage         : 'split_by_sample',
        artifact_type : 'sample_manifest'
    ]

    summary_meta = meta + [
        stage         : 'split_by_sample',
        artifact_type : 'split_summary'
    ]

    """
    Rscript ${split_script} \
        --input ${input_object} \
        --output-dir sample_objects \
        --manifest 01f_split_by_sample_manifest.tsv \
        --summary 01f_split_by_sample_summary.csv
    """

    stub:

    sample_objects_meta = meta + [
        stage         : 'split_by_sample',
        artifact_type : 'sample_seurat_objects'
    ]

    manifest_meta = meta + [
        stage         : 'split_by_sample',
        artifact_type : 'sample_manifest'
    ]

    summary_meta = meta + [
        stage         : 'split_by_sample',
        artifact_type : 'split_summary'
    ]

    """
    mkdir -p sample_objects/stub_sample
    touch sample_objects/stub_sample/seu_qc_filtered.stub_sample.rds

    cat > 01f_split_by_sample_manifest.tsv <<'TSV'
    sample	sample_safe_id	sample_column	cells	features	object_path
    stub_sample	stub_sample	sample	0	0	sample_objects/stub_sample/seu_qc_filtered.stub_sample.rds
    TSV

    cat > 01f_split_by_sample_summary.csv <<'CSV'
    metric,value
    cells_input,0
    features_input,0
    samples,1
    cells_output,0
    features_output,0
    sample_column,sample
    CSV

    cat > 01f_sessionInfo.txt <<'TXT'
    stub=true
    process=SPLIT_SEURAT_BY_SAMPLE
    TXT
    """
}

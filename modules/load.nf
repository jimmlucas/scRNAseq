process LOAD_DATA {

    tag "load_data"
    label "high_memory"

    publishDir "${params.outdir}/00_load_data", mode: 'copy'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.r_image.docker}" :
        params.r_image}"

    input:
    val meta
    path config_file
    path load_script

    output:
    tuple val(output_meta), path("seu.rds"), emit: seurat_object
    path "00_load_data_summary.csv", emit: summary
    path "00_sessionInfo.txt", emit: session_info

    script:
    output_meta = meta + [
        stage         : 'load_data',
        artifact_type : 'loaded_seurat'
    ]

    """
    Rscript ${load_script} \
        --config ${config_file} \
        --project_dir ${projectDir}
    """
}
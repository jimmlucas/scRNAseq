/*
 * Escape a value so it can be safely written inside a
 * double-quoted YAML scalar.
 */
def yaml_escape(value) {

    return (value ?: "")
        .toString()
        .replace("\\", "\\\\")
        .replace("\"", "\\\"")
        .replace("\r", "\\r")
        .replace("\n", "\\n")
}


process RUN_MANIFEST {

    tag "run_manifest"
    label "lightweight"

    publishDir(
        params.pipeline_info_dir
            ? "${params.pipeline_info_dir}/run_metadata"
            : "${params.outdir}/run_metadata",
        mode: 'copy',
        overwrite: true
    )

    output:
    path "run_manifest.yaml", emit: manifest

    script:

    def project_id =
        params.project_id ?:
        params.project_name ?:
        "default_project"

    def project_name =
        params.project_name ?:
        project_id

    def biological_run_id =
        params.biological_run_id ?: ""

    def execution_id =
        params.execution_id ?: "manual_execution"

    def parent_execution_id =
        params.parent_execution_id ?: ""

    def branch_id =
        params.branch_id ?: ""

    def objective =
        params.objective ?: ""

    def execution_mode =
        params.mode ?: ""

    def start_from =
        params.start_from ?: execution_mode

    def input_artifact_stage =
        params.input_artifact_stage ?: ""

    def input_artifact_type =
        params.input_artifact_type ?: ""

    def input_artifact_path =
        params.input_artifact ?:
        params.input_object ?:
        ""

    def expected_output_stage =
        params.expected_output_stage ?: ""

    def expected_output_type =
        params.expected_output_type ?: ""

    def expected_output_path =
        params.expected_output_path ?: ""

    if (execution_mode == "load_phase") {

        expected_output_stage =
            expected_output_stage ?: "load_data"

        expected_output_type =
            expected_output_type ?: "loaded_seurat"

        expected_output_path =
            expected_output_path ?:
            "${params.outdir}/00_load_data/seu.rds"

    } else if (execution_mode == "qc_filter_phase") {

        input_artifact_stage =
            input_artifact_stage ?: "load_data"

        input_artifact_type =
            input_artifact_type ?: "loaded_seurat"

        expected_output_stage =
            expected_output_stage ?: "qc_filter"

        expected_output_type =
            expected_output_type ?: "filtered_seurat"

        expected_output_path =
            expected_output_path ?:
            "${params.outdir}/01_qc_filter/seu_qc_filtered.rds"

    } else if (execution_mode == "qc_visualization_phase") {

        input_artifact_stage =
            input_artifact_stage ?: "qc_filter"

        input_artifact_type =
            input_artifact_type ?: "filtered_seurat"

        expected_output_stage =
            expected_output_stage ?: "qc_visualization"

        expected_output_type =
            expected_output_type ?: "qc_report_bundle"

        expected_output_path =
            expected_output_path ?:
            "${params.outdir}/01b_qc_visualization"

    } else if (execution_mode == "qc_phase_1") {

        expected_output_stage =
            expected_output_stage ?: "qc_filter"

        expected_output_type =
            expected_output_type ?: "filtered_seurat"

        expected_output_path =
            expected_output_path ?:
            "${params.outdir}/01_qc_filter/seu_qc_filtered.rds"

    } else if (execution_mode == "doublet_phase") {

        input_artifact_stage =
            input_artifact_stage ?: "qc_filter"

        input_artifact_type =
            input_artifact_type ?: "filtered_seurat"

        expected_output_stage =
            expected_output_stage ?: "doublet_detection"

        expected_output_type =
            expected_output_type ?: "singlet_seurat"

        expected_output_path =
            expected_output_path ?:
            "${params.outdir}/01c_doublet_detection/seu_singlets.rds"
    
    } else if (execution_mode == "qc_scores_phase") {

        input_artifact_stage =
            input_artifact_stage ?: "doublet_detection"

        input_artifact_type =
            input_artifact_type ?: "singlet_seurat"

        expected_output_stage =
            expected_output_stage ?: "qc_scores"

        expected_output_type =
            expected_output_type ?: "scored_seurat"

        expected_output_path =
            expected_output_path ?:
            "${params.outdir}/01d_qc_scores/seu_qc_scored.rds"
        
    } else if (execution_mode == "qc_phase_2") {

        input_artifact_stage =
            input_artifact_stage ?: "qc_filter"

        input_artifact_type =
            input_artifact_type ?: "filtered_seurat"

        expected_output_stage =
            expected_output_stage ?: "qc_scores"

        expected_output_type =
            expected_output_type ?: "scored_seurat"

        expected_output_path =
            expected_output_path ?:
            "${params.outdir}/01d_qc_scores/seu_qc_scored.rds"
    } else if (execution_mode == "qc_process") {

        expected_output_stage =
            expected_output_stage ?: "qc_scores"

        expected_output_type =
            expected_output_type ?: "scored_seurat"

        expected_output_path =
            expected_output_path ?:
            "${params.outdir}/01d_qc_scores/seu_qc_scored.rds"
    }

    def legacy_config_file =
        params.config_file ?: ""

    def load_config =
        params.load_config ?: ""

    def qc_filter_config =
        params.qc_filter_config ?: ""

    def qc_visualization_config =
        params.qc_visualization_config ?: ""

    def doublet_config =
        params.doublet_config ?: ""

    def qc_scores_config =
        params.qc_scores_config ?: ""

    def pipeline_info_dir =
        params.pipeline_info_dir ?: ""

    def yaml_project_id =
        yaml_escape(project_id)

    def yaml_project_name =
        yaml_escape(project_name)

    def yaml_biological_run_id =
        yaml_escape(biological_run_id)

    def yaml_execution_id =
        yaml_escape(execution_id)

    def yaml_parent_execution_id =
        yaml_escape(parent_execution_id)

    def yaml_execution_mode =
        yaml_escape(execution_mode)

    def yaml_branch_id =
        yaml_escape(branch_id)

    def yaml_objective =
        yaml_escape(objective)

    def yaml_start_from =
        yaml_escape(start_from)

    def yaml_input_artifact_stage =
        yaml_escape(input_artifact_stage)

    def yaml_input_artifact_type =
        yaml_escape(input_artifact_type)

    def yaml_input_artifact_path =
        yaml_escape(input_artifact_path)

    def yaml_expected_output_stage =
        yaml_escape(expected_output_stage)

    def yaml_expected_output_type =
        yaml_escape(expected_output_type)

    def yaml_expected_output_path =
        yaml_escape(expected_output_path)

    def yaml_pipeline_root =
        yaml_escape(projectDir)

    def yaml_outdir =
        yaml_escape(params.outdir)

    def yaml_pipeline_info_dir =
        yaml_escape(pipeline_info_dir)

    def yaml_container_image =
        yaml_escape(params.r_image)

    def yaml_profile =
        yaml_escape(workflow.profile)

    def yaml_work_dir =
        yaml_escape(workflow.workDir)

    def yaml_launch_dir =
        yaml_escape(workflow.launchDir)

    def yaml_run_name =
        yaml_escape(workflow.runName)

    def yaml_nextflow_version =
        yaml_escape(workflow.nextflow.version)

    def yaml_legacy_config_file =
        yaml_escape(legacy_config_file)

    def yaml_load_config =
        yaml_escape(load_config)

    def yaml_qc_filter_config =
        yaml_escape(qc_filter_config)

    def yaml_qc_visualization_config =
        yaml_escape(qc_visualization_config)

    def yaml_doublet_config =
        yaml_escape(doublet_config)

    def yaml_qc_scores_config =
        yaml_escape(qc_scores_config)

    """
    cat > run_manifest.yaml <<'YAML'
    schema_version: "1.2"

    manifest:
      kind: "execution_plan"
      status: "planned"

    scientific_context:
      project_id: "${yaml_project_id}"
      project_name: "${yaml_project_name}"
      biological_run_id: "${yaml_biological_run_id}"
      scope: "run"

    execution:
      id: "${yaml_execution_id}"
      parent_execution_id: "${yaml_parent_execution_id}"
      mode: "${yaml_execution_mode}"
      branch_id: "${yaml_branch_id}"
      objective: "${yaml_objective}"
      start_from: "${yaml_start_from}"

    input_artifact:
      stage: "${yaml_input_artifact_stage}"
      type: "${yaml_input_artifact_type}"
      path: "${yaml_input_artifact_path}"

    expected_output_artifact:
      stage: "${yaml_expected_output_stage}"
      type: "${yaml_expected_output_type}"
      path: "${yaml_expected_output_path}"

    runtime:
      pipeline_root: "${yaml_pipeline_root}"
      outdir: "${yaml_outdir}"
      pipeline_info_dir: "${yaml_pipeline_info_dir}"
      container_image: "${yaml_container_image}"
      profile: "${yaml_profile}"
      work_dir: "${yaml_work_dir}"
      launch_dir: "${yaml_launch_dir}"
      run_name: "${yaml_run_name}"
      nextflow_version: "${yaml_nextflow_version}"

    runtime_configs:
      legacy_config_file: "${yaml_legacy_config_file}"
      load_config: "${yaml_load_config}"
      qc_filter_config: "${yaml_qc_filter_config}"
      qc_visualization_config: "${yaml_qc_visualization_config}"
      doublet_config: "${yaml_doublet_config}"
      qc_scores_config: "${yaml_qc_scores_config}"

    notes:
      generated_by: "Nextflow RUN_MANIFEST process"
      description: "Execution intent generated before scientific processes complete."
    YAML
    """
}
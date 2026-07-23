nextflow.enable.dsl=2

include { QC_PHASE_1 }             from './subworkflow/qc_phase_1'
include { QC_PHASE_2 }             from './subworkflow/qc_phase_2'
include { qc_process }             from './subworkflow/qc_process'
include { LOAD_PHASE }             from './subworkflow/load_phase'
include { BBUDY_QC_PROCESS }       from './subworkflow/bbudy/bbudy_qc_process'
include { QC_FILTER_PHASE }        from './subworkflow/qc_filter_phase'
include { QC_VISUALIZATION_PHASE } from './subworkflow/qc_visualization_phase'
include { DOUBLET_PHASE }          from './subworkflow/doublet_phase'
include { QC_SCORES_PHASE }        from './subworkflow/qc_scores_phase'
include { RUN_MANIFEST }           from './modules/run_manifest'


workflow {

    /*
     * Write the execution manifest independently from the
     * scientific processes.
     */
    RUN_MANIFEST()

    /*
     * Stable scientific identity.
     *
     * execution_id, parent_execution_id and branch_id are deliberately
     * excluded so that administrative retries do not invalidate the
     * cache of scientific processes.
     */
    scientific_meta = [
        project_id        : params.project_id ?: '',
        project_name      : params.project_name ?: '',
        biological_run_id : params.biological_run_id ?: '',
        scope             : 'run'
    ]

    scientific_meta_ch = Channel.value(scientific_meta)

    /*
     * MODE: bbudy_qc_process
     *
     * New plan-driven QC orchestration.
     *
     * The analysis plan defines the experiment structure,
     * enabled modules and decisions requiring approval.
     */
    
    if (params.mode == 'bbudy_qc_process') {

        if (!params.analysis_plan) {
            error """
            The bbudy_qc_process mode requires:

            --analysis_plan <path_to_analysis_plan.yaml>
            """.stripIndent()
        }

        analysis_plan_ch = Channel.value(
            file(
                params.analysis_plan,
                checkIfExists: true
            )
        )

        BBUDY_QC_PROCESS(
            scientific_meta_ch,
            analysis_plan_ch
        )


    /*
     * MODE: load_phase
     */
    } else if (params.mode == 'load_phase') {
        
        if (!params.load_config) {
            error """
            The load_phase mode requires:

            --load_config <path_to_load_config.yaml>
            """.stripIndent()
        }

        load_config_ch = Channel.value(
            file(
                params.load_config,
                checkIfExists: true
            )
        )

        LOAD_PHASE(
            scientific_meta_ch,
            load_config_ch
        )

    /*
     * MODE: qc_filter_phase
     *
     * Start from an existing loaded Seurat object and execute only:
     *
     * QC_FILTER
     */
    } else if (params.mode == 'qc_filter_phase') {

        if (!params.input_object) {
            error """
            The qc_filter_phase mode requires:

            --input_object <path_to_loaded_object.rds>
            """.stripIndent()
        }

        if (!params.qc_filter_config) {
            error """
            The qc_filter_phase mode requires:

            --qc_filter_config <path_to_qc_filter_config.yaml>
            """.stripIndent()
        }

        loaded_object_ch = Channel.value(
            file(
                params.input_object,
                checkIfExists: true
            )
        )

        qc_filter_config_ch = Channel.value(
            file(
                params.qc_filter_config,
                checkIfExists: true
            )
        )

        QC_FILTER_PHASE(
            scientific_meta_ch,
            loaded_object_ch,
            qc_filter_config_ch
        )


    /*
     * MODE: qc_phase_1
     *
     * Execute:
     *
     * LOAD_DATA
     *     ↓
     * QC_FILTER
     *
     * LOAD_DATA receives the stable scientific metadata.
     * Each phase retains its own configuration file.
     */
    } else if (params.mode == 'qc_phase_1') {

        if (!params.load_config) {
            error """
            The qc_phase_1 mode requires:

            --load_config <path_to_load_config.yaml>
            """.stripIndent()
        }

        if (!params.qc_filter_config) {
            error """
            The qc_phase_1 mode requires:

            --qc_filter_config <path_to_qc_filter_config.yaml>
            """.stripIndent()
        }

        load_config_ch = Channel.value(
            file(
                params.load_config,
                checkIfExists: true
            )
        )

        qc_filter_config_ch = Channel.value(
            file(
                params.qc_filter_config,
                checkIfExists: true
            )
        )

        QC_PHASE_1(
            scientific_meta_ch,
            load_config_ch,
            qc_filter_config_ch
        )


    /*
    * MODE: qc_phase_2
    *
    * Start from an existing filtered Seurat artifact and execute:
    *
    * QC_VISUALIZATION
    * QC_DOUBLET_DETECTION
    * QC_SCORES
    *
    * All scientific processes use tuple(meta, path) artifact contracts.
    */
    
    } else if (params.mode == 'qc_phase_2') {

        if (!params.input_object) {
            error """
            The qc_phase_2 mode requires:

            --input_object <path_to_filtered_object.rds>
            """.stripIndent()
        }

        if (!params.config_file) {
            error """
            The qc_phase_2 mode requires:

            --config_file <path_to_config.yaml>
            """.stripIndent()
        }

        filtered_object_ch = Channel.value(
            file(
                params.input_object,
                checkIfExists: true
            )
        )

        config_ch = Channel.value(
            file(
                params.config_file,
                checkIfExists: true
            )
        )

        QC_PHASE_2(
            scientific_meta_ch,
            filtered_object_ch,
            config_ch
        )


    /*
     * MODE: qc_visualization_phase
     *
     * Start from an existing filtered object and execute only:
     *
     * QC_VISUALIZATION
     */
    } else if (params.mode == 'qc_visualization_phase') {

        if (!params.input_object) {
            error """
            The qc_visualization_phase mode requires:

            --input_object <path_to_filtered_object.rds>
            """.stripIndent()
        }

        phase_config = (
            params.qc_visualization_config
                ?: params.config_file
        )

        if (!phase_config) {
            error """
            The qc_visualization_phase mode requires one of:

            --qc_visualization_config <path_to_visualization_config.yaml>

            or:

            --config_file <path_to_config.yaml>
            """.stripIndent()
        }

        filtered_object_ch = Channel.value(
            file(
                params.input_object,
                checkIfExists: true
            )
        )

        visualization_config_ch = Channel.value(
            file(
                phase_config,
                checkIfExists: true
            )
        )

        QC_VISUALIZATION_PHASE(
            scientific_meta_ch,
            filtered_object_ch,
            visualization_config_ch
        )


    /*
     * MODE: doublet_phase
     *
     * Start from an existing filtered object and execute only:
     *
     * QC_DOUBLET_DETECTION
     */
    } else if (params.mode == 'doublet_phase') {

        if (!params.input_object) {
            error """
            The doublet_phase mode requires:

            --input_object <path_to_filtered_object.rds>
            """.stripIndent()
        }

        phase_config = (
            params.doublet_config
                ?: params.config_file
        )

        if (!phase_config) {
            error """
            The doublet_phase mode requires one of:

            --doublet_config <path_to_doublet_config.yaml>

            or:

            --config_file <path_to_config.yaml>
            """.stripIndent()
        }

        filtered_object_ch = Channel.value(
            file(
                params.input_object,
                checkIfExists: true
            )
        )

        doublet_config_ch = Channel.value(
            file(
                phase_config,
                checkIfExists: true
            )
        )

        DOUBLET_PHASE(
            scientific_meta_ch,
            filtered_object_ch,
            doublet_config_ch
        )


    /*
     * MODE: qc_scores_phase
     *
     * Start from an existing singlet object and execute only:
     *
     * QC_SCORES
     */
    } else if (params.mode == 'qc_scores_phase') {

        if (!params.input_object) {
            error """
            The qc_scores_phase mode requires:

            --input_object <path_to_singlet_object.rds>
            """.stripIndent()
        }

        phase_config = (
            params.qc_scores_config
                ?: params.config_file
        )

        if (!phase_config) {
            error """
            The qc_scores_phase mode requires one of:

            --qc_scores_config <path_to_scores_config.yaml>

            or:

            --config_file <path_to_config.yaml>
            """.stripIndent()
        }

        singlet_object_ch = Channel.value(
            file(
                params.input_object,
                checkIfExists: true
            )
        )

        scores_config_ch = Channel.value(
            file(
                phase_config,
                checkIfExists: true
            )
        )

        QC_SCORES_PHASE(
            scientific_meta_ch,
            singlet_object_ch,
            scores_config_ch
        )


    /*
     * MODE: qc_process
     *
     * Execute the complete legacy workflow.
     *
     * This mode temporarily retains its internal config_file handling
     * and remains available for backward compatibility.
     */
    } else if (params.mode == 'qc_process') {

        if (!params.config_file) {
            error """
            The qc_process mode requires:

            --config_file <path_to_config.yaml>
            """.stripIndent()
        }

        config_ch = Channel.value(
            file(
                params.config_file,
                checkIfExists: true
            )
        )

        qc_process(
            scientific_meta_ch,
            config_ch
        )


        /*
        * Invalid execution mode.
        */
        } else {

            error """
            Invalid mode: ${params.mode}

            Available modes:

            - bbudy_qc_process
            - qc_process
            - load_phase
            - qc_filter_phase
            - qc_phase_1
            - qc_phase_2
            - qc_visualization_phase
            - doublet_phase
            - qc_scores_phase
            """.stripIndent()
        }
}
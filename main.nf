nextflow.enable.dsl=2

include { QC_PHASE_1 }             from './subworkflow/qc_phase_1'
include { QC_PHASE_2 }             from './subworkflow/qc_phase_2'
include { qc_process }             from './subworkflow/qc_process'
include { LOAD_PHASE }             from './subworkflow/load_phase'
include { ANALYSIS_PLAN_QC_PROCESS } from './subworkflow/control/analysis_plan_qc_process'
include { QC_THRESHOLD_APPROVAL_PHASE } from './subworkflow/control/qc_threshold_approval_phase'
include { QC_FILTER_POLICY_PHASE } from './subworkflow/control/qc_filter_policy_phase'
include { QC_VISUALIZATION_POLICY_PHASE } from './subworkflow/control/qc_visualization_policy_phase'
include { DOUBLET_DETECTION_SCORING_PHASE } from './subworkflow/control/doublet_detection_scoring_phase'
include { DOUBLET_POLICY_RESOLUTION_PHASE } from './subworkflow/control/doublet_policy_resolution_phase'
include { DOUBLET_APPROVAL_PHASE } from './subworkflow/control/doublet_approval_phase'
include { DOUBLET_FILTER_PHASE } from './subworkflow/control/doublet_filter_phase'
include { POST_DOUBLET_QC_SUMMARY_PHASE } from './subworkflow/control/post_doublet_qc_summary_phase'
include { NORMALIZATION_POLICY_RESOLUTION_PHASE } from './subworkflow/control/normalization_policy_resolution_phase'
include { NORMALIZATION_APPROVAL_PHASE } from './subworkflow/control/normalization_approval_phase'
include { NORMALIZATION_PHASE } from './subworkflow/control/normalization_phase'
include { DIMENSIONALITY_PHASE } from './subworkflow/control/dimensionality_phase'
include { PC_SELECTION_POLICY_RESOLUTION_PHASE } from './subworkflow/control/pc_selection_policy_resolution_phase'
include { PC_SELECTION_APPROVAL_PHASE } from './subworkflow/control/pc_selection_approval_phase'
include { INTEGRATION_POLICY_RESOLUTION_PHASE } from './subworkflow/control/integration_policy_resolution_phase'
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
     * MODE: analysis_plan_qc_process
     *
     * Plan-driven QC orchestration.
     *
     * The analysis plan defines the experiment structure,
     * enabled modules and decisions requiring approval.
     */

    if (params.mode == 'analysis_plan_qc_process') {

        if (!params.analysis_plan) {
            error """
            The analysis_plan_qc_process mode requires:

            --analysis_plan <path_to_analysis_plan.yaml>
            """.stripIndent()
        }

        analysis_plan_ch = Channel.value(
            file(
                params.analysis_plan,
                checkIfExists: true
            )
        )

        ANALYSIS_PLAN_QC_PROCESS(
            scientific_meta_ch,
            analysis_plan_ch
        )


    /*
     * MODE: qc_threshold_approval_phase
     *
     * Validate a human-reviewed QC threshold policy.
     */
    } else if (
        params.mode == 'qc_threshold_approval_phase'
    ) {

        if (!params.analysis_plan) {
            error """
            The qc_threshold_approval_phase mode requires:

            --analysis_plan <path_to_analysis_plan.yaml>
            """.stripIndent()
        }

        if (!params.approved_qc_threshold_policy) {
            error """
            The qc_threshold_approval_phase mode requires:

            --approved_qc_threshold_policy <path_to_reviewed_policy.yaml>
            """.stripIndent()
        }

        reviewed_plan = file(
            params.analysis_plan,
            checkIfExists: true
        )

        reviewed_policy = file(
            params.approved_qc_threshold_policy,
            checkIfExists: true
        )

        plan_artifact_ch = Channel.value(
            tuple(
                scientific_meta + [
                    stage         : 'analysis_plan',
                    artifact_type : 'analysis_plan_yaml'
                ],
                reviewed_plan
            )
        )

        reviewed_policy_artifact_ch = Channel.value(
            tuple(
                scientific_meta + [
                    stage           : 'qc_threshold_review',
                    artifact_type   : 'qc_threshold_policy',
                    approval_status : 'human_reviewed'
                ],
                reviewed_policy
            )
        )

        QC_THRESHOLD_APPROVAL_PHASE(
            plan_artifact_ch,
            reviewed_policy_artifact_ch
        )


    /*
     * MODE: qc_filter_policy_phase
     *
     * Apply a human-approved and validated QC threshold policy
     * to an existing QC-metrics Seurat artifact.
     */
    } else if (
        params.mode == 'qc_filter_policy_phase'
    ) {

        if (!params.input_object) {
            error """
            The qc_filter_policy_phase mode requires:

            --input_object <path_to_qc_metrics_object.rds>
            """.stripIndent()
        }

        if (!params.approved_qc_threshold_policy) {
            error """
            The qc_filter_policy_phase mode requires:

            --approved_qc_threshold_policy <path_to_validated_policy.yaml>
            """.stripIndent()
        }

        qc_metrics_object_ch = Channel.value(
            file(
                params.input_object,
                checkIfExists: true
            )
        )

        approved_policy_ch = Channel.value(
            file(
                params.approved_qc_threshold_policy,
                checkIfExists: true
            )
        )

        QC_FILTER_POLICY_PHASE(
            scientific_meta_ch,
            qc_metrics_object_ch,
            approved_policy_ch
        )


    /*
     * MODE: qc_visualization_policy_phase
     *
     * Generate traceable before/after QC visualizations from
     * previously published control artifacts.
     */
    } else if (
        params.mode == 'qc_visualization_policy_phase'
    ) {

        if (!params.before_object) {
            error """
            The qc_visualization_policy_phase mode requires:

            --before_object <path_to_qc_metrics_object.rds>
            """.stripIndent()
        }

        if (!params.after_object) {
            error """
            The qc_visualization_policy_phase mode requires:

            --after_object <path_to_filtered_object.rds>
            """.stripIndent()
        }

        if (!params.cell_decisions) {
            error """
            The qc_visualization_policy_phase mode requires:

            --cell_decisions <path_to_cell_decisions.csv>
            """.stripIndent()
        }

        if (!params.approved_qc_threshold_policy) {
            error """
            The qc_visualization_policy_phase mode requires:

            --approved_qc_threshold_policy <path_to_validated_policy.yaml>
            """.stripIndent()
        }

        qc_metrics_object_ch = Channel.value(
            file(
                params.before_object,
                checkIfExists: true
            )
        )

        filtered_object_ch = Channel.value(
            file(
                params.after_object,
                checkIfExists: true
            )
        )

        cell_decisions_ch = Channel.value(
            file(
                params.cell_decisions,
                checkIfExists: true
            )
        )

        approved_policy_ch = Channel.value(
            file(
                params.approved_qc_threshold_policy,
                checkIfExists: true
            )
        )

        QC_VISUALIZATION_POLICY_PHASE(
            scientific_meta_ch,
            qc_metrics_object_ch,
            filtered_object_ch,
            cell_decisions_ch,
            approved_policy_ch
        )


    /*
     * MODE: doublet_detection_scoring_phase
     *
     * Score doublet evidence from a QC-filtered object without removing
     * cells. Scientific filtering remains a separate approved policy step.
     */
    } else if (
        params.mode == 'doublet_detection_scoring_phase'
    ) {

        if (!params.input_object) {
            error """
            The doublet_detection_scoring_phase mode requires:

            --input_object <path_to_qc_filtered_object.rds>
            """.stripIndent()
        }

        if (!params.doublet_config) {
            error """
            The doublet_detection_scoring_phase mode requires:

            --doublet_config <path_to_doublet_config.yaml>
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
                params.doublet_config,
                checkIfExists: true
            )
        )

        DOUBLET_DETECTION_SCORING_PHASE(
            scientific_meta_ch,
            filtered_object_ch,
            doublet_config_ch
        )


    /*
     * MODE: doublet_policy_resolution_phase
     *
     * Resolve a proposed human-reviewable doublet policy from
     * technical doublet scoring outputs. This mode does not modify
     * Seurat objects and does not approve filtering decisions.
     */
    } else if (
        params.mode == 'doublet_policy_resolution_phase'
    ) {

        if (!params.doublet_summary) {
            error """
            The doublet_policy_resolution_phase mode requires:

            --doublet_summary <path_to_01g_doublet_detection_summary.csv>
            """.stripIndent()
        }

        if (!params.doublet_by_sample) {
            error """
            The doublet_policy_resolution_phase mode requires:

            --doublet_by_sample <path_to_01g_doublet_detection_by_sample.csv>
            """.stripIndent()
        }

        if (!params.doublet_cell_scores) {
            error """
            The doublet_policy_resolution_phase mode requires:

            --doublet_cell_scores <path_to_01g_doublet_detection_cell_scores.csv>
            """.stripIndent()
        }

        doublet_summary_ch = Channel.value(
            file(
                params.doublet_summary,
                checkIfExists: true
            )
        )

        doublet_by_sample_ch = Channel.value(
            file(
                params.doublet_by_sample,
                checkIfExists: true
            )
        )

        doublet_cell_scores_ch = Channel.value(
            file(
                params.doublet_cell_scores,
                checkIfExists: true
            )
        )

        DOUBLET_POLICY_RESOLUTION_PHASE(
            scientific_meta_ch,
            doublet_summary_ch,
            doublet_by_sample_ch,
            doublet_cell_scores_ch
        )


    /*
     * MODE: doublet_approval_phase
     *
     * Validate a human-reviewed doublet policy before it can be
     * consumed by a future DOUBLET_FILTER step.
     */
    } else if (
        params.mode == 'doublet_approval_phase'
    ) {

        if (!params.analysis_plan) {
            error """
            The doublet_approval_phase mode requires:

            --analysis_plan <path_to_analysis_plan.yaml>
            """.stripIndent()
        }

        if (!params.approved_doublet_policy) {
            error """
            The doublet_approval_phase mode requires:

            --approved_doublet_policy <path_to_reviewed_doublet_policy.yaml>
            """.stripIndent()
        }

        reviewed_plan = file(
            params.analysis_plan,
            checkIfExists: true
        )

        reviewed_policy = file(
            params.approved_doublet_policy,
            checkIfExists: true
        )

        plan_artifact_ch = Channel.value(
            tuple(
                scientific_meta + [
                    stage         : 'analysis_plan',
                    artifact_type : 'analysis_plan_yaml'
                ],
                reviewed_plan
            )
        )

        reviewed_policy_artifact_ch = Channel.value(
            tuple(
                scientific_meta + [
                    stage           : 'doublet_policy_review',
                    artifact_type   : 'doublet_detection_policy',
                    approval_status : 'human_reviewed'
                ],
                reviewed_policy
            )
        )

        DOUBLET_APPROVAL_PHASE(
            plan_artifact_ch,
            reviewed_policy_artifact_ch
        )


    /*
     * MODE: doublet_filter_phase
     *
     * Apply a human-approved and technically validated doublet
     * policy to an existing doublet-scored Seurat artifact.
     */
    } else if (
        params.mode == 'doublet_filter_phase'
    ) {

        if (!params.input_object) {
            error """
            The doublet_filter_phase mode requires:

            --input_object <path_to_seu_doublet_scored.rds>
            """.stripIndent()
        }

        if (!params.approved_doublet_policy) {
            error """
            The doublet_filter_phase mode requires:

            --approved_doublet_policy <path_to_approved_doublet_policy.yaml>
            """.stripIndent()
        }

        scored_object_ch = Channel.value(
            file(
                params.input_object,
                checkIfExists: true
            )
        )

        approved_policy_ch = Channel.value(
            file(
                params.approved_doublet_policy,
                checkIfExists: true
            )
        )

        DOUBLET_FILTER_PHASE(
            scientific_meta_ch,
            scored_object_ch,
            approved_policy_ch
        )


    /*
     * MODE: post_doublet_qc_summary_phase
     *
     * Produce a non-destructive readiness summary from the canonical
     * post-QC, post-doublet Seurat artifact.
     */
    } else if (
        params.mode == 'post_doublet_qc_summary_phase'
    ) {

        if (!params.input_object) {
            error """
            The post_doublet_qc_summary_phase mode requires:

            --input_object <path_to_seu_doublet_filtered.rds>
            """.stripIndent()
        }

        if (!params.applied_doublet_policy) {
            error """
            The post_doublet_qc_summary_phase mode requires:

            --applied_doublet_policy <path_to_applied_doublet_policy.yaml>
            """.stripIndent()
        }

        filtered_object_ch = Channel.value(
            file(
                params.input_object,
                checkIfExists: true
            )
        )

        applied_policy_ch = Channel.value(
            file(
                params.applied_doublet_policy,
                checkIfExists: true
            )
        )

        POST_DOUBLET_QC_SUMMARY_PHASE(
            scientific_meta_ch,
            filtered_object_ch,
            applied_policy_ch
        )


    /*
     * MODE: normalization_policy_resolution_phase
     *
     * Resolve a proposed normalization policy from the
     * post-doublet analysis readiness contract. This mode does not
     * normalize or modify the Seurat object.
     */
    } else if (
        params.mode == 'normalization_policy_resolution_phase'
    ) {

        if (!params.analysis_readiness) {
            error """
            The normalization_policy_resolution_phase mode requires:

            --analysis_readiness <path_to_03_analysis_readiness.yaml>
            """.stripIndent()
        }

        if (!params.post_doublet_qc_summary) {
            error """
            The normalization_policy_resolution_phase mode requires:

            --post_doublet_qc_summary <path_to_03_post_doublet_qc_summary.csv>
            """.stripIndent()
        }

        if (!params.post_doublet_qc_by_sample) {
            error """
            The normalization_policy_resolution_phase mode requires:

            --post_doublet_qc_by_sample <path_to_03_post_doublet_qc_by_sample.csv>
            """.stripIndent()
        }

        if (!params.post_doublet_qc_warnings) {
            error """
            The normalization_policy_resolution_phase mode requires:

            --post_doublet_qc_warnings <path_to_03_post_doublet_qc_warnings.tsv>
            """.stripIndent()
        }

        readiness_ch = Channel.value(
            file(
                params.analysis_readiness,
                checkIfExists: true
            )
        )

        post_doublet_summary_ch = Channel.value(
            file(
                params.post_doublet_qc_summary,
                checkIfExists: true
            )
        )

        post_doublet_by_sample_ch = Channel.value(
            file(
                params.post_doublet_qc_by_sample,
                checkIfExists: true
            )
        )

        post_doublet_warnings_ch = Channel.value(
            file(
                params.post_doublet_qc_warnings,
                checkIfExists: true
            )
        )

        NORMALIZATION_POLICY_RESOLUTION_PHASE(
            scientific_meta_ch,
            readiness_ch,
            post_doublet_summary_ch,
            post_doublet_by_sample_ch,
            post_doublet_warnings_ch
        )


    /*
     * MODE: normalization_approval_phase
     *
     * Validate a human-reviewed normalization policy before any
     * normalization process can transform the object.
     */
    } else if (
        params.mode == 'normalization_approval_phase'
    ) {

        if (!params.analysis_plan) {
            error """
            The normalization_approval_phase mode requires:

            --analysis_plan <path_to_analysis_plan.yaml>
            """.stripIndent()
        }

        if (!params.approved_normalization_policy) {
            error """
            The normalization_approval_phase mode requires:

            --approved_normalization_policy <path_to_reviewed_normalization_policy.yaml>
            """.stripIndent()
        }

        reviewed_plan = file(
            params.analysis_plan,
            checkIfExists: true
        )

        reviewed_policy = file(
            params.approved_normalization_policy,
            checkIfExists: true
        )

        plan_artifact_ch = Channel.value(
            tuple(
                scientific_meta + [
                    stage         : 'analysis_plan',
                    artifact_type : 'analysis_plan_yaml'
                ],
                reviewed_plan
            )
        )

        reviewed_policy_artifact_ch = Channel.value(
            tuple(
                scientific_meta + [
                    stage           : 'normalization_policy_review',
                    artifact_type   : 'normalization_policy',
                    approval_status : 'human_reviewed'
                ],
                reviewed_policy
            )
        )

        NORMALIZATION_APPROVAL_PHASE(
            plan_artifact_ch,
            reviewed_policy_artifact_ch
        )


    /*
     * MODE: normalization_phase
     *
     * Apply a human-approved and validated normalization policy
     * to the post-doublet Seurat artifact.
     */
    } else if (
        params.mode == 'normalization_phase'
    ) {

        if (!params.input_object) {
            error """
            The normalization_phase mode requires:

            --input_object <path_to_seu_doublet_filtered.rds>
            """.stripIndent()
        }

        if (!params.approved_normalization_policy) {
            error """
            The normalization_phase mode requires:

            --approved_normalization_policy <path_to_approved_normalization_policy.yaml>
            """.stripIndent()
        }

        filtered_object_ch = Channel.value(
            file(
                params.input_object,
                checkIfExists: true
            )
        )

        approved_policy_ch = Channel.value(
            file(
                params.approved_normalization_policy,
                checkIfExists: true
            )
        )

        NORMALIZATION_PHASE(
            scientific_meta_ch,
            filtered_object_ch,
            approved_policy_ch
        )


    /*
     * MODE: dimensionality_phase
     *
     * Apply deterministic technical dimensionality preparation:
     * ScaleData over variable features and RunPCA from the normalized
     * Seurat artifact. This phase publishes evidence for later
     * integration decisions but does not decide integration.
     */
    } else if (
        params.mode == 'dimensionality_phase'
    ) {

        if (!params.input_object) {
            error """
            The dimensionality_phase mode requires:

            --input_object <path_to_seu_normalized.rds>
            """.stripIndent()
        }

        if (!params.applied_normalization_policy) {
            error """
            The dimensionality_phase mode requires:

            --applied_normalization_policy <path_to_applied_normalization_policy.yaml>
            """.stripIndent()
        }

        normalized_object_ch = Channel.value(
            file(
                params.input_object,
                checkIfExists: true
            )
        )

        applied_policy_ch = Channel.value(
            file(
                params.applied_normalization_policy,
                checkIfExists: true
            )
        )

        DIMENSIONALITY_PHASE(
            scientific_meta_ch,
            normalized_object_ch,
            applied_policy_ch
        )


    /*
     * MODE: pc_selection_policy_resolution_phase
     *
     * Resolve a proposed human-reviewable PC-selection policy from
     * published PCA metrics. This phase does not run neighbors, UMAP,
     * clustering or integration.
     */
    } else if (
        params.mode == 'pc_selection_policy_resolution_phase'
    ) {

        if (!params.dimensionality_summary) {
            error """
            The pc_selection_policy_resolution_phase mode requires:

            --dimensionality_summary <path_to_05_dimensionality_summary.csv>
            """.stripIndent()
        }

        if (!params.variable_features_summary) {
            error """
            The pc_selection_policy_resolution_phase mode requires:

            --variable_features_summary <path_to_05_variable_features_summary.csv>
            """.stripIndent()
        }

        if (!params.pca_elbow_metrics) {
            error """
            The pc_selection_policy_resolution_phase mode requires:

            --pca_elbow_metrics <path_to_05_pca_elbow_metrics.csv>
            """.stripIndent()
        }

        dimensionality_summary_ch = Channel.value(
            file(
                params.dimensionality_summary,
                checkIfExists: true
            )
        )

        variable_features_ch = Channel.value(
            file(
                params.variable_features_summary,
                checkIfExists: true
            )
        )

        elbow_metrics_ch = Channel.value(
            file(
                params.pca_elbow_metrics,
                checkIfExists: true
            )
        )

        PC_SELECTION_POLICY_RESOLUTION_PHASE(
            scientific_meta_ch,
            dimensionality_summary_ch,
            variable_features_ch,
            elbow_metrics_ch
        )


    /*
     * MODE: pc_selection_approval_phase
     *
     * Validate a human-reviewed PC-selection policy. This phase does
     * not run neighbors, UMAP, clustering or integration.
     */
    } else if (
        params.mode == 'pc_selection_approval_phase'
    ) {

        if (!params.analysis_plan) {
            error """
            The pc_selection_approval_phase mode requires:

            --analysis_plan <path_to_analysis_plan.yaml>
            """.stripIndent()
        }

        if (!params.approved_pc_selection_policy) {
            error """
            The pc_selection_approval_phase mode requires:

            --approved_pc_selection_policy <path_to_reviewed_pc_selection_policy.yaml>
            """.stripIndent()
        }

        reviewed_plan = file(
            params.analysis_plan,
            checkIfExists: true
        )

        reviewed_policy = file(
            params.approved_pc_selection_policy,
            checkIfExists: true
        )

        plan_artifact_ch = Channel.value(
            tuple(
                scientific_meta + [
                    stage         : 'analysis_plan',
                    artifact_type : 'analysis_plan_yaml'
                ],
                reviewed_plan
            )
        )

        reviewed_policy_artifact_ch = Channel.value(
            tuple(
                scientific_meta + [
                    stage           : 'pc_selection_policy_review',
                    artifact_type   : 'pc_selection_policy',
                    approval_status : 'human_reviewed'
                ],
                reviewed_policy
            )
        )

        PC_SELECTION_APPROVAL_PHASE(
            plan_artifact_ch,
            reviewed_policy_artifact_ch
        )


    /*
     * MODE: integration_policy_resolution_phase
     *
     * Resolve a proposed human-reviewable integration policy from
     * PCA evidence, approved PC selection, readiness and applied
     * normalization. This phase does not execute Harmony, Seurat
     * integration, neighbors, UMAP or clustering.
     */
    } else if (
        params.mode == 'integration_policy_resolution_phase'
    ) {

        if (!params.input_object) {
            error """
            The integration_policy_resolution_phase mode requires:

            --input_object <path_to_seu_pca.rds>
            """.stripIndent()
        }

        if (!params.dimensionality_summary) {
            error """
            The integration_policy_resolution_phase mode requires:

            --dimensionality_summary <path_to_05_dimensionality_summary.csv>
            """.stripIndent()
        }

        if (!params.approved_pc_selection_policy) {
            error """
            The integration_policy_resolution_phase mode requires:

            --approved_pc_selection_policy <path_to_approved_pc_selection_policy.yaml>
            """.stripIndent()
        }

        if (!params.analysis_readiness) {
            error """
            The integration_policy_resolution_phase mode requires:

            --analysis_readiness <path_to_03_analysis_readiness.yaml>
            """.stripIndent()
        }

        if (!params.applied_normalization_policy) {
            error """
            The integration_policy_resolution_phase mode requires:

            --applied_normalization_policy <path_to_applied_normalization_policy.yaml>
            """.stripIndent()
        }

        pca_object_ch = Channel.value(
            file(
                params.input_object,
                checkIfExists: true
            )
        )

        dimensionality_summary_ch = Channel.value(
            file(
                params.dimensionality_summary,
                checkIfExists: true
            )
        )

        approved_pc_selection_policy_ch = Channel.value(
            file(
                params.approved_pc_selection_policy,
                checkIfExists: true
            )
        )

        analysis_readiness_ch = Channel.value(
            file(
                params.analysis_readiness,
                checkIfExists: true
            )
        )

        applied_normalization_policy_ch = Channel.value(
            file(
                params.applied_normalization_policy,
                checkIfExists: true
            )
        )

        INTEGRATION_POLICY_RESOLUTION_PHASE(
            scientific_meta_ch,
            pca_object_ch,
            dimensionality_summary_ch,
            approved_pc_selection_policy_ch,
            analysis_readiness_ch,
            applied_normalization_policy_ch
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

            - analysis_plan_qc_process
            - qc_threshold_approval_phase
            - qc_filter_policy_phase
            - qc_visualization_policy_phase
            - doublet_detection_scoring_phase
            - doublet_policy_resolution_phase
            - doublet_approval_phase
            - doublet_filter_phase
            - post_doublet_qc_summary_phase
            - normalization_policy_resolution_phase
            - normalization_approval_phase
            - normalization_phase
            - dimensionality_phase
            - pc_selection_policy_resolution_phase
            - pc_selection_approval_phase
            - integration_policy_resolution_phase
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

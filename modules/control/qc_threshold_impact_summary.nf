nextflow.enable.dsl=2


process QC_THRESHOLD_IMPACT_SUMMARY {

    tag "${policy_meta.project_id ?: 'project'}:${policy_meta.biological_run_id ?: 'run'}"

    label 'lightweight'

    publishDir(
        "${params.outdir}/control/qc_threshold_impact",
        mode: 'copy',
        overwrite: true
    )

    input:

    tuple val(policy_meta), path(policy_yaml)
    path cell_metadata
    path by_sample
    path impact_script

    output:

    tuple val(impact_meta), path("qc_threshold_impact.yaml"),
        emit: impact_yaml

    tuple val(impact_meta), path("qc_threshold_impact_by_sample.csv"),
        emit: by_sample

    tuple val(impact_meta), path("qc_threshold_impact_by_metric.csv"),
        emit: by_metric

    tuple val(impact_meta), path("qc_threshold_overlap_summary.csv"),
        emit: overlap

    tuple val(impact_meta), path("qc_lower_tail_by_sample.csv"),
        emit: lower_tail_by_sample

    tuple val(impact_meta), path("qc_lower_tail_summary.csv"),
        emit: lower_tail_summary

    tuple val(summary_meta), path("qc_threshold_impact.summary.tsv"),
        emit: summary

    path "qc_threshold_impact.session.txt",
        emit: session_info

    script:

    impact_meta = policy_meta + [
        stage         : 'qc_threshold_impact_summary',
        artifact_type : 'qc_threshold_impact',
        application_status: 'not_applied'
    ]

    summary_meta = policy_meta + [
        stage         : 'qc_threshold_impact_summary',
        artifact_type : 'qc_threshold_impact_summary'
    ]

    """
    python3 ${impact_script} \
        --policy-yaml ${policy_yaml} \
        --cell-metadata ${cell_metadata} \
        --by-sample ${by_sample} \
        --impact-yaml qc_threshold_impact.yaml \
        --by-sample-output qc_threshold_impact_by_sample.csv \
        --by-metric-output qc_threshold_impact_by_metric.csv \
        --overlap-output qc_threshold_overlap_summary.csv \
        --lower-tail-by-sample-output qc_lower_tail_by_sample.csv \
        --lower-tail-summary-output qc_lower_tail_summary.csv \
        --summary qc_threshold_impact.summary.tsv \
        --session-info qc_threshold_impact.session.txt \
        --project-id '${policy_meta.project_id ?: ""}' \
        --biological-run-id '${policy_meta.biological_run_id ?: ""}'
    """

    stub:

    impact_meta = policy_meta + [
        stage         : 'qc_threshold_impact_summary',
        artifact_type : 'qc_threshold_impact',
        application_status: 'not_applied'
    ]

    summary_meta = policy_meta + [
        stage         : 'qc_threshold_impact_summary',
        artifact_type : 'qc_threshold_impact_summary'
    ]

    """
    cat > qc_threshold_impact.yaml <<'YAML'
    qc_threshold_impact:
      schema_version: "1.0"
      plan_id: "stub_plan"
      threshold_mode: "data_driven"
      approval_status: "proposed"
      application_status: "not_applied"
      cell_count: 0
      flagged_any: 0
      retained: 0
      flagged_fraction: 0.0
      enabled_metrics: []
    YAML

    cat > qc_threshold_impact_by_sample.csv <<'CSV'
    sample,Group,cells_before,flagged_any,retained,flagged_fraction
    CSV

    cat > qc_threshold_impact_by_metric.csv <<'CSV'
    sample,Group,metric,cells_before,flagged_cells,flagged_fraction,unique_flagged_cells,unique_flagged_fraction,overlapped_flagged_cells
    CSV

    cat > qc_threshold_overlap_summary.csv <<'CSV'
    metric_a,metric_b,metric_a_flagged,metric_b_flagged,flagged_both,flagged_union,overlap_fraction_of_metric_a,overlap_fraction_of_metric_b
    CSV

    cat > qc_lower_tail_by_sample.csv <<'CSV'
    sample,Group,metric,threshold,cells_evaluated,cells_below_threshold,fraction_below_threshold
    CSV

    cat > qc_lower_tail_summary.csv <<'CSV'
    metric,threshold,cells_evaluated,cells_below_threshold,fraction_below_threshold
    CSV

    cat > qc_threshold_impact.summary.tsv <<'TSV'
    plan_id	project_id	biological_run_id	threshold_mode	approval_status	cell_count	flagged_any	retained	flagged_fraction	enabled_metric_count
    stub_plan	stub_project	stub_run	data_driven	proposed	0	0	0	0.0	0
    TSV

    cat > qc_threshold_impact.session.txt <<'TXT'
    stub=true
    process=QC_THRESHOLD_IMPACT_SUMMARY
    application_status=not_applied
    TXT
    """
}

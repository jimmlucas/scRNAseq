nextflow.enable.dsl=2


process NORMALIZATION_POLICY_RESOLUTION {

    tag "${meta.project_id ?: 'project'}:${meta.biological_run_id ?: 'run'}"

    label 'lightweight'

    publishDir(
        "${params.outdir}/control/normalization_policy",
        mode: 'copy',
        overwrite: true
    )

    input:

    val meta
    path readiness
    path post_doublet_summary
    path post_doublet_by_sample
    path post_doublet_warnings
    path resolver_script

    output:

    tuple val(policy_meta), path("normalization_policy.yaml"),
        emit: policy_yaml

    tuple val(policy_meta), path("normalization_policy.tsv"),
        emit: policy_tsv

    tuple val(summary_meta), path("normalization_policy_resolution.summary.tsv"),
        emit: summary

    path "normalization_policy_resolution.session.txt",
        emit: session_info

    script:

    policy_meta = meta + [
        stage          : 'normalization_policy_resolution',
        artifact_type  : 'normalization_policy',
        approval_status: 'proposed'
    ]

    summary_meta = meta + [
        stage         : 'normalization_policy_resolution',
        artifact_type : 'normalization_policy_resolution'
    ]

    """
    python3 ${resolver_script} \
        --readiness ${readiness} \
        --summary ${post_doublet_summary} \
        --by-sample ${post_doublet_by_sample} \
        --warnings ${post_doublet_warnings} \
        --policy-yaml normalization_policy.yaml \
        --policy-tsv normalization_policy.tsv \
        --resolution-summary normalization_policy_resolution.summary.tsv \
        --session-info normalization_policy_resolution.session.txt \
        --project-id '${meta.project_id ?: ""}' \
        --project-name '${meta.project_name ?: ""}' \
        --biological-run-id '${meta.biological_run_id ?: ""}'
    """

    stub:

    policy_meta = meta + [
        stage          : 'normalization_policy_resolution',
        artifact_type  : 'normalization_policy',
        approval_status: 'proposed'
    ]

    summary_meta = meta + [
        stage         : 'normalization_policy_resolution',
        artifact_type : 'normalization_policy_resolution'
    ]

    """
    cat > normalization_policy.yaml <<'YAML'
    normalization_policy:
      schema_version: "1.0"
      status: "proposed"
      strategy:
        method: "LogNormalize"
        variable_features_method: "vst"
        nfeatures: 2000
        scale_factor: 10000
        scope: "per_sample"
        assay: "RNA"
      approval:
        required: true
        status: "proposed"
      application:
        status: "not_applied"
        normalization_applied: false
    YAML

    cat > normalization_policy.tsv <<'TSV'
    project_id	biological_run_id	readiness	cells	features	samples	method	variable_features_method	nfeatures	scale_factor	scope	assay	integration_required	approval_status	application_status
    stub_project	stub_run	ready	0	0	0	LogNormalize	vst	2000	10000	per_sample	RNA	false	proposed	not_applied
    TSV

    cat > normalization_policy_resolution.summary.tsv <<'TSV'
    project_id	biological_run_id	readiness	cells	features	samples	method	scope	integration_required	approval_status
    stub_project	stub_run	ready	0	0	0	LogNormalize	per_sample	false	proposed
    TSV

    cat > normalization_policy_resolution.session.txt <<'TXT'
    stub=true
    process=NORMALIZATION_POLICY_RESOLUTION
    approval_status=proposed
    TXT
    """
}

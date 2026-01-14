# Prepare Mime-compatible datasets for machine learning analysis
# Mime: https://github.com/l-magnificence/Mime

# rule prepare_mime_dataset_survival:
#     """Prepare Mime-compatible dataset for survival analysis"""
#     input:
#         counts="{project}/quantification/{tool}/{project}_count_matrix_corrected.txt",
#         samples=get_info,
#     output:
#         mime_rds="{project}/mime/{tool}/survival_dataset.rds",
#     log:
#         "logs/{project}/prepare_mime_survival_{tool}.log",
#     container:
#         (
#             "docker://btrspg/rlan:20260104"
#             if config["container"].get("r", None) is None
#             else config["container"].get("r", None)
#         )
#     threads: 1
#     resources:
#         mem_mb=config["resources"]["mem_mb"].get("prepare_mime", 8192),
#     params:
#         project=config["project"],
#         analysis_type="survival",
#         time_col=config.get("mime_time_col", "OS.time"),
#         status_col=config.get("mime_status_col", "OS"),
#         discovery=config.get("discovery_dataset", ""),
#     script:
#         "../../scripts/utilities/prepare_mime_dataset.R"


rule prepare_mime_dataset_response:
    """Prepare Mime-compatible dataset for response prediction"""
    input:
        expression="{project}/quantification/{tool}/{dataset}_TPM_matrix_corrected.txt",
        samples=get_info,
    output:
        mime_rds="{project}/mime/{tool}/{dataset}_response_dataset.rds",
    log:
        "logs/{project}/prepare_mime_response_{dataset}_{tool}.log",
    container:
        (
            "docker://btrspg/rlan:20260104"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    threads: 1
    resources:
        mem_mb=config["resources"]["mem_mb"].get("prepare_mime", 8192),
    params:
        project=get_project,
        dataset=get_dataset,
        response_col=config.get("mime", "response_col"),
    script:
        "../../scripts/utilities/prepare_mime_dataset.R"


rule combine_mime_datasets:
    """Combine individual Mime dataset RDS files into one list for Mime"""
    input:
        rds_files=get_mime_rds_files,
    output:
        combined_rds="{project}/mime/{tool}/combined_response_datasets.rds",
    log:
        "logs/{project}/combine_mime_datasets_{tool}.log",
    container:
        (
            "docker://btrspg/rlan:20260114"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    threads: 1
    resources:
        mem_mb=config["resources"]["mem_mb"].get("prepare_mime", 8192),
    params:
        samples=config["samples"],
        project=get_project,
    script:
        "../../scripts/utilities/combine_mime_datasets.R"


rule mime_response_analysis:
    """
    Run Mime response prediction analysis with model saving.

    Trains ML models on discovery dataset, validates on all datasets,
    and saves trained models for future implementation.
    """
    input:
        combined_rds="{project}/mime/{tool}/combined_response_datasets.rds",
        genelist="{project}/DEG/{tool}_{dataset}_combined_degs.tsv",  # Or use a specific gene list
    output:
        directory="{project}/mime/{tool}/response_analysis/",
    params:
        project=config["project"],
        tool=get_tool,
        methods=config.get("mime", {}).get("response_methods", ["nb", "svmRadialWeights", "rf", "kknn", "adaboost", "LogitBoost", "cancerclass"]),
        seed=config.get("mime", {}).get("seed", 5201314),
    container:
        (
            "docker://btrspg/rlan:20260114"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    threads: config["threads"].get("mime", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("mime", 16384),
    log:
        "logs/{project}/mime_response_analysis_{tool}.log",
    script:
        "../../scripts/utilities/mime_response_analysis.R"


rule mime_apply_model:
    """
    Apply a trained Mime model to new datasets.

    Use this rule to implement a saved model on new data
    without retraining or benchmarking.
    """
    input:
        model="{project}/mime/{tool}/response_analysis/trained_models.rds",
        new_data="{project}/new_data/{dataset}_expression.tsv",
    output:
        predictions="{project}/mime/{tool}/predictions/{dataset}_predictions.rds",
    params:
        project=config["project"],
        tool=get_tool,
        dataset=get_dataset,
    container:
        (
            "docker://btrspg/rlan:20260114"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    threads: 1
    resources:
        mem_mb=config["resources"]["mem_mb"].get("mime", 8192),
    log:
        "logs/{project}/mime_apply_{tool}_{dataset}.log",
    script:
        "../../scripts/utilities/mime_apply_model.R"

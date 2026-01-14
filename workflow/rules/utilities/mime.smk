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
        counts="{project}/quantification/{tool}/{dataset}_TPM_matrix_corrected.txt",
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

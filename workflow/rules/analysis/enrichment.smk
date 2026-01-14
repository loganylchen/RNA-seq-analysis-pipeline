rule clusterprofiler_enrichment:
    input:
        combined_deg_tsv="{project}/DEG/STAR_FC_{dataset}_combined_degs.tsv",
    output:
        enrichment="{project}/enrichment/clusterprofiler/STAR_FC_{dataset}_enrichment.rds",
    log:
        "logs/{project}/enrichment_{dataset}_clusterprofiler.log",
    container:
        (
            "docker://btrspg/clusterprofiler:4.14.0"
            if config["container"].get("clusterprofiler", None) is None
            else config["container"].get("clusterprofiler", None)
        )
    params:
        species=config["clusterprofiler"]["species"],
        padj_threshold=config["clusterprofiler"]["padj_threshold"],
        deg_tool_n_threshold=config["clusterprofiler"]["deg_tool_n_threshold"],
    threads: config["threads"].get("default", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("clusterprofiler", 8192),
    priority: 10
    script:
        "../../scripts/analysis/enrichment/clusterprofiler.R"


# rule parse_clusterprofiler_enrichment:
#     input:
#         enrichment="{project}/enrichment/clusterprofiler/all_enrichment.rds",
#     output:
#         discovery_go="{project}/enrichment/clusterprofiler/discovery_go_enrichment.tsv",
#         discovery_kegg="{project}/enrichment/clusterprofiler/discovery_kegg_enrichment.tsv",
#         discovery_others="{project}/enrichment/clusterprofiler/discovery_others_enrichment.tsv",
#         discovery_gsea="{project}/enrichment/clusterprofiler/discovery_gsea_enrichment.tsv",
#         validation_go="{project}/enrichment/clusterprofiler/validation_go_enrichment.tsv",
#         validation_kegg="{project}/enrichment/clusterprofiler/validation_kegg_enrichment.tsv",
#         validation_others="{project}/enrichment/clusterprofiler/validation_others_enrichment.tsv",
#         validation_gsea="{project}/enrichment/clusterprofiler/validation_gsea_enrichment.tsv",
#     log:
#         "logs/{project}/parse_enrichment_clusterprofiler.log",
#     container:
#         (
#             "docker://btrspg/clusterprofiler:4.14.0"
#             if config["container"].get("clusterprofiler", None) is None
#             else config["container"].get("clusterprofiler", None)
#         )
#     threads: config["threads"].get("default", 1)
#     resources:
#         mem_mb=config["resources"]["mem_mb"].get("clusterprofiler", 8192),
#     priority: 10
#     script:
#         "../../../scripts/analysis/enrichment/parse_clusterprofiler.R"

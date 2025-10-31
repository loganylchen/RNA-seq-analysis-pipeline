rule clusterprofiler_enrichment:
    input:
        discovery_deg_tsv="{project}/deseq2/discovery_deg.tsv",
        validation_deg_tsv="{project}/deseq2/validation_deg.tsv",
    output:
        enrichment="{project}/enrichment/all_enrichment.rds",
    log:
        "logs/{project}/enrichment_clusterprofiler.log",
    container:
        (
            "docker://btrspg/clusterprofiler:4.14.0"
            if config["container"].get("clusterprofiler", None) is None
            else config["container"].get("clusterprofiler", None)
        )
    params:
        species=config["clusterprofiler"]["species"],
        padj_threshold=config["deg"].get("padj", 0.05),
        log2fc_threshold=config["deg"].get("log2fc", 2),
    threads: 1
    priority: 10
    script:
        "../scripts/clusterprofiler.R"


rule parse_clusterprofiler_enrichment:
    input:
        enrichment="{project}/enrichment/all_enrichment.rds",
    output:
        discovery_go="{project}/enrichment/discovery_go_enrichment.tsv",
        discovery_kegg="{project}/enrichment/discovery_kegg_enrichment.tsv",
        discovery_others="{project}/enrichment/discovery_others_enrichment.tsv",
        discovery_gsea="{project}/enrichment/discovery_gsea_enrichment.tsv",
        validation_go="{project}/enrichment/validation_go_enrichment.tsv",
        validation_kegg="{project}/enrichment/validation_kegg_enrichment.tsv",
        validation_others="{project}/enrichment/validation_others_enrichment.tsv",
        validation_gsea="{project}/enrichment/validation_gsea_enrichment.tsv",
    log:
        "logs/{project}/parse_enrichment_clusterprofiler.log",
    container:
        (
            "docker://btrspg/clusterprofiler:4.14.0"
            if config["container"].get("clusterprofiler", None) is None
            else config["container"].get("clusterprofiler", None)
        )
    threads: 1
    priority: 10
    script:
        "../scripts/parse_clusterprofiler.R"

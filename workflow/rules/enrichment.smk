# rule genekitr_enrichment:
#     input:
#         "results/diffexp/{contrast}/{subclass}.diffexp.tsv",
#     output:
#         directory("results/enrichment/{contrast}_{subclass}"),
#     params:
#         log2foldchange_threshold=config['thresholds']['deseq2']['log2foldchange'],
#         padj_threshold=config['thresholds']['deseq2']['padj'],
#         go_name=config['ref']['go_name'],
#         kegg_name=config['ref']['kegg_name']
#     container:
#         "docker://btrspg/genekitr:latest"
#     benchmark:
#         "benchmarks/{contrast}_{subclass}.enrichment.benchmark.txt"
#     log:
#         "logs/enrichment/{contrast}_{subclass}.log",
#     benchmark:
#         "benchmarks/{contrast}_{subclass}.enrichment.genekitr.benchmark.txt",
#     script:
#         "../scripts/genekitr.R"


rule clusterprofiler_enrichment:
    input:
        discovery_deg_tsv="{project}/deseq2/discovery_deg.tsv",
        validation_deg_tsv="{project}/deseq2/validation_deg.tsv",
    output:
        discovery_enrichment="{project}/enrichment/discovery.rds",
        validation_enrichment="{project}/enrichment/validation.rds",
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
    threads: 1
    script:
        "../scripts/clusterprofiler.R"

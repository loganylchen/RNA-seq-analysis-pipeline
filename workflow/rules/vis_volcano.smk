rule volcano_vis:
    input:
        discovery_deg_rds="{project}/DEG/deseq2/discovery_deg.rds",
        validation_deg_rds="{project}/DEG/deseq2/validation_deg.rds",
        geneid_to_genename="resources/gene_id_to_gene_name.tsv",
    output:
        discovery_png="{project}/visualization/Volcano_discovery.png",
        discovery_pdf="{project}/visualization/Volcano_discovery.pdf",
        validation_png="{project}/visualization/Volcano_validation.png",
        validation_pdf="{project}/visualization/Volcano_validation.pdf",
    container:
        (
            "docker://btrspg/rlan:20251027"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    log:
        "logs/{project}/visualization_Volcano.diffexp.log",
    threads: 1
    benchmark:
        "benchmarks/{project}/Volcano.benchmark.txt"
    script:
        "../scripts/volcano.R"

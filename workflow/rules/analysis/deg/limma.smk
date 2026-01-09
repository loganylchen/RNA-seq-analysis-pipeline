# tool could be STAR_FC or salmon or kallisto
rule limma_trend:
    input:
        counts="{project}/quantification/{tool}/{dataset}_count_matrix.txt",
    output:
        deg_rds="{project}/DEG/limma_trend/{tool}/{dataset}_deg.rds",
        deg_tsv="{project}/DEG/limma_trend/{tool}/{dataset}_deg.tsv",
    params:
        samples=config["samples"],
        dataset=get_dataset,
        project=get_project,
        case_condition=get_case_condition,
        control_condition=get_control_condition,
    container:
        (
            "docker://btrspg/limma:3.62.1"
            if config["container"].get("limma", None) is None
            else config["container"].get("limma", None)
        )
    log:
        "logs/{project}/limma_trend_{tool}_{dataset}.log",
    threads: config["threads"].get("limma", 4)
    benchmark:
        "benchmarks/{project}/limma_trend_{tool}_{dataset}.benchmark.txt"
    resources:
        mem_mb=config["resources"]["mem_mb"].get("limma", 8192),
    script:
        "../../../scripts/analysis/deg/limma_trend.R"


rule limma_voom:
    input:
        counts="{project}/quantification/{tool}/{dataset}_count_matrix.txt",
    output:
        deg_rds="{project}/DEG/limma_voom/{tool}/{dataset}_deg.rds",
        deg_tsv="{project}/DEG/limma_voom/{tool}/{dataset}_deg.tsv",
    params:
        samples=config["samples"],
        dataset=get_dataset,
        project=get_project,
        case_condition=get_case_condition,
        control_condition=get_control_condition,
    container:
        (
            "docker://btrspg/limma:3.62.1"
            if config["container"].get("limma", None) is None
            else config["container"].get("limma", None)
        )
    log:
        "logs/{project}/limma_voom_{tool}_{dataset}.log",
    threads: config["threads"].get("limma", 4)
    benchmark:
        "benchmarks/{project}/limma_voom_{tool}_{dataset}.benchmark.txt"
    resources:
        mem_mb=config["resources"]["mem_mb"].get("limma", 8192),
    script:
        "../../../scripts/analysis/deg/limma_voom.R"

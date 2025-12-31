# tool could be STAR_FC or salmon or kallisto
rule limma_trend:
    input:
        counts="{project}/quantification/{tool}/count_matrix.txt",
    output:
        discovery_deg_rds="{project}/DEG/limma_trend/{tool}/discovery_deg.rds",
        validation_deg_rds="{project}/DEG/limma_trend/{tool}/validation_deg.rds",
        discovery_deg_tsv="{project}/DEG/limma_trend/{tool}/discovery_deg.tsv",
        validation_deg_tsv="{project}/DEG/limma_trend/{tool}/validation_deg.tsv",
    params:
        samples=config["samples"],
        project=project,
        case_condition=case_condition,
        control_condition=control_condition,
        discovery_sample_type=discovery_sample_type,
    container:
        (
            "docker://btrspg/limma:3.62.1"
            if config["container"].get("limma", None) is None
            else config["container"].get("limma", None)
        )
    log:
        "logs/{project}/limma_trend_{tool}.log",
    threads: config["threads"].get("limma", 4)
    benchmark:
        "benchmarks/{project}/limma_trend_{tool}.benchmark.txt"
    resources:
        mem_mb=config["resources"]["mem_mb"].get("limma", 8192),
    script:
        "../../../scripts/analysis/deg/limma_trend.R"


rule limma_voom:
    input:
        counts="{project}/quantification/{tool}/count_matrix.txt",
    output:
        discovery_deg_rds="{project}/DEG/limma_voom/{tool}/discovery_deg.rds",
        validation_deg_rds="{project}/DEG/limma_voom/{tool}/validation_deg.rds",
        discovery_deg_tsv="{project}/DEG/limma_voom/{tool}/discovery_deg.tsv",
        validation_deg_tsv="{project}/DEG/limma_voom/{tool}/validation_deg.tsv",
    params:
        samples=config["samples"],
        project=project,
        case_condition=case_condition,
        control_condition=control_condition,
        discovery_sample_type=discovery_sample_type,
    container:
        (
            "docker://btrspg/limma:3.62.1"
            if config["container"].get("limma", None) is None
            else config["container"].get("limma", None)
        )
    log:
        "logs/{project}/limma_voom_{tool}.log",
    threads: config["threads"].get("limma", 4)
    benchmark:
        "benchmarks/{project}/limma_voom_{tool}.benchmark.txt"
    resources:
        mem_mb=config["resources"]["mem_mb"].get("limma", 8192),
    script:
        "../../../scripts/analysis/deg/limma_voom.R"

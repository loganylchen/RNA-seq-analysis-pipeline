# tool could be STAR_FC or salmon or kallisto
rule edger:
    input:
        counts="{project}/quantification/{tool}/{dataset}_count_matrix.txt",
    output:
        deg_rds="{project}/DEG/edger/{tool}/{dataset}_deg.rds",
        deg_tsv="{project}/DEG/edger/{tool}/{dataset}_deg.tsv",
    params:
        samples=config["samples"],
        dataset=get_dataset,
        project=get_project,
        design=get_edeger_design,
        case_condition=get_case_condition,
        control_condition=get_control_condition,
    container:
        (
            "docker://btrspg/edger:4.4.0"
            if config["container"].get("edger", None) is None
            else config["container"].get("edger", None)
        )
    log:
        "logs/{project}/edger_{tool}_{dataset}.log",
    threads: config["threads"].get("edger", 4)
    benchmark:
        "benchmarks/{project}/edger_{tool}_{dataset}.benchmark.txt"
    resources:
        mem_mb=config["resources"]["mem_mb"].get("edger", 8192),
    script:
        "../../../scripts/analysis/deg/edger.R"

# tool could be STAR_FC or salmon or kallisto
rule edger:
    input:
        counts="{project}/quantification/{tool}/count_matrix.txt",
    output:
        discovery_deg_rds="{project}/DEG/edger/{tool}/discovery_deg.rds",
        validation_deg_rds="{project}/DEG/edger/{tool}/validation_deg.rds",
        discovery_deg_tsv="{project}/DEG/edger/{tool}/discovery_deg.tsv",
        validation_deg_tsv="{project}/DEG/edger/{tool}/validation_deg.tsv",
    params:
        samples=config["samples"],
        project=project,
        case_condition=case_condition,
        control_condition=control_condition,
        discovery_sample_type=discovery_sample_type,
    container:
        (
            "docker://btrspg/edger:4.4.0"
            if config["container"].get("edger", None) is None
            else config["container"].get("edger", None)
        )
    log:
        "logs/{project}/edger_{tool}.log",
    threads: config["threads"].get("edger", 4)
    benchmark:
        "benchmarks/{project}/edger_{tool}.benchmark.txt"
    resources:
        mem_mb=config["resources"]["mem_mb"].get("edger", 8192),
    script:
        "../../../scripts/analysis/deg/edger.R"

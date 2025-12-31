# tool could be STAR_FC or salmon or kallisto
rule deseq2:
    input:
        counts="{project}/quantification/{tool}/count_matrix.txt",
    output:
        discovery_count_rds="{project}/DEG/deseq2/{tool}/discovery_count_matrix.rds",
        validation_count_rds="{project}/DEG/deseq2/{tool}/validation_count_matrix.rds",
        discovery_vst_rds="{project}/DEG/deseq2/{tool}/discovery_vst_matrix.rds",
        validation_vst_rds="{project}/DEG/deseq2/{tool}/validation_vst_matrix.rds",
        discovery_deg_rds="{project}/DEG/deseq2/{tool}/discovery_deg.rds",
        validation_deg_rds="{project}/DEG/deseq2/{tool}/validation_deg.rds",
        discovery_deg_tsv="{project}/DEG/deseq2/{tool}/discovery_deg.tsv",
        validation_deg_tsv="{project}/DEG/deseq2/{tool}/validation_deg.tsv",
    params:
        samples=config["samples"],
        project=project,
        case_condition=case_condition,
        control_condition=control_condition,
        discovery_sample_type=discovery_sample_type,
    container:
        (
            "docker://btrspg/deseq2:1.46.0"
            if config["container"].get("deseq2", None) is None
            else config["container"].get("deseq2", None)
        )
    log:
        "logs/{project}/deseq2_{tool}.log",
    threads: config["threads"].get("deseq2", 4)
    benchmark:
        "benchmarks/{project}/deseq2_{tool}.benchmark.txt"
    resources:
        mem_mb=config["resources"]["mem_mb"].get("deseq2", 8192),
    script:
        "../../scripts/analysis/deg/deseq2.R"

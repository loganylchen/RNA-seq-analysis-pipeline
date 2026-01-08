# tool could be STAR_FC or salmon or kallisto
rule deseq2:
    input:
        counts="{project}/quantification/{tool}/{dataset}_count_matrix.txt",
    output:
        count_rds="{project}/DEG/deseq2/{tool}/{dataset}_count_matrix.rds",
        vst_rds="{project}/DEG/deseq2/{tool}/{dataset}_vst_matrix.rds",
        deg_rds="{project}/DEG/deseq2/{tool}/{dataset}_deg.rds",
        deg_tsv="{project}/DEG/deseq2/{tool}/{dataset}_deg.tsv",
    params:
        samples=config["samples"],
        dataset="{dataset}",
        project="{project}",
        case_condition=config["datasets"]["{dataset}"]["case_condition"],
        control_condition=config["datasets"]["{dataset}"]["control_condition"],
    container:
        (
            "docker://btrspg/deseq2:1.46.0"
            if config["container"].get("deseq2", None) is None
            else config["container"].get("deseq2", None)
        )
    log:
        "logs/{project}/deseq2_{tool}_{dataset}.log",
    threads: config["threads"].get("deseq2", 4)
    benchmark:
        "benchmarks/{project}/deseq2_{tool}_{dataset}.benchmark.txt"
    resources:
        mem_mb=config["resources"]["mem_mb"].get("deseq2", 8192),
    script:
        "../../../scripts/analysis/deg/deseq2.R"

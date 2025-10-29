rule deseq2_init:
    input:
        counts="{project}/quantification/STAR_fc_count_matrix.txt",
    output:
        discovery_count_rds="{project}/deseq2/discovery_count_matrix.rds",
        validation_count_rds="{project}/deseq2/validation_count_matrix.rds",
        discovery_vst_rds="{project}/deseq2/discovery_vst_matrix.rds",
        validation_vst_rds="{project}/deseq2/validation_vst_matrix.rds",
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
        "logs/{project}/deseq2_init.log",
    threads: config["threads"]["deseq2"]
    script:
        "../scripts/deseq2-init.R"


rule deseq2:
    input:
        discovery_count_rds="{project}/deseq2/discovery_count_matrix.rds",
        validation_count_rds="{project}/deseq2/validation_count_matrix.rds",
    output:
        discovery_deg_rds="{project}/deseq2/discovery_deg.rds",
        validation_deg_rds="{project}/deseq2/validation_deg.rds",
        discovery_deg_tsv="{project}/deseq2/discovery_deg.tsv",
        validation_deg_tsv="{project}/deseq2/validation_deg.tsv",
    params:
        case_condition=case_condition,
        control_condition=control_condition,
    container:
        (
            "docker://btrspg/deseq2:1.46.0"
            if config["container"].get("deseq2", None) is None
            else config["container"].get("deseq2", None)
        )
    benchmark:
        "benchmarks/{project}/deseq2.benchmark.txt"
    log:
        "logs/{project}/DESeq2_diffexp.log",
    threads: config["threads"]["deseq2"]
    script:
        "../scripts/deseq2.R"

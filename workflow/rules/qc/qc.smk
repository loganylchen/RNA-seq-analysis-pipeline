rule multiqc_qc:
    input:
        get_qc_files(),
    output:
        outdir=directory(
            "{project}/qc/multiqc/",
        ),
    log:
        "logs/{project}/qc.log",
    container:
        (
            "docker://btrspg/multiqc:1.32"
            if config["container"].get("multiqc") is None
            else config["container"].get("multiqc", None)
        )
    threads: config["threads"].get("default", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("multiqc", 4096),
    priority: 10
    shell:
        "multiqc -f "
        "--outdir {output.outdir} {input} &>{log}"


rule qc_summary:
    input:
        samples=config["samples"],
        multiqc_dir="{project}/qc/multiqc/",
    output:
        summary="{project}/qc/qc_summary.tsv",
    log:
        "logs/{project}/qc_summary.log",
    container:
        (
            "docker://btrspg/rlan:20260104"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    threads: config["threads"].get("default", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("qc_summary", 16384),
    priority: 10
    params:
        project=config["project"],
        samples=config["samples"],
    script:
        "../../scripts/qc/qc_summary.R"


rule qc_summary_visualization:
    input:
        summary="{project}/qc/qc_summary.tsv",
    output:
        figures_dir=directory("{project}/qc/qc_summary_figures/"),
    log:
        "logs/{project}/qc_summary_viz.log",
    container:
        (
            "docker://btrspg/rlan:20260104"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    threads: config["threads"].get("default", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("qc_summary_viz", 16384),
    priority: 10

    script:
        "../../scripts/qc/qc_summary_visualization.R"

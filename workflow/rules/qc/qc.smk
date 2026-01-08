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
        qc_files=get_qc_files(),
    output:
        summary="{project}/qc/qc_summary.tsv",
    log:
        "logs/{project}/qc_summary.log",
    container:
        (
            "docker://btrspg/rlan:20251120"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    threads: config["threads"].get("default", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("qc_summary", 8192),
    priority: 10
    params:
        project=config["project"],
        samples=config["samples"],
    script:
        "../../scripts/qc/qc_summary.R"

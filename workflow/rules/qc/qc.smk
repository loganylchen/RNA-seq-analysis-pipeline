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

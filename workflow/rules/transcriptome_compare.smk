rule transcript_compare_gffcompare:
    input:
        predicted_gtf="{project}/assembly/stringtie/merged.gtf",
        reference_gtf="resources/genome.gtf",
    output:
        cmp="{project}/assembly/stringtie/gffcompare.annotated.gtf",
    params:
        out_prefix="{project}/assembly/stringtie/gffcompare",
        extra=config["gffcompare"]["extra"],
    container:
        (
            "docker://btrspg/gffcompare:0.12.10"
            if config["container"].get("gffcompare", None) is None
            else config["container"].get("gffcompare", None)
        )
    threads: 1
    log:
        "logs/{project}/gffcompare.log",
    shell:
        "gffcompare -r {input.reference_gtf} "
        "-o {params.out_prefix} "
        "{input.predicted_gtf} &>{log}"

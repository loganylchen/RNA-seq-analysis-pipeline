rule transcript_compare_gffcompare:
    input:
        predicted_gtf="{project}/assembly/merged.gtf",
        reference_gtf="resources/genome.gtf",
    output:
        cmp="{project}/assembly/gffcompare.annotated.gtf",
    params:
        out_prefix="{project}/assembly/gffcompare",
        extra=config["gffcompare"]["extra"],
    threads: 1
    log:
        "logs/{project}/gffcompare.log",
    shell:
        "gffcompare -r {input.reference_gtf} "
        "-o {params.out_prefix} "
        "{input.predicted_gtf} &>{log}"

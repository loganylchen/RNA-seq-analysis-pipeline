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
    threads: config["threads"].get("default", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("gffcompare", 4096),
    log:
        "logs/{project}/gffcompare.log",
    shell:
        "gffcompare -r {input.reference_gtf} "
        "-o {params.out_prefix} "
        "{input.predicted_gtf} &>{log}"


rule gtf_to_bed12:
    input:
        gtf="{project}/assembly/stringtie/gffcompare.annotated.gtf",
    output:
        bed="{project}/assembly/stringtie/gffcompare.annotated.bed",
    log:
        "logs/{project}/gtf2bed12.log",
    benchmark:
        "benchmarks/{project}/gtf2bed12.benchmark.txt"
    container:
        (
            "docker://btrspg/python3:20251024"
            if config["container"].get("python3", None) is None
            else config["container"].get("python3", None)
        )
    resources:
        mem_mb=config["resources"]["mem_mb"].get("default", 4096),
    script:
        "../scripts/gtf2bed12.py"


rule sort_bed12:
    input:
        bed="{project}/assembly/stringtie/gffcompare.annotated.bed",
    output:
        bed="{project}/assembly/stringtie/gffcompare.sorted.bed",
    log:
        "logs/{project}/sort_bed.log",
    benchmark:
        "benchmarks/{project}/sort_bed.benchmark.txt"
    container:
        (
            "docker://btrspg/bedtools:2.31.1"
            if config["container"].get("bedtools", None) is None
            else config["container"].get("bedtools", None)
        )
    resources:
        mem_mb=config["resources"]["mem_mb"].get("default", 4096),
    shell:
        "bedtools sort -i {input.bed} > {output.bed} 2>{log}"

rule qualimap_rnaseq_qc:
    input:
        bam="{project}/alignment/{sample}/{sample}.star.bam",
        gtf="resources/genome.gtf",
    output:
        outdir="{project}/{sample}/qc/qualimap",
    params:
        pe="-pe" if is_pe else "",
    container:
        (
            "docker://btrspg/qualimap:2.3"
            if config["container"].get("qualimap", None) is None
            else config["container"].get("qualimap", None)
        )
    threads: 1
    resources:
        mem_mb=1024 * 20,
    log:
        "logs/{project}/{sample}_star.log",
    shell:
        "qualimap rnaseq "
        "-bam {input.bam} "
        "-gtf {input.gtf} "
        "-outdir {output.outdir} "
        "{params.pe} "
        "--java-mem-size=20G &>{log}"

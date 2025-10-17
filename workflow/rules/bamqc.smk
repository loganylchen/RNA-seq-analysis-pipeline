rule qualimap_rnaseq_qc:
    input:
        bam="{project}/alignment/{sample}/{sample}.star.bam",
        gtf="resources/genome.gtf",
    output:
        rnaseq_qc="{project}/qc/{sample}/rnaseq_qc_results.txt",
        temp_dir=temp(directory("{project}/qc/{sample}/temp")),
    params:
        pe="-pe" if is_pe else "",
        outdir=lambda wc, output: os.path.dirname(output.rnaseq_qc),
    container:
        (
            "docker://btrspg/qualimap:2.3"
            if config["container"].get("qualimap", None) is None
            else config["container"].get("qualimap", None)
        )
    threads: 1
    resources:
        mem_mb=1024 * 20,
        tmpdir="{project}/qc/{sample}/temp",
    log:
        "logs/{project}/{sample}_star.log",
    shell:
        "JAVA_OPTS='-Djava.io.tmpdir={resources.tmpdir}' qualimap rnaseq "
        "-bam {input.bam} "
        "-gtf {input.gtf} "
        "-outdir {params.outdir} "
        "{params.pe} "
        "--java-mem-size=20G &>{log}"

rule qualimap_rnaseq_qc:
    input:
        bam="{project}/alignment/STAR/{sample}/{sample}.bam",
        gtf="resources/genome.gtf",
    output:
        rnaseq_qc="{project}/qc/qualimap-rnaseq/{sample}/rnaseq_qc_results.txt",
        temp_dir=temp(directory("{project}/qc/temp/{sample}")),
    params:
        pe="-pe" if is_pe else "",
        outdir=lambda wc, output: os.path.dirname(output.rnaseq_qc),
    container:
        (
            "docker://btrspg/qualimap:2.3"
            if config["container"].get("qualimap", None) is None
            else config["container"].get("qualimap", None)
        )
    threads: config["threads"].get("qualimap", 10)
    resources:
        mem_mb=1024 * 20,
    log:
        "logs/{project}/{sample}_qualimap-rnaseq-qc.log",
    shell:
        "mkdir -p {output.temp_dir};"
        "JAVA_OPTS='-Djava.io.tmpdir={output.temp_dir}' qualimap rnaseq "
        "-bam {input.bam} "
        "-gtf {input.gtf} "
        "-outdir {params.outdir} "
        "{params.pe} "
        "--java-mem-size=20G &>{log}"

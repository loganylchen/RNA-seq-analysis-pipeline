rule modtect_mod:
    input:
        bam="{project}/alignment/STAR/{sample}/{sample}.split.bam",
        genome="resources/genome.fasta",
        bed="resources/transcriptome.bed",
    output:
        output="{project}/modification/modtect/{sample}/{sample}.modtect.combined.txt",
        temp=temp(
            directory("{project}/modification/modtect/{sample}/{sample}.modtect_tmp")
        ),
    log:
        log="logs/{project}/modtect_{sample}.log",
    params:
        label=lambda wc, output: output.output.replace(".combined.txt", ""),
    benchmark:
        "benchmarks/{project}/{sample}.modtect.benchmark.txt"
    container:
        (
            "docker://btrspg/modtect:20251029"
            if config["container"].get("modtect", None) is None
            else config["container"].get("modtect", None)
        )
    threads: config["threads"]["modtect"]
    shell:
        "modtect.py {input.bam} {input.genome} 0 0 0 "
        "--threads {threads} --regionFile {input.bed} "
        "--label {params.label}"
        "&>{log}"

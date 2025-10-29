rule modtect_mod:
    input:
        bam="{project}/alignment/{sample}/{sample}.star.bam",
        genome="resources/genome.fasta",
        bed="resources/transcriptome.bed",
    output:
        output="{project}/modification/{sample}/modtect/{sample}.modtect.combined.txt",
    log:
        log="logs/{project}/modtect_{sample}.log",
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
        "&>{log}"

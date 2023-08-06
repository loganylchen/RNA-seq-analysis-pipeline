rule modtect:
    input:
        bam='results/star/{sample}/split.bam',
        genome="resources/genome.fa",
        bed="resources/genome.sorted.bed",
    output:
        output='results/modtect/{sample}/modtect.combined.txt'
    log:
        log='logs/modtect/{sample}.log',
        error='logs/modtect/{sample}.err',
    params:
        label='results/modtect/{sample}/modtect',
    container:
        "docker://btrspg/modtect:latest"
    threads: config["threads"]["modtect"]
    shell:
        "modtect.py {input.bam} {input.genome} 0 0 0 "
        "--threads {threads} --regionFile {input.bed} "
        "--label {params.label} 1>{log.log} 2>{log.error}"


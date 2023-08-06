rule modtect:
    input:
        bam='results/star/{sample}/split.bam',
        genome="resources/genome.fa",
    output:
        output='results/modtect/{sample}/modtect.combined.txt'
    log:
        log='logs/modtect/{sample}.log',
        error='logs/modtect/{sample}.err',
    params:

        bed="resources/genome.sorted.bed",
        label='results/modtect/{sample}/modtect',
    threads: config["threads"]["modtect"]
    shell:
        "modtect.py {input.bam} {input.genome} 0 0 0 --threads {threads} --label {params.label} 1>{log.log} 2>{log.error}"


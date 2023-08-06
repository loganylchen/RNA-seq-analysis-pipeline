rule hamr:
    input:
        bam='results/star/{sample}/split.bam',
    output:
        outdir=directory('results/hamr/{sample}/'),
        output="results/hamr/{sample}/hamr.mods.txt"
    log:
        log='logs/hamr/{sample}.log',
        error='logs/hamr/{sample}.err',
    params:
        genome="resources/genome.fa",
        bed="resources/genome.sorted.bed",
        extra=config['params']['hamr']
    threads: config["threads"]["star"]
    container:
        "docker://btrspg/hamr:latest"
    shell:
        "hamr.py {input.bam} {params.genome} {output.outdir} -n {params.bed} {params.extra} 1>{log.log} 2>{log.error}"


rule hamr:
    input:
        bam='results/star/{sample}/split.bam',
        genome="resources/genome.fa",
        bed="resources/genome.merged.bed",
    output:
        outdir=directory('results/hamr/{sample}/'),
        output="results/hamr/{sample}/hamr.mods.txt"
    log:
        log='logs/hamr/{sample}.log',
        error='logs/hamr/{sample}.err',
    params:
        extra=config['params']['hamr']
    threads: config["threads"]["star"]
    container:
        "docker://btrspg/hamr:latest"
    shell:
        "hamr.py {input.bam} {input.genome} {output.outdir} -n {input.bed} {params.extra} 1>{log.log} 2>{log.error}"


rule modtect:
    input:
        bam='results/star/{sample}/split.bam',
    output:
        outdir=directory('results/hamr/{sample}/'),
        output="results/hamr/{sample}/hamr.mod.txt"
    log:
        log='logs/hamr/{sample}.log',
        error='logs/hamr/{sample}.err',
    params:
        genome="resources/genome.fa",
        extra=config['params']['hamr']
    threads: config["threads"]["star"]
    shell:
        "hamr.py {input.bam} {params.genome} {output.outdir} {params.extra} 1>{log.log} 2>{log.error}"


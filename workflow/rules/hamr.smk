rule hamr:
    input:
        bam='results/star/{sample}/split.bam',
        genome="resources/genome.fasta",
        bed="resources/genome.merged.bed",
    output:
        outdir=directory('results/hamr/{sample}/'),
        output="results/hamr/{sample}/hamr.mods.txt"
    log:
        log='logs/hamr/{sample}.log',
        error='logs/hamr/{sample}.err',
    benchmark:
        "benchmarks/{sample}.hamr.benchmark.txt"
    params:
        extra=config['params']['hamr']
    container:
        "docker://btrspg/hamr:latest"
    shell:
        "hamr.py {input.bam} {input.genome} {output.outdir} -n {input.bed} {params.extra} 1>{log.log} 2>{log.error}"


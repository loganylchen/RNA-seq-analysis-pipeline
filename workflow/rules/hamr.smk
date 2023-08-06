rule hamr:
    input:
        bam='results/star/{sample}/split.bam',
    output:
        outdir=directory('results/hamr/{sample}/')
    log:
        log='logs/hamr/{sample}.log',
        error='logs/hamr/{sample}.err',
    params:
        genome="resources/genome.fa",
        is_pe=' --pe ' if lambda w: samples.loc[w.sample].loc['seq_type']  == 'pe' else ' ',
    threads: config["threads"]["star"]
    shell:
        "hamr.py {input.bam} {params.genome} {output.outdir} {params.is_pe} 1>{log.log} 2>{log.error}"

clusters =
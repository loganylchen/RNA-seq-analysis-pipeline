rule align:
    input:
        unpack(get_clean_data),
    output:
        aln="results/star/{sample}/Aligned.sortedByCoord.out.bam",
        reads_per_gene="results/star/{sample}/ReadsPerGene.out.tab",
    log:
        "logs/star/{sample}.log",
    params:
        idx=lambda wc, input: input.index,
        extra="--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile {} {}".format(
            "resources/genome.gtf", config["params"]["star"]
        ),
    threads: config["threads"]["star"]
    wrapper:
        "v1.21.4/bio/star/align"

rule split_bam:
    input:
        aln = "results/star/{sample}/Aligned.sortedByCoord.out.bam",
        genome_dict="resources/genome.dict",
        genome='resources/genome.fasta'
    output:
        split_bam = "results/star/{sample}/split.bam",
    params:
        extra=config['params']['gatk']
    log:
        "logs/split_bam/{sample}.log",
    conda:
        "../envs/gatk4.yaml"
    threads: config["threads"]["gatk4"]
    shell:
        "gatk SplitNCigarReads --create-output-bam-index -R {input.genome} -I {input.aln} -o {output.split_bam} -U ALLOW_N_CIGAR_READS {params.extra} 2>{log}"


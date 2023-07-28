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
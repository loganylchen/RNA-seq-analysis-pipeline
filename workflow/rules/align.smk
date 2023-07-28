rule align:
    input:
        *get_clean_data(),
        fq1="results/clean_fastq/{sample}/{sample}_1.fastq.gz",
        fq2="results/clean_fastq/{sample}/{sample}_2.fastq.gz",
        index="resources/star_genome",
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
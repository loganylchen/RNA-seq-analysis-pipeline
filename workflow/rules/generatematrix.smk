rule count_matrix:
    input:
        expand(
            "results/star/{sample.sample_name}/ReadsPerGene.out.tab",
            sample=samples.itertuples(),
        ),
    output:
        "results/counts/count_matrix.tsv",
        "results/counts/count_matrix_unstrandedness.tsv",
        "results/counts/count_matrix_strandedness.tsv",
        "results/counts/count_matrix_reverse.tsv",
        "results/counts/strandedness_check.pdf"
    log:
        "logs/count-matrix.log",
    params:
        samples=samples["sample_name"].tolist(),
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"


rule feature_counts:
    input:

        samples=expand(
            "results/star/{sample.sample_name}/Aligned.sortedByCoord.out.bam",
            sample=samples.itertuples(),
        ),
        annotation="annotation.gtf",

    output:
        multiext(
            "results/counts/count_matrix",
            ".featureCounts",
            ".featureCounts.summary",
            ".featureCounts.jcounts",
        ),
    threads: 10
    params:
        strand=0,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        r_path="",  # implicitly sets the --Rpath flag
        extra=config['params']['featurecounts']
    log:
        "logs/count-matrix_featurecounts.log",
    wrapper:
        "v2.2.1/bio/subread/featurecounts"
rule wgcna:
    input:
        normalized_matrix="results/deseq2/vst_count_matrix.tsv",
    output:
        "results/deseq2/count_matrix.rds",
        "results/deseq2/normalized_count_matrix.tsv",
        "results/deseq2/vst_count_matrix.tsv"
    params:
        samples=config["samples"],
        model=config["diffexp"]["model"],
        count_threshold=config['thresholds']['deseq2']['count']
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log",
    threads: config['threads']['deseq2']
    script:
        "../scripts/deseq2-init.R"

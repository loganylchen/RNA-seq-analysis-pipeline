rule pca_vis:
    input:
        counts="results/counts/count_matrix.tsv",
    output:
        png="results/visualization/PCA.png",
        pdf="results/visualization/PCA.pdf",
    params:
        samples=config["samples"],
        model=config["diffexp"]["model"],
        count_threshold=config['thresholds']['deseq2']['count'],
        color_by=config['vis']['pcatools']['color_by'],
        shape_by=config['vis']['pcatools']['shape_by'],
    conda:
        "../envs/pcatools.yaml"
    log:
        "logs/pcatools/pcatools.log",
    threads: config['threads']['deseq2']
    script:
        "../scripts/pcatools.R"
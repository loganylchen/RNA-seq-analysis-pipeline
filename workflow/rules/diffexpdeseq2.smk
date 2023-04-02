rule deseq2_init:
    input:
        counts="results/counts/count_matrix.tsv",
    output:
        "results/deseq2/count_matrix.rds",
        "results/deseq2/normalized_count_matrix.tsv",
    params:
        samples=config["samples"],
        model=config["diffexp"]["model"],
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log",
    threads: config['threads']['deseq2']
    script:
        "../scripts/deseq2-init.R"


# rule pca:
#     input:
#         "results/deseq2/all.rds",
#     output:
#         report("results/pca.svg", "../report/pca.rst"),
#     params:
#         pca_labels=config["pca"]["labels"],
#     conda:
#         "../envs/deseq2.yaml"
#     log:
#         "logs/pca.log",
#     script:
#         "../scripts/plot-pca.R"
#
#
# rule deseq2:
#     input:
#         "results/deseq2/all.rds",
#     output:
#         table=report("results/diffexp/{contrast}.diffexp.tsv", "../report/diffexp.rst"),
#         ma_plot=report("results/diffexp/{contrast}.ma-plot.svg", "../report/ma.rst"),
#     params:
#         contrast=get_contrast,
#     conda:
#         "../envs/deseq2.yaml"
#     log:
#         "logs/deseq2/{contrast}.diffexp.log",
#     threads: get_deseq2_threads
#     script:
#         "../scripts/deseq2.R"
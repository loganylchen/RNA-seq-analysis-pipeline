rule wgcna:
    input:
        normalized_matrix="results/deseq2/vst_count_matrix.tsv",
    output:
        bwnet_rds = "results/wgcna/wgcna.rds",
        gene_module_key = "results/wgcna/gene_module_key.tsv",
    params:
        samples=config["samples"],
        phenotype=config['params']['wgcna_phenotype'],
        fig_outdir="results/wgcna/figures/"
    conda:
        "../envs/wgcna.yaml"
    benchmark:
        "benchmarks/wgcna.benchmark.txt"
    log:
        "logs/wgcna/init.log",
    threads: config['threads']['wgcna']
    script:
        "../scripts/wgcna.R"

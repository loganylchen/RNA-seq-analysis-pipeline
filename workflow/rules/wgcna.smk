rule wgcna:
    input:
        normalized_matrix="results/deseq2/vst_count_matrix.tsv",
    output:
        bwnet_rds="results/wgcna/wgcna.rds",
        gene_module_key="results/wgcna/gene_module_key.tsv",
    params:
        samples=config["samples"],
        phenotype=config["params"]["wgcna_phenotype"],
        fig_outdir="results/wgcna/figures/",
    conda:
        "../envs/wgcna.yaml"
    benchmark:
        "benchmarks/wgcna.benchmark.txt"
    log:
        "logs/wgcna/init.log",
    threads: config["threads"].get("wgcna", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("wgcna", 16384),
    script:
        "../scripts/wgcna.R"

rule wgcna:
    input:
        normalized_matrix="{project}/DEG/deseq2/discovery_vst_matrix.rds",
    output:
        bwnet_rds="{project}/wgcna/wgcna.rds",
        gene_module_key="{project}/wgcna/gene_module_key.tsv",
    params:
        samples=config["samples"],
        project=project,
        discovery_sample_type=discovery_sample_type,
        phenotype=config.get("params", {}).get("wgcna_phenotype", "condition"),
        fig_outdir="{project}/wgcna/figures/",
    conda:
        "../envs/wgcna.yaml"
    benchmark:
        "benchmarks/{project}/wgcna.benchmark.txt"
    log:
        "logs/{project}/wgcna_init.log",
    threads: config["threads"].get("wgcna", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("wgcna", 16384),
    script:
        "../scripts/wgcna.R"

rule volcano_vis:
    input:
        deg_exp = "results/diffexp/{contrast}/{subclass}.diffexp.tsv",
        geneid_to_genename = "resources/gene_id_to_gene_name.tsv"

    output:
        png="results/visualization/Volcano.{contrast}_{subclass}.diffexp.png",
        pdf="results/visualization/Volcano.{contrast}_{subclass}.diffexp.pdf",
    params:
        contrast='{contrast}',
        fc_threshold=config['thresholds']['volcano']['log2foldchange'],
        p_threshold=config['thresholds']['volcano']['padj']

    conda:
        "../envs/enhancedvolcano.yaml"
    log:
        "logs/visualization/Volcano.{contrast}_{subclass}.diffexp.log",
    benchmark:
        "benchmarks/Volcano.{contrast}_{subclass}.benchmark.txt"
    script:
        "../scripts/volcano.R"
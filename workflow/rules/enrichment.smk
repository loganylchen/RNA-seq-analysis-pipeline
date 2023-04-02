rule genekitr_enrichment:
    input:
        "results/diffexp/{contrast}.diffexp.tsv",
    output:
        directory("results/enrichment/{contrast}"),
    params:
        log2foldchange_threshold=config['thresholds']['deseq2']['log2foldchange'],
        padj_threshold=config['thresholds']['deseq2']['padj'],
        go_name=config['ref']['go_name'],
        kegg_name=config['ref']['kegg_name']
    container:
        "docker://btrspg/genekitr:main"
    log:
        "logs/enrichment/{contrast}.log",
    script:
        "../scripts/enrichment.R"



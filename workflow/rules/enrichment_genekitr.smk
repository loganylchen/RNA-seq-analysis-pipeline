rule genekitr_enrichment:
    input:
        "results/diffexp/{contrast}/{subclass}.diffexp.tsv",
    output:
        directory("results/enrichment/{contrast}_{subclass}"),
    params:
        log2foldchange_threshold=config['thresholds']['deseq2']['log2foldchange'],
        padj_threshold=config['thresholds']['deseq2']['padj'],
        go_name=config['ref']['go_name'],
        kegg_name=config['ref']['kegg_name']
    container:
        "docker://btrspg/genekitr:latest"
    benchmark:
        "benchmarks/{contrast}_{subclass}.enrichment.benchmark.txt"
    log:
        "logs/enrichment/{contrast}_{subclass}.log",
    benchmark:
        "benchmarks/{contrast}_{subclass}.enrichment.genekitr.benchmark.txt",
    script:
        "../scripts/genekitr.R"



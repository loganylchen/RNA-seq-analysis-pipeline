rule volcano_vis:
    input:
        "results/diffexp/{contrast}/{subclass}.diffexp.tsv",
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
    script:
        "../scripts/volcano.R"
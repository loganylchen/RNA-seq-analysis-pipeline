rule volcano_vis:
    input:
        "results/diffexp/{contrast}.diffexp.tsv",
    output:
        png="results/visualization/Volcano.{contrast}.diffexp.png",
        pdf="results/visualization/Volcano.{contrast}.diffexp.pdf",
    params:
        contrast=get_contrast,
        fc_threshold=config['thresholds']['volcano']['log2foldchange'],
        p_threshold=config['thresholds']['volcano']['padj']

    conda:
        "../envs/enhancedvolcano.yaml"
    log:
        "logs/visualization/Volcano.{contrast}.diffexp.log",
    script:
        "../scripts/volcano.R"
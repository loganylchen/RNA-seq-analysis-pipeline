# DEG Visualization Rules
# Create heatmaps, scatter plots, and boxplots for DEG results


rule deg_visualization:
    """
    Generate visualizations for DEG results including:
    - Heatmap of DEGs using ComplexHeatmap
    - Scatter plot comparing discovery vs validation log2FC
    - Boxplots for genes significant in both cohorts
    """
    input:
        discovery_deg_rds="{project}/DEG/deseq2/discovery_deg.rds",
        validation_deg_rds="{project}/DEG/deseq2/validation_deg.rds",
        discovery_vst_rds="{project}/DEG/deseq2/discovery_vst_matrix.rds",
        validation_vst_rds="{project}/DEG/deseq2/validation_vst_matrix.rds",
        samples=config["samples"],
    output:
        heatmap="{project}/DEG/visualization/deg_heatmap.png",
        scatter="{project}/DEG/visualization/discovery_vs_validation_scatter.png",
        boxplot="{project}/DEG/visualization/significant_genes_boxplot.png",
        discovery_deg_list="{project}/DEG/visualization/discovery_DEG_list.tsv",
        validation_deg_list="{project}/DEG/visualization/validation_DEG_list.tsv",
        comparison_table="{project}/DEG/visualization/discovery_validation_comparison.tsv",
    params:
        project=project,
        case_condition=case_condition,
        control_condition=control_condition,
        discovery_sample_type=discovery_sample_type,
        log2fc_threshold=config.get("deg", {}).get("log2fc", 2),
        padj_threshold=config.get("deg", {}).get("padj", 0.05),
        heatmap_top_n=config.get("deg_vis", {}).get("heatmap_top_n", 100),
    container:
        (
            "docker://btrspg/rlan:20251229"
            if config["container"].get("deg_vis", None) is None
            else config["container"].get("deg_vis", None)
        )
    benchmark:
        "benchmarks/{project}/deg_visualization.benchmark.txt"
    log:
        "logs/{project}/deg_visualization.log",
    threads: config["threads"].get("deg_vis", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("deg_vis", 16384),
    script:
        "../scripts/deg_visualization.R"


rule deg_visualization_salmon:
    """
    DEG visualization using Salmon quantification results.
    """
    input:
        discovery_deg_rds="{project}/DEG/deseq2/salmon_discovery_deg.rds",
        validation_deg_rds="{project}/DEG/deseq2/salmon_validation_deg.rds",
        discovery_vst_rds="{project}/DEG/deseq2/salmon_discovery_vst_matrix.rds",
        validation_vst_rds="{project}/DEG/deseq2/salmon_validation_vst_matrix.rds",
        samples=config["samples"],
    output:
        heatmap="{project}/DEG/visualization/salmon_deg_heatmap.png",
        scatter="{project}/DEG/visualization/salmon_discovery_vs_validation_scatter.png",
        boxplot="{project}/DEG/visualization/salmon_significant_genes_boxplot.png",
        discovery_deg_list="{project}/DEG/visualization/salmon_discovery_DEG_list.tsv",
        validation_deg_list="{project}/DEG/visualization/salmon_validation_DEG_list.tsv",
        comparison_table="{project}/DEG/visualization/salmon_discovery_validation_comparison.tsv",
    params:
        project=project,
        case_condition=case_condition,
        control_condition=control_condition,
        discovery_sample_type=discovery_sample_type,
        log2fc_threshold=config.get("deg", {}).get("log2fc", 2),
        padj_threshold=config.get("deg", {}).get("padj", 0.05),
        heatmap_top_n=config.get("deg_vis", {}).get("heatmap_top_n", 100),
    container:
        (
            "docker://btrspg/rlan:20251229"
            if config["container"].get("deg_vis", None) is None
            else config["container"].get("deg_vis", None)
        )
    benchmark:
        "benchmarks/{project}/deg_visualization_salmon.benchmark.txt"
    log:
        "logs/{project}/deg_visualization_salmon.log",
    threads: config["threads"].get("deg_vis", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("deg_vis", 16384),
    script:
        "../scripts/deg_visualization.R"


rule deg_visualization_kallisto:
    """
    DEG visualization using Kallisto quantification results.
    """
    input:
        discovery_deg_rds="{project}/DEG/deseq2/kallisto_discovery_deg.rds",
        validation_deg_rds="{project}/DEG/deseq2/kallisto_validation_deg.rds",
        discovery_vst_rds="{project}/DEG/deseq2/kallisto_discovery_vst_matrix.rds",
        validation_vst_rds="{project}/DEG/deseq2/kallisto_validation_vst_matrix.rds",
        samples=config["samples"],
    output:
        heatmap="{project}/DEG/visualization/kallisto_deg_heatmap.png",
        scatter="{project}/DEG/visualization/kallisto_discovery_vs_validation_scatter.png",
        boxplot="{project}/DEG/visualization/kallisto_significant_genes_boxplot.png",
        discovery_deg_list="{project}/DEG/visualization/kallisto_discovery_DEG_list.tsv",
        validation_deg_list="{project}/DEG/visualization/kallisto_validation_DEG_list.tsv",
        comparison_table="{project}/DEG/visualization/kallisto_discovery_validation_comparison.tsv",
    params:
        project=project,
        case_condition=case_condition,
        control_condition=control_condition,
        discovery_sample_type=discovery_sample_type,
        log2fc_threshold=config.get("deg", {}).get("log2fc", 2),
        padj_threshold=config.get("deg", {}).get("padj", 0.05),
        heatmap_top_n=config.get("deg_vis", {}).get("heatmap_top_n", 100),
    container:
        (
            "docker://btrspg/rlan:20251229"
            if config["container"].get("deg_vis", None) is None
            else config["container"].get("deg_vis", None)
        )
    benchmark:
        "benchmarks/{project}/deg_visualization_kallisto.benchmark.txt"
    log:
        "logs/{project}/deg_visualization_kallisto.log",
    threads: config["threads"].get("deg_vis", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("deg_vis", 16384),
    script:
        "../scripts/deg_visualization.R"


rule deg_lasso_classifier:
    """
    Build LASSO logistic regression classifier using DEGs from discovery dataset
    and test performance on validation dataset.
    """
    input:
        discovery_deg_rds="{project}/DEG/deseq2/discovery_deg.rds",
        validation_deg_rds="{project}/DEG/deseq2/validation_deg.rds",
        discovery_vst_rds="{project}/DEG/deseq2/discovery_vst_matrix.rds",
        validation_vst_rds="{project}/DEG/deseq2/validation_vst_matrix.rds",
        samples=config["samples"],
    output:
        signature="{project}/classifier/DEG/lasso_signature_genes.tsv",
        coefficients="{project}/classifier/DEG/lasso_coefficients.tsv",
        discovery_predictions="{project}/classifier/DEG/discovery_predictions.tsv",
        validation_predictions="{project}/classifier/DEG/validation_predictions.tsv",
        roc_plot="{project}/classifier/DEG/lasso_roc_curve.png",
        summary="{project}/classifier/DEG/lasso_summary.txt",
    params:
        project=project,
        case_condition=case_condition,
        control_condition=control_condition,
        discovery_sample_type=discovery_sample_type,
        log2fc_threshold=config.get("deg", {}).get("log2fc", 2),
        padj_threshold=config.get("deg", {}).get("padj", 0.05),
    container:
        (
            "docker://btrspg/glmnet:4.1_10"
            if config["container"].get("lasso", None) is None
            else config["container"].get("lasso", None)
        )
    benchmark:
        "benchmarks/{project}/deg_lasso_classifier.benchmark.txt"
    log:
        "logs/{project}/deg_lasso_classifier.log",
    threads: config["threads"].get("deg_vis", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("deg_vis", 16384),
    script:
        "../scripts/deg_lasso_classifier.R"


rule deg_lasso_classifier_salmon:
    """
    LASSO classifier using Salmon quantification results.
    """
    input:
        discovery_deg_rds="{project}/DEG/deseq2/salmon_discovery_deg.rds",
        validation_deg_rds="{project}/DEG/deseq2/salmon_validation_deg.rds",
        discovery_vst_rds="{project}/DEG/deseq2/salmon_discovery_vst_matrix.rds",
        validation_vst_rds="{project}/DEG/deseq2/salmon_validation_vst_matrix.rds",
        samples=config["samples"],
    output:
        signature="{project}/DEG/classifier/salmon_lasso_signature_genes.tsv",
        coefficients="{project}/DEG/classifier/salmon_lasso_coefficients.tsv",
        discovery_predictions="{project}/DEG/classifier/salmon_discovery_predictions.tsv",
        validation_predictions="{project}/DEG/classifier/salmon_validation_predictions.tsv",
        roc_plot="{project}/DEG/classifier/salmon_lasso_roc_curve.png",
        summary="{project}/DEG/classifier/salmon_lasso_summary.txt",
    params:
        project=project,
        case_condition=case_condition,
        control_condition=control_condition,
        discovery_sample_type=discovery_sample_type,
        log2fc_threshold=config.get("deg", {}).get("log2fc", 2),
        padj_threshold=config.get("deg", {}).get("padj", 0.05),
    container:
        (
            "docker://btrspg/rlan:20251230"
            if config["container"].get("deg_vis", None) is None
            else config["container"].get("deg_vis", None)
        )
    benchmark:
        "benchmarks/{project}/deg_lasso_classifier_salmon.benchmark.txt"
    log:
        "logs/{project}/deg_lasso_classifier_salmon.log",
    threads: config["threads"].get("deg_vis", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("deg_vis", 16384),
    script:
        "../scripts/deg_lasso_classifier.R"


rule deg_lasso_classifier_kallisto:
    """
    LASSO classifier using Kallisto quantification results.
    """
    input:
        discovery_deg_rds="{project}/DEG/deseq2/kallisto_discovery_deg.rds",
        validation_deg_rds="{project}/DEG/deseq2/kallisto_validation_deg.rds",
        discovery_vst_rds="{project}/DEG/deseq2/kallisto_discovery_vst_matrix.rds",
        validation_vst_rds="{project}/DEG/deseq2/kallisto_validation_vst_matrix.rds",
        samples=config["samples"],
    output:
        signature="{project}/DEG/classifier/kallisto_lasso_signature_genes.tsv",
        coefficients="{project}/DEG/classifier/kallisto_lasso_coefficients.tsv",
        discovery_predictions="{project}/DEG/classifier/kallisto_discovery_predictions.tsv",
        validation_predictions="{project}/DEG/classifier/kallisto_validation_predictions.tsv",
        roc_plot="{project}/DEG/classifier/kallisto_lasso_roc_curve.png",
        summary="{project}/DEG/classifier/kallisto_lasso_summary.txt",
    params:
        project=project,
        case_condition=case_condition,
        control_condition=control_condition,
        discovery_sample_type=discovery_sample_type,
        log2fc_threshold=config.get("deg", {}).get("log2fc", 2),
        padj_threshold=config.get("deg", {}).get("padj", 0.05),
    container:
        (
            "docker://btrspg/rlan:20251230"
            if config["container"].get("deg_vis", None) is None
            else config["container"].get("deg_vis", None)
        )
    benchmark:
        "benchmarks/{project}/deg_lasso_classifier_kallisto.benchmark.txt"
    log:
        "logs/{project}/deg_lasso_classifier_kallisto.log",
    threads: config["threads"].get("deg_vis", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("deg_vis", 16384),
    script:
        "../scripts/deg_lasso_classifier.R"

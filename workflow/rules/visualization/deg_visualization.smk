# DEG Visualization Rules
# Create heatmaps, scatter plots, and boxplots for DEG results

# DEG Upset Plot
# Visualize overlap of DEGs from different tools (DESeq2, edgeR, limma-trend, limma-voom)
# for each quantification method (STAR_FC, salmon, kallisto)
# Creates a single upset plot with 8 sets (4 tools x 2 directions: up/down)


rule deg_upset_plot:
    input:
        deseq2="{project}/DEG/deseq2/{tool}/discovery_deg.tsv",
        edger="{project}/DEG/edger/{tool}/discovery_deg.tsv",
        limma_trend="{project}/DEG/limma_trend/{tool}/discovery_deg.tsv",
        limma_voom="{project}/DEG/limma_voom/{tool}/discovery_deg.tsv",
    output:
        upset_plot="{project}/visualization/DEG_{tool}_upset.pdf",
        upset_data="{project}/DEG/{tool}_upset_data.tsv",
        summary="{project}/DEG/{tool}_upset_summary.txt",
    params:
        project=project,
        tool="{tool}",
        log2fc_threshold=config.get("deg", {}).get("log2fc", 1),
        padj_threshold=config.get("deg", {}).get("padj", 0.05),
    container:
        (
            "docker://btrspg/rlan:20251229"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    log:
        "logs/{project}/deg_upset_{tool}.log",
    threads: config["threads"].get("deg_vis", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("deg_vis", 16384),
    script:
        "../../scripts/visualization/upset_plot.R"


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
        "../../scripts/visualization/deg_visualization.R"


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
        "../../scripts/visualization/deg_visualization.R"


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
        "../../scripts/visualization/deg_visualization.R"


rule deg_lasso_classifier:
    """
    Build LASSO logistic regression classifier using DEGs from discovery dataset
    and test performance on validation dataset.
    """
    input:
        discovery_deg_tsv="{project}/DEG/deseq2/discovery_deg.tsv",
        validation_deg_tsv="{project}/DEG/deseq2/validation_deg.tsv",
        expression_tpm="{project}/quantification/STAR_FC/TPM_matrix.txt",
    output:
        signature="{project}/DEG/classifier/lasso_signature_genes.tsv",
        coefficients="{project}/DEG/classifier/lasso_coefficients.tsv",
        discovery_predictions="{project}/DEG/classifier/discovery_predictions.tsv",
        validation_predictions="{project}/DEG/classifier/validation_predictions.tsv",
        roc_plot="{project}/DEG/classifier/lasso_roc_curve.png",
        summary="{project}/DEG/classifier/lasso_summary.txt",
    params:
        project=project,
        samples=config["samples"],
        case_condition=case_condition,
        control_condition=control_condition,
        discovery_sample_type=discovery_sample_type,
        log2fc_threshold=config.get("deg", {}).get("log2fc", 1),
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
        "../../scripts/visualization/classifier/lasso.R"


rule deg_lasso_classifier_salmon:
    """
    LASSO classifier using Salmon quantification results.
    """
    input:
        discovery_deg_tsv="{project}/DEG/deseq2/salmon_discovery_deg.tsv",
        validation_deg_tsv="{project}/DEG/deseq2/salmon_validation_deg.tsv",
        discovery_tpm="{project}/quantification/salmon/TPM_matrix.txt",
        validation_tpm="{project}/quantification/salmon/TPM_matrix.txt",
    output:
        signature="{project}/DEG/classifier/salmon_lasso_signature_genes.tsv",
        coefficients="{project}/DEG/classifier/salmon_lasso_coefficients.tsv",
        discovery_predictions="{project}/DEG/classifier/salmon_discovery_predictions.tsv",
        validation_predictions="{project}/DEG/classifier/salmon_validation_predictions.tsv",
        roc_plot="{project}/DEG/classifier/salmon_lasso_roc_curve.png",
        summary="{project}/DEG/classifier/salmon_lasso_summary.txt",
    params:
        project=project,
        samples=config["samples"],
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
        "../../scripts/visualization/classifier/lasso.R"


rule deg_lasso_classifier_kallisto:
    """
    LASSO classifier using Kallisto quantification results.
    """
    input:
        discovery_deg_tsv="{project}/DEG/deseq2/kallisto_discovery_deg.tsv",
        validation_deg_tsv="{project}/DEG/deseq2/kallisto_validation_deg.tsv",
        discovery_tpm="{project}/quantification/kallisto/TPM_matrix.txt",
        validation_tpm="{project}/quantification/kallisto/TPM_matrix.txt",
    output:
        signature="{project}/DEG/classifier/kallisto_lasso_signature_genes.tsv",
        coefficients="{project}/DEG/classifier/kallisto_lasso_coefficients.tsv",
        discovery_predictions="{project}/DEG/classifier/kallisto_discovery_predictions.tsv",
        validation_predictions="{project}/DEG/classifier/kallisto_validation_predictions.tsv",
        roc_plot="{project}/DEG/classifier/kallisto_lasso_roc_curve.png",
        summary="{project}/DEG/classifier/kallisto_lasso_summary.txt",
    params:
        project=project,
        samples=config["samples"],
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
        "../../scripts/visualization/classifier/lasso.R"


rule common_deg_heatmap:
    """
    ComplexHeatmap visualization of common DEGs across all 4 tools (DESeq2, edgeR,
    limma-trend, limma-voom). Shows expression heatmap using log10(TPM+1) with
    sample annotations from clinical data and gene annotations from binned
    log2FC and padj values. Only discovery samples are visualized.
    """
    input:
        deseq2="{project}/DEG/deseq2/{tool}/discovery_deg.tsv",
        edger="{project}/DEG/edger/{tool}/discovery_deg.tsv",
        limma_trend="{project}/DEG/limma_trend/{tool}/discovery_deg.tsv",
        limma_voom="{project}/DEG/limma_voom/{tool}/discovery_deg.tsv",
        tpm="{project}/quantification/{tool}/TPM_matrix.txt",
        samples=config["samples"],
        gene_name_map="resources/gene_id_to_gene_name.tsv",
    output:
        heatmap="{project}/visualization/common_DEGs_{tool}_heatmap.pdf",
        gene_list="{project}/visualization/common_DEGs_{tool}_gene_list.tsv",
        annotation_data="{project}/visualization/common_DEGs_{tool}_annotations.tsv",
    params:
        project=project,
        discovery_sample_type=discovery_sample_type,
        log2fc_threshold=config.get("deg", {}).get("log2fc", 1),
        padj_threshold=config.get("deg", {}).get("padj", 0.05),
        top_n=config.get("deg_vis", {}).get("common_deg_top_n", 100),
    container:
        (
            "docker://btrspg/rlan:20251229"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    log:
        "logs/{project}/common_deg_heatmap_{tool}.log",
    threads: config["threads"].get("deg_vis", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("deg_vis", 32768),
    script:
        "../../scripts/visualization/common_deg_heatmap.R"


rule pca_visualization:
    """
    PCA visualization using log10(TPM+1) transformed data.
    Performs PCA analysis on discovery samples and generates comprehensive plots
    including scree plot, pairs plot, biplot, loadings plot, and eigencorplot.
    """
    input:
        tpm="{project}/quantification/{tool}/TPM_matrix.txt",
        samples=config["samples"],
    output:
        png="{project}/visualization/PCA_{tool}_pca.png",
        pdf="{project}/visualization/PCA_{tool}_pca.pdf",
        pca_data="{project}/visualization/PCA_{tool}_pca_data.tsv",
        variance="{project}/visualization/PCA_{tool}_variance.tsv",
    params:
        project=project,
        discovery_sample_type=discovery_sample_type,
        color_by=config.get("deg_vis", {}).get("pca_color_by", "condition"),
        shape_by=config.get("deg_vis", {}).get("pca_shape_by", "sample_type"),
        removeVar=config.get("deg_vis", {}).get("pca_removeVar", 0.1),
    container:
        (
            "docker://btrspg/rlan:20251229"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    log:
        "logs/{project}/pca_visualization_{tool}.log",
    threads: config["threads"].get("deg_vis", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("deg_vis", 32768),
    script:
        "../../scripts/visualization/pca_visualization.R"

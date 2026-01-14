rule deg_upset_plot:
    input:
        deseq2="{project}/DEG/deseq2/{tool}/{dataset}_deg.tsv",
        edger="{project}/DEG/edger/{tool}/{dataset}_deg.tsv",
        limma_trend="{project}/DEG/limma_trend/{tool}/{dataset}_deg.tsv",
        limma_voom="{project}/DEG/limma_voom/{tool}/{dataset}_deg.tsv",
    output:
        upset_plot="{project}/visualization/{tool}_{dataset}_upset.pdf",
        upset_data="{project}/DEG/{tool}_{dataset}_upset_data.tsv",
        summary="{project}/DEG/{tool}_{dataset}_upset_summary.txt",
    params:
        project=get_project,
        log2fc_threshold=config.get("deg", {}).get("log2fc", 1),
        padj_threshold=config.get("deg", {}).get("padj", 0.05),
    container:
        (
            "docker://btrspg/rlan:20251229"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    log:
        "logs/{project}/deg_upset_{dataset}_{tool}.log",
    threads: config["threads"].get("deg_vis", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("deg_vis", 16384),
    script:
        "../../scripts/visualization/upset_plot.R"


# rule volcano_vis:
#     input:
#         discovery_deg_rds="{project}/DEG/deseq2/discovery_deg.rds",
#         validation_deg_rds="{project}/DEG/deseq2/validation_deg.rds",
#         geneid_to_genename="resources/gene_id_to_gene_name.tsv",
#     output:
#         discovery_png="{project}/visualization/Volcano_discovery.png",
#         discovery_pdf="{project}/visualization/Volcano_discovery.pdf",
#         validation_png="{project}/visualization/Volcano_validation.png",
#         validation_pdf="{project}/visualization/Volcano_validation.pdf",
#     container:
#         (
#             "docker://btrspg/rlan:20251027"
#             if config["container"].get("r", None) is None
#             else config["container"].get("r", None)
#         )
#     log:
#         "logs/{project}/visualization_Volcano.diffexp.log",
#     threads: config["threads"].get("default", 1)
#     resources:
#         mem_mb=config["resources"]["mem_mb"].get("default", 4096),
#     benchmark:
#         "benchmarks/{project}/Volcano.benchmark.txt"
#     script:
#         "../../../scripts/visualization/volcano.R"


rule pcatools_vis_database:
    input:
        counts="{project}/quantification/{tool}/{dataset}_count_matrix.txt",
        qc_files=expand(
            "{project}/qc/qualimap-rnaseq/{sample}/rnaseq_qc_results.txt",
            project=project,
            sample=samples.index.tolist(),
        ),
    output:
        pdf="{project}/visualization/{tool}_{dataset}_pca.pdf",
        png="{project}/visualization/{tool}_{dataset}_pca.png",
        clinical_info="{project}/visualization/{tool}_{dataset}_pca_clinical_info.tsv",
    params:
        samples=config["samples"],
        strandness=lambda w, input: get_samples_strandness(input.qc_files),
        dataset=get_dataset,
        project=get_project,
        case_condition=get_case_condition,
        control_condition=get_control_condition,
    container:
        (
            "docker://btrspg/rlan:20251027"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    log:
        "logs/{project}/{tool}_{dataset}_vis_pca.log",
    threads: config["threads"].get("deseq2", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("deseq2", 8192),
    script:
        "../../scripts/visualization/pca.R"


rule pcatools_vis_database_batchcorrected:
    input:
        counts="{project}/quantification/{tool}/{dataset}_count_matrix_corrected.txt",
        qc_files=expand(
            "{project}/qc/qualimap-rnaseq/{sample}/rnaseq_qc_results.txt",
            project=project,
            sample=samples.index.tolist(),
        ),
    output:
        pdf="{project}/visualization/{tool}_{dataset}_pca_batchcorrected.pdf",
        png="{project}/visualization/{tool}_{dataset}_pca_batchcorrected.png",
        clinical_info="{project}/visualization/{tool}_{dataset}_pca_batchcorrected_clinical_info.tsv",
    params:
        samples=config["samples"],
        strandness=lambda w, input: get_samples_strandness(input.qc_files),
        dataset=get_dataset,
        project=get_project,
        case_condition=get_case_condition,
        control_condition=get_control_condition,
    container:
        (
            "docker://btrspg/rlan:20251027"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    log:
        "logs/{project}/{tool}_{dataset}_vis_pca_batchcorrected.log",
    threads: config["threads"].get("deseq2", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("deseq2", 8192),
    script:
        "../../scripts/visualization/pca.R"


rule extract_strandness_to_samples:
    """
    Extract strandness information from QualiMap RNA-seq QC results and
    add it as a new column to the samples.tsv file.

    This rule reads the rnaseq_qc_results.txt files from QualiMap RNA-seq,
    extracts the strand specificity percentage for each sample, and adds
    it as a 'strandness' column to an updated samples file.
    """
    input:
        qualimap=expand(
            "{project}/qc/qualimap-rnaseq/{sample}/rnaseq_qc_results.txt",
            project=project,
            sample=samples.index.tolist(),
        ),
        samples=config["samples"],
    output:
        samples_with_strandness="{project}/config/samples_with_strandness.tsv",
    params:
        project=project,
    container:
        (
            "docker://btrspg/rlan:20251027"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    log:
        "logs/{project}/extract_strandness_to_samples.log",
    threads: config["threads"].get("default", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("default", 4096),
    script:
        "../../scripts/qc/extract_strandness_to_samples.R"


rule deg_summary_thresholds:
    """
    Generate comprehensive DEG summaries across multiple thresholds.

    Creates:
    1. Summary table with DEG counts for all threshold combinations
    2. Upset-style data for overlap analysis
    3. Heatmap of DEG counts across thresholds
    4. Bar plot comparing methods across thresholds

    Thresholds tested:
    - padj: 0.05, 0.01, 0.001, 0.0001
    - log2FC: 1, 1.2, 1.5, 2
    """
    input:
        deseq2="{project}/DEG/deseq2/{tool}/{dataset}_deg.tsv",
        edger="{project}/DEG/edger/{tool}/{dataset}_deg.tsv",
        limma_trend="{project}/DEG/limma_trend/{tool}/{dataset}_deg.tsv",
        limma_voom="{project}/DEG/limma_voom/{tool}/{dataset}_deg.tsv",
    output:
        summary_table="{project}/DEG/{tool}_{dataset}_summary_thresholds.tsv",
        upset_data="{project}/DEG/{tool}_{dataset}_upset_data.tsv",
        heatmap="{project}/visualization/DEG_{tool}_{dataset}_thresholds_heatmap.pdf",
        comparison_plot="{project}/visualization/DEG_{tool}_{dataset}_thresholds_comparison.pdf",
    params:
        project=project,
        dataset="{dataset}",
        tool="{tool}",
    container:
        (
            "docker://btrspg/rlan:20251229"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    log:
        "logs/{project}/deg_summary_thresholds_{tool}_{dataset}.log",
    threads: config["threads"].get("deg_vis", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("deg_vis", 16384),
    script:
        "../../scripts/analysis/deg/deg_summary_thresholds.R"

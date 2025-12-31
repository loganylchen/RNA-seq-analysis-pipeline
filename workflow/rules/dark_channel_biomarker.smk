# Dark Channel Biomarker (DCB) Analysis
# DCB genes are defined as genes that are:
# 1. Undetected in control/normal samples (low-noise regions)
# 2. Recurrently detected in case/cancer samples
# Based on cfRNA cancer detection methodology (Circulating Cell-free Genome Atlas)


rule dcb_analysis:
    """
    Detect Dark Channel Biomarker (DCB) genes from RNA-seq count/expression data.
    DCB genes are tissue- and cancer-specific genes detected in cancer samples
    but absent in non-cancer individuals within low-noise regions.
    """
    input:
        # Use TPM matrix for detection rate analysis
        tpm="{project}/quantification/STAR_FC/TPM_matrix.txt",
        discovery_deg_tsv="{project}/DEG/deseq2/discovery_deg.tsv",
        validation_deg_tsv="{project}/DEG/deseq2/validation_deg.tsv",
    output:
        discovery_dcb_tsv="{project}/dcb/discovery_dcb.tsv",
        discovery_dcb_rds="{project}/dcb/discovery_dcb.rds",
        validation_dcb_tsv="{project}/dcb/validation_dcb.tsv",
        validation_dcb_rds="{project}/dcb/validation_dcb.rds",
        discovery_summary="{project}/dcb/discovery_summary.txt",
        validation_summary="{project}/dcb/validation_summary.txt",
        discovery_plot="{project}/dcb/discovery_dcb_distribution.png",
        validation_plot="{project}/dcb/validation_dcb_distribution.png",
    params:
        samples=config["samples"],
        project=project,
        case_condition=case_condition,
        control_condition=control_condition,
        discovery_sample_type=discovery_sample_type,
        log2fc=config.get("dcb", {}).get("log2fc", 1),
        padj=config.get("dcb", {}).get("padj", 0.05),
        control_tpm_threshold=config.get("dcb", {}).get("control_tpm_threshold", 1),
        case_tpm_threshold=config.get("dcb", {}).get("case_tpm_threshold", 1),
        control_detection_rate=config.get("dcb", {}).get("control_detection_rate", 0.1),
        case_detection_rate=config.get("dcb", {}).get("case_detection_rate", 0.3),
    container:
        (
            "docker://btrspg/rlan:20251027"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    benchmark:
        "benchmarks/{project}/dcb_analysis.benchmark.txt"
    log:
        "logs/{project}/dcb_analysis.log",
    threads: config["threads"].get("dcb", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("dcb", 8192),
    script:
        "../scripts/dcb_analysis.R"


rule dcb_analysis_salmon:
    """
    DCB analysis using Salmon quantification results.
    """
    input:
        tpm="{project}/quantification/salmon/TPM_matrix.txt",
        count_matrix="{project}/quantification/salmon/count_matrix.txt",
    output:
        discovery_dcb_tsv="{project}/dcb/salmon_discovery_dcb.tsv",
        discovery_dcb_rds="{project}/dcb/salmon_discovery_dcb.rds",
        validation_dcb_tsv="{project}/dcb/salmon_validation_dcb.tsv",
        validation_dcb_rds="{project}/dcb/salmon_validation_dcb.rds",
        discovery_summary="{project}/dcb/salmon_discovery_summary.txt",
        validation_summary="{project}/dcb/salmon_validation_summary.txt",
        discovery_plot="{project}/dcb/salmon_discovery_dcb_distribution.png",
        validation_plot="{project}/dcb/salmon_validation_dcb_distribution.png",
    params:
        samples=config["samples"],
        project=project,
        case_condition=case_condition,
        control_condition=control_condition,
        discovery_sample_type=discovery_sample_type,
        log2fc=config.get("dcb", {}).get("log2fc", 1),
        control_tpm_threshold=config.get("dcb", {}).get("control_tpm_threshold", 1),
        case_tpm_threshold=config.get("dcb", {}).get("case_tpm_threshold", 1),
        control_detection_rate=config.get("dcb", {}).get("control_detection_rate", 0.1),
        case_detection_rate=config.get("dcb", {}).get("case_detection_rate", 0.3),
    container:
        (
            "docker://btrspg/rlan:20251027"
            if config["container"].get("dcb", None) is None
            else config["container"].get("dcb", None)
        )
    benchmark:
        "benchmarks/{project}/dcb_analysis_salmon.benchmark.txt"
    log:
        "logs/{project}/dcb_analysis_salmon.log",
    threads: config["threads"].get("dcb", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("dcb", 8192),
    script:
        "../scripts/dcb_analysis.R"


rule dcb_analysis_kallisto:
    """
    DCB analysis using Kallisto quantification results.
    """
    input:
        tpm="{project}/quantification/kallisto/TPM_matrix.txt",
        count_matrix="{project}/quantification/kallisto/count_matrix.txt",
    output:
        discovery_dcb_tsv="{project}/dcb/kallisto_discovery_dcb.tsv",
        discovery_dcb_rds="{project}/dcb/kallisto_discovery_dcb.rds",
        validation_dcb_tsv="{project}/dcb/kallisto_validation_dcb.tsv",
        validation_dcb_rds="{project}/dcb/kallisto_validation_dcb.rds",
        discovery_summary="{project}/dcb/kallisto_discovery_summary.txt",
        validation_summary="{project}/dcb/kallisto_validation_summary.txt",
        discovery_plot="{project}/dcb/kallisto_discovery_dcb_distribution.png",
        validation_plot="{project}/dcb/kallisto_validation_dcb_distribution.png",
    params:
        samples=config["samples"],
        project=project,
        case_condition=case_condition,
        control_condition=control_condition,
        discovery_sample_type=discovery_sample_type,
        log2fc=config.get("dcb", {}).get("log2fc", 1),
        control_tpm_threshold=config.get("dcb", {}).get("control_tpm_threshold", 1),
        case_tpm_threshold=config.get("dcb", {}).get("case_tpm_threshold", 1),
        control_detection_rate=config.get("dcb", {}).get("control_detection_rate", 0.1),
        case_detection_rate=config.get("dcb", {}).get("case_detection_rate", 0.3),
    container:
        (
            "docker://btrspg/rlan:20251027"
            if config["container"].get("dcb", None) is None
            else config["container"].get("dcb", None)
        )
    benchmark:
        "benchmarks/{project}/dcb_analysis_kallisto.benchmark.txt"
    log:
        "logs/{project}/dcb_analysis_kallisto.log",
    threads: config["threads"].get("dcb", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("dcb", 8192),
    script:
        "../scripts/dcb_analysis.R"

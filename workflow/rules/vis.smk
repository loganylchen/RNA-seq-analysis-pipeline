rule volcano_vis:
    input:
        discovery_deg_rds="{project}/DEG/deseq2/discovery_deg.rds",
        validation_deg_rds="{project}/DEG/deseq2/validation_deg.rds",
        geneid_to_genename="resources/gene_id_to_gene_name.tsv",
    output:
        discovery_png="{project}/visualization/Volcano_discovery.png",
        discovery_pdf="{project}/visualization/Volcano_discovery.pdf",
        validation_png="{project}/visualization/Volcano_validation.png",
        validation_pdf="{project}/visualization/Volcano_validation.pdf",
    container:
        (
            "docker://btrspg/rlan:20251027"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    log:
        "logs/{project}/visualization_Volcano.diffexp.log",
    threads: config["threads"].get("default", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("default", 4096),
    benchmark:
        "benchmarks/{project}/Volcano.benchmark.txt"
    script:
        "../scripts/volcano.R"


rule pcatools_vis:
    input:
        counts="{project}/quantification/STAR_FC/count_matrix.txt",
    output:
        pdf="{project}/visualization/pca.pdf",
        png="{project}/visualization/pca.png",
        discovery_pdf="{project}/visualization/discovery_pca.pdf",
        discovery_png="{project}/visualization/discovery_pca.png",
        validation_pdf="{project}/visualization/validation_pca.pdf",
        validation_png="{project}/visualization/validation_pca.png",
    params:
        samples=config["samples"],
        project=project,
        case_condition=case_condition,
        control_condition=control_condition,
        discovery_sample_type=discovery_sample_type,
    container:
        (
            "docker://btrspg/rlan:20251027"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    log:
        "logs/{project}/vis_pca.log",
    threads: config["threads"].get("deseq2", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("deseq2", 8192),
    script:
        "../scripts/PCA.R"


rule pcatools_vis_salmon:
    input:
        counts="{project}/quantification/salmon/count_matrix.txt",
    output:
        pdf="{project}/visualization/salmon_pca.pdf",
        png="{project}/visualization/salmon_pca.png",
        discovery_pdf="{project}/visualization/salmon_discovery_pca.pdf",
        discovery_png="{project}/visualization/salmon_discovery_pca.png",
        validation_pdf="{project}/visualization/salmon_validation_pca.pdf",
        validation_png="{project}/visualization/salmon_validation_pca.png",
    params:
        samples=config["samples"],
        project=project,
        case_condition=case_condition,
        control_condition=control_condition,
        discovery_sample_type=discovery_sample_type,
    container:
        (
            "docker://btrspg/rlan:20251027"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    log:
        "logs/{project}/vis_pca_salmon.log",
    threads: config["threads"].get("deseq2", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("deseq2", 8192),
    script:
        "../scripts/PCA.R"


rule pcatools_vis_kallisto:
    input:
        counts="{project}/quantification/kallisto/count_matrix.txt",
    output:
        pdf="{project}/visualization/kallisto_pca.pdf",
        png="{project}/visualization/kallisto_pca.png",
        discovery_pdf="{project}/visualization/kallisto_discovery_pca.pdf",
        discovery_png="{project}/visualization/kallisto_discovery_pca.png",
        validation_pdf="{project}/visualization/kallisto_validation_pca.pdf",
        validation_png="{project}/visualization/kallisto_validation_pca.png",
    params:
        samples=config["samples"],
        project=project,
        case_condition=case_condition,
        control_condition=control_condition,
        discovery_sample_type=discovery_sample_type,
    container:
        (
            "docker://btrspg/rlan:20251027"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    log:
        "logs/{project}/vis_pca_kallisto.log",
    threads: config["threads"].get("deseq2", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("deseq2", 8192),
    script:
        "../scripts/PCA.R"

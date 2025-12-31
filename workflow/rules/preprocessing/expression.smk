rule count_matrix:
    input:
        expand(
            "{project}/quantification/STAR/{sample}/{sample}.ReadsPerGene.out.tab",
            project=project,
            sample=samples.index.tolist(),
        ),
    output:
        "{project}/quantification/STAR_count_matrix.txt",
        "{project}/quantification/STAR_count_matrix_un.txt",
        "{project}/quantification/STAR_count_matrix_strand.txt",
        "{project}/quantification/STAR_count_matrix_reverse.txt",
        "{project}/quantification/STAR_strandness.pdf",
    log:
        "logs/{project}/count-matrix.log",
    params:
        samples=samples.index.tolist(),
    container:
        (
            "docker://btrspg/python3:20251024"
            if config["container"].get("python3", None) is None
            else config["container"].get("python3", None)
        )
    script:
        "../scripts/count-matrix.py"


rule count_matrix_star_FC:
    input:
        expand(
            "{project}/quantification/featurecounts/{sample}.txt",
            project=project,
            sample=samples.index.tolist(),
        ),
    output:
        count_matrix="{project}/quantification/STAR_FC/count_matrix.txt",
        puree_count_matrix="{project}/quantification/STAR_FC/count_matrix_PUREE.txt",
    log:
        "logs/{project}/count-matrix_star2fc.log",
    params:
        samples=samples.index.tolist(),
    container:
        (
            "docker://btrspg/rlan:20251110"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    script:
        "../scripts/count-matrix.R"


rule TPM_matrix_star_FC_RAW:
    input:
        expand(
            "{project}/quantification/featurecounts/{sample}.txt",
            project=project,
            sample=samples.index.tolist(),
        ),
    output:
        tpm_matrix="{project}/quantification/STAR_FC/TPM_matrix.txt",
    log:
        "logs/{project}/tpm-matrix_star2fc.log",
    params:
        samples=samples.index.tolist(),
    container:
        (
            "docker://btrspg/rlan:20251110"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    script:
        "../scripts/TPM-matrix.R"


rule TPM_matrix_star_FC:
    input:
        expand(
            "{project}/quantification/featurecounts_novel/{sample}.txt",
            project=project,
            sample=samples.index.tolist(),
        ),
    output:
        tpm_matrix="{project}/quantification/STAR_FC4splicetool/TPM_matrix.txt",
    log:
        "logs/{project}/tpm-matrix_star2fc4splicetool.log",
    params:
        samples=samples.index.tolist(),
    container:
        (
            "docker://btrspg/rlan:20251110"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    script:
        "../scripts/TPM-matrix.R"


rule TPM_separate_matrix:
    input:
        tpm_matrix="{project}/quantification/STAR_FC4splicetool/TPM_matrix.txt",
    output:
        discovery_tpm_matrix="{project}/quantification/STAR_FC4splicetool/Discovery_TPM_matrix.txt",
    log:
        "logs/{project}/tpm-matrix_star2fc4splicetool_sep.log",
    params:
        case_samples=discovery_case_samples.index.tolist(),
        control_samples=discovery_control_samples.index.tolist(),
    container:
        (
            "docker://btrspg/rlan:20251110"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    script:
        "../scripts/TPM-matrix_sep.R"


rule count_matrix_salmon:
    input:
        expand(
            "{project}/quantification/salmon/{sample}/quant.sf",
            project=project,
            sample=samples.index.tolist(),
        ),
    output:
        count_matrix="{project}/quantification/salmon/count_matrix.txt",
    log:
        "logs/{project}/count-matrix_salmon.log",
    params:
        samples=samples.index.tolist(),
    container:
        (
            "docker://btrspg/rlan:20251110"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    script:
        "../scripts/count-matrix_salmon.R"


rule count_matrix_kallisto:
    input:
        expand(
            "{project}/quantification/kallisto/{sample}/abundance.tsv",
            project=project,
            sample=samples.index.tolist(),
        ),
    output:
        count_matrix="{project}/quantification/kallisto/count_matrix.txt",
    log:
        "logs/{project}/count-matrix_kallisto.log",
    params:
        samples=samples.index.tolist(),
    container:
        (
            "docker://btrspg/rlan:20251110"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    script:
        "../scripts/count-matrix_kallisto.R"

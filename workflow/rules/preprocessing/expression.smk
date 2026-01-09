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
        "../../scripts/quantification/count_matrix.py"


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
        "../../scripts/quantification/count_matrix.R"


rule count_matrix_star_dataset:
    input:
        expand(
            "{project}/quantification/featurecounts/{sample}.txt",
            project=project,
            sample=samples.index.tolist(),
        ),
    output:
        count_matrix="{project}/quantification/STAR_FC/{dataset}_count_matrix.txt",
        puree_count_matrix="{project}/quantification/STAR_FC/{dataset}_count_matrix_PUREE.txt",
    log:
        "logs/{project}/count-matrix_star2fc_{dataset}.log",
    params:
        samples=get_dataset_samples,
    container:
        (
            "docker://btrspg/rlan:20251110"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    script:
        "../../scripts/quantification/count_matrix.R"


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
        "../../scripts/quantification/tpm_matrix.R"


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
        "../../scripts/quantification/tpm_matrix.R"


rule count_matrix_salmon_dataset:
    input:
        expand(
            "{project}/quantification/salmon/{sample}/quant.sf",
            project=project,
            sample=samples.index.tolist(),
        ),
    output:
        count_matrix="{project}/quantification/salmon/{dataset}_count_matrix.txt",
    log:
        "logs/{project}/count-matrix_salmon_{dataset}.log",
    params:
        samples=get_dataset_samples,
    container:
        (
            "docker://btrspg/rlan:20251110"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    script:
        "../../scripts/quantification/count_matrix_salmon.R"


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
        "../../scripts/quantification/count_matrix_salmon.R"


rule TPM_matrix_salmon:
    input:
        expand(
            "{project}/quantification/salmon/{sample}/quant.sf",
            project=project,
            sample=samples.index.tolist(),
        ),
    output:
        tpm_matrix="{project}/quantification/salmon/TPM_matrix.txt",
    log:
        "logs/{project}/tpm-matrix_salmon.log",
    params:
        samples=samples.index.tolist(),
    container:
        (
            "docker://btrspg/rlan:20251110"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    script:
        "../../scripts/quantification/tpm_matrix_salmon.R"


rule count_matrix_kallisto_dataset:
    input:
        expand(
            "{project}/quantification/kallisto/{sample}/abundance.tsv",
            project=project,
            sample=samples.index.tolist(),
        ),
    output:
        count_matrix="{project}/quantification/kallisto/{dataset}_count_matrix.txt",
    log:
        "logs/{project}/count-matrix_kallisto_{dataset}.log",
    params:
        samples=get_dataset_samples,
    container:
        (
            "docker://btrspg/rlan:20251110"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    script:
        "../../scripts/quantification/count_matrix_kallisto.R"


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
        "../../scripts/quantification/count_matrix_kallisto.R"


rule TPM_matrix_kallisto:
    input:
        expand(
            "{project}/quantification/kallisto/{sample}/abundance.tsv",
            project=project,
            sample=samples.index.tolist(),
        ),
    output:
        tpm_matrix="{project}/quantification/kallisto/TPM_matrix.txt",
    log:
        "logs/{project}/tpm-matrix_kallisto.log",
    params:
        samples=samples.index.tolist(),
    container:
        (
            "docker://btrspg/rlan:20251110"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    script:
        "../../scripts/quantification/tpm_matrix_kallisto.R"


# rule sva_remove_batch_effect:
#     input:
#         count_matrix="{project}/quantification/{tool}/{dataset}_count_matrix.txt",
#     output:
#         corrected_matrix="{project}/quantification/{tool}/{dataset}_corrected_count_matrix.txt",
#     log:
#         "logs/{project}/sva_remove_batch_effect_{tool}_{dataset}.log",
#     params:
#         samples=get_dataset_samples,
#     container:
#         (
#             "docker://btrspg/rlan:20251110"
#             if config["container"].get("r", None) is None
#             else config["container"].get("r", None)
#         )
#     script:
#         "../../scripts/quantification/sva_remove_batch_effect.R"

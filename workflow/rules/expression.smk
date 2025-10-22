rule count_matrix:
    input:
        expand(
            "{project}/quantification/{sample}/{sample}.ReadsPerGene.out.tab",
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
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"


rule count_matrix_star_FC:
    input:
        expand(
            "{project}/quantification/{sample}.star_counts.txt",
            project=project,
            sample=samples.index.tolist(),
        ),
    output:
        "{project}/quantification/STAR_fc_count_matrix.txt",
    log:
        "logs/{project}/count-matrix_star2fc.log",
    params:
        samples=samples.index.tolist(),
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix_fc.py"

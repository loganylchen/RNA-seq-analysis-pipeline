rule fastp:
    input:
        fq1=(
            "{project}/data/{sample}/{sample}_1.fastq.gz"
            if is_pe
            else "{project}/data/{sample}/{sample}.fastq.gz"
        ),
        fq2="{project}/data/{sample}/{sample}_2.fastq.gz",
    output:
        fq1=(
            "{project}/clean_data/{sample}_1.fastq.gz"
            if is_pe
            else "{project}/clean_data/{sample}.fastq.gz"
        ),
        fq2="{project}/clean_data/{sample}_2.fastq.gz",
        html="{project}/report/fastp/{sample}.fastp.html",
        json="{project}/qc/fastp/{sample}/{sample}.fastp.json",
    log:
        "logs/{project}/fastp_{sample}.log",
    container:
        (
            "docker://btrspg/fastp:0.24.0"
            if config["container"].get("fastp", None) is None
            else config["container"].get("fastp", None)
        )
    params:
        pe=is_pe,
        extra=config["fastp"]["extra"],
    threads: config["threads"].get("fastp", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("fastp", 8192),
    script:
        "../scripts/fastp.sh"

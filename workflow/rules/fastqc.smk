rule fastp:
    input:
        unpack(get_raw_fq),
    output:
        fq1=(
            "{project}/clean_data/{sample}_1.fastq.gz"
            if is_pe
            else "{project}/clean_data/{sample}.fastq.gz"
        ),
        fq2="{project}/clean_data/{sample}_2.fastq.gz",
        html="{project}/report/{sample}.fastp.html",
        json="{project}/qc/{sample}/{sample}.fastp.json",
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
    threads: config["threads"]["fastp"]
    script:
        "../scripts/fastp.sh"

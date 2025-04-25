rule fastp:
    input:
        sample=["{project}/data/{sample}/{sample}_1.fastq.gz", "{project}/data/{sample}/{sample}_2.fastq.gz"] if is_pe else ["{project}/data/{sample}/{sample}.fastq.gz"]
    output:
        trimmed=["{project}/clean_data/{sample}_1.fastq.gz","{project}/clean_data/{sample}_2.fastq.gz"] if is_pe else "{project}/clean_data/{sample}.fastq.gz",
        html="{project}/report/{sample}.fastp.html",
        json="{project}/report/{sample}.fastp.json"
    log:
        "logs/{project}/fastp_{sample}.log"
    container:
        "docker://btrspg/fastp:0.24.0"
    params:
        extra=config['fastp']['extra']
    threads: config['threads']['fastp']
    script:
        "../scripts/fastp.sh"

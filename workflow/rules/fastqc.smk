rule fastp:
    input:
        sample=["{project}/data/{sample}_1.fastq.gz", "{project}/data/{sample}_2.fastq.gz"] if is_pe else ["{project}/data/{sample}.fastq.gz"]
    output:
        trimmed=["{project}/clean_data/{sample}_1.fastq.gz","{project}/clean_data/{sample}_2.fastq.gz"] if is_pe else "{project}/clean_data/{sample}.fastq.gz",
        html="{project}/report/{sample}.fastp.html",
        json="{project}/report/{sample}.fastp.json"
    log:
        "logs/{project}/fastp_{sample}.log"
    params:
        extra=config['fastp']['pe_extra']
    threads: config['threads']['fastp']
    wrapper:
        "v6.0.0/bio/fastp"



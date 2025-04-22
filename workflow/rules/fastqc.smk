rule fastp_se:
    input:
        unpack(get_raw_fq_4qc)
    output:
        unpack(get_trimmed_data),
        html="{project}/report/{sample}.fastp.html",
        json="{project}/report/{sample}.fastp.json",
    log:
        "logs/{project}/fastp_{sample}.log"
    params:
        extra=config['fastp']['se_extra']
    threads: config['threads']['fastp']
    wrapper:
        "v6.0.0/bio/fastp"

rule fastp_pe:
    input:
        unpack(get_raw_fq_4qc)
    output:
        unpack(get_trimmed_data),
        html="{project}/report/{sample}.fastp.html",
        json="{project}/report/{sample}.fastp.json"
    log:
        "logs/{project}/fastp_{sample}.log"
    params:
        extra=config['fastp']['pe_extra']
    threads: config['threads']['fastp']
    wrapper:
        "v6.0.0/bio/fastp"



rule preparing_rmats:
    input:
        case_bams=expand(
            "{project}/alignment/STAR/{sample}/{sample}.bam",
            project=project,
            sample=discovery_case_samples.index.tolist(),
        ),
        control_bams=expand(
            "{project}/alignment/STAR/{sample}/{sample}.bam",
            project=project,
            sample=discovery_control_samples.index.tolist(),
        ),
    output:
        case_bam_list_f="{project}/transcript_splicing/rmats-temp/case.list",
        control_bam_list_f="{project}/transcript_splicing/rmats-temp/control.list",
    log:
        "logs/{project}/rmats_sample_list.log",
    threads: config["threads"].get("default", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("default", 4096),
    shell:
        "echo {input.case_bams} | sed 's/ /,/g' > {output.case_bam_list_f}; "
        "echo {input.control_bams} | sed 's/ /,/g' > {output.control_bam_list_f}; "
        "echo `date` > {log}"


rule splicing_rmats:
    input:
        case_bam_list_f="{project}/transcript_splicing/rmats-temp/case.list",
        control_bam_list_f="{project}/transcript_splicing/rmats-temp/control.list",
        gtf="{project}/assembly/stringtie/gffcompare.annotated.gtf",
    output:
        outdir=directory("{project}/transcript_splicing/rmats/"),
        ri_jcec="{project}/transcript_splicing/rmats/RI.MATS.JCEC.txt",
        se_jcec="{project}/transcript_splicing/rmats/SE.MATS.JCEC.txt",
        temp_dir=temp(directory("{project}/transcript_splicing/rmats_temp")),
    log:
        "logs/{project}/splicing_rmats.log",
    container:
        (
            "docker://btrspg/rmatsturbo:4.3.0"
            if config["container"].get("rmatsturbo", None) is None
            else config["container"].get("rmatsturbo", None)
        )
    params:
        extra=config["rmats"]["extra"],
    threads: config["threads"].get("rmats", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("rmats", 16384),
    shell:
        "rmats.py "
        "--b1 {input.case_bam_list_f} "
        "--b2 {input.control_bam_list_f} "
        "--gtf {input.gtf} "
        "{params.extra} "
        "--nthread {threads} "
        "--od {output.outdir} "
        "--tmp {output.temp_dir} &>{log}"


rule splicetools:
    input:
        ri_jcec="{project}/transcript_splicing/rmats/RI.MATS.JCEC.txt",
        se_jcec="{project}/transcript_splicing/rmats/SE.MATS.JCEC.txt",
        input_dir="{project}/transcript_splicing/rmats/",
        annotation_bed12="{project}/assembly/stringtie/gffcompare.sorted.bed",
        genome_fasta="resources/genome.fasta",
        tpm_file="{project}/quantification/STAR_FC4splicetool/Discovery_TPM_matrix.txt",
    output:
        output_dir=directory("{project}/transcript_splicing/splicetools/"),
    log:
        "logs/{project}/splicetools.log",
    container:
        (
            "docker://btrspg/splicetools:c9fd38227fcdf43d1e08e919480372751f2ee5a4"
            if config["container"].get("splicetools", None) is None
            else config["container"].get("splicetools", None)
        )
    params:
        control_tpm_threshold=config["splicetools"].get("control_tpm_threshold", 1),
        case_tpm_threshold=config["splicetools"].get("case_tpm_threshold", 1),
        control_n=len(discovery_control_samples.index),
        case_n=len(discovery_case_samples.index),
        fdr=config["splicetools"].get("fdr", 0.05),
    threads: config["threads"].get("splicetools", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("splicetools", 8192),
    script:
        "../scripts/splicetools.sh"

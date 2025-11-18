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
    threads: 1
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
    threads: config["threads"]["rmats"]
    shell:
        "rmats.py "
        "--b1 {input.case_bam_list_f} "
        "--b2 {input.control_bam_list_f} "
        "--gtf {input.gtf} "
        "{params.extra} "
        "--nthread {threads} "
        "--od {output.outdir} "
        "--tmp {output.temp_dir} &>{log}"

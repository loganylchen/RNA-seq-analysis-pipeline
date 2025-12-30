rule modtect_mod:
    input:
        bam="{project}/alignment/STAR/{sample}/{sample}.split.bam",
        genome="resources/genome.fasta",
        bed="resources/transcriptome.bed",
    output:
        output="{project}/modification/modtect/{sample}/{sample}.modtect.combined.txt",
        temp=temp(
            directory("{project}/modification/modtect/{sample}/{sample}.modtect_tmp")
        ),
    log:
        log="logs/{project}/modtect_{sample}.log",
    params:
        label=lambda wc, output: output.output.replace(".combined.txt", ""),
    benchmark:
        "benchmarks/{project}/{sample}.modtect.benchmark.txt"
    container:
        (
            "docker://btrspg/modtect:20251029"
            if config["container"].get("modtect", None) is None
            else config["container"].get("modtect", None)
        )
    threads: config["threads"].get("modtect", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("modtect", 8192),
    shell:
        "modtect.py {input.bam} {input.genome} 0 0 0 "
        "--threads {threads} --regionFile {input.bed} "
        "--label {params.label}"
        "&>{log}"


rule modtect_mod_merge:
    input:
        bam=expand(
            "{project}/modification/modtect/{sample}/{sample}.modtect.combined.txt",
            project=project,
            sample=samples.index.tolist(),
        ),
    output:
        output="{project}/modification/modtect/merged.modtect.txt",
    log:
        log="logs/{project}/modtect_merged.log",
    benchmark:
        "benchmarks/{project}/modtect_merged.benchmark.txt"
    container:
        (
            "docker://btrspg/rlan:20251110"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    threads: config["threads"].get("default", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("default", 4096),
    script:
        "../scripts/modtect_merge.R"


rule modtect_statistical_analysis:
    input:
        merged="{project}/modification/modtect/merged.modtect.txt",
    output:
        results="{project}/modification/modtect/modtect_statistical_results.tsv",
        significant="{project}/modification/modtect/modtect_significant_modifications.tsv",
        summary="{project}/modification/modtect/modtect_analysis_summary.txt",
    params:
        samples=config["samples"],
        project=project,
        case_condition=case_condition,
        control_condition=control_condition,
        discovery_sample_type=discovery_sample_type,
    log:
        log="logs/{project}/modtect_statistical_analysis.log",
    container:
        (
            "docker://btrspg/rlan:20251110"
            if config["container"].get("r", None) is None
            else config["container"].get("r", None)
        )
    threads: config["threads"].get("default", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("default", 8192),
    script:
        "../scripts/modtect_statistical_analysis.R"

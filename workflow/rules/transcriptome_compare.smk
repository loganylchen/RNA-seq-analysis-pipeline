rule transcript_compare_gffcompare:
    input:
        predicted_gtf="{project}/assembly/merged.gtf",
        reference_gtf="resources/genome.gtf",
    output:
        cmp="results/gffcompare/gffcompare.combined.gtf",
        stats="results/gffcompare/gffcompare.stats",
    params:
        out_prefix="results/gffcompare/gffcompare",
    shell:
        """
        mkdir -p results/gffcompare
        gffcompare -r {input.reference_gtf} -o {params.out_prefix} {input.predicted_gtf}
        """

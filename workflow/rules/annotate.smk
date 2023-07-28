rule annotate_variants:
    input:
        calls="variants.bcf",  # .vcf, .vcf.gz or .bcf
        cache="resources/vep/cache",  # can be omitted if fasta and gff are specified
        plugins="resources/vep/plugins",
    output:
        calls="variants.annotated.bcf",  # .vcf, .vcf.gz or .bcf
        stats="variants.html",
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=["LoFtool"],
        extra="--everything",  # optional: extra arguments
    log:
        "logs/vep/annotate.log",
    threads: 4
    wrapper:
        "v1.25.0/bio/vep/annotate"

rule geneid_to_genename:
    input:
        gtf="resources/genome.gtf",
    output:
        tsv="resources/gene_id_to_gene_name.tsv"
    log:
        "logs/geneid2genename.log"
    script:
        "../scripts/get_genename.py"
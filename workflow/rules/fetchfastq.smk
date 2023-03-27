rule sratools_fetchfastq:
	output:
		fq1="results/raw_fastq/{sample}/{sample}_1.fastq.gz",
		fq2="results/raw_fastq/{sample}/{sample}_2.fastq.gz",
		outdir="results/raw_fastq/{sample}"
	log:
	    "logs/raw_fastq/{sample}_fetchfastq.log"
	params:
		extra = config['params']['sratools_fetchfastq'],
	benchmark:
	    "benchmarks/{sample}.sratools_fetchfastq.benchmark.txt"
	conda:
		"../envs/sratools.yaml"
	shell:
		"fastq-dump "
		"{params.extra} "
		"--outdir {output.outdir} 2>{log}"



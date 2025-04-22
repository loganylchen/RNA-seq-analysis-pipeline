import sys
from Bio import SeqIO

def filter_fasta(input_fasta, output_fasta, contigs):
    with open(input_fasta) as infile, open(output_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id in contigs:
                SeqIO.write(record, outfile, "fasta")

def filter_gtf(input_gtf, output_gtf, contigs):
    with open(input_gtf) as infile, open(output_gtf, "w") as outfile:
        for line in infile:
            if line.startswith("#") or line.split("\t")[0] in contigs:
                outfile.write(line)


input_fasta = snakemake.input.fasta
input_gtf = snakemake.input.gtf
output_fasta = snakemake.output.fasta
output_gtf = snakemake.output.gtf
contigs = snakemake.params.contigs

filter_fasta(input_fasta, output_fasta, contigs)
filter_gtf(input_gtf, output_gtf, contigs)
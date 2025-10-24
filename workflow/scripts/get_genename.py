import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

import re


ATTRIBUTE=re.compile('(\w+?) "(.+?)";')

def _parse_attributes(attributes):
    att_dict = dict()
    for (key,value) in ATTRIBUTE.findall(attributes):
        att_dict[key] = value
    return att_dict['gene_id'],att_dict.get('gene_name',att_dict['gene_id']),att_dict.get('transcript_id',None)




with open(snakemake.input.gtf, 'r') as gtf, open(snakemake.output.tsv,'w') as f:
    f.write('gene_id\tgene_name\ttranscript_id\n')
    for line in gtf:
        if line.startswith("#"):
            continue
        else:
            elements = line.strip().split('\t')
            line_type,attributes = elements[2], elements[8]
            if line_type == 'transcript':
                gene_id, gene_name,transcript_id = _parse_attributes(attributes)
                f.write(f'{gene_id}\t{gene_name}\t{transcript_id}\n')








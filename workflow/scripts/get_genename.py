import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

import json




def _parse_attributes(attributes):
    atts = attributes.split(';')
    att_dict = dict()
    for att in atts:
        if att != '':
            key,value = att.strip().split(' ')
            value = value.replace('"','')
            att_dict[key] = value
    return att_dict['gene_id'],att_dict['gene_name']




with open(snakemake.input.gtf, 'r') as gtf, open(snakemake.output.tsv,'w') as f:
    f.write('gene_id\tgene_name\n')
    for line in gtf:
        if line.startswith("#"):
            continue
        else:
            elements = line.strip().split('\t')
            line_type,attributes = elements[2], elements[8]
            if line_type == 'gene':
                gene_id, gene_name = _parse_attributes(attributes)
                f.write(f'{gene_id}\t{gene_name}\n')








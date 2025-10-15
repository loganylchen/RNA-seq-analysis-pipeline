import sys
import os
# logging
sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import matplotlib.pyplot as plt

count_table = pd.read_csv(snakemake.input[0],index_col=0,comment='#',sep='\t')

colnames = count_table.columns

def _get_sample(bamfile):
    if bamfile.endswith('Aligned.sortedByCoord.out.bam'):
        return os.path.basename(os.path.dirname(bamfile))
    else:
        return bamfile


rename_colnames_dict = {coln:_get_sample(coln) for coln in colnames}

count_table.rename(columns=rename_colnames_dict,inplace=True)

count_table.to_csv(snakemake.output[0],sep='\t')



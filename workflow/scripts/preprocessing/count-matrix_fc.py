import sys

# logging
sys.stderr = open(snakemake.log[0], "w")
sys.stdout = sys.stderr
import pandas as pd




files = snakemake.input
_tmp_df = pd.read_csv(files[0],sep='\t',comment='#',index_col=0).iloc[:,5].to_frame()
for f in files[1:]:
    _tmp_df = _tmp_df.merge(pd.read_csv(f,sep='\t',comment='#',index_col=0).iloc[:,5].to_frame())
_tmp_df.columns = snakemake.params.samples
_tmp_df.to_csv(snakemake.output.count_matrix,sep='\t')
_tmp_df.T.to_csv(snakemake.output.puree_count_matrix,sep='\t')

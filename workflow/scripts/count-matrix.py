import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import matplotlib.pyplot as plt

counts_unstrandedness = [
    pd.read_table(
        f, index_col=0, usecols=[0, 1], header=None, skiprows=4
    )
    for f in snakemake.input
]

counts_strandedness = [
    pd.read_table(
        f, index_col=0, usecols=[0, 2], header=None, skiprows=4
    )
    for f in snakemake.input
]

counts_reverse = [
    pd.read_table(
        f, index_col=0, usecols=[0, 3], header=None, skiprows=4
    )
    for f in snakemake.input
]

for t, sample in zip(counts_unstrandedness, snakemake.params.samples):
    t.columns = [sample]

for t, sample in zip(counts_strandedness, snakemake.params.samples):
    t.columns = [sample]

for t, sample in zip(counts_reverse, snakemake.params.samples):
    t.columns = [sample]


unstrandedness_matrix = pd.concat(counts_unstrandedness, axis=1)
unstrandedness_matrix.index.name = "gene"

strandedness_matrix = pd.concat(counts_strandedness, axis=1)
strandedness_matrix.index.name = "gene"

reverse_matrix = pd.concat(counts_reverse, axis=1)
reverse_matrix.index.name = "gene"

unstrandedness_matrix.to_csv(snakemake.output[1], sep="\t")
strandedness_matrix.to_csv(snakemake.output[2], sep="\t")
reverse_matrix.to_csv(snakemake.output[3], sep="\t")

unstrandedness_sum = unstrandedness_matrix.sum(axis=0)
strandedness_sum = strandedness_matrix.sum(axis=0)
reverse_sum = reverse_matrix.sum(axis=0)


fig,ax = plt.subplots(1,1,figsize=(15,7))
ax.plot(unstrandedness_sum,label='unstrandedness')
ax.plot(strandedness_sum,label='strandedness')
ax.plot(reverse_sum,label='reverse')
ax.set_title('RNA-seq strandedness check')
ax.legend()
ax.set_xticklabels(ax.get_xticklabels(),rotation=75)
fig.savefig(snakemake.output[4],dpi=300)


u_s = unstrandedness_sum - strandedness_sum
s_u = -u_s
u_r = unstrandedness_sum - reverse_sum
r_u = -u_r
r_s = strandedness_sum - reverse_sum
s_r = -r_s




if ((u_s>0).count() > (s_u > 0).count()) and ((u_r>0).count() > (r_u > 0).count()):
    unstrandedness_matrix.to_csv(snakemake.output[0], sep="\t")
elif ((s_u>0).count() > (u_s > 0).count()) and ((s_r>0).count() > (r_s < 0).count()):
    strandedness_matrix.to_csv(snakemake.output[0], sep="\t")
elif ((r_u>0).count() > (u_r > 0).count()) and ((r_s>0).count() > (s_r < 0).count()):
    reverse_matrix.to_csv(snakemake.output[0], sep="\t")
else:
    raise ValueError("Can't decide the strandedness of the RNA-seq, please check by yourself")
# collapse technical replicates
# matrix = matrix.groupby(matrix.columns, axis=1, sort=False).sum()


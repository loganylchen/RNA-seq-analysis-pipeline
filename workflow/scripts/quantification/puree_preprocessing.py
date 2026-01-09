#!/usr/bin/env python3
"""
PUREE Tumor Purity Estimation
Estimates tumor purity (cancer cell fraction) from bulk gene expression data
"""

import sys
import os

import pandas as pd


# Logging
log_file = snakemake.log[0]
sys.stderr = open(log_file, "w")
sys.stdout = sys.stderr


df = pd.read_csv(
    snakemake.input.counts,
    sep="\t",
    index_col=0,
)

df.T.to_csv(
    snakemake.output.counts,
    sep="\t",
    header=True,
    index=True,
)

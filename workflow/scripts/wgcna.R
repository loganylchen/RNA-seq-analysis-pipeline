log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(WGCNA)
library(dplyr)


normalized_counts <- t(read.csv(snakemake@input[['normalized_matrix']],check.names=FALSE,row.names=1,header=TRUE,sep='\t'))

sft <- pickSoftThreshold(normalized_counts,
  dataIsExpr = TRUE,
  corFnc = cor,
  networkType = "signed"
)

sft_df <- data.frame(sft$fitIndices) %>%
  mutate(model_fit = -sign(slope) * SFT.R.sq)
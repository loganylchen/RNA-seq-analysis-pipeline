log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(geneset)
library(genekitr)
library(dplyr)

mf_set <- geneset::getGO(org = snakemake@params[["go_name"]], ont = "mf")
cc_set <- geneset::getGO(org = snakemake@params[["go_name"]], ont = "cc")
bp_set <- geneset::getGO(org = snakemake@params[["go_name"]], ont = "bp")
kegg_set <- getKEGG(org = snakemake@params[["kegg_name"]], category = "pathway")



df <- read.table(snakemake@input[[1]], header = T, row.names = 1,sep='\t')

geneList <- df$log2FoldChange
names(geneList) <- rownames(df)
geneList <- sort(geneList, decreasing = TRUE)


deg_df <- df %>%
    filter(abs(log2FoldChange) > snakemake@params[["log2foldchange_threshold"]]) %>%
    filter(padj < snakemake@params[["padj_threshold"]])

outdir <- snakemake@output[[1]]

up_enrich_mf <- genORA(rownames(deg_df %>% filter(log2FoldChange > 0)), geneset = mf_set)
write.table(up_enrich_mf, paste0(outdir,"/up_enrich_mf.tsv"), sep = "\t", quote = F, row.names = F)
down_enrich_mf <- genORA(rownames(deg_df %>% filter(log2FoldChange < 0)), geneset = mf_set)
write.table(up_enrich_mf, paste0(outdir,"/down_enrich_mf.tsv"), sep = "\t", quote = F, row.names = F)

up_enrich_cc <- genORA(rownames(deg_df %>% filter(log2FoldChange > 0)), geneset = cc_set)
write.table(up_enrich_cc, paste0(outdir,"/up_enrich_cc.tsv"), sep = "\t", quote = F, row.names = F)
down_enrich_cc <- genORA(rownames(deg_df %>% filter(log2FoldChange < 0)), geneset = cc_set)
write.table(down_enrich_cc, paste0(outdir,"/down_enrich_cc.tsv"), sep = "\t", quote = F, row.names = F)
up_enrich_bp <- genORA(rownames(deg_df %>% filter(log2FoldChange > 0)), geneset = bp_set)
write.table(up_enrich_bp, paste0(outdir,"/up_enrich_bp.tsv"), sep = "\t", quote = F, row.names = F)
down_enrich_bp <- genORA(rownames(deg_df %>% filter(log2FoldChange < 0)), geneset = bp_set)
write.table(down_enrich_bp, paste0(outdir,"/down_enrich_bp.tsv"), sep = "\t", quote = F, row.names = F)
gse <- genGSEA(genelist = geneList, geneset = kegg_set)
expoSheet(
    data_list = gse,
    data_name = names(gse),
    filename = "gsea_result.xlsx",
    dir = outdir
)















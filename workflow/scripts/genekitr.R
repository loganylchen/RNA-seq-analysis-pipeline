log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(geneset)
library(genekitr)
library(dplyr)
library(ggplot2)

source(file.path(snakemake@scriptdir, "defaultfunction.R"))

df <- read.table(snakemake@input[[1]], header = T, row.names = 1,sep='\t')

geneList <- df$log2FoldChange
names(geneList) <- rownames(df)
geneList <- sort(geneList, decreasing = TRUE)


deg_df <- df %>%
    filter(abs(log2FoldChange) > snakemake@params[["log2foldchange_threshold"]]) %>%
    filter(padj < snakemake@params[["padj_threshold"]])

outdir <- snakemake@output[[1]]

if (!dir.exists(outdir)) {dir.create(outdir)}


if(snakemake@params[["go_name"]] != "UNKNOWN"){

    mf_set <- geneset::getGO(org = snakemake@params[["go_name"]], ont = "mf")
    cc_set <- geneset::getGO(org = snakemake@params[["go_name"]], ont = "cc")
    bp_set <- geneset::getGO(org = snakemake@params[["go_name"]], ont = "bp")


    up_regulated_gene <- rownames(deg_df %>% filter(log2FoldChange > 0))
    down_regulated_gene <- rownames(deg_df %>% filter(log2FoldChange < 0))

    down_go<-FALSE
    up_go<-FALSE
    if(length(up_regulated_gene)>0){
        up_go<-TRUE
        up_enrich_mf <- genORA(up_regulated_gene, geneset = mf_set,q_cutoff=1,p_cutoff=1)
        if(dim(up_enrich_mf)[1]>0){
            write.table(up_enrich_mf, paste0(outdir,"/up_enrich_mf.tsv"), sep = "\t", quote = F, row.names = F)
        }
        up_enrich_cc <- genORA(up_regulated_gene, geneset = cc_set,q_cutoff=1,p_cutoff=1)
        if(dim(up_enrich_cc)[1]>0){
            write.table(up_enrich_cc, paste0(outdir,"/up_enrich_cc.tsv"), sep = "\t", quote = F, row.names = F)
        }
        up_enrich_bp <- genORA(up_regulated_gene, geneset = bp_set,q_cutoff=1,p_cutoff=1)
        if(dim(up_enrich_bp)[1]>0){
            write.table(up_enrich_bp, paste0(outdir,"/up_enrich_bp.tsv"), sep = "\t", quote = F, row.names = F)
        }
    }
    if(length(down_regulated_gene)>0){
        down_go<-TRUE
        down_enrich_mf <- genORA(down_regulated_gene, geneset = mf_set,q_cutoff=1,p_cutoff=1)
        if(dim(down_enrich_mf)[1]>0){
            write.table(down_enrich_mf, paste0(outdir,"/down_enrich_mf.tsv"), sep = "\t", quote = F, row.names = F)
        }
        down_enrich_cc <- genORA(down_regulated_gene, geneset = cc_set,q_cutoff=1,p_cutoff=1)
        if(dim(down_enrich_cc)[1]>0){
            write.table(down_enrich_cc, paste0(outdir,"/down_enrich_cc.tsv"), sep = "\t", quote = F, row.names = F)
        }
        down_enrich_bp <- genORA(down_regulated_gene, geneset = bp_set,q_cutoff=1,p_cutoff=1)
        if(dim(down_enrich_bp)[1]>0){
            write.table(down_enrich_bp, paste0(outdir,"/down_enrich_bp.tsv"), sep = "\t", quote = F, row.names = F)
        }
    }
    if(down_go & up_go){
        if((dim(up_enrich_mf)[1]>0)&(dim(down_enrich_mf)[1]>0)){
            go_mf_plot<-plotEnrichAdv(head(up_enrich_mf,10), head(down_enrich_mf,10),
                          plot_type = "one",
                          term_metric = "FoldEnrich",
                          stats_metric = "qvalue",
                          xlim_left = 25, xlim_right = 15) +
              theme(legend.position = c(0.15, 0.9))
            ggsave(paste0(outdir,"/enrich_mf.pdf"),go_mf_plot,width=10,height=4)
            ggsave(paste0(outdir,"/enrich_mf.png"),go_mf_plot,width=10,height=4)

        }
        if((dim(up_enrich_bp)[1]>0)&(dim(down_enrich_bp)[1]>0)){
            go_bp_plot<-plotEnrichAdv(head(up_enrich_bp,10), head(down_enrich_bp,10),
                          plot_type = "one",
                          term_metric = "FoldEnrich",
                          stats_metric = "qvalue",
                          xlim_left = 25, xlim_right = 15) +
              theme(legend.position = c(0.15, 0.9))
            ggsave(paste0(outdir,"/enrich_bp.pdf"),go_bp_plot,width=10,height=4)
            ggsave(paste0(outdir,"/enrich_bp.png"),go_bp_plot,width=10,height=4)
        }
        if((dim(up_enrich_cc)[1]>0)&(dim(down_enrich_cc)[1]>0)){
            go_cc_plot<-plotEnrichAdv(head(up_enrich_cc,10), head(down_enrich_cc,10),
                          plot_type = "one",
                          term_metric = "FoldEnrich",
                          stats_metric = "qvalue",
                          xlim_left = 25, xlim_right = 15) +
              theme(legend.position = c(0.15, 0.9))
            ggsave(paste0(outdir,"/enrich_cc.pdf"),go_cc_plot,width=10,height=4)
            ggsave(paste0(outdir,"/enrich_cc.png"),go_cc_plot,width=10,height=4)
        }





    }
}


if(snakemake@params[["kegg_name"]] != "UNKNOWN"){
    kegg_set <- getKEGG(org = snakemake@params[["kegg_name"]], category = "pathway")


    gse <- genGSEA(genelist = geneList, geneset = kegg_set,p_cutoff = 1, q_cutoff = 1)
    expoSheet(
        data_list = gse,
        data_name = names(gse),
        filename = "gsea_result.xlsx",
        dir = outdir
    )


    pathways <- rownames(gse$gsea_df %>% arrange(desc(abs(NES))))[1:3]
    gsea_plot <- plotGSEA(gse, plot_type = "classic", show_pathway = pathways)

    ggsave(paste0(outdir,"/enrich_gsea_kegg.png"),gsea_plot,width=10,height=10)
    ggsave(paste0(outdir,"/enrich_gsea_kegg.pdf"),gsea_plot,width=10,height=10)

    gsea_volcano_plot <- plotGSEA(gse, plot_type = "volcano", show_pathway = 5)

    ggsave(paste0(outdir,"/enrich_gsea_kegg_volcano.png"),gsea_volcano_plot,width=10,height=10)
    ggsave(paste0(outdir,"/enrich_gsea_kegg_volcano.pdf"),gsea_volcano_plot,width=10,height=10)

    gsea_bar_plot <- plotGSEA(gse, plot_type = "bar", colour = c("navyblue", "orange"),show_pathway = 5)

    ggsave(paste0(outdir,"/enrich_gsea_kegg_bar.png"),gsea_bar_plot,width=10,height=4)
    ggsave(paste0(outdir,"/enrich_gsea_kegg_bar.pdf"),gsea_bar_plot,width=10,height=4)
}
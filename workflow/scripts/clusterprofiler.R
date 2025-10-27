#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


suppressPackageStartupMessages({
    library(clusterProfiler)
    library(DOSE)
    library(dplyr)
    library(org.Mm.eg.db)
    library(org.Hs.eg.db)
})

species <- snakemake@params[['species']]

message(paste0('The species is ',species))
if(species == 'human'){
    kegg_org <- 'hsa'
    wp_org <- 'Homo sapiens'
    org.eg.db<-org.Hs.eg.db
}else if (species == 'mouse'){
    kegg_org <- 'mmu'
    wp_org <- 'Mus musculus'
    org.eg.db<-org.Mm.eg.db
}else{
    message('Only support human and mouse now')
    Sys.exit(1)
}


log2foldchange_threshold <- 1.5
padj_threshold <- 0.05




loading_data <- function(deg_tsv){

    message(paste0('Loading:',deg_tsv))
    DEG_list <- read.table(deg_tsv) %>%
            dplyr::filter(!is.na(baseMean)) %>% 
            dplyr::mutate(Ensembl_ID=rownames(.))
            
    ID_CONV <- bitr(DEG_list$Ensembl_ID, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"),OrgDb=org.eg.db)

    DEG_list <- DEG_list %>% 
                dplyr::full_join(ID_CONV, by=c('Ensembl_ID'='ENSEMBL'))
                

    sig_deg_list <- DEG_list %>% 
                    dplyr::filter(padj<padj_threshold,abs(log2FoldChange)>log2foldchange_threshold)

    up_regulated_deg_list <- sig_deg_list %>% dplyr::filter(log2FoldChange>0)
    down_regulated_deg_list <- sig_deg_list %>% dplyr::filter(log2FoldChange<0)
    return(list(
        deg_list=DEG_list,
        up_deg_list=up_regulated_deg_list,
        down_deg_list=down_regulated_deg_list
    ))
}



gsea_enrichment <- function(full_deg_list){
    message(paste0('full_deg_list: ', dim(full_deg_list)[1]))
    sorted_gene_list  <- full_deg_list %>% 
    dplyr::filter(!is.na(ENTREZID)) %>% 
    distinct(ENTREZID,.keep_all=TRUE) %>%
    arrange(-log2FoldChange)
    gene_list <- sorted_gene_list$log2FoldChange
    names(gene_list) <- sorted_gene_list$ENTREZID
    message(paste0('gene_list: ', length(gene_list)))
    message("GSEA on KEGG")
    gsea_kegg <- gseKEGG(geneList     =  gene_list,
                    organism     = kegg_org,
                    pvalueCutoff = 0.05,
                    verbose      = FALSE) %>% setReadable(.,org.eg.db)
    message("GSEA on WP")
    gsea_wp <- gseWP(geneList     =  gene_list,
                    organism     = wp_org,
                    pvalueCutoff = 0.05,
                    verbose      = FALSE) %>% setReadable(.,org.eg.db)
    message("GSEA on DO")
    gsea_do <- gseDO(gene_list,
           pAdjustMethod = "BH",
           verbose       = FALSE) %>% setReadable(.,org.eg.db)

    message("GSEA on NCG")
    gsea_ncg <-  gseNCG(gene_list,
              pAdjustMethod = "BH",
              verbose       = FALSE) %>% setReadable(.,org.eg.db)
    message("GSEA on DGN")
    gsea_dgn <- gseDGN(geneList,
              pAdjustMethod = "BH",
              verbose       = FALSE) %>% setReadable(.,org.eg.db)
    return(list(
        gsea=gsea_kegg,
        wp=gsea_wp,
        do=gsea_do,
        ncg=gsea_ncg,
        dgn=gsea_dgn,
    ))
}

ora_enrichment <- function(deg_list){
    message(paste0('The shape of the ORA deg_list:',dim(deg_list)[1]))
    if(dim(deg_list)[1]>=5){

   
    message("GO MF enrichment")
    go_mf <- enrichGO(gene= unique(deg_list$Ensembl_ID),
                    OrgDb         = org.eg.db,
                    ont='MF',
                    readable=TRUE,
                    keyType       = 'ENSEMBL',
                    pAdjustMethod = "BH",
                    qvalueCutoff  = 0.05) 
    message("GO CC enrichment")
    go_cc <- enrichGO(gene=  unique(deg_list$Ensembl_ID),
                    OrgDb         = org.eg.db,
                    ont='CC',
                    readable=TRUE,
                    keyType       = 'ENSEMBL',
                    pAdjustMethod = "BH",
                    qvalueCutoff  = 0.05) 

    message("GO BP enrichment")
    go_bp <- enrichGO(gene=  unique(deg_list$Ensembl_ID),
                    OrgDb         = org.eg.db,
                    ont='BP',
                    readable=TRUE,
                    keyType       = 'ENSEMBL',
                    pAdjustMethod = "BH",
                    qvalueCutoff  = 0.05) 

    message("KEGG enrichment")
    kegg_id <- enrichKEGG(gene=  unique(deg_list$ENTREZID),
                 organism     = kegg_org,pAdjustMethod = "BH",
                 pvalueCutoff = 0.05) 
    message('Readable on KEGG')
    # message(kegg_id%>% as.data.frame() %>% head())
    kegg_dim <- kegg_id%>% as.data.frame() %>% dim()
    if(kegg_dim[1]>0){
        kegg_res <- setReadable(kegg_id,org.eg.db,keyType='ENTREZID')
    }else{
        message('KEGG is NULL')
    }
    

    message("WIKIPATHWAY enrichment")
    wp_res<- enrichWP(gene=  unique(deg_list$ENTREZID), organism = wp_org,
                 pvalueCutoff = 0.05) 
    # message(wp_res%>% as.data.frame() %>% head())
    wp_dim <- wp_res%>% as.data.frame() %>% dim()
    if(wp_dim[1]>0){
        wp_res <- setReadable(wp_res,org.eg.db,keyType='ENTREZID')
    }else{
        message('WP is NULL')
    }


    message("DO enrichment")
    do_res <- enrichDO(gene  = unique(deg_list$ENTREZID),
              ont           = "HDO",
              pAdjustMethod = "BH",
              qvalueCutoff  = 0.05,
              readable      = FALSE)
    # message(do_res%>% as.data.frame() %>% head())
    do_dim <- do_res%>% as.data.frame() %>% dim()
    if(do_dim[1]>0){
        do_res <- setReadable(do_res,org.eg.db,keyType='ENTREZID')
    }else{
        message('DO is NULL')
    }
    message("NCG enrichment")
    ncg_res <- enrichNCG(gene  = unique(deg_list$ENTREZID),pAdjustMethod = "BH",
              readable      = FALSE) 
    # message(ncg_res %>% as.data.frame() %>% head())
    ncg_dim <- ncg_res%>% as.data.frame() %>% dim()
    if(ncg_dim[1]>0){
        ncg_res <- setReadable(ncg_res,org.eg.db,keyType='ENTREZID')
    }else{
        message('NCG is NULL')
    }
    
    message("DGN enrichment")
    dgn_res <- enrichDGN(gene  = unique(deg_list$ENTREZID),pAdjustMethod = "BH",
              readable      = FALSE) 
    # message(dgn_res%>% as.data.frame() %>% head())
    dgn_dim <- dgn_res%>% as.data.frame() %>% dim()
    if(dgn_dim[1]>0){
        dgn_res <- setReadable(dgn_res,org.eg.db,keyType='ENTREZID')
    }else{
        message('DGN is NULL')
    }
    return(list(
        go_cc=go_cc,
        go_bp=go_bp,
        go_mf=go_mf,
        kegg=kegg_res,
        wp=wp_res,
        do=do_res,
        ncg=ncg_res,
        dgn=dgn_res
    ))
     }else{
        return(list(
        go_cc=NULL,
        go_bp=NULL,
        go_mf=NULL,
        kegg=NULL,
        wp=NULL,
        do=NULL,
        ncg=NULL,
        dgn=NULL
    ))
     }
}




discovery_data_list <- loading_data(snakemake@input[['discovery_deg_tsv']])
validation_data_list <- loading_data(snakemake@input[['validation_deg_tsv']])
message('Discovery')
message('Up ORA')
up_ora <- ora_enrichment(discovery_data_list$up_deg_list)
message('Down ORA')
down_ora<-ora_enrichment(discovery_data_list$down_deg_list)
message('Running GSEA')
gsea_res <- gsea_enrichment(discovery_data_list$deg_list)
discovery=list(
        up_ora=up_ora,
        down_ora=down_ora,
        gsea=gsea_res)
message('Validation')

message('Up ORA')
up_ora <- ora_enrichment(validation_data_list$up_deg_list)
message('Down ORA')
down_ora<-ora_enrichment(validation_data_list$down_deg_list)

message('Running GSEA')
gsea_res <- gsea_enrichment(validation_data_list$deg_list)


validation=list(
        up_ora=up_ora,
        down_ora=down_ora,
        gsea=gsea_res)

final_res <- list(
    discovery=discovery,
    validation=validation
)







output<-snakemake@output[['enrichment']]
message("Wrting results into file")
saveRDS(final_res,output)

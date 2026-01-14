#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


cat("==============================================================\n")
cat("ClusterProfiler Enrichment Analysis\n")
cat("==============================================================\n\n")

cat("Session Info:\n")
print(sessionInfo())
cat("\n")

cat("Snakemake Parameters:\n")
cat("  params:", paste(names(snakemake@params), collapse=", "), "\n")
for (param_name in names(snakemake@params)) {
    cat("  ", param_name, ":", as.character(snakemake@params[[param_name]]), "\n")
}
cat("\n")

cat("Snakemake Input:\n")
cat("  input:", paste(names(snakemake@input), collapse=", "), "\n")
for (input_name in names(snakemake@input)) {
    cat("  ", input_name, ":", as.character(snakemake@input[[input_name]]), "\n")
}
cat("\n")

cat("Snakemake Output:\n")
cat("  output:", paste(names(snakemake@output), collapse=", "), "\n")
for (output_name in names(snakemake@output)) {
    cat("  ", output_name, ":", as.character(snakemake@output[[output_name]]), "\n")
}
cat("\n")

suppressPackageStartupMessages({
    cat("Loading required packages...\n")
    library(clusterProfiler)
    library(DOSE)
    library(dplyr)
    library(org.Mm.eg.db)
    library(org.Hs.eg.db)
    cat("All packages loaded successfully\n\n")
})

species <- snakemake@params[['species']]
cat("Species parameter:", species, "\n")

if(species == 'human'){
    kegg_org <- 'hsa'
    wp_org <- 'Homo sapiens'
    org.eg.db<-org.Hs.eg.db
    cat("Using human organism databases\n")
    cat("  KEGG organism:", kegg_org, "\n")
    cat("  WikiPathway organism:", wp_org, "\n")
    cat("  OrgDb:", class(org.eg.db), "\n\n")
}else if (species == 'mouse'){
    kegg_org <- 'mmu'
    wp_org <- 'Mus musculus'
    org.eg.db<-org.Mm.eg.db
    cat("Using mouse organism databases\n")
    cat("  KEGG organism:", kegg_org, "\n")
    cat("  WikiPathway organism:", wp_org, "\n")
    cat("  OrgDb:", class(org.eg.db), "\n\n")
}else{
    cat('ERROR: Only support human and mouse now, got:', species, "\n")
    Sys.exit(1)
}

padj_threshold <- as.numeric(snakemake@params[["padj_threshold"]])
deg_tool_n_threshold <- as.numeric(snakemake@params[["deg_tool_n_threshold"]])

cat("Analysis thresholds:\n")
cat("  padj_threshold:", padj_threshold, "\n")
cat("  deg_tool_n_threshold:", deg_tool_n_threshold, "\n\n")

loading_data <- function(deg_tsv, deg_tool_n_threshold){
    cat("\n==============================================================\n")
    cat("FUNCTION: loading_data\n")
    cat("==============================================================\n")

    cat("Input file:", deg_tsv, "\n")
    cat("deg_tool_n_threshold:", deg_tool_n_threshold, "\n")

    # Check if file exists
    if (!file.exists(deg_tsv)) {
        stop("ERROR: Input file does not exist: ", deg_tsv)
    }

    cat("Reading DEG TSV file...\n")
    DEG_df <- read.table(deg_tsv, header=TRUE, row.names=1, check.names=FALSE, sep='\t')
    cat("  Dimensions:", nrow(DEG_df), "genes x", ncol(DEG_df), "columns\n")
    cat("  Column names:", paste(colnames(DEG_df), collapse=", "), "\n")

    # Add Ensembl_ID column
    DEG_df <- DEG_df %>%
        dplyr::mutate(Ensembl_ID=rownames(.)) %>% as.data.frame()

    cat("  First 5 rows of DEG data:\n")
    print(head(DEG_df))

    cat("\nFiltering DEGs with up_regulated_count >= ", deg_tool_n_threshold,
        " OR down_regulated_count >= ", deg_tool_n_threshold, "\n")

    # Check if columns exist
    if (!"up_regulated_count" %in% colnames(DEG_df)) {
        stop("ERROR: 'up_regulated_count' column not found in input file")
    }
    if (!"down_regulated_count" %in% colnames(DEG_df)) {
        stop("ERROR: 'down_regulated_count' column not found in input file")
    }

    # Show distribution of counts
    cat("  up_regulated_count distribution:\n")
    print(table(DEG_df$up_regulated_count))
    cat("  down_regulated_count distribution:\n")
    print(table(DEG_df$down_regulated_count))

    # Filter DEGs (genes that are consistently up OR down regulated)
    DEG_list <- DEG_df %>%
        dplyr::filter(up_regulated_count >= deg_tool_n_threshold | down_regulated_count >= deg_tool_n_threshold)

    cat("  After filtering:", nrow(DEG_list), "genes remaining\n")

    if(dim(DEG_list)[1]==0){
        cat("WARNING: No DEGs passed the threshold. Relaxing to threshold=1\n")
        deg_tool_n_threshold <- 1
        DEG_list <- DEG_df %>%
            dplyr::filter(up_regulated_count >= deg_tool_n_threshold | down_regulated_count >= deg_tool_n_threshold)
        cat("  After relaxing threshold:", nrow(DEG_list), "genes\n")
    }

    if (nrow(DEG_list) == 0) {
        stop("ERROR: No DEGs found even with relaxed threshold")
    }

    cat("\nConverting Ensembl IDs to ENTREZID and SYMBOL...\n")
    cat("  Input Ensembl IDs:", length(unique(DEG_list$Ensembl_ID)), "\n")
    cat("  First 10 Ensembl IDs:", paste(head(DEG_list$Ensembl_ID, 10), collapse=", "), "\n")

    ID_CONV <- bitr(DEG_list$Ensembl_ID, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.eg.db)
    cat("  Conversion successful\n")
    cat("  Converted", nrow(ID_CONV), "IDs\n")
    cat("  First 5 conversions:\n")
    print(head(ID_CONV))

    DEG_list <- DEG_list %>%
        dplyr::full_join(ID_CONV, by=c('Ensembl_ID'='ENSEMBL'))

    cat("  After join:", nrow(DEG_list), "genes\n")
    cat("  Genes with ENTREZID:", sum(!is.na(DEG_list$ENTREZID)), "\n")
    cat("  Genes with SYMBOL:", sum(!is.na(DEG_list$SYMBOL)), "\n")

    # Split into up and down regulated
    cat("\nSplitting into up-regulated and down-regulated genes...\n")

    up_regulated_deg_list <- DEG_list %>% dplyr::filter(up_regulated_count >= deg_tool_n_threshold)
    down_regulated_deg_list <- DEG_list %>% dplyr::filter(down_regulated_count >= deg_tool_n_threshold)

    cat("  Up-regulated genes:", nrow(up_regulated_deg_list), "\n")
    cat("  Down-regulated genes:", nrow(down_regulated_deg_list), "\n")

    if (nrow(up_regulated_deg_list) > 0) {
        cat("  First 5 up-regulated genes:\n")
        print(head(up_regulated_deg_list[, c("Ensembl_ID", "SYMBOL", "up_regulated_count", "down_regulated_count")]))
    }

    if (nrow(down_regulated_deg_list) > 0) {
        cat("  First 5 down-regulated genes:\n")
        print(head(down_regulated_deg_list[, c("Ensembl_ID", "SYMBOL", "up_regulated_count", "down_regulated_count")]))
    }

    cat("\n==============================================================\n")
    cat("FUNCTION: loading_data - COMPLETE\n")
    cat("==============================================================\n\n")

    return(list(
        deg_list=DEG_list,
        up_deg_list=up_regulated_deg_list,
        down_deg_list=down_regulated_deg_list
    ))
}


ora_enrichment <- function(deg_list, padj_threshold){
    cat("\n==============================================================\n")
    cat("FUNCTION: ora_enrichment\n")
    cat("==============================================================\n")

    cat("Input deg_list dimensions:", nrow(deg_list), "x", ncol(deg_list), "\n")
    cat("padj_threshold:", padj_threshold, "\n")

    if(dim(deg_list)[1] < 5){
        cat("WARNING: Less than 5 genes, skipping enrichment\n")
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

    # Show genes being used for enrichment
    cat("Genes for enrichment analysis:\n")
    cat("  Total genes:", nrow(deg_list), "\n")
    cat("  Genes with Ensembl_ID:", sum(!is.na(deg_list$Ensembl_ID)), "\n")
    cat("  Genes with ENTREZID:", sum(!is.na(deg_list$ENTREZID)), "\n")

    # Filter to only genes with ENTREZID for databases that need it
    deg_list_entrez <- deg_list[!is.na(deg_list$ENTREZID), ]
    cat("  Genes with valid ENTREZID for enrichment:", nrow(deg_list_entrez), "\n")

    if (nrow(deg_list_entrez) < 3) {
        cat("WARNING: Less than 3 genes with ENTREZID, skipping enrichment\n")
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

    cat("\n--- Starting GO MF enrichment ---\n")
    go_mf <- tryCatch({
        result <- enrichGO(gene = unique(deg_list$Ensembl_ID),
                          OrgDb = org.eg.db,
                          ont = 'MF',
                          readable = TRUE,
                          keyType = 'ENSEMBL',
                          pAdjustMethod = "BH",
                          pvalueCutoff = padj_threshold,
                          qvalueCutoff = padj_threshold)
        if (is.null(result) || nrow(as.data.frame(result)) == 0) {
            cat("  No significant GO MF terms found\n")
        } else {
            cat("  Found", nrow(as.data.frame(result)), "significant GO MF terms\n")
        }
        result
    }, error = function(e) {
        cat("  ERROR in GO MF:", e$message, "\n")
        NULL
    })

    cat("\n--- Starting GO CC enrichment ---\n")
    go_cc <- tryCatch({
        result <- enrichGO(gene = unique(deg_list$Ensembl_ID),
                          OrgDb = org.eg.db,
                          ont = 'CC',
                          readable = TRUE,
                          keyType = 'ENSEMBL',
                          pAdjustMethod = "BH",
                          pvalueCutoff = padj_threshold,
                          qvalueCutoff = padj_threshold)
        if (is.null(result) || nrow(as.data.frame(result)) == 0) {
            cat("  No significant GO CC terms found\n")
        } else {
            cat("  Found", nrow(as.data.frame(result)), "significant GO CC terms\n")
        }
        result
    }, error = function(e) {
        cat("  ERROR in GO CC:", e$message, "\n")
        NULL
    })

    cat("\n--- Starting GO BP enrichment ---\n")
    go_bp <- tryCatch({
        result <- enrichGO(gene = unique(deg_list$Ensembl_ID),
                          OrgDb = org.eg.db,
                          ont = 'BP',
                          readable = TRUE,
                          keyType = 'ENSEMBL',
                          pAdjustMethod = "BH",
                          pvalueCutoff = padj_threshold,
                          qvalueCutoff = padj_threshold)
        if (is.null(result) || nrow(as.data.frame(result)) == 0) {
            cat("  No significant GO BP terms found\n")
        } else {
            cat("  Found", nrow(as.data.frame(result)), "significant GO BP terms\n")
        }
        result
    }, error = function(e) {
        cat("  ERROR in GO BP:", e$message, "\n")
        NULL
    })

    cat("\n--- Starting KEGG enrichment ---\n")
    cat("  Using", length(unique(deg_list_entrez$ENTREZID)), "ENTREZIDs\n")
    kegg_id <- tryCatch({
        result <- enrichKEGG(gene = unique(deg_list_entrez$ENTREZID),
                            organism = kegg_org,
                            pAdjustMethod = "BH",
                            pvalueCutoff = padj_threshold,
                            qvalueCutoff = padj_threshold)
        if (is.null(result) || nrow(as.data.frame(result)) == 0) {
            cat("  No significant KEGG pathways found\n")
            kegg_res <- NULL
        } else {
            cat("  Found", nrow(as.data.frame(result)), "significant KEGG pathways\n")
            kegg_res <- setReadable(result, org.eg.db, keyType='ENTREZID')
        }
        kegg_res
    }, error = function(e) {
        cat("  ERROR in KEGG:", e$message, "\n")
        NULL
    })

    cat("\n--- Starting WikiPathway enrichment ---\n")
    wp_res <- tryCatch({
        result <- enrichWP(gene = unique(deg_list_entrez$ENTREZID),
                          organism = wp_org,
                          pvalueCutoff = padj_threshold,
                          qvalueCutoff = padj_threshold)
        if (is.null(result) || nrow(as.data.frame(result)) == 0) {
            cat("  No significant WikiPathways found\n")
            wp_result <- NULL
        } else {
            cat("  Found", nrow(as.data.frame(result)), "significant WikiPathways\n")
            wp_result <- setReadable(result, org.eg.db, keyType='ENTREZID')
        }
        wp_result
    }, error = function(e) {
        cat("  ERROR in WikiPathway:", e$message, "\n")
        NULL
    })

    cat("\n--- Starting DO enrichment ---\n")
    do_res <- tryCatch({
        result <- enrichDO(gene = unique(deg_list_entrez$ENTREZID),
                          ont = "HDO",
                          pAdjustMethod = "BH",
                          pvalueCutoff = padj_threshold,
                          qvalueCutoff = padj_threshold,
                          readable = FALSE)
        if (is.null(result) || nrow(as.data.frame(result)) == 0) {
            cat("  No significant DO terms found\n")
            do_result <- NULL
        } else {
            cat("  Found", nrow(as.data.frame(result)), "significant DO terms\n")
            do_result <- setReadable(result, org.eg.db, keyType='ENTREZID')
        }
        do_result
    }, error = function(e) {
        cat("  ERROR in DO:", e$message, "\n")
        NULL
    })

    cat("\n--- Starting NCG enrichment ---\n")
    ncg_res <- tryCatch({
        result <- enrichNCG(gene = unique(deg_list_entrez$ENTREZID),
                           pAdjustMethod = "BH",
                           pvalueCutoff = padj_threshold,
                           qvalueCutoff = padj_threshold,
                           readable = FALSE)
        if (is.null(result) || nrow(as.data.frame(result)) == 0) {
            cat("  No significant NCG genes found\n")
            ncg_result <- NULL
        } else {
            cat("  Found", nrow(as.data.frame(result)), "significant NCG genes\n")
            ncg_result <- setReadable(result, org.eg.db, keyType='ENTREZID')
        }
        ncg_result
    }, error = function(e) {
        cat("  ERROR in NCG:", e$message, "\n")
        NULL
    })

    cat("\n--- Starting DGN enrichment ---\n")
    dgn_res <- tryCatch({
        result <- enrichDGN(gene = unique(deg_list_entrez$ENTREZID),
                           pAdjustMethod = "BH",
                           pvalueCutoff = padj_threshold,
                           qvalueCutoff = padj_threshold,
                           readable = FALSE)
        if (is.null(result) || nrow(as.data.frame(result)) == 0) {
            cat("  No significant DGN genes found\n")
            dgn_result <- NULL
        } else {
            cat("  Found", nrow(as.data.frame(result)), "significant DGN genes\n")
            dgn_result <- setReadable(result, org.eg.db, keyType='ENTREZID')
        }
        dgn_result
    }, error = function(e) {
        cat("  ERROR in DGN:", e$message, "\n")
        NULL
    })

    cat("\n==============================================================\n")
    cat("FUNCTION: ora_enrichment - COMPLETE\n")
    cat("==============================================================\n")

    return(list(
        go_cc=go_cc,
        go_bp=go_bp,
        go_mf=go_mf,
        kegg=kegg_id,
        wp=wp_res,
        do=do_res,
        ncg=ncg_res,
        dgn=dgn_res
    ))
}

cat("\n==============================================================\n")
cat("MAIN ANALYSIS STARTING\n")
cat("==============================================================\n\n")

cat("Loading combined DEG data...\n")
cat("Input file:", snakemake@input[['combined_deg_tsv']], "\n")
data_list <- loading_data(snakemake@input[['combined_deg_tsv']], deg_tool_n_threshold)

cat("\nData summary:\n")
cat("  Total DEGs:", nrow(data_list$deg_list), "\n")
cat("  Up-regulated DEGs:", nrow(data_list$up_deg_list), "\n")
cat("  Down-regulated DEGs:", nrow(data_list$down_deg_list), "\n")

message('\n--- Processing Up-regulated Genes ---')
up_ora <- ora_enrichment(data_list$up_deg_list, padj_threshold)

message('\n--- Processing Down-regulated Genes ---')
down_ora <- ora_enrichment(data_list$down_deg_list, padj_threshold)

discovery <- list(
    up_ora=up_ora,
    down_ora=down_ora
)

cat("\n==============================================================\n")
cat("Saving results\n")
cat("==============================================================\n")

output <- snakemake@output[['enrichment']]
cat("Output file:", output, "\n")

saveRDS(discovery, file=output)
cat("Results saved successfully\n")

cat("\n==============================================================\n")
cat("ClusterProfiler Analysis Complete!\n")
cat("==============================================================\n\n")

cat("Output summary:\n")
cat("  Up-regulated enrichment results:\n")
cat("    GO MF:", if(!is.null(discovery$up_ora$go_mf)) nrow(as.data.frame(discovery$up_ora$go_mf)) else 0, "terms\n")
cat("    GO CC:", if(!is.null(discovery$up_ora$go_cc)) nrow(as.data.frame(discovery$up_ora$go_cc)) else 0, "terms\n")
cat("    GO BP:", if(!is.null(discovery$up_ora$go_bp)) nrow(as.data.frame(discovery$up_ora$go_bp)) else 0, "terms\n")
cat("    KEGG:", if(!is.null(discovery$up_ora$kegg)) nrow(as.data.frame(discovery$up_ora$kegg)) else 0, "pathways\n")
cat("\n")
cat("  Down-regulated enrichment results:\n")
cat("    GO MF:", if(!is.null(discovery$down_ora$go_mf)) nrow(as.data.frame(discovery$down_ora$go_mf)) else 0, "terms\n")
cat("    GO CC:", if(!is.null(discovery$down_ora$go_cc)) nrow(as.data.frame(discovery$down_ora$go_cc)) else 0, "terms\n")
cat("    GO BP:", if(!is.null(discovery$down_ora$go_bp)) nrow(as.data.frame(discovery$down_ora$go_bp)) else 0, "terms\n")
cat("    KEGG:", if(!is.null(discovery$down_ora$kegg)) nrow(as.data.frame(discovery$down_ora$kegg)) else 0, "pathways\n")
cat("\n")

# Close logging
sink()
sink(type="message")
close(log)
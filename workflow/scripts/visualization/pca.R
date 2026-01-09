log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(dplyr)
library(DESeq2)
library(PCAtools)
library(cowplot)
library(ggplotify)
library(ggsci)



parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}


project<- snakemake@params[["project"]]
case_condition<-snakemake@params[["case_condition"]]
control_condition<-snakemake@params[["control_condition"]]
samples<-snakemake@params[["samples"]]
counts <- snakemake@input[["counts"]]
strandness <- snakemake@params[["strandness"]]
database <-snakemake@params[["database"]]

# output
output_pdf=snakemake@output[['pdf']]
output_png=snakemake@output[['png']]
output_clinical_info=snakemake@output[['clinical_info']]

# Extract output directory and prefix for individual subplots
output_dir <- dirname(output_png)
output_prefix <- tools::file_path_sans_ext(basename(output_png))

# Define individual subplot paths
scree_pdf <- file.path(output_dir, paste0(output_prefix, "_scree.pdf"))
scree_png <- file.path(output_dir, paste0(output_prefix, "_scree.png"))
pairs_pdf <- file.path(output_dir, paste0(output_prefix, "_pairs.pdf"))
pairs_png <- file.path(output_dir, paste0(output_prefix, "_pairs.png"))
biplot_pdf <- file.path(output_dir, paste0(output_prefix, "_biplot.pdf"))
biplot_png <- file.path(output_dir, paste0(output_prefix, "_biplot.png"))
loadings_pdf <- file.path(output_dir, paste0(output_prefix, "_loadings.pdf"))
loadings_png <- file.path(output_dir, paste0(output_prefix, "_loadings.png"))
eigencor_pdf <- file.path(output_dir, paste0(output_prefix, "_eigencor.pdf"))
eigencor_png <- file.path(output_dir, paste0(output_prefix, "_eigencor.png"))

strandness_df <- data.frame(
    sample_name = names(strandness),
    strandness = unlist(strandness)
)

coldata <- read.table(samples, header=TRUE, row.names="sample_name", check.names=FALSE,sep='\t')%>%
            dplyr::mutate(strandness= strandness_df$strandness[match(rownames(.),strandness_df$sample_name)]) %>%
            dplyr::mutate(plot_condition=paste0(strandness,':',condition))

tryCatch(
    {
      cts <- read.table(counts, header=TRUE, row.names="Geneid", check.names=FALSE,sep='\t')
    },
    error = function(e) {
        message("Error reading counts file: ", e$message)
        cts <- read.table(counts, header=TRUE, check.names=FALSE,sep='\t')
    },
    finally = {
        print(head(cts))
        message("Finished attempting to read counts file.")
    }
)


target_samples <- intersect(rownames(coldata), colnames(cts))
coldata <- coldata[target_samples, ]
cts <- cts[,target_samples]


dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=~condition)




identify_columns <- function(df) {
  n_rows <- nrow(df)
  selected_cols <- character()

  for (col_name in names(df)) {
    # Skip non-numeric or non-factor columns if needed
    unique_vals <- unique(df[[col_name]])
    n_unique <- length(unique_vals)

    # Check conditions
    if (n_unique > 1 && n_unique < n_rows) {
      selected_cols <- c(selected_cols, col_name)
    }
  }

  return(selected_cols)
}


get_pcatools_correlation <- function(pcaobj,
                                     components = NULL,
                                     metavars = NULL,
                                     corFUN = "pearson",
                                     corUSE = "pairwise.complete.obs",
                                     corMultipleTestCorrection = "BH") {

  # Use the code from PCAtools directly
  data <- pcaobj$rotated
  metadata <- pcaobj$metadata

  # If components not specified, use all
  if (is.null(components)) {
    components <- paste0("PC", 1:min(10, ncol(data)))  # PCAtools default
  }

  # If metavars not specified, use all
  if (is.null(metavars)) {
    metavars <- colnames(metadata)
  }

  # Convert data to matrix
  xvals <- data.matrix(data[, which(colnames(data) %in% components), drop = FALSE])
  yvals <- metadata[, which(colnames(metadata) %in% metavars), drop = FALSE]

  # Convert character columns to numeric (same as PCAtools)
  character_columns <- !unlist(lapply(yvals, is.numeric))
  character_columns <- names(which(character_columns))

  for (c in character_columns) {
    yvals[, c] <- as.numeric(as.factor(yvals[, c]))
  }

  yvals <- data.matrix(yvals)

  # Create correlation table
  corvals <- cor(xvals, yvals, use = corUSE, method = corFUN)

  # Calculate p-values
  N <- ncol(xvals) * ncol(yvals)
  pvals <- data.frame(
    pval = numeric(N),
    i = numeric(N),
    j = numeric(N)
  )

  k <- 0
  for (i in seq_len(ncol(xvals))) {
    for (j in seq_len(ncol(yvals))) {
      k <- k + 1
      pvals[k, 'pval'] <- cor.test(
        xvals[, i],
        yvals[, j],
        use = corUSE,
        method = corFUN
      )$p.value
      pvals[k, "i"] <- colnames(xvals)[i]
      pvals[k, "j"] <- colnames(yvals)[j]
    }
  }

  # Adjust for multiple testing
  if (corMultipleTestCorrection != "none") {
    pvals$pval <- p.adjust(pvals$pval, method = corMultipleTestCorrection)
  }

  # Reshape p-values to match correlation matrix
  pvals_wide <- reshape2::dcast(pvals, i ~ j, value.var = "pval")
  rownames(pvals_wide) <- pvals_wide$i
  pvals_wide$i <- NULL
  pvals_wide <- pvals_wide[match(rownames(corvals), rownames(pvals_wide)), ]
  pvals_wide <- pvals_wide[colnames(corvals)]

  # Convert to matrix
  pvals_matrix <- as.matrix(pvals_wide)

  return(list(
    correlation_matrix = corvals,
    p_value_matrix = pvals_matrix,
    xvals = xvals,
    yvals = yvals
  ))
}





draw_pca <- function(dds,coldata,output_pdf,output_png,output_clinical_info){
    message('DESeq')
    dds<- DESeq(dds)
    message('DESeq:vst')
    vst <- assay(vst(dds))
    message('DESeq:PCA')
    p <- pca(vst, metadata = colData(dds), removeVar = 0.5)
    message('DESeq:scree')
    pscree <- screeplot(p, components = getComponents(p, 1:30),
            hline = 80, axisLabSize = 14, titleLabSize = 20,
            returnPlot = FALSE)
    message('DESeq:parisplot')
    ppairs <- pairsplot(p, components = getComponents(p, c(1:5)),
            triangle = TRUE, trianglelabSize = 12,
            hline = 0, vline = 0,
            pointSize = 0.8, gridlines.major = FALSE, gridlines.minor = FALSE,
            colby = 'condition',  shape='strandness',
            title = '', plotaxes = FALSE,
            returnPlot = FALSE)
    message('DESeq:biplot')
    pbiplot <- biplot(p,
        # loadings parameters
            lab=rownames(coldata),
            showLoadings = FALSE,
            colby = 'condition', shape='strandness',
            hline = 0, vline = 0,
            gridlines.major = FALSE, gridlines.minor = FALSE,
            pointSize = 5,
            legendLabSize = 16, legendIconSize = 8.0,
            drawConnectors = TRUE,
            encircle = TRUE,
            encircleFill = TRUE,
            title = 'PCA bi-plot',
            subtitle = 'PC1 versus PC2',
            returnPlot = FALSE,legendPosition = 'top') + scale_color_bmj()
    message('DESeq:plotloadings')
    ploadings <- plotloadings(p, rangeRetain = 0.01, labSize = 4,
    title = 'Loadings plot', axisLabSize = 12,
    subtitle = 'PC1, PC2, PC3, PC4, PC5',
    caption = 'Top 1% variables',
    shape = 24, shapeSizeRange = c(4, 8),
    col = c('limegreen', 'black', 'red3'),
    legendPosition = 'top',
    drawConnectors = FALSE,
    returnPlot = FALSE)


    metavars<- identify_columns(as.data.frame(colData(dds)))
    message(metavars)
    message('DESeq:epigencorplot')
    write.table(as.data.frame(colData(dds))[ ,metavars],file=output_clinical_info,sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)
    cor_data <- get_pcatools_correlation(p,metavars=metavars,components = getComponents(p, 1:10),corFUN = 'pearson',corUSE = 'pairwise.complete.obs')
    print(cor_data)
    peigencor <- eigencorplot(p,
    components = getComponents(p, 1:10),
    metavars = metavars,
    cexCorval = 1.0,
    fontCorval = 2,
    posLab = 'all',
    rotLabX = 45,
    scale = TRUE,
    main = "PC clinical correlates",
    col = c('white', 'cornsilk1', 'gold', 'forestgreen', 'darkgreen'),
    cexMain = 1.5,
    plotRsquared = FALSE,
    corFUN = 'pearson',
    corUSE = 'pairwise.complete.obs',
    signifSymbols = c('****', '***', '**', '*', ''),
    signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
    returnPlot = FALSE)
    cor_matrix <- attributes(peigencor)$cor_matrix
    print(cor_matrix)


    top_row <- plot_grid(pscree, ppairs, pbiplot,
        ncol = 3,
        labels = c('A', 'B  Pairs plot', 'C'),
        label_fontfamily = 'serif',
        label_fontface = 'bold',
        label_size = 22,
        align = 'h',
        rel_widths = c(1.10, 0.80, 1.10))

    bottom_row <- plot_grid(ploadings,
        as.grob(peigencor),
        ncol = 2,
        labels = c('D', 'E'),
        label_fontfamily = 'serif',
        label_fontface = 'bold',
        label_size = 22,
        align = 'h',
        rel_widths = c(0.8, 1.2))

    fig<- plot_grid(top_row, bottom_row, ncol = 1,
        rel_heights = c(1.1, 0.9))

    # Save individual subplots
    message('Saving individual subplots...')
    message('  Saving scree plot...')
    ggsave(scree_pdf, pscree, width=10, height=8)
    ggsave(scree_png, pscree, width=10, height=8)
    message('    Saved: ', scree_pdf)

    message('  Saving pairs plot...')
    ggsave(pairs_pdf, ppairs, width=12, height=10)
    ggsave(pairs_png, ppairs, width=12, height=10)
    message('    Saved: ', pairs_pdf)

    message('  Saving biplot...')
    ggsave(biplot_pdf, pbiplot, width=10, height=8)
    ggsave(biplot_png, pbiplot, width=10, height=8)
    message('    Saved: ', biplot_pdf)

    message('  Saving loadings plot...')
    ggsave(loadings_pdf, ploadings, width=10, height=8)
    ggsave(loadings_png, ploadings, width=10, height=8)
    message('    Saved: ', loadings_pdf)

    message('  Saving eigencorplot...')
    ggsave(eigencor_pdf, as.grob(peigencor), width=12, height=10)
    ggsave(eigencor_png, as.grob(peigencor), width=12, height=10)
    message('    Saved: ', eigencor_pdf)

    # Save combined plot
    message('Saving combined PCA plot...')
    message('  Saving PDF: ', output_pdf)
    message('  Saving PNG: ', output_png)
    ggsave(output_pdf,fig,width=20,height=13)
    ggsave(output_png,fig,width=20,height=13)

    message('Output files:')
    message('  Combined PCA plot:')
    message('    PDF: ', output_pdf)
    message('    PNG: ', output_png)
    message('  Individual subplots:')
    message('    Scree plot: ', scree_pdf, ' / ', scree_png)
    message('    Pairs plot: ', pairs_pdf, ' / ', pairs_png)
    message('    Biplot: ', biplot_pdf, ' / ', biplot_png)
    message('    Loadings plot: ', loadings_pdf, ' / ', loadings_png)
    message('    Eigencorplot: ', eigencor_pdf, ' / ', eigencor_png)
}
message('ALL')
draw_pca(dds,coldata,output_pdf,output_png,output_clinical_info)
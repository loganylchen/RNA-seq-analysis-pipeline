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

strandness_df <- data.frame(
    sample_name = names(strandness),
    strandness = unlist(strandness)
)

coldata <- read.table(samples, header=TRUE, row.names="sample_name", check.names=FALSE,sep='\t')%>% 
            dplyr::mutate(strandness= strandness_df$strandness[match(rownames(.),strandness_df$sample_name)]) %>%
            dplyr::mutate(plot_condition=paste0(strandness,':',condition)) 

cts <- read.table(counts, header=TRUE, row.names="Geneid", check.names=FALSE,sep='\t')

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


draw_pca <- function(dds,coldata,output_pdf,output_png){
    message('DESeq')
    dds<- DESeq(dds)
    message('DESeq:vst')
    vst <- assay(vst(dds))
    message('DESeq:PCA')
    p <- pca(vst, metadata = colData(dds), removeVar = 0.1)
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
            colby = 'plot_condition', 
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

    ggsave(output_pdf,fig,width=20,height=13)
    ggsave(output_png,fig,width=20,height=13)
}
message('ALL')
draw_pca(dds,coldata,output_pdf,output_png)

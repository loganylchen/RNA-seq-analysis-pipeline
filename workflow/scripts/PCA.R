log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


suppressPackageStartupMessages({
    library(dplyr)
    library(DESeq2)
    library(PCAtools)
    library(cowplot)
    library(ggplotify)
    library(ggsci)
})


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
discovery_sample_type<-snakemake@params[["discovery_sample_type"]]
samples<-snakemake@params[["samples"]]
counts <- snakemake@input[["counts"]]


# output
output_pdf=snakemake@output[['pdf']]
output_png=snakemake@output[['png']]
output_discovery_pdf=snakemake@output[['discovery_pdf']]
output_discovery_png=snakemake@output[['discovery_png']]
output_validation_pdf=snakemake@output[['validation_pdf']]
output_validation_png=snakemake@output[['validation_png']]

coldata <- read.table(samples, header=TRUE, row.names="sample_name", check.names=FALSE,sep='\t')
cts <- read.table(counts, header=TRUE, row.names="Geneid", check.names=FALSE,sep='\t')
cts <- cts[,rownames(coldata)]
dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=~sample_type+condition)
coldata_discovery <- coldata %>% 
                    dplyr::filter(sample_type==discovery_sample_type)
coldata_validation <- coldata %>% 
                    dplyr::filter(sample_type != discovery_sample_type)
cts_discovery <- cts[,rownames(coldata_discovery)]
cts_validation <- cts[,rownames(coldata_validation)]
dds_discovery <- DESeqDataSetFromMatrix(countData=cts_discovery,
                              colData=coldata_discovery,
                              design=~condition)
dds_validation <- DESeqDataSetFromMatrix(countData=cts_validation,
                              colData=coldata_validation,
                              design=~condition)

draw_pca <- function(dds,coldata,output_pdf,output_png){
    dds<- DESeq(dds, parallel=parallel)
    vst <- assay(vst(dds))
    p <- pca(vst, metadata = colData(dds), removeVar = 0.1)
    pscree <- screeplot(p, components = getComponents(p, 1:30),
            hline = 80, axisLabSize = 14, titleLabSize = 20,
            returnPlot = FALSE) 

    ppairs <- pairsplot(p, components = getComponents(p, c(1:5)),
            triangle = TRUE, trianglelabSize = 12,
            hline = 0, vline = 0,
            pointSize = 0.8, gridlines.major = FALSE, gridlines.minor = FALSE,
            colby = 'sample_type',
            title = '', plotaxes = FALSE,
            returnPlot = FALSE)

    pbiplot <- biplot(p,
        # loadings parameters
            lab=rownames(coldata),
            showLoadings = TRUE,
            lengthLoadingsArrowsFactor = 1.5,
            sizeLoadingsNames = 4,
            colLoadingsNames = 'red4',
    # other parameters
            colby = 'sample_type', 
            hline = 0, vline = 0,
            gridlines.major = FALSE, gridlines.minor = FALSE,
            pointSize = 5,
            legendPosition = 'none', legendLabSize = 16, legendIconSize = 8.0,
            shape = 'condition', 
            drawConnectors = TRUE,
            title = 'PCA bi-plot',
            subtitle = 'PC1 versus PC2',
            returnPlot = FALSE) + scale_color_png()

    ploadings <- plotloadings(p, rangeRetain = 0.01, labSize = 4,
    title = 'Loadings plot', axisLabSize = 12,
    subtitle = 'PC1, PC2, PC3, PC4, PC5',
    caption = 'Top 1% variables',
    shape = 24, shapeSizeRange = c(4, 8),
    col = c('limegreen', 'black', 'red3'),
    legendPosition = 'none',
    drawConnectors = FALSE,
    returnPlot = FALSE)

    peigencor <- eigencorplot(p,
    components = getComponents(p, 1:10),
    metavars = c('condition','sample_type'),
    cexCorval = 1.0,
    fontCorval = 2,
    posLab = 'all', 
    rotLabX = 45,
    scale = TRUE,
    main = "PC clinical correlates",
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

    ggsave(output_pdf,fig,width=20,height=10)
    ggsave(output_png,fig,width=20,height=10)
}


draw_pca(dds,coldata,output_pdf,output_png)
draw_pca(dds_discovery,coldata_discovery,output_discovery_pdf,output_discovery_png)
draw_pca(dds_validation,coldata_validation,output_validation_pdf,output_validation_png)
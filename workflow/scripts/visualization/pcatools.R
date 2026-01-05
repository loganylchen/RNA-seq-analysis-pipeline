log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(DESeq2)
library(PCAtools)
library(ggplot2)
library(cowplot)
library(ggplotify)



parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}


# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene", check.names=FALSE,sep='\t')
cts <- cts[ , order(names(cts))]

coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample_name", check.names=FALSE,sep='\t')
coldata <- coldata[order(row.names(coldata)), , drop=F]

dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=as.formula(snakemake@params[["model"]]))

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > snakemake@params[["count_threshold"]], ]
# normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)

vst_data <- assay(vst(dds))
p <- pca(vst_data, metadata = coldata, removeVar = 0.1)
horn <- parallelPCA(vst_data)
elbow <- findElbowPoint(p$variance)
color_by <- snakemake@params[["color_by"]]
shape_by <- snakemake@params[["shape_by"]]

# Different plots

vline <- as.numeric(which(cumsum(p$variance) > 80)[1])




pscree <- screeplot(p,
    components = getComponents(p, 1:min(20,length(p$components))),
    vline = c(horn$n, elbow)) +
    geom_label(aes(x = horn$n + 1, y = 50,
      label = 'Horn\'s', vjust = -1, size = 8)) +
    geom_label(aes(x = elbow + 1, y = 60,
      label = 'Elbow method', vjust = -1, size = 8))

ppairs <- pairsplot(p, components = getComponents(p,  1:min(5,length(p$components))),
triangle = TRUE, trianglelabSize = 12,
hline = 0, vline = 0,
pointSize = 0.8, gridlines.major = FALSE, gridlines.minor = FALSE,
colby = color_by,
title = '', plotaxes = FALSE,
margingaps = unit(c(0.01, 0.01, 0.01, 0.01), 'cm'),
returnPlot = FALSE)

pbiplot <- biplot(p,
# loadings parameters
  showLoadings = TRUE,
  lengthLoadingsArrowsFactor = 1.5,
  sizeLoadingsNames = 4,
  colLoadingsNames = 'red4',
# other parameters
  lab = NULL,
  colby = color_by,
  hline = 0, vline = 0,
  gridlines.major = FALSE, gridlines.minor = FALSE,
  pointSize = 5,
  legendLabSize = 16, legendIconSize = 8.0,
  shape = shape_by,
  drawConnectors = FALSE,
  title = 'PCA',
  subtitle = 'PC1 versus PC2',
  returnPlot = FALSE)

ploadings <- plotloadings(p, rangeRetain = 0.01, labSize = 4,
title = 'Loadings plot', axisLabSize = 12,
subtitle = 'PC1, PC2, PC3, PC4, PC5',
caption = 'Top 1% variables',
shape = 24, shapeSizeRange = c(4, 8),
col = c('limegreen', 'black', 'red3'),
drawConnectors = FALSE,
returnPlot = FALSE)

peigencor <- eigencorplot(p,
components = getComponents(p, 1:max(horn$n, elbow)),
metavars = colnames(coldata),
col = c('white', 'cornsilk1', 'gold', 'forestgreen', 'darkgreen'),
cexCorval = 1.0,
fontCorval = 2,
posLab = 'all',
rotLabX = 45,
scale = FALSE,
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
  labels = c('A', 'B ', 'C'),
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

final_plot <- plot_grid(top_row, bottom_row, ncol = 1,
  rel_heights = c(1.1, 0.9))

png_out <- snakemake@output[['png']]
pdf_out <- snakemake@output[['pdf']]

ggsave(png_out,final_plot,width=20,height=15)
ggsave(pdf_out,final_plot,width=20,height=15)
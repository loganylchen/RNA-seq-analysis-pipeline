log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(EnhancedVolcano)
library(ggplot2)
library(dplyr)
source('defaultfunction.R')

df <- read.table(snakemake@input[["deg_exp"]], header = T, row.names = 1,sep='\t', check.names=FALSE)
geneid_to_genename <- read.table(snakemake@input[["geneid_to_genename"]], header = T, row.names = 1,sep='\t', check.names=FALSE)



g <- EnhancedVolcano(df,
        lab = geneid_to_genename[rownames(df),'gene_name'],
        title=snakemake@params[["contrast"]],
        x = "log2FoldChange",
        y = "padj",
        pCutoff = as.numeric(snakemake@params[["p_threshold"]]),
        FCcutoff = as.numeric(snakemake@params[["fc_threshold"]]),
        pointSize = 3.0,
        labSize = 6.0,
        legendLabels = c(
            "Not sig.", "Log2FC", "q-value",
            "q-value & Log2FC"
        ),
        gridlines.major = FALSE,
        gridlines.minor = FALSE,
        raster=TRUE,
        drawConnectors = TRUE,
        maxoverlapsConnectors = 20
    )

png_out = snakemake@output[['png']]
pdf_out = snakemake@output[['pdf']]

ggsave(png_out,g+annotation_custom(xmin=-Inf, ymin=-Inf, xmax=Inf, ymax=Inf, watermarkGrob()),width=10,height=10)
ggsave(pdf_out,g,width=10,height=10)



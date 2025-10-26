log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(EnhancedVolcano)
library(ggplot2)
library(dplyr)


discovery_deg <- readRDS(snakemake@input[['discovery_deg_rds']])
validation_deg <- readRDS(snakemake@input[['validation_deg_rds']])

geneid_to_genename <- read.table(snakemake@input[["geneid_to_genename"]], header = T, row.names = F,sep='\t', check.names=FALSE)


discovery_png <- snakemake@output[['discovery_png']]
discovery_pdf <- snakemake@output[['discovery_pdf']]
validation_png <- snakemake@output[['validation_png']]
validation_pdf <- snakemake@output[['validation_pdf']]

g1 <- EnhancedVolcano(discovery_deg,
        lab = geneid_to_genename[rownames(discovery_deg),'gene_name'],
        title="Discovery",
        x = "log2FoldChange",
        y = "padj",
        pCutoff = 0.05,
        FCcutoff = 1.5,
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

g2 <- EnhancedVolcano(validation_deg,
        lab = geneid_to_genename[rownames(validation_deg),'gene_name'],
        title="Discovery",
        x = "log2FoldChange",
        y = "padj",
        pCutoff = 0.05,
        FCcutoff = 1.5,
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



ggsave(discovery_png,g1,width=10,height=10)
ggsave(discovery_pdf,g1,width=10,height=10)
ggsave(validation_png,g2,width=10,height=10)
ggsave(validation_pdf,g2,width=10,height=10)


log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(WGCNA)
library(dplyr)
library(ggplot2)
library(limma)
library(tibble)
library(ComplexHeatmap)

allowWGCNAThreads(snakemake@threads)


source(file.path(snakemake@scriptdir, "defaultfunction.R"))

#### Function
make_module_heatmap <- function(module_name,
                                expression_mat,
                                meta_and_moduleeigengenes,
                                gene_module_key_df,
                                treatment
                                ) {



  # Create the ComplexHeatmap column annotation object
  col_annot <- HeatmapAnnotation(
    # Supply treatment labels
    Treatment=meta_and_moduleeigengenes[,treatment],
    module_eigengene = anno_barplot(dplyr::select(meta_and_moduleeigengenes, module_name)),
  )

  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    filter(module == module_name) %>%
    dplyr::pull(gene)

  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    t() %>%
    as.data.frame() %>%
    # Only keep genes from this module
    dplyr::filter(rownames(.) %in% module_genes) %>%
    # Order the samples to match col_annot_df
    dplyr::select(rownames(meta_and_moduleeigengenes)) %>%
    # Data needs to be a matrix
    as.matrix()

  # Normalize the gene expression values
  mod_mat <- mod_mat %>%
    # Scale can work on matrices, but it does it by column so we will need to
    # transpose first
    t() %>%
    scale() %>%
    # And now we need to transpose back
    t()

  # Create a color function based on standardized scale
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#67a9cf", "#f7f7f7", "#ef8a62")
  )

  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
    name = module_name,
    # Supply color function
    col = color_func,
    # Supply column annotation
    bottom_annotation = col_annot,
    # We don't want to cluster samples
    cluster_columns = FALSE,
    # We don't need to show sample or gene labels
    show_row_names = FALSE,
    show_column_names = FALSE
  )

  # Return heatmap
  return(heatmap)
}

#######################################################################################################################








outdir <- dirname(snakemake@output[["bwnet_rds"]])
if (!dir.exists(outdir)) {dir.create(outdir)}

fig_outdir <- snakemake@params[["fig_outdir"]]
if (!dir.exists(fig_outdir)) {dir.create(fig_outdir)}

normalized_counts <- t(read.csv(snakemake@input[['normalized_matrix']],check.names=FALSE,row.names=1,header=TRUE,sep='\t'))
metadata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample_name", check.names=FALSE,sep='\t')
phenotype <- snakemake@params[['phenotype']]

gsg <- goodSamplesGenes(normalized_counts, verbose = 3);
print('================================')
print(gsg$allOK)
print('================================')

if (!gsg$allOK) {
# Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0){
    printFlush(paste("Removing genes:", paste(names(normalized_counts)[!gsg$goodGenes], collapse = ", ")));
  }
  if (sum(!gsg$goodSamples)>0){
    printFlush(paste("Removing samples:", paste(rownames(normalized_counts)[!gsg$goodSamples], collapse = ", ")));
  }
# Remove the offending genes and samples from the data:
  normalized_counts <- normalized_counts[gsg$goodSamples, gsg$goodGenes]
}

sampleTree <- hclust(dist(normalized_counts), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.

pdf(file.path(fig_outdir,'sampleTree.pdf'), width = 12, height = 9);
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
dev.off()



sft <- pickSoftThreshold(normalized_counts,
  dataIsExpr = TRUE,
  corFnc = cor,
  networkType = "signed"
)

sft_df <- data.frame(sft$fitIndices) %>%
  mutate(model_fit = -sign(slope) * SFT.R.sq)


soft_threshold <- sft_df$Power[find_peaks(sft_df$model_fit,length(sft_df$model_fit))]

bwnet <- blockwiseModules(normalized_counts,
  maxBlockSize = 15000, # What size chunks (how many genes) the calculations should be run in
  TOMType = "signed", # topological overlap matrix
  power = soft_threshold, # soft threshold for network construction
  numericLabels = TRUE, # Let's use numbers instead of colors for module labels
  randomSeed = 1234, # there's some randomness associated with this calculation
  # so we should set a seed
)

saveRDS(bwnet,snakemake@output[["bwnet_rds"]])


module_eigengenes <- bwnet$MEs



# Convert labels to colors for plotting
mergedColors = labels2colors(bwnet$colors)
# Plot the dendrogram and the module colors underneath
pdf(file.path(fig_outdir,'ClusterDendrogram.pdf'), width = 12, height = 9);
plotDendroAndColors(bwnet$dendrograms[[1]], mergedColors[bwnet$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()






metadata <- metadata[rownames(module_eigengenes),]


des_mat <- model.matrix(~ metadata[,phenotype])

fit <- lmFit(t(module_eigengenes), design = des_mat)
fit <- eBayes(fit)

stats_df <- topTable(fit, number = ncol(module_eigengenes))

merge_df <- merge(metadata,module_eigengenes,by='row.names')
rownames(merge_df) <- merge_df[,'Row.names']
merge_df <- merge_df[,-1]




gene_module_key <- enframe(bwnet$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))


write.table(gene_module_key,snakemake@output[['gene_module_key']],sep='\t',quote=FALSE)



pdf(file.path(fig_outdir,'module_heatmap.pdf'))
for(moduel_name in unique(gene_module_key$module)){
  make_module_heatmap(module_name = module_name,expression_mat=normalized_counts,
                                meta_and_moduleeigengenes=merge_df,
                                gene_module_key_df=gene_module_key,
                                treatment=phenotype)
}
dev.off()



## vis

scale_threshold_plot <- ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  geom_vline(xintercept=soft_threshold, col='green') +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(sft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()

ggsave(file.path(fig_outdir,'scale_threshold.pdf'),width=10,height=10)
ggsave(file.path(fig_outdir,'scale_threshold.png'),width=10,height=10)


for(moduel_name in unique(gene_module_key$module)){
  g <- ggplot(merge_df, aes_string(x = phenotype ,y = module_name, color = phenotype )) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  theme_classic()
  ggsave(file.path(fig_outdir,paste0('ME_',moduel_name,'.pdf')),width=10,height=10)
  ggsave(file.path(fig_outdir,paste0('ME_',moduel_name,'.png')),width=10,height=10)
}


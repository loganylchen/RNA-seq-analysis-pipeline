#!/usr/bin/env Rscript
"""
Analyze and visualize rMATS alternative splicing results.
Creates summary statistics, volcano plots, and event-specific visualizations.
"""

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

# Parse Snakemake inputs
rmats_dir <- snakemake@input[["rmats_dir"]]
output_dir <- snakemake@output[["output_dir"]]
fdr_threshold <- snakemake@params[["fdr_threshold"]]
dpsi_threshold <- snakemake@params[["dpsi_threshold"]]

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Event types
event_types <- c("SE", "RI", "MXE", "A5SS", "A3SS")
event_names <- c(
  SE = "Skipped Exon",
  RI = "Retained Intron",
  MXE = "Mutually Exclusive Exons",
  A5SS = "Alternative 5' Splice Site",
  A3SS = "Alternative 3' Splice Site"
)

# Function to read rMATS file
read_rmats <- function(event_type, dir) {
  file_path <- file.path(dir, paste0(event_type, ".MATS.JCEC.txt"))
  if (!file.exists(file_path)) {
    return(NULL)
  }
  
  df <- read.table(file_path, header = TRUE, sep = "\t", quote = "", 
                   stringsAsFactors = FALSE, fill = TRUE)
  
  if (nrow(df) == 0) {
    return(NULL)
  }
  
  df$EventType <- event_type
  return(df)
}

# Read all event types
all_events <- lapply(event_types, read_rmats, dir = rmats_dir)
names(all_events) <- event_types
all_events <- all_events[!sapply(all_events, is.null)]

if (length(all_events) == 0) {
  stop("No rMATS results found!")
}

# Combine all events
combined_events <- bind_rows(all_events)

# Add significance flag
combined_events$Significant <- combined_events$FDR < fdr_threshold & 
  abs(combined_events$IncLevelDifference) > dpsi_threshold

# Summary statistics
summary_stats <- combined_events %>%
  group_by(EventType) %>%
  summarise(
    Total = n(),
    Significant = sum(Significant, na.rm = TRUE),
    Increased = sum(Significant & IncLevelDifference > 0, na.rm = TRUE),
    Decreased = sum(Significant & IncLevelDifference < 0, na.rm = TRUE),
    Mean_dPSI = mean(IncLevelDifference, na.rm = TRUE),
    Median_dPSI = median(IncLevelDifference, na.rm = TRUE)
  )

write.csv(summary_stats, file.path(output_dir, "summary_statistics.csv"), 
          row.names = FALSE)

cat("Summary Statistics:\n")
print(summary_stats)

# 1. Bar plot of event counts
p1 <- ggplot(summary_stats, aes(x = EventType, y = Significant)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Significant), vjust = -0.5) +
  labs(title = "Significant Alternative Splicing Events",
       subtitle = paste0("FDR < ", fdr_threshold, ", |dPSI| > ", dpsi_threshold),
       x = "Event Type", y = "Number of Significant Events") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(output_dir, "event_counts.pdf"), p1, width = 8, height = 6)

# 2. Stacked bar plot (increased vs decreased)
summary_long <- summary_stats %>%
  select(EventType, Increased, Decreased) %>%
  pivot_longer(cols = c(Increased, Decreased), 
               names_to = "Direction", values_to = "Count")

p2 <- ggplot(summary_long, aes(x = EventType, y = Count, fill = Direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Increased" = "#D55E00", "Decreased" = "#0072B2")) +
  labs(title = "Direction of Splicing Changes",
       x = "Event Type", y = "Number of Events",
       fill = "PSI Change") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "direction_stacked.pdf"), p2, width = 8, height = 6)

# 3. Volcano plots for each event type
volcano_plots <- lapply(names(all_events), function(event_type) {
  df <- all_events[[event_type]]
  
  # Add color column
  df$Color <- "Not Significant"
  df$Color[df$FDR < fdr_threshold & df$IncLevelDifference > dpsi_threshold] <- "Increased"
  df$Color[df$FDR < fdr_threshold & df$IncLevelDifference < -dpsi_threshold] <- "Decreased"
  
  # Create volcano plot
  p <- ggplot(df, aes(x = IncLevelDifference, y = -log10(FDR), color = Color)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Increased" = "#D55E00", 
                                   "Decreased" = "#0072B2",
                                   "Not Significant" = "grey70")) +
    geom_vline(xintercept = c(-dpsi_threshold, dpsi_threshold), 
               linetype = "dashed", color = "grey30") +
    geom_hline(yintercept = -log10(fdr_threshold), 
               linetype = "dashed", color = "grey30") +
    labs(title = event_names[event_type],
         x = "IncLevel Difference (dPSI)",
         y = "-log10(FDR)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "bottom")
  
  return(p)
})

# Combine volcano plots
combined_volcano <- wrap_plots(volcano_plots, ncol = 2)
ggsave(file.path(output_dir, "volcano_plots.pdf"), combined_volcano, 
       width = 12, height = 15)

# 4. Distribution of dPSI values
p4 <- ggplot(combined_events, aes(x = IncLevelDifference, fill = EventType)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Distribution of dPSI Values",
       x = "IncLevel Difference (dPSI)", y = "Density",
       fill = "Event Type") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "dpsi_distribution.pdf"), p4, width = 10, height = 6)

# 5. Box plot of dPSI by event type
p5 <- ggplot(combined_events, aes(x = EventType, y = IncLevelDifference, fill = EventType)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "dPSI Distribution by Event Type",
       x = "Event Type", y = "IncLevel Difference (dPSI)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")

ggsave(file.path(output_dir, "dpsi_boxplot.pdf"), p5, width = 8, height = 6)

# 6. Top significant events table
top_events <- combined_events %>%
  filter(Significant) %>%
  arrange(FDR, desc(abs(IncLevelDifference))) %>%
  select(GeneID, EventType, chr, strand, IncLevelDifference, FDR) %>%
  head(50)

write.csv(top_events, file.path(output_dir, "top_50_significant_events.csv"), 
          row.names = FALSE)

# 7. Gene-level summary (genes with multiple events)
gene_summary <- combined_events %>%
  filter(Significant) %>%
  group_by(GeneID) %>%
  summarise(
    NumEvents = n(),
    EventTypes = paste(unique(EventType), collapse = ","),
    MeanDPSI = mean(IncLevelDifference),
    MinFDR = min(FDR)
  ) %>%
  arrange(desc(NumEvents), MinFDR)

write.csv(gene_summary, file.path(output_dir, "gene_level_summary.csv"), 
          row.names = FALSE)

# 8. Heatmap of top genes with multiple events (if any)
top_genes <- gene_summary %>% filter(NumEvents >= 3) %>% head(20)

if (nrow(top_genes) > 0) {
  # Create matrix for heatmap
  heatmap_data <- combined_events %>%
    filter(GeneID %in% top_genes$GeneID, Significant) %>%
    select(GeneID, EventType, IncLevelDifference) %>%
    mutate(EventID = paste(EventType, 1:n(), sep = "_")) %>%
    pivot_wider(id_cols = GeneID, names_from = EventType, 
                values_from = IncLevelDifference, 
                values_fn = list(IncLevelDifference = mean))
  
  mat <- as.matrix(heatmap_data[, -1])
  rownames(mat) <- heatmap_data$GeneID
  mat[is.na(mat)] <- 0
  
  pdf(file.path(output_dir, "multi_event_genes_heatmap.pdf"), width = 8, height = 10)
  Heatmap(mat,
          name = "dPSI",
          col = colorRamp2(c(-0.5, 0, 0.5), c("#0072B2", "white", "#D55E00")),
          cluster_rows = TRUE,
          cluster_columns = FALSE,
          row_names_gp = gpar(fontsize = 8),
          column_title = "Genes with Multiple Significant Splicing Events",
          heatmap_legend_param = list(title = "dPSI"))
  dev.off()
}

# Save significant events for downstream analysis
significant_events <- combined_events %>% filter(Significant)
write.csv(significant_events, file.path(output_dir, "significant_events_all.csv"), 
          row.names = FALSE)

cat("\nAnalysis complete! Results saved to:", output_dir, "\n")
cat("Total significant events:", nrow(significant_events), "\n")

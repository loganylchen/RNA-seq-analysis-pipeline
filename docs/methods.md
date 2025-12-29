# Methods

## 2.1 RNA-Sequencing Data Processing

RNA-seq data processing and analysis were performed using a comprehensive Snakemake workflow (version ≥6.4.1) with containerized software tools to ensure reproducibility.

## 2.2 Quality Control and Preprocessing

Raw sequencing data in FASTQ format were subjected to quality control and preprocessing using **fastp** (version 0.24.0). This step included automatic adapter trimming, quality filtering with a Phred score threshold of Q20, polyG tail trimming (for Illumina NovaSeq data), and removal of low-quality reads. Per-sample quality metrics including read quality scores, GC content distribution, adapter contamination rates, and sequence duplication levels were collected. **MultiQC** was used to generate aggregated quality control reports across all samples.

## 2.3 Read Alignment

Quality-controlled reads were aligned to the human reference genome (**GRCh38**, Ensembl release 104) using **STAR aligner** (version 2.7.11b). STAR was run with the following parameters:

- `--twopassMode Basic` to improve alignment accuracy by performing a two-pass mapping for novel splice junction detection
- `--quantMode GeneCounts` to generate gene-level counts during alignment
- `--outSAMtype BAM SortedByCoordinate` to produce coordinate-sorted BAM files
- Strand-specific parameters were automatically applied based on library strandedness

As an alternative alignment strategy, **HISAT2** (version 2.2.1) was employed with transcriptome-aware indexing and the `--dta` flag for transcript assembly compatibility. Alignment quality metrics were assessed using **Qualimap** (version 2.3), including mapping rates, insert size distributions, and coverage uniformity. **Picard** (version 3.4.0) was used to collect additional RNA-seq metrics including ribosomal RNA content and 5'/3' bias.

## 2.4 Strandedness Assessment

Library strandedness was automatically assessed from STAR alignment output using SSP (Strand Specificity Prediction). Forward/reverse read assignment ratios were calculated from splice junction data, and libraries were classified as forward-stranded (FWD), reverse-stranded (REV), or unstranded (UNSTRAND) based on SSP scores > 0.8. The detected strandedness was applied to all downstream quantification steps.

## 2.5 Gene and Transcript Quantification

Gene-level read counts were generated using **featureCounts** from the **Subread** package (version 2.1.1). The aligned BAM files were processed with the following parameters:

- `-p --countReadPairs` for paired-end read counting
- `-t exon -g gene_id` for exon-level aggregation by gene identifier
- Strand-specific parameter (`-s 1` for forward-stranded, `-s 2` for reverse-stranded, `-s 0` for unstranded) automatically applied based on strandedness assessment

Transcript-level quantification was additionally performed using alignment-free methods. **Salmon** (version 1.10.3) performed quasi-mapping with GC bias correction and bootstrap confidence intervals using the library type parameter automatically detected (`-l ISF` for forward-stranded, `-l ISR` for reverse-stranded, `-l IU` for unstranded). **Kallisto** (version 0.51.1) provided fast pseudo-alignment for transcript abundance estimation with analogous strand-specific parameters. Transcripts per million (TPM) values were calculated for all samples to enable cross-sample comparison.

## 2.6 Transcript Assembly and Annotation

Reference-guided transcript assembly was performed using **StringTie** (version 2.2.3) with strand-aware transcript construction. Assembled transcripts from all samples were merged using StringTie merge operation with parameters `-c 0 -F 0 -T 0 -i` to combine novel isoforms and low-abundance transcripts. Transcript comparison against reference annotation (Ensembl release 104) was conducted using **gffcompare** (version 0.12.10), classifying transcripts as known (=), novel isoforms in known loci (j, c), or intergenic transcripts (u, x, i).

## 2.7 Differential Gene Expression Analysis

Differential gene expression (DEG) analysis was performed using **DESeq2** (version 1.46.0) in R. The analysis implemented a negative binomial generalized linear model with Wald tests for significance testing. Independent filtering based on mean normalized count and Benjamini-Hochberg multiple testing correction were applied to optimize statistical power. The pipeline supported a discovery-validation cohort design, with separate analyses performed on tissue samples (discovery cohort) and independent sample sets such as cell-free RNA (validation cohort). Genes with an adjusted p-value (FDR) < 0.00001 and absolute log2 fold change > 2 were considered significantly differentially expressed. Variance-stabilized transformation (VST) was applied to count data for visualization and clustering analyses.

## 2.8 Alternative Splicing Analysis

Differential alternative splicing analysis was conducted using **rMATS** (version 4.3.0, rMATS-turbo). Five splicing event types were analyzed: skipped exon (SE), retained intron (RI), mutually exclusive exons (MXE), alternative 5' splice site (A5SS), and alternative 3' splice site (A3SS). rMATS was executed with parameters `--variable-read-length --readLength 100 --novelSS --individual-counts` to enable variable read length support, novel splice site detection, and individual sample quantification. Splicing events with false discovery rate (FDR) < 0.05 and absolute inclusion level difference (dPSI) > 0.1 were considered statistically significant. Splicing results were further validated using **SpliceTools** with TPM-based filtering (control TPM threshold: 1, case TPM threshold: 1) for robust retained intron and skipped exon medley analysis.

## 2.9 RNA Modification Detection

Detection of RNA modifications, including 5-methylcytosine (m5C), was performed using **modtect**. The analysis utilized split BAM files with junction-aware alignments and incorporated statistical confidence scoring for modification site predictions.

## 2.10 Dark Channel Biomarker Analysis

A novel dark channel biomarker (DCB) analysis was implemented to identify cancer-specific biomarkers detectable in cell-free RNA. The DCB workflow comprised: (1) Discovery phase in tumor tissue to identify upregulated genes (TPM ≥ 1, |log2FC| ≥ 2 relative to normal tissue); (2) Validation phase in cell-free RNA to confirm tissue-specific detection. DCB genes were defined by: (i) low detection rate in normal cfRNA (< 10% of samples with TPM > 1), (ii) high detection rate in cancer cfRNA (≥ 30% of samples with TPM > 1), and (iii) significant upregulation in tumor tissue. DCB candidates were ranked by tissue fold-change and cfRNA detection statistics.

## 2.11 Functional Enrichment Analysis

Functional enrichment analysis of differentially expressed genes was conducted using **clusterProfiler** (version 4.14.0) in R. Over-representation analysis (ORA) and gene set enrichment analysis (GSEA) were performed across multiple ontologies and databases:

- **Gene Ontology (GO)**: Biological Process (BP), Cellular Component (CC), and Molecular Function (MF)
- **KEGG** (Kyoto Encyclopedia of Genes and Genomes) for metabolic and signaling pathways
- **Reactome** for curated biological pathways
- **WikiPathways** for community-maintained pathway knowledge
- **Disease Ontology (DO)** for disease associations
- **NCG** (Network of Cancer Genes) for cancer gene annotations

Enriched terms with adjusted p-value < 0.05 were considered statistically significant. Visualization of enrichment results was performed using dot plots, bar plots, enrichment maps, and cnetplots via clusterProfiler and custom R scripts.

## 2.12 Co-expression Network Analysis

Weighted gene co-expression network analysis (**WGCNA**) was performed to identify modules of co-expressed genes. The analysis included automatic power selection for scale-free topology, hierarchical clustering with dynamic tree cutting for module detection, and correlation of module eigengenes with clinical traits.

## 2.13 Data Visualization

Multiple visualization approaches were employed for result interpretation:

- **Principal Component Analysis**: **PCAtools** was used to assess sample clustering and identify batch effects with generation of scree plots, biplots, and variance decomposition
- **Volcano Plots**: **EnhancedVolcano** displayed differential expression results with gene annotation for top significant genes
- **Heatmaps**: **ComplexHeatmap** illustrated expression patterns of differentially expressed genes across samples with clinical annotation and Z-score scaling
- **Scatter Plots**: Comparison of log2 fold changes between discovery and validation cohorts with significance categorization
- **Boxplots**: Expression distributions of significant genes across conditions with faceted display by significance category

All visualizations were generated in publication-ready formats (PNG at 300 DPI, PDF) with customizable styling.

## 2.14 Software and Resource Availability

All analyses were performed using containerized software tools via Docker/Singularity for reproducibility. The pipeline was implemented in Snakemake and is publicly available at https://github.com/loganylchen/RNA-seq-analysis-pipeline. Key software versions included: STAR (2.7.11b), HISAT2 (2.2.1), fastp (0.24.0), Subread (2.1.1), Salmon (1.10.3), Kallisto (0.51.1), StringTie (2.2.3), rMATS (4.3.0), Qualimap (2.3), Picard (3.4.0), DESeq2 (1.46.0), clusterProfiler (4.14.0), and SpliceTools (custom build). Computational resources included multi-core processing (6-10 threads per job) and memory allocation ranging from 4GB to 100GB depending on analysis step.

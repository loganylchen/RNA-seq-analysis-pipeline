# RNA-seq Analysis Pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.4.1-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/loganylchen/RNA-seq-analysis-pipeline/workflows/Tests/badge.svg?branch=main)](https://github.com/loganylchen/RNA-seq-analysis-pipeline/actions?query=branch%3Amain+workflow%3ATests)

A comprehensive Snakemake workflow for RNA-seq data analysis, from raw reads to differential expression and functional enrichment analysis.

## Overview

This pipeline provides an end-to-end solution for RNA-seq data analysis with support for:
- Multiple aligners (STAR, HISAT2)
- Multiple quantification tools (featureCounts, Salmon, Kallisto)
- Comprehensive quality control
- Differential expression analysis (DESeq2)
- Alternative splicing analysis (rMATS, SpliceTools)
- RNA modification detection (HAMR, modtect)
- Functional enrichment analysis (clusterProfiler)
- Advanced visualization (PCA, volcano plots, WGCNA)

## Features

### Data Processing
- ✅ **Quality Control**: FastP for adapter trimming and quality filtering
- ✅ **Alignment**: STAR (2-pass mode) and HISAT2 with automatic strand detection
- ✅ **Quantification**: featureCounts, Salmon, and Kallisto with TPM normalization
- ✅ **Assembly**: StringTie for transcript assembly and quantification

### Analysis
- ✅ **Differential Expression**: DESeq2 with discovery/validation cohort support
- ✅ **Alternative Splicing**: rMATS for splicing event detection and SpliceTools for medley analysis
- ✅ **RNA Modifications**: HAMR and modtect for m5C and other modifications
- ✅ **Functional Enrichment**: GO, KEGG, WikiPathways, Disease Ontology, and more
- ✅ **Co-expression Networks**: WGCNA module detection and analysis

### Visualization
- ✅ **PCA Analysis**: Comprehensive PCA with PCAtools including scree plots, biplots, and clinical correlates
- ✅ **Volcano Plots**: Enhanced volcano plots with gene labels for discovery and validation sets
- ✅ **Quality Reports**: MultiQC reports integrating all QC metrics

## Requirements

### Software
- Snakemake ≥ 6.4.1
- Singularity/Docker (for containerized execution)
- Conda/Mamba (optional, for environment management)

### System Resources
- Minimum 32 GB RAM (recommended: 64+ GB for human genome)
- 100+ GB storage for reference files and results
- Multi-core processor (8+ cores recommended)

## Installation

### 1. Clone the repository
```bash
git clone https://github.com/loganylchen/RNA-seq-analysis-pipeline.git
cd RNA-seq-analysis-pipeline
```

### 2. Install Snakemake
```bash
# Using conda
conda create -n snakemake snakemake -c conda-forge -c bioconda
conda activate snakemake

# Or using mamba (faster)
mamba create -n snakemake snakemake -c conda-forge -c bioconda
mamba activate snakemake
```

### 3. Install Singularity (optional but recommended)
Follow instructions at: https://sylabs.io/guides/latest/user-guide/

## Configuration

### 1. Prepare Sample Sheet

Create or edit `config/samples.tsv` with your samples:

```tsv
sample_name	condition	patient	raw_data	fq1	fq2	project	batch	seq_type	sample_type	raw_project	raw_condition
Sample1	Tumor	P1	SRR123456	/path/to/R1.fq.gz	/path/to/R2.fq.gz	MyProject	batch1	pe	tissue	MyProject	Tumor
Sample2	Normal	P1	SRR123457	/path/to/R1.fq.gz	/path/to/R2.fq.gz	MyProject	batch1	pe	tissue	MyProject	Normal
```

**Required columns:**
- `sample_name`: Unique sample identifier
- `condition`: Experimental condition (e.g., Tumor, Normal, Treatment, Control)
- `project`: Project name (alphanumeric only)
- `seq_type`: Sequencing type (`pe` for paired-end, `se` for single-end)
- `sample_type`: Sample type for discovery/validation split

**Optional columns:**
- `raw_data`: SRA accession number (for automatic download)
- `fq1`, `fq2`: Paths to local FASTQ files
- `patient`: Patient/subject identifier
- `batch`: Batch identifier for batch effect analysis

### 2. Edit Configuration File

Edit `config/config.yaml`:

```yaml
# Sample information
samples: config/samples.tsv
project: MyProject
case_condition: Tumor
control_condition: Normal
discovery_sample_type: tissue

# Reference genome
reference:
  species: homo_sapiens  # or mus_musculus
  release: 104           # Ensembl release
  build: GRCh38          # or GRCh37, GRCm39, etc.
  select_contigs:        # Chromosomes to include
    - "1"
    - "2"
    # ... add all chromosomes
    - "X"
    - "Y"

# Differential expression thresholds
deg:
  padj: 0.05
  log2fc: 1.5

# Threading (adjust based on your system)
threads:
  star: 8
  salmon: 8
  kallisto: 8
  deseq2: 4
```

### 3. Configure Tool Parameters

Adjust tool-specific parameters in `config/config.yaml`:

```yaml
# STAR alignment
star:
  extra: " --twopassMode Basic "
  index_extra: ""

# FastP trimming
fastp:
  extra: " -c "

# featureCounts
featurecounts:
  extra: " -p --countReadPairs -t exon -g gene_id "

# rMATS splicing analysis
rmats:
  extra: " --variable-read-length --readLength 100 --novelSS --individual-counts "
```

## Usage

### Quick Start

```bash
# Dry run to check the workflow
snakemake --use-singularity -n

# Run with 8 cores
snakemake --use-singularity --cores 8

# Run on cluster (SLURM example)
snakemake --use-singularity --cluster "sbatch -p normal -n {threads} -t 24:00:00" --jobs 100
```

### Step-by-Step Execution

#### 1. Download Reference Genome (if needed)
```bash
snakemake --use-singularity --cores 8 resources/star_genome
```

#### 2. Quality Control and Alignment
```bash
snakemake --use-singularity --cores 8 {project}/qc/multiqc/
```

#### 3. Quantification
```bash
# STAR + featureCounts
snakemake --use-singularity --cores 8 {project}/quantification/STAR_FC/

# Salmon
snakemake --use-singularity --cores 8 {project}/quantification/salmon/

# Kallisto
snakemake --use-singularity --cores 8 {project}/quantification/kallisto/
```

#### 4. Differential Expression
```bash
snakemake --use-singularity --cores 8 {project}/diffexp/deseq2/
```

#### 5. Visualization
```bash
snakemake --use-singularity --cores 8 {project}/visualization/
```

#### 6. Run Complete Pipeline
```bash
snakemake --use-singularity --cores 16 all
```

## Output Structure

```
{project}/
├── data/                          # Raw data (downloaded or linked)
├── clean_data/                    # Trimmed reads
├── qc/
│   ├── fastp/                    # FastP QC reports
│   ├── STAR/                     # STAR alignment logs
│   ├── hisat2/                   # HISAT2 alignment logs
│   ├── qualimap-rnaseq/         # QualiMap reports
│   ├── kallisto/                 # Kallisto logs
│   └── multiqc/                  # Aggregated QC report
├── alignment/
│   ├── STAR/                     # STAR BAM files
│   └── hisat2/                   # HISAT2 BAM files
├── quantification/
│   ├── STAR_FC/                  # featureCounts results
│   │   ├── count_matrix.txt
│   │   └── TPM_matrix.txt
│   ├── salmon/                   # Salmon quantification
│   └── kallisto/                 # Kallisto quantification
├── diffexp/
│   ├── deseq2/
│   │   ├── discovery_deg.tsv    # Discovery cohort DEGs
│   │   ├── validation_deg.tsv   # Validation cohort DEGs
│   │   ├── discovery_deg.rds    # DESeq2 results object
│   │   └── validation_deg.rds
├── enrichment/
│   └── clusterprofiler/
│       ├── discovery_gsea_enrichment.tsv
│       ├── validation_gsea_enrichment.tsv
│       └── plots/
├── transcript_splicing/
│   ├── rmats/                    # rMATS splicing events
│   └── splicetools/              # SpliceTools analysis
├── modification/
│   └── modtect/                  # RNA modification calls
├── visualization/
│   ├── pca.png                   # PCA plots (STAR)
│   ├── salmon_pca.png           # PCA plots (Salmon)
│   ├── kallisto_pca.png         # PCA plots (Kallisto)
│   ├── Volcano_discovery.pdf    # Volcano plots
│   └── Volcano_validation.pdf
└── assembly/
    └── stringtie/                # Transcript assemblies
```

## Pipeline Components

### 1. Quality Control & Preprocessing
- **FastP**: Adapter trimming, quality filtering, and per-sample QC
- **MultiQC**: Aggregated QC report across all samples
- **Automatic strand detection**: Infers library strandedness from alignment

### 2. Alignment
- **STAR**: Two-pass mode for splice-aware alignment
- **HISAT2**: Fast alignment with splice-aware mapping
- **QualiMap**: RNA-seq specific QC metrics

### 3. Quantification
- **featureCounts**: Gene-level read counting from BAM files
- **Salmon**: Transcript-level quasi-mapping and quantification
- **Kallisto**: Fast transcript quantification
- **Automatic strandedness**: Detects and applies correct strandedness

### 4. Differential Expression
- **DESeq2**: Statistical analysis with discovery/validation cohort design
- **Normalization**: Automatic normalization and variance stabilization
- **Multiple testing correction**: Benjamini-Hochberg FDR correction

### 5. Alternative Splicing
- **rMATS**: Detects differential splicing events (SE, MXE, A5SS, A3SS, RI)
- **SpliceTools**: Retained intron and skipped exon medley analysis

### 6. Functional Analysis
- **Gene Ontology**: Biological Process, Molecular Function, Cellular Component
- **KEGG Pathways**: Metabolic and signaling pathways
- **WikiPathways**: Curated biological pathways
- **Disease Ontology**: Disease associations
- **GSEA**: Gene Set Enrichment Analysis

### 7. Visualization
- **PCAtools**: Comprehensive PCA with multiple plot types
- **EnhancedVolcano**: Publication-ready volcano plots
- **WGCNA**: Co-expression network modules (optional)

## Advanced Usage

### Custom Genome/Annotation

To use custom genome and annotation files:

```yaml
reference:
  custom: true
  genome_fasta: /path/to/genome.fa
  genome_gtf: /path/to/annotation.gtf
```

### Batch Effect Correction

Include batch information in `samples.tsv` and update the design formula:

```yaml
model: "~ batch + condition"
```

### Multiple Comparisons

Edit `config/config.yaml` to specify multiple condition comparisons:

```yaml
comparisons:
  - case: Treatment1
    control: Control
  - case: Treatment2
    control: Control
```

### Running Specific Modules

```bash
# Only QC and alignment
snakemake --use-singularity --cores 8 qc_all

# Only quantification
snakemake --use-singularity --cores 8 quantification_all

# Only differential expression
snakemake --use-singularity --cores 8 diffexp_all
```

## Troubleshooting

### Common Issues

1. **Out of memory errors**
   - Reduce the number of parallel jobs
   - Increase `--resources mem_mb=64000`
   - Use HISAT2 instead of STAR (lower memory)

2. **Alignment fails**
   - Check FASTQ file integrity
   - Verify reference genome matches your data
   - Check for sufficient disk space

3. **Empty quantification files**
   - Verify correct strandedness in config
   - Check GTF file format and compatibility
   - Ensure reads are mapping to genes (check MultiQC report)

4. **DESeq2 errors**
   - Ensure at least 2 replicates per condition
   - Check for outlier samples (remove from samples.tsv)
   - Verify count matrix has sufficient genes

### Getting Help

- Check the [Snakemake documentation](https://snakemake.readthedocs.io/)
- Review tool-specific documentation in the workflow
- Open an issue on GitHub

## Citation

If you use this pipeline in your research, please cite:

```
RNA-seq Analysis Pipeline
https://github.com/loganylchen/RNA-seq-analysis-pipeline
```

And cite the individual tools used in your analysis:
- **Snakemake**: Mölder et al., F1000Research 2021
- **STAR**: Dobin et al., Bioinformatics 2013
- **DESeq2**: Love et al., Genome Biology 2014
- **Salmon**: Patro et al., Nature Methods 2017
- **MultiQC**: Ewels et al., Bioinformatics 2016
- **clusterProfiler**: Yu et al., OMICS 2012
- (and others based on your specific analysis)

## License

This project is licensed under the terms specified in the LICENSE file.

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## Authors

- Logan Chen ([@loganylchen](https://github.com/loganylchen))

## Changelog

### Current Version (DEG branch)
- Comprehensive differential expression analysis
- Discovery/validation cohort design
- Multiple quantification methods
- RNA modification detection
- Alternative splicing analysis

### Planned Features
- [ ] WGCNA co-expression analysis (multiple modes)
- [ ] De novo annotation for non-model organisms
- [ ] Additional visualization options (heatmaps, etc.)
- [ ] Enhanced QC reporting
- [ ] Automated report generation
- [ ] Support for single-cell RNA-seq
- [ ] Integration with pathway databases
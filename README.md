# Snakemake workflow: `RNA-seq-analysis-pipeline`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/loganylchen/RNA-seq-analysis-pipeline/workflows/Tests/badge.svg?branch=main)](https://github.com/loganylchen/RNA-seq-analysis-pipeline/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for `General RNA-seq data analysis`


## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=loganylchen%2FRNA-seq-analysis-pipeline).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) <repo>sitory and its DOI (see above).

# TODO

* [ ] add WGCNA in the analysis, including such settings:
  * [ ] Using all the gene
  * [ ] Using proportion highly various gene
  * [ ] Using DEGs
* [ ] add ANNOTATION for un-annotated species.
  * [ ] Enrichment analysis.
* [ ] Visualization:
  * [ ] Heatmap for DEGs
  * [ ] PCA plots
* [ ] QC
* [ ] Reports
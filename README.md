# HSPC_omics

This repository contains the scripts, meta data files used in the manuscripts:  **Epigenetic fate determinants of human hematopoietic stem and progenitor cells revealed by single-cell multi-omics analysis**.

## Dependencies

- python(2/3)
- R(v4.0.3)


## Scripts

Scripts and snakemake files for the reproduce figures and analysis.

#### Transcriptome_cluster.R

Single-cell RNA-seq data dimensionality reduction and clustering.

#### DMR_identification.R

Differential methylated region identify.

#### 5mc_DR.R

DNA methylation data dimensionality reduction and unsupervised clustering.

#### 5mc_enrichment.R

DNA methylation data genomic enrichment and motif enrichment analysis of specific regions(eg. DMRs). 

#### 5mc_Monocle.R

DNA methylation data pseudotime analysis.

#### mtDNA

Optimized mtDNA mutation calling strategy by combining both DNA and RNA sequences in single-cell multi-omics datasets.

#### eDMR

Enhancer-DMR identification by using REPTILE.

#### Multi-omics_integration.R

Multi-omics integration analysis.

#### Plot.R

Ploting codes for figures reproduce.

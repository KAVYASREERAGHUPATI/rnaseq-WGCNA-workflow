# RNA-seq WGCNA Workflow
A complete R workflow for Weighted Gene Co-expression Network Analysis (WGCNA) using normalized RNA-seq expression data, including quality control, co-expression module detection, module–trait association, hub gene identification, TOM-based Cytoscape network export and publication-ready tables and figures.


A comprehensive and reproducible **Weighted Gene Co-expression Network Analysis (WGCNA)** workflow implemented in **R** for normalized RNA-seq expression data.
This pipeline performs an end-to-end WGCNA analysis starting from a normalized gene expression matrix and a sample trait file. It automatically generates publication-ready figures, comprehensive result tables, hub gene statistics, TOM-based Cytoscape network files and reusable R objects.

## Input Files

Only **two input files** are required.

### 1. Normalized Expression Matrix

- Genes in rows
- Samples in columns
- First column contains Gene IDs

Example

| GeneID | Sample_1 | Sample_2 | Sample_3 |
|---------|----------|----------|----------|
| Gene_1 | 8.42 | 8.13 | 7.84 |
| Gene_2 | 5.67 | 6.04 | 6.31 |
**Note:** Raw read counts should not be used directly for WGCNA.

### 2. Sample Trait File

| sample | condition |
|---------|-----------|
| Sample_1 | Control |
| Sample_2 | Control |
| Sample_3 | Treatment |

Sample names must exactly match the column names in the expression matrix.
## Workflow


     Normalized Expression Matrix
              │
              ▼
     Expression Filtering
              │
              ▼
     Sample Quality Control
              │
              ▼
     Sample Clustering
              │
              ▼
     Soft-threshold Selection
              │
              ▼
     Co-expression Network Construction
              │
              ▼
     Module Detection
              │
              ▼
     Module Eigengene Analysis
              │
              ▼
     Module–Trait Correlation
              │
              ▼
     Gene Significance & Module Membership
              │
              ▼
    Hub Gene Identification
              │
              ▼
    TOM-based Cytoscape Network Export
              │
              ▼
    Publication-ready Results

## Output Files

### Quality Control

- Filtered normalized expression matrix
- Matched sample and trait table
- Sample clustering plot

### Network Construction

- Soft-threshold plots
- Gene dendrogram with module colors
- Module eigengene clustering
- Module eigengene correlation heatmap

### Module Analysis

- Module size table
- Gene-module assignment table
- Module–trait correlation table
- Module–trait P-value table

### Gene Statistics

- Gene significance (GS)
- Module membership (MM)
- Module membership P-values
- Intramodular connectivity (kWithin, kTotal, kOut and kDiff)

### Hub Genes

For every module:

- Complete gene list
- Top 50 hub genes

### Network Matrices

- Gene correlation matrix
- Adjacency matrix
- Topological Overlap Matrix (TOM)

### Cytoscape Files

- Edge table
- Node table

### Figures

- Sample clustering
- Soft-threshold plots
- Gene dendrogram with module colors
- Module eigengene clustering
- Module eigengene correlation heatmap
- Module–trait relationship heatmap
- GS vs MM scatter plots
- TOM heatmap

### Additional Outputs

- Analysis summary
- Session information
- Reusable R objects

## Requirements

- R (version 4.2 or later recommended)

Required R packages:

- WGCNA
- dynamicTreeCut
- dplyr

The script automatically installs any missing packages before running the analysis.

## Repository Structure

rnaseq-wgcna-workflow/
│
├── scripts
       |-- Complete_General_WGCNA_Pipeline.R
├── README.md
└── LICENSE

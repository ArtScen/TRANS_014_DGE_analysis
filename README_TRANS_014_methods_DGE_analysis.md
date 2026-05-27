# TRANS_014 : RNA-seq DGE Analysis | Mus musculus BMDMs Spic KO × Heme RNAseq

**Author:** Hugo Lainé  
**Affiliation:** Gulbenkian Institute for Molecular Medicine, Oeiras, Portugal  
**Project:** *Mus musculus* BMDMs | DGE - Spic KO × Heme treatment (2×2 factorial design)  
**Contact:** hugo.laine@gimm.pt

> **For the self-contained publishable methods script** that reproduces all pairwise comparisons
> DGE analyses results and figures described here, see `TRANS_014_methods_DGE_analysis.R`.  
> Edit **SECTION 1 only** to adapt it to your project paths and identifiers.

---

## Overview

This repository contains the code for the differential gene expression (DGE)
analysis of RNA-seq data from bone marrow-derived macrophages (BMDMs) with a
2×2 factorial design:

| Factor | Levels |
|---|---|
| **Genotype** | Wild Type (WT) · Spic knockout (KO) |
| **Treatment** | Control (Ctrl) · Heme |

**n = 3** biological replicates per group · **12 samples total**  
**Reference genome:** *Mus musculus* GRCm39 (Ensembl release 111)

---

## Script

`TRANS_014_methods_DGE_analysis.R` is the single self-contained script covering the
complete workflow:

| Section | Role |
|---------|------|
| SECTION 0 | Dependencies - R packages loading |
| SECTION 1 | **User configuration : ⚠️⚠️⚠️ only section requiring edits ⚠️⚠️⚠️** |
| SECTION 2 | DGE Analysis parameters : thresholds for mean average counts, log2FC and padj  |
| SECTION 3 | GTF annotation utility : mapping Ensembl gene IDs (i.e. `ENSEMBL_ID`) to HGNC/MGI gene symbols (i.e. `gene_name`) by parsing the GTF file |
| SECTION 4 | Load and prepare count data |
| SECTION 5 | DESeq2 dataset construction |
| SECTION 6 | Principal Component Analysis |
| SECTION 7 | Heatmaps |
    | SECTION 7 - Heatmap A | Top 500 most variable genes |
    | SECTION 7 - Heatmap B | Gene panel I: iron metabolism & oxidative stress markers |
    | SECTION 7 - Heatmap C | Gene panel II: Nrf2/heme response pathway |
| SECTION 8 | Volcano plot function |
| SECTION 9 | Differential Gene Expression: 4 pairwise comparisons |
| SECTION 10 | Session information: R environment used for the workflow : OS, R version, packages attached & versions |

---

## Dependencies

**R version:** 4.4.x or later recommended (script developed and tested under R 4.4.2)  
**Bioconductor version:** 3.20 (automatically matched to your R version by BiocManager)  

Install required R packages before running:

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2"))

install.packages(c("ggplot2", "ggrepel", "dplyr", "pheatmap",
                   "data.table", "stringr", "matrixStats"))
```

---

## Reproducible environment (renv)

To ensure full reproducibility, this repository includes an `renv.lock` file that records
the exact version of every R package used in the original analysis.

**To restore the exact environment on your machine (run once, from the R console):**

```r
# Install renv if not already installed
install.packages("renv")

# Restore the locked environment (downloads and installs the exact package versions)
renv::restore()
```

Once `renv::restore()` completes, run `TRANS_014_methods_DGE_analysis.R` as usual.
No need to run `renv::init()` or `renv::snapshot()`, those were used by the
author to generate the `renv.lock` file and are not part of the analysis workflow.

**If you encounter a Bioconductor version conflict error** during `renv::restore()`, 
it means your machine has a stale `BiocVersion` package incompatible with your R version.
Fix it with:

```r
remove.packages("BiocVersion")
install.packages("BiocManager")
BiocManager::install(version = "3.20")   # use the version matching R 4.4.x
renv::restore()
```

> **Note:** `renv` manages package versions but not the R version itself.
> This workflow was developed under **R 4.4.2 / Bioconductor 3.20**.
> Using a substantially different R version may require re-running `renv::init()`

---

## Usage

1. Open `TRANS_014_methods_DGE_analysis.R`
2. Edit the three path variables at the top of **Section 1**:

```r
COUNTS_FILE <- "path/to/featurecounts_output.fc"        # Raw reads per gene counts table from 
GTF_FILE    <- "path/to/Mus_musculus.GRCm39.111.gtf"    # Annotation file used with the reference genome to align the fastq file and count the features (here gene-level raw counts) 
                                                        # Used here to Map Ensembl gene IDs (i.e. `ENSEMBL_ID`) to HGNC/MGI gene symbols (i.e. `gene_name`) by parsing the GTF file
OUTPUT_DIR  <- "path/to/output/directory"               # Directory for `TRANS_014_methods_DGE_analysis.R` results output 
```

3. Run the script in full (`Rscript TRANS_014_methods_DGE_analysis.R` or source in RStudio)

---

## Input data

| File | Description |
|---|---|
| `*.fc` | featureCounts output (gene-level raw counts, 57 180 genes) |
| `Mus_musculus.GRCm39.111.gtf` | Ensembl 111 annotation (gene ID → symbol mapping) |

> **Note on sample labels:** Raw featureCounts columns were delivered with
> incorrect labels. The script renames columns to match the true experimental
> design, as documented in Section 4 of `TRANS_014_methods_DGE_analysis.R`.
 
> **Raw data availability:** RNA-seq raw data (FASTQ files) and the featureCounts
> count matrix will be deposited in NCBI GEO upon publication (accession: GSExxxxxxx - 
> placeholder; update before making this repository public).

---

## Outputs

For each of the four pairwise comparisons:

| Output | Description |
|---|---|
| `DGE_<comparison>.csv` | Full DESeq2 result table (all tested genes) |
| `DEGs_<comparison>.csv` | Filtered list: padj < 0.05, \|log2FC\| > 0.585 |
| `Volcano_<comparison>.pdf` | Volcano plot labelled with gene symbols |
| `PCA.pdf` | Principal Component Analysis plot : unsupervised ~ no prior on groups |
| `Heatmap_top_500_variable_genes.pdf` | Top 500 variable genes heatmap |
| `Heatmap_GOI_gene_panel_I.pdf` | Gene panel I  heatmap: iron metabolism & oxidative stress markers |
| `Heatmap_GOI_gene_panel_II.pdf` | Gene panel II  heatmap: Nrf2/heme response pathway |


> **Suggestion about 1 item deliberately not included:**
> This workflow does not directly extract and export of the raw interaction term (`genotypeKO.treatmentHeme` coefficient from the *linear model design*) and respective results. 
> This a scientific suggestion,  not a lack or a publication issue, whether to include it depends on whether the paper makes a claim about genotype × treatment synergy. 
> If it does, add `results(dds, name = "genotypeKO.treatmentHeme")` as a fifth output in Section 9. If the paper does not directly test that question, omit it.

---

## Statistical methods

- **Framework:** DESeq2 (Love, Anders & Huber, 2014)
- **Model:** `~ genotype + treatment + genotype:treatment`
- **Normalisation:** Median-of-ratios (DESeq2 size factors)
- **Visualisation:** Variance-stabilising transformation (`vst`):
   - `blind = TRUE` for PCA unsupervised | Dispersion estimates are not informed by group structure
   - `blind = FALSE` for heatmaps | Dispersion estimates use the fitted model  --> Appropriate when comparing known experimental groups
- **Significance thresholds:** padj < 0.05 · |log2FC| > 0.585 (~1.5-fold)
- **Multiple testing correction:** Benjamini-Hochberg FDR

---

## Citation

If you use this code or the analysis approach, please cite:

**DESeq2** (primary statistical framework):
> Love MI, Huber W, Anders S (2014). Moderated estimation of fold change
> and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15:550.
> https://doi.org/10.1186/s13059-014-0550-8

**featureCounts** (input data generation):
> Liao Y, Smyth GK, Shi W (2014). featureCounts: an efficient general purpose
> program for assigning sequence reads to genomic features. *Bioinformatics*, 30(7):923–930.
> https://doi.org/10.1093/bioinformatics/btt656

**ggplot2** (visualisation):
> Wickham H (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag New York.
> https://ggplot2.tidyverse.org

**pheatmap** (heatmap visualisation):
> Kolde R (2019). pheatmap: Pretty Heatmaps. R package version 1.0.12.
> https://CRAN.R-project.org/package=pheatmap
```

---

## Acknowledgement

`TRANS_014_methods_DGE_analysis.R` and this README were written by **Hugo Lainé**
(Gulbenkian Institute for Molecular Medicine, 2025). If you use or adapt this
workflow, please cite the associated publication and credit the original author.
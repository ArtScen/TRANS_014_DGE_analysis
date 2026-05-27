# =============================================================================
# TRANS_014_methods_DGE_analysis.R
# RNA-seq Differential Gene Expression Analysis
#
# Author      : Hugo Lainé
# Affiliation : Gulbenkian Institute for Molecular Medicine, Oeiras, PT
# Contact     : hugo.laine@gimm.pt
#
# Project     : Transcriptomic characterisation of BMDMs from Spic KO and
#               Wild Type mice under Control and Heme treatment conditions.
#
# Organism    : Mus musculus (strain C57BL/6J)
# Cell type   : Bone marrow-derived macrophages (BMDMs)
# Reference   : GRCm39 (Ensembl release 111)
#
# Description :
#   This script performs the complete differential gene expression (DGE)
#   analysis from raw featureCounts output through to annotated result tables
#   and publication oriented figures.
#
#   Experimental design: 2×2 factorial
#     Factor 1 - Genotype  | Wild Type (WT) vs Spic knockout (KO)
#     Factor 2 - Treatment | Control (Ctrl) vs Heme
#     n = 3 biological replicates per group (12 samples total)
#
#   Four pairwise comparisons are performed --> four respective "effects tested":
#     1. KO_Ctrl vs WT_Ctrl  -->  Genotype effect at baseline
#     2. KO_Heme vs WT_Heme  -->  Genotype effect under Heme treatment
#     3. WT_Heme vs WT_Ctrl  -->  Heme treatment effect in WT macrophages
#     4. KO_Heme vs KO_Ctrl  -->  Heme treatment effect in Spic KO macrophages
#
# Reproducibility :
#   Set the three path variables in Section 1 to point to your local copies
#   of the input files. All other parameters are set as named constants to
#   facilitate reporting and reproduction of results.
#
# =============================================================================


# =============================================================================
# SECTION 0 | Dependencies
# =============================================================================

library("DESeq2")         # DGE analysis framework
library("ggplot2")        # Plotting
library("ggrepel")        # Non-overlapping labels on plots
library("dplyr")          # Data manipulation
library("pheatmap")       # Clustered heatmaps
library("data.table")     # Fast GTF parsing
library("stringr")        # String operations (GTF attribute extraction)
library("grid")           # Required by save_pheatmap_pdf()
library("matrixStats")    # rowVars() for top-variable gene selection

set.seed(42)              # Ensures reproducibility of stochastic label placement in ggrepel volcano plots
                          # DESeq2 results are deterministic and unaffected by this seed


# =============================================================================
# SECTION 1 - Input paths  (edit these three lines)
# =============================================================================

COUNTS_FILE <- "path/to/featurecounts_output.fc"
GTF_FILE    <- "path/to/Mus_musculus.GRCm39.111.gtf"
OUTPUT_DIR  <- "path/to/output/directory"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# SECTION 2 | Analysis parameters
# =============================================================================

# Minimum average raw count threshold for gene inclusion
MIN_AVG_COUNT    <- 5

# Differential expression significance thresholds
PADJ_THRESHOLD   <- 0.05    # Benjamini-Hochberg adjusted p-value

LOG2FC_THRESHOLD <- 0.585   # ~1.5-fold change (2^0.585 ≈ 1.5).
                            # Rationale: a 1.5-fold threshold is appropriate for detecting
                            # moderate regulatory changes relevant to transcription factor
                            # perturbation (Spic KO); a 2-fold cutoff would be overly
                            # conservative given the expected phenotype magnitude

# =============================================================================
# SECTION 3 | GTF annotation utility
#
#   Maps Ensembl gene IDs to HGNC/MGI gene symbols by parsing the GTF file.
#   Called once; the resulting gene_map is reused for all downstream tables.
# =============================================================================

create_gene_map <- function(gtf_path) {
  gtf         <- fread(gtf_path, header = FALSE, sep = "\t", skip = "#",
                       showProgress = FALSE)
  gene_entries <- gtf[V3 == "gene"]

  extract_attr <- function(attr_string, key) {
    pattern <- paste0(key, ' "([^"]*)"')
    match   <- str_match(attr_string, pattern)
    ifelse(is.na(match[, 2]), NA_character_, match[, 2])
  }

  gene_entries[, gene_id   := extract_attr(V9, "gene_id")]
  gene_entries[, gene_name := extract_attr(V9, "gene_name")]
  gene_map    <- unique(gene_entries[!is.na(gene_id), .(gene_id, gene_name)])
  gene_map[is.na(gene_name), gene_name := gene_id]
  return(gene_map)
}

add_gene_names <- function(de_table, gene_map) {
  idx        <- match(rownames(de_table), gene_map$gene_id)
  gene_names <- ifelse(is.na(idx), rownames(de_table), gene_map$gene_name[idx])
  cbind(gene_name = gene_names, de_table)
}

gene_map <- create_gene_map(GTF_FILE)


# =============================================================================
# SECTION 4 | Load and prepare count data
# =============================================================================

DF <- read.table(COUNTS_FILE, check.names = FALSE, header = TRUE, sep = "\t")
rownames(DF) <- DF[, 1]
DF[, 1]      <- NULL

# Note: samples were delivered with incorrect labels (see project documentation).
# Column names are corrected here to reflect the true experimental design.
colnames(DF) <- c(
  "10_WT_Heme", "12_WT_Heme", "13_KO_Ctrl", "14_KO_Ctrl", "15_KO_Ctrl",
  "19_KO_Heme",  "1_WT_Ctrl",  "20_KO_Heme", "23_KO_Heme",  "2_WT_Ctrl",
   "3_WT_Ctrl",  "7_WT_Heme"
)

# Reorder: WT before KO, Control before Heme within each genotype block.
# This ordering ensures WT and Control are reference levels in DESeq2.
DF <- DF[, c("1_WT_Ctrl",  "2_WT_Ctrl",  "3_WT_Ctrl",
             "7_WT_Heme",  "10_WT_Heme", "12_WT_Heme",
             "13_KO_Ctrl", "14_KO_Ctrl", "15_KO_Ctrl",
             "19_KO_Heme", "20_KO_Heme", "23_KO_Heme")]

# Remove genes with insufficient expression evidence
DF <- DF[rowMeans(DF) > MIN_AVG_COUNT, ]
message("Genes retained after low-count filtering:", nrow(DF), "\n")

# Annotate filtered count table
DF_GN <- add_gene_names(DF, gene_map)


# =============================================================================
# SECTION 5 | DESeq2 dataset construction
# =============================================================================

# Sample metadata
meta_DF <- DataFrame(
  row.names = colnames(DF),
  genotype  = factor(rep(c("WT", "KO"), times = c(6, 6)),
                     levels = c("WT", "KO")),
  treatment = factor(rep(c("Control", "Heme", "Control", "Heme"),
                         times = c(3, 3, 3, 3)),
                     levels = c("Control", "Heme"))
)

# DESeq2 object with interaction term to model condition-specific effects
# Note: gene_name column is excluded explicitly by name to avoid fragile positional indexing
count_cols <- setdiff(colnames(DF_GN), "gene_name")
dds <- DESeqDataSetFromMatrix(
  countData = DF_GN[, count_cols],
  colData   = meta_DF,
  design    = ~ genotype + treatment + genotype:treatment
)

# Explicit relevel: confirms reference levels (already set above; kept for clarity)
# Relevel: WT = reference genotype, Control = reference treatment
dds$genotype  <- relevel(dds$genotype,  ref = "WT")
dds$treatment <- relevel(dds$treatment, ref = "Control")

dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

message("\nDESeq2 model coefficients:\n")
print(resultsNames(dds))

# Normalised counts
norm_counts <- as.data.frame(counts(dds, normalized = TRUE))

# Variance-stabilising transformation for visualisation
# For unsupervised PCA plot = no prior on groups: VST computed with blind = TRUE 
vsd_pca     <- vst(dds, blind = TRUE)
vsd_pca_mat <- assay(vsd_pca)

# For Heatmaps showing gene expression profiles between known conditions (~ groups): VST computed with blind = FALSE. 
# Uses fitted dispersion estimates; appropriate for Heatmaps with prior on groups
vsd_hm     <- vst(dds, blind = FALSE)
vsd_hm_mat <- assay(vsd_hm)


# =============================================================================
# SECTION 6 | Principal Component Analysis
# =============================================================================

PCA  <- plotPCA(vsd_pca, intgroup = c("genotype", "treatment"), returnData = TRUE)
pct  <- round(100 * attr(PCA, "percentVar"))

p_pca <- ggplot(PCA, aes(x = PC1, y = PC2, fill = group)) +
  geom_point(size = 4, pch = 21, colour = "black") +
  geom_text_repel(aes(label = name), size = 4, colour = "#030303") +
  labs(title = "Principal Components Analysis",
       subtitle = "Unsupervised",
       x     = paste0("PC1: ", pct[1], "% variance"),
       y     = paste0("PC2: ", pct[2], "% variance"),
       fill  = "Genotype : Treatment") +
  coord_fixed() +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, face = "italic"))

pdf(file.path(OUTPUT_DIR, "PCA.pdf"), width = 12, height = 12)
print(p_pca)
dev.off()


# =============================================================================
# SECTION 7 | Heatmaps
# =============================================================================

#   Heatmap A | Top 500 most variably expressed genes (genome-wide overview)
#   Heatmap B | Gene panel I  - 9 iron metabolism / oxidative stress markers
#   Heatmap C | Gene panel II - 23 Nrf2/heme pathway genes


# -------------- Helper function for Heatmaps -----------------------------------
#
#' Build standard column annotation data frame for pheatmap
#'
#' Constructs the annotation_col data.frame expected by pheatmap from a
#' DESeqTransform object (output of vst() or rlog()).
#'
#' @param vsd_obj    DESeqTransform object with colData containing genotype
#'                   and treatment factors.
#' @param ann_colors Named list of colour vectors (see ANN_COLORS below).
#' @return           List with elements: annotation_col, annotation_colors.

build_heatmap_annotation <- function(vsd_obj, ann_colors) {

  cd <- colData(vsd_obj)

  annotation_col <- data.frame(
    Genotype  = cd$genotype,
    Treatment = cd$treatment,
    Group     = paste(cd$genotype, cd$treatment, sep = "_"),
    row.names = colnames(assay(vsd_obj))
  )

  return(list(annotation_col   = annotation_col,
              annotation_colors = ann_colors))
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Save heatmap to PDF
save_pheatmap_pdf <- function(x, filename, width = 15, height = 10) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Plot annotation colours  (shared across heatmaps and PCA)


ANN_COLORS <- list(
  Genotype  = c(WT = "#008B00",  KO = "#EEEE00"),
  Treatment = c(Control = "blue", Heme = "#8B1A1A"),
  Group     = c(
    WT_Control = "#00FF00",
    WT_Heme    = "#CD2626",
    KO_Control = "pink",
    KO_Heme    = "#FF6A6A"
  )
)

# Build column annotation (shared across all three heatmaps)
ann_info <- build_heatmap_annotation(vsd_hm, ANN_COLORS)

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Heatmap A | Top 500 most variable genes
# -----------------------------------------------------------------------------

message("Generating Heatmap A: top 500 variable genes...")

top_500_idx <- head(order(rowVars(vsd_hm_mat), decreasing = TRUE), 500)
mat       <- vsd_hm_mat[top_500_idx, ]

# Z-score scaling (scale="row") applied by pheatmap automatically
hm_top500 <- pheatmap(
  mat,
  scale            = "row",
  border_color     = "grey15",
  annotation_col   = ann_info$annotation_col,
  annotation_colors = ann_info$annotation_colors,
  annotation_legend = TRUE,
  show_rownames     = FALSE,
  main              = "Top 500 most variably expressed genes",
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  legend_breaks     = c(-2, 0, 2),
  legend_labels     = c("Low (-2)", "Medium", "High (+2)"),
  angle_col        = 315,
  cutree_cols      = 2,
  fontsize         = 12,
  fontsize_col     = 14,
  silent            = TRUE
)

save_pheatmap_pdf(hm_top500,
                  file.path(OUTPUT_DIR, "Heatmap_top_500_variable_genes.pdf"),
                  width = 14, height = 10)


# -----------------------------------------------------------------------------
# Heatmap B | Gene panel I: iron metabolism & oxidative stress markers
# -----------------------------------------------------------------------------

# Panel I - Iron metabolism & general oxidative stress markers
GOI_PANEL_I <- c(
  "ENSMUSG00000004359" = "Spic",
  "ENSMUSG00000015839" = "Nfe2l2",   # Nrf2
  "ENSMUSG00000001999" = "Blvra",
  "ENSMUSG00000024661" = "Fth1",
  "ENSMUSG00000006818" = "Sod2",
  "ENSMUSG00000022982" = "Sod1",
  "ENSMUSG00000025993" = "Slc40a1",
  "ENSMUSG00000027737" = "Slc7a11",
  "ENSMUSG00000022797" = "Tfrc"
)

message("Generating Heatmap B: Gene panel I (", length(GOI_PANEL_I), " genes)...")

# Checking if some GOI from Panel I have been filtered out for data curation at the begining of the analysis
# Remove any GOI from Panel I that did not pass the low-count filter
missing_B <- setdiff(names(GOI_PANEL_I), rownames(vsd_hm_mat))
if (length(missing_B) > 0) {
  warning("The following GOI Panel I genes were filtered out and will be skipped: ",
          paste(GOI_PANEL_I[missing_B], collapse = ", "))
  GOI_PANEL_I <- GOI_PANEL_I[!names(GOI_PANEL_I) %in% missing_B]
}

mat_B           <- vsd_hm_mat[names(GOI_PANEL_I), ]
rownames(mat_B) <- GOI_PANEL_I   # replace ENSEMBL IDs with gene symbols

# Heatmap B - Gene panel I
# cluster_rows = FALSE: rows are ordered as defined in GOI_PANEL_I to preserve
# biological pathway grouping; hierarchical clustering is not applied to curated panels.

heatmap_goi1 <- pheatmap(
  mat_B,
  scale             = "row",
  border_color      = "grey15",
  annotation_col    = ann_info$annotation_col,
  annotation_colors = ann_info$annotation_colors,
  annotation_legend = TRUE,
  show_rownames     = TRUE,
  main              = "Iron metabolism & oxidative stress markers",
  cluster_rows      = FALSE,
  cluster_cols      = TRUE,
  cellwidth         = 80,
  cellheight        = 30,
  legend_breaks     = c(-1.5, 0, 1.5),   # Narrower scale than Heatmap A/C: Panel I genes span a smaller dynamic range
  legend_labels     = c("Low (-1.5)", "Medium (0)", "High (+1.5)"),
  angle_col         = 315,
  fontsize          = 12,
  fontsize_row      = 13,
  fontsize_col      = 13,
  silent            = TRUE
)

out_B <- file.path(OUTPUT_DIR, "Heatmap_GOI_gene_panel_I.pdf")
save_pheatmap_pdf(heatmap_goi1, out_B, width = 20, height = 15)


# -----------------------------------------------------------------------------
# Heatmap C | Gene panel II: Nrf2/heme response pathway
# -----------------------------------------------------------------------------

# Panel II - Nrf2/heme response pathway
GOI_PANEL_II <- c(
  "ENSMUSG00000004359" = "Spic",
  "ENSMUSG00000005413" = "Hmox1",
  "ENSMUSG00000027962" = "Vcam1",
  "ENSMUSG00000051682" = "Treml4",
  "ENSMUSG00000025993" = "Slc40a1",
  "ENSMUSG00000003849" = "Nqo1",
  "ENSMUSG00000031584" = "Gsr",
  "ENSMUSG00000028124" = "Gclm",
  "ENSMUSG00000032350" = "Gclc",
  "ENSMUSG00000027610" = "Gss",
  "ENSMUSG00000060803" = "Gstp1",
  "ENSMUSG00000020250" = "Txnrd1",
  "ENSMUSG00000028691" = "Prdx1",
  "ENSMUSG00000027187" = "Cat",
  "ENSMUSG00000027737" = "Slc7a11",
  "ENSMUSG00000031400" = "G6pdx",
  "ENSMUSG00000028961" = "Pgd",
  "ENSMUSG00000032418" = "Me1",
  "ENSMUSG00000025950" = "Idh1",
  "ENSMUSG00000001999" = "Blvra",
  "ENSMUSG00000021696" = "Elovl7",
  "ENSMUSG00000031257" = "Nox1",
  "ENSMUSG00000027313" = "Chac1"
)

message("Generating Heatmap C: Gene panel II (", length(GOI_PANEL_II), " genes)...")

# Checking if some GOI from Panel II have been filtered out for data curation at the begining of the analysis
# Remove any GOI from Panel II that did not pass the low-count filter
missing_C <- setdiff(names(GOI_PANEL_II), rownames(vsd_hm_mat))
if (length(missing_C) > 0) {
  warning("The following GOI Panel II genes were filtered out and will be skipped: ",
          paste(GOI_PANEL_II[missing_C], collapse = ", "))
  GOI_PANEL_II <- GOI_PANEL_II[!names(GOI_PANEL_II) %in% missing_C]
}

mat_C           <- vsd_hm_mat[names(GOI_PANEL_II), ]
rownames(mat_C) <- GOI_PANEL_II

# Heatmap C - Gene panel II
# cluster_rows = FALSE: rows are ordered as defined in GOI_PANEL_II to preserve
# biological pathway grouping; hierarchical clustering is not applied to curated panels.

heatmap_goi2 <- pheatmap(
  mat_C,
  scale             = "row",
  border_color      = "grey15",
  annotation_col    = ann_info$annotation_col,
  annotation_colors = ann_info$annotation_colors,
  annotation_legend = TRUE,
  show_rownames     = TRUE,
  main              = "Nrf2/heme response pathway",
  cluster_rows      = FALSE,
  cluster_cols      = TRUE,
  cellwidth         = 60,
  cellheight        = 25,
  legend_breaks     = c(-2, 0, 2),
  legend_labels     = c("Low (-2)", "Medium (0)", "High (+2)"),
  angle_col         = 315,
  fontsize          = 12,
  fontsize_row      = 12,
  fontsize_col      = 12,
  silent            = TRUE
)

out_C <- file.path(OUTPUT_DIR, "Heatmap_GOI_gene_panel_II.pdf")
save_pheatmap_pdf(heatmap_goi2, out_C, width = 20, height = 20)


message("All heatmaps complete.")


# =============================================================================
# SECTION 8 | Volcano plot function
# =============================================================================

volcano_plot <- function(dge_df, log2FC, padj_cutoff, top = 10) {

  dge_df$Expression <- "Non-DE"
  dge_df$Expression[dge_df$log2FoldChange >  log2FC & dge_df$padj < padj_cutoff] <- "Up"
  dge_df$Expression[dge_df$log2FoldChange < -log2FC & dge_df$padj < padj_cutoff] <- "Down"
  dge_df$Expression <- factor(dge_df$Expression, levels = c("Up", "Down", "Non-DE"))

  top_genes <- bind_rows(
    dge_df %>% filter(Expression == "Up")   %>% slice_max(log2FoldChange,  n = top),
    dge_df %>% filter(Expression == "Down")  %>% slice_min(log2FoldChange, n = top)
  )

  ggplot(dge_df, aes(x = log2FoldChange, y = -log10(padj),
                     colour = Expression, size = Mean_Gene_Exp.)) +
    geom_point(alpha = 0.6) +
    geom_vline(xintercept = c(-log2FC, log2FC),
               linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    scale_colour_manual(values = c(Up = "#CD2626", Down = "#00688B",
                                   `Non-DE` = "grey70")) +
    geom_text_repel(data = top_genes, aes(label = GeneID),
                    size = 3.5, colour = "#1A1A1A",
                    max.overlaps = 1000) +
    xlab(bquote("log"[2] ~ "FC")) +
    ylab(bquote("-log"[10] ~ "(padj)")) +
    theme_bw(base_size = 13) +
    theme(plot.title   = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, face = "italic"))
}


# =============================================================================
# SECTION 9 | Differential Gene Expression: four pairwise comparisons
#
#   Results are extracted using the contrast= argument of results() in DESeq2.
#   The interaction model allows all four biologically meaningful comparisons
#   to be computed from a single fitted model.
#
#   For each comparison:
#     - Full result table saved as CSV (all tested genes, no NAs)
#     - Filtered DEG list saved as CSV (padj < 0.05, |log2FC| > 0.585)
#     - Volcano plot saved as PDF (labelled with gene symbols)
# =============================================================================

run_comparison <- function(label, contrast, samples_vec) {

  COMP_DIR <- file.path(OUTPUT_DIR, label)
  dir.create(COMP_DIR, showWarnings = FALSE, recursive = TRUE)
    
  res      <- results(dds, contrast = contrast)
  res_noNA <- res[complete.cases(as.data.frame(res)), ]

  # Build annotated result table
  dge <- cbind(norm_counts[rownames(res_noNA), samples_vec],
               as.data.frame(res_noNA))
  dge <- add_gene_names(dge, gene_map)

  # Save full result table
  write.csv(dge, file.path(COMP_DIR, paste0("DGE_", label, ".csv")),
            row.names = TRUE)

  # Filter DEGs
  degs <- dge[!is.na(dge$padj) &
                dge$padj < PADJ_THRESHOLD &
                abs(dge$log2FoldChange) > LOG2FC_THRESHOLD, ]
  degs$Direction <- ifelse(degs$log2FoldChange > 0, "Up", "Down")

  write.csv(degs, file.path(COMP_DIR, paste0("DEGs_", label, ".csv")),
            row.names = TRUE)

  message(sprintf("\n%s - Up: %d | Down: %d | Total DEGs: %d\n",
              label,
              sum(degs$Direction == "Up"),
              sum(degs$Direction == "Down"),
              nrow(degs)))

  # Volcano plot (gene symbol labels)
  vdf        <- dge
  colnames(vdf)[colnames(vdf) == "baseMean"] <- "Mean_Gene_Exp."
  vdf$GeneID <- vdf$gene_name

  p <- volcano_plot(vdf, LOG2FC_THRESHOLD, PADJ_THRESHOLD, top = 10) +
    ggtitle("Volcano plot", subtitle = label) +
    labs(caption = paste0(nrow(res_noNA), " genes  |  padj < ", PADJ_THRESHOLD))

  pdf(file.path(COMP_DIR, paste0("Volcano_", label, ".pdf")),
      width = 12, height = 18)
  print(p)
  dev.off()

  return(invisible(list(full = dge, DEGs = degs)))
}


# --- Comparison 1: Genotype effect at baseline (Control treatment) -----------
res_KO_Ctrl_vs_WT_Ctrl <- run_comparison(
  label       = "KO_Ctrl_vs_WT_Ctrl",
  contrast    = list("genotype_KO_vs_WT"),
  samples_vec = c("1_WT_Ctrl", "2_WT_Ctrl", "3_WT_Ctrl",
                  "13_KO_Ctrl", "14_KO_Ctrl", "15_KO_Ctrl")
)

# --- Comparison 2: Genotype effect under Heme treatment ---------------------
res_KO_Heme_vs_WT_Heme <- run_comparison(
  label       = "KO_Heme_vs_WT_Heme",
  contrast    = list(c("genotype_KO_vs_WT", "genotypeKO.treatmentHeme")),
  samples_vec = c("7_WT_Heme", "10_WT_Heme", "12_WT_Heme",
                  "19_KO_Heme", "20_KO_Heme", "23_KO_Heme")
)

# --- Comparison 3: Heme treatment effect in WT macrophages ------------------
res_WT_Heme_vs_WT_Ctrl <- run_comparison(
  label       = "WT_Heme_vs_WT_Ctrl",
  contrast    = list("treatment_Heme_vs_Control"),
  samples_vec = c("1_WT_Ctrl", "2_WT_Ctrl", "3_WT_Ctrl",
                  "7_WT_Heme", "10_WT_Heme", "12_WT_Heme")
)

# --- Comparison 4: Heme treatment effect in Spic KO macrophages -------------
res_KO_Heme_vs_KO_Ctrl <- run_comparison(
  label       = "KO_Heme_vs_KO_Ctrl",
  contrast    = list(c("treatment_Heme_vs_Control", "genotypeKO.treatmentHeme")),
  samples_vec = c("13_KO_Ctrl", "14_KO_Ctrl", "15_KO_Ctrl",
                  "19_KO_Heme", "20_KO_Heme", "23_KO_Heme")
)


# =============================================================================
# SECTION 10 | Session information
# =============================================================================

session_log <- file.path(OUTPUT_DIR, "session_info.txt")
writeLines(capture.output(sessionInfo()), session_log)
message("Session information written to:", session_log, "\n")

message("\n── Session information ───────────────────────────────────────────────\n")
sessionInfo()
################################################################################
# Seurat Clustering Analysis and Visualization
# 
# Description: Downstream analysis of integrated single-cell RNA-seq data
# including cluster identification, marker gene analysis, and visualization
# for manuscript figures
#
# Requirements: 
# - Integrated Seurat object (Cond.combineds_integrated.rds)
# - R packages: Seurat, dplyr, tidyverse, ggplot2, pheatmap, 
#   RColorBrewer, SeuratObject, scales
################################################################################

# Load required libraries
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(SeuratObject)

################################################################################
# USER-DEFINED PATHS - MODIFY THESE FOR YOUR SYSTEM
################################################################################

# Path to integrated Seurat object
integrated_object_path <- "Data_Files/Files/Cond.combineds_integrated.rds"

# Path to save Schwann Cell subset (optional, for later use)
schwann_cells_path <- "Data_Files/Files/SchwannCells.rds"

# Path to save marker genes table
markers_output_path <- "Data_Files/Files/Table1.csv"

################################################################################
# LOAD DATA
################################################################################

cat("=== Loading Integrated Seurat Object ===\n")

# Read in already saved object (can do this after generating Cond.combined_integrated.rds)
Cond.combineds_integrated <- readRDS(integrated_object_path)

# SchwannCells can be read in later after this code has been run:
# SchwannCells <- readRDS(schwann_cells_path)

################################################################################
# FIGURE 1G - PART 1: UMAP WITH CLUSTER LABELS
################################################################################

cat("\n=== Generating Figure 1G - Part 1 ===\n")

Idents(Cond.combineds_integrated) <- "seurat_clusters"
DimPlot(Cond.combineds_integrated, reduction = "umap", label = TRUE, pt.size = 0.01)

# Validate that clustering is not happening on amount of genes detected in each cell
# FeaturePlot(Cond.combineds_integrated, features = c('nFeature_RNA'), pt.size = 0.5)
# VlnPlot(Cond.combineds_integrated, features = c("nFeature_RNA"))

################################################################################
# FIGURE 1G - PART 2: UMAP SPLIT BY CONDITION
################################################################################

cat("\n=== Generating Figure 1G - Part 2 ===\n")

# Create a new column for plotting groups, preserving the original condition names
Cond.combineds_integrated$plot_groups <- as.character(Cond.combineds_integrated$condition)

# Reassign both Sham and Naive to the same value to combine them on one plot
Cond.combineds_integrated$plot_groups[Cond.combineds_integrated$plot_groups %in% c("Sham", "Naive")] <- "Sham_and_Naive"

# Convert to factor with specific level order to control plot arrangement
Cond.combineds_integrated$plot_groups <- factor(
  Cond.combineds_integrated$plot_groups,
  levels = c("Sham_and_Naive", "2hpi", "1dpi")
)

# Generate split UMAP
DimPlot(Cond.combineds_integrated, 
        reduction = "umap", 
        label = FALSE, 
        split.by = "plot_groups",
        ncol = 2)

################################################################################
# FIND CLUSTER MARKERS
################################################################################

cat("\n=== Finding Cluster Marker Genes ===\n")

# Find top significant genes enriched in each cluster
Idents(Cond.combineds_integrated) <- "seurat_clusters"
combined.integrated.markers <- FindAllMarkers(Cond.combineds_integrated, 
                                              only.pos = TRUE, 
                                              min.pct = 0.25, 
                                              logfc.threshold = 0.25)

cat("Total markers found:", nrow(combined.integrated.markers), "\n")

# Save marker genes table (Table 1)
write.csv(combined.integrated.markers, file = markers_output_path, row.names = FALSE)
cat("Marker genes saved to:", markers_output_path, "\n")

# Get top 5 markers per cluster
top5 <- combined.integrated.markers %>% 
  filter(p_val_adj < 0.05) %>%          # Filter for significant p-values first
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)         # Then select top 5 by fold change

cat("Top 5 markers per cluster identified\n")

################################################################################
# FIGURE 1H: STACKED VIOLIN PLOT OF TOP MARKERS
################################################################################

cat("\n=== Generating Figure 1H ===\n")

DefaultAssay(Cond.combineds_integrated) <- "RNA"
VlnPlot(Cond.combineds_integrated, 
        features = top5$gene,
        add.noise = FALSE,
        stack = TRUE, 
        flip = TRUE, 
        fill.by = "ident") +
  theme(legend.position = "none",
        strip.text.x = element_text(angle = 90, 
                                    hjust = 0, 
                                    vjust = 0.5,
                                    face = "plain",
                                    size = 8))

################################################################################
# FIGURE S2B: HEATMAP OF TOP MARKERS
################################################################################

cat("\n=== Generating Figure S2B ===\n")

# Set RNA assay and get average expression
DefaultAssay(Cond.combineds_integrated) <- "RNA"
cluster.averages <- AverageExpression(Cond.combineds_integrated, 
                                      group.by = c("seurat_clusters"),
                                      assays = "RNA", 
                                      features = top5$gene)

# Convert the average expression data to a matrix
expr_mat <- as.matrix(cluster.averages$RNA)

# Order genes by cluster
ordered_genes <- top5 %>%
  arrange(cluster) %>%
  pull(gene)

# Reorder the matrix rows according to the ordered genes
expr_mat <- expr_mat[ordered_genes, ]

# Generate heatmap
pheatmap(expr_mat, 
         border_color = "NA", 
         scale = "row", 
         color = cm.colors(256), 
         show_rownames = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         angle_col = "45",
         main = "Cluster average expressoin of top 5 expressed markers for each cluster")

################################################################################
# PREPARE DATA FOR STACKED BAR PLOTS
################################################################################

cat("\n=== Preparing Data for Composition Plots ===\n")

# Filter the data to include only "Sham" or "Naive" in Hashtags
filtered_data_shams_and_naive <- Cond.combineds_integrated@meta.data %>%
  filter(grepl("Sham|Naive", Hashtags))

# Filter the data to include only injury samples
filtered_data_injury <- Cond.combineds_integrated@meta.data %>%
  filter(!grepl("Sham|Naive", Hashtags) & !grepl("Naive", Hashtags))

# Create color palettes
color_palette_paired_sham <- brewer.pal(nlevels(factor(filtered_data_shams_and_naive$Hashtags)), 
                                        "Paired")
color_palette_paired_injury <- brewer.pal(nlevels(factor(filtered_data_injury$Hashtags)), 
                                          "Paired")

# Define label order
label_order <- c(
  "Naive",
  "Sham-for-2hpi-WT",
  "Sham-for-1dpi-WT",
  "Sham-for-1dpi-Sarm1-KO",
  "Sham-for-3dpi-WT",
  "2hpi-WT",
  "1dpi-WT",
  "1dpi-Sarm1-KO",
  "3dpi-WT",
  "3dpi-Sarm1-KO"
)

# Convert label to factor with specific order in both datasets
filtered_data_shams_and_naive$label <- factor(filtered_data_shams_and_naive$label, 
                                              levels = label_order)
filtered_data_injury$label <- factor(filtered_data_injury$label, 
                                     levels = label_order)

################################################################################
# FIGURE 1I: COMPOSITION OF SHAM/NAIVE SAMPLES
################################################################################

cat("\n=== Generating Figure 1I ===\n")

ggplot(filtered_data_shams_and_naive, aes(x = label, fill = seurat_clusters)) + 
  geom_bar(position = "fill") +
  theme(panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

################################################################################
# FIGURE 1J: COMPOSITION OF INJURY SAMPLES
################################################################################

cat("\n=== Generating Figure 1J ===\n")

ggplot(filtered_data_injury, aes(x = label, fill = seurat_clusters)) + 
  geom_bar(position = "fill") +
  theme(panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

################################################################################
# FIGURE S2J: CLUSTER COMPOSITION FOR SHAM/NAIVE (SPLIT)
################################################################################

cat("\n=== Generating Figure S2J ===\n")

# Calculate total count from original dataset
total_count <- nrow(filtered_data_shams_and_naive)

# Figure S2J - Part 1: Clusters 0-6 for sham cells
cluster0_6 <- filtered_data_shams_and_naive[filtered_data_shams_and_naive$seurat_clusters %in% 0:6, ]
ggplot(cluster0_6, aes(x = seurat_clusters, fill = label)) + 
  geom_bar(aes(y = after_stat(count)/total_count), position = "stack") +
  theme(panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank()) +
  scale_fill_manual(values = color_palette_paired_injury) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

# Figure S2J - Part 2: Clusters 7-17 for sham cells
cluster7_17 <- filtered_data_shams_and_naive[filtered_data_shams_and_naive$seurat_clusters %in% 7:17, ]
ggplot(cluster7_17, aes(x = seurat_clusters, fill = label)) + 
  geom_bar(aes(y = after_stat(count)/total_count), position = "stack") +
  theme(panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank()) +
  scale_fill_manual(values = color_palette_paired_injury) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

################################################################################
# FIGURE S2K: CLUSTER COMPOSITION FOR INJURY (SPLIT)
################################################################################

cat("\n=== Generating Figure S2K ===\n")

# Calculate total count from original dataset
total_count <- nrow(filtered_data_injury)

# Figure S2K - Part 1: Clusters 0-7 for injured SCs
cluster0_7 <- filtered_data_injury[filtered_data_injury$seurat_clusters %in% 0:7, ]
ggplot(cluster0_7, aes(x = seurat_clusters, fill = label)) + 
  geom_bar(aes(y = after_stat(count)/total_count), position = "stack") +
  theme(panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank()) +
  scale_fill_manual(values = color_palette_paired_injury) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

# Figure S2K - Part 2: Clusters 8-17 for injured SCs
cluster8_17 <- filtered_data_injury[filtered_data_injury$seurat_clusters %in% 8:17, ]
ggplot(cluster8_17, aes(x = seurat_clusters, fill = label)) + 
  geom_bar(aes(y = after_stat(count)/total_count), position = "stack") +
  theme(panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank()) +
  scale_fill_manual(values = color_palette_paired_injury) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

################################################################################
# FIGURE S2I: SCHWANN CELL MARKER EXPRESSION
################################################################################

cat("\n=== Generating Figure S2I ===\n")

# Which clusters are Schwann Cells based on commonly used myelinating 
# and non-myelinating Schwann Cell markers
# Expected: clusters 2, 5, 6, 7, 8
DotPlot(Cond.combineds_integrated, features = c("Mbp", "Ncam1", "Cadm2"))

################################################################################
# FIGURE S2C-E: QC METRICS BY CLUSTER
################################################################################

cat("\n=== Generating Figures S2C-E (QC by Cluster) ===\n")

# Figure S2C: UMI counts by cluster
VlnPlot(Cond.combineds_integrated, 
        features = "nCount_RNA", 
        group.by = "seurat_clusters", 
        pt.size = 0)

# Figure S2D: Gene counts by cluster
VlnPlot(Cond.combineds_integrated, 
        features = "nFeature_RNA", 
        group.by = "seurat_clusters", 
        pt.size = 0)

# Figure S2E: Mitochondrial percentage by cluster
VlnPlot(Cond.combineds_integrated, 
        features = "percent.mt", 
        group.by = "seurat_clusters", 
        pt.size = 0)

################################################################################
# FIGURE S2F-H: QC METRICS BY BATCH
################################################################################

cat("\n=== Generating Figures S2F-H (QC by Batch) ===\n")

# Figure S2F: UMI counts by batch
VlnPlot(Cond.combineds_integrated, 
        features = "nCount_RNA", 
        group.by = "batch.orig", 
        pt.size = 0)

# Figure S2G: Gene counts by batch
VlnPlot(Cond.combineds_integrated, 
        features = "nFeature_RNA", 
        group.by = "batch.orig", 
        pt.size = 0)

# Figure S2H: Mitochondrial percentage by batch
VlnPlot(Cond.combineds_integrated, 
        features = "percent.mt", 
        group.by = "batch.orig", 
        pt.size = 0)

cat("\n=== Analysis Complete ===\n")
cat("All figures and tables have been generated\n")

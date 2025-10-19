################################################################################
# Injured Schwann Cell Analysis
# 
# Description: Analysis of injured Schwann cells including subsetting,
# reclustering, and marker identification for 1dpi and 2hpi samples
#
# Requirements: 
# - Schwann Cells Seurat object (SchwannCells.rds)
# - R packages: Seurat, dplyr, tidyverse, ggplot2, pheatmap, RColorBrewer,
#   SeuratObject, gridExtra
################################################################################

# Load required libraries
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(SeuratObject)
library(gridExtra)

################################################################################
# USER-DEFINED PATHS - MODIFY THESE FOR YOUR SYSTEM
################################################################################

# Path to Schwann Cells Seurat object
schwann_cells_path <- "Data_Files/Files/SchwannCells.rds"

# Path to save injured SC marker genes table
injured_markers_output_path <- "Data_Files/Files/Table5.csv"

################################################################################
# LOAD DATA
################################################################################

cat("=== Loading Schwann Cells Object ===\n")

# SchwannCells <- readRDS(schwann_cells_path)

################################################################################
# FIGURE S4B: IDENTIFY INJURED SCHWANN CELLS
################################################################################

cat("\n=== Generating Figure S4B ===\n")

cat("Current layers:\n")
print(Layers(SchwannCells))

# Create a data frame with the necessary information to identify injured SCs
df <- data.frame(UMAP1 = SchwannCells@reductions$umap@cell.embeddings[, 1],
                 UMAP2 = SchwannCells@reductions$umap@cell.embeddings[, 2],
                 label = SchwannCells@meta.data$label)

# Create a new column to identify Injured samples
cat("Unique labels:\n")
print(unique(df$label))

df$is_injury <- !grepl("Sham|Naive", df$label)

# Figure S4B
ggplot(data = df) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = is_injury),
             size = 0.5) +
  labs(title = "Injury Schwann Cell") +
  scale_color_manual(values = c("FALSE" = "lightgrey", "TRUE" = "red"),
                     labels = c("Sham", "Injury"),
                     name = "Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_line())

################################################################################
# SUBSET INJURED SCHWANN CELLS
################################################################################

cat("\n=== Subsetting Injured Schwann Cells ===\n")

# Read in Schwann Cells object if not already loaded
SchwannCells <- readRDS(schwann_cells_path)

DefaultAssay(SchwannCells) <- "RNA"
Idents(SchwannCells) <- "label"

cat("Total Schwann cells:", length(WhichCells(SchwannCells)), "\n")
cat("Unique labels:\n")
print(unique(SchwannCells@meta.data$label))

# Check if layers are split by batch
cat("Current layers:\n")
print(Layers(SchwannCells))

# Run only if the layers are not split by batch already
#SchwannCells[["RNA"]] <- split(SchwannCells[["RNA"]], f = SchwannCells$batch.orig)

cat("Layers after splitting:\n")
print(Layers(SchwannCells))

# Subset only injured samples
InjuredSCs <- subset(SchwannCells, 
                     cells = WhichCells(SchwannCells, 
                                        idents = c("1dpi-Sarm1-KO", "1dpi-WT", "2hpi-WT")))

cat("Injured Schwann cells:", length(WhichCells(InjuredSCs)), "\n")
cat("Unique batches:\n")
print(unique(InjuredSCs@meta.data$batch.orig))
cat("Unique labels:\n")
print(unique(InjuredSCs@meta.data$label))

################################################################################
# RECLUSTER INJURED SCHWANN CELLS
################################################################################

cat("\n=== Reclustering Injured Schwann Cells ===\n")

# Use your existing InjuredSCs object
DefaultAssay(InjuredSCs) <- "RNA"
InjuredSCs <- NormalizeData(InjuredSCs)
InjuredSCs <- FindVariableFeatures(InjuredSCs)
InjuredSCs <- ScaleData(InjuredSCs)
InjuredSCs <- RunPCA(InjuredSCs)

# Visualize elbow plot
ElbowPlot(InjuredSCs)

# Number of principal components
x <- 8

# Continue with downstream analysis
InjuredSCs <- FindNeighbors(InjuredSCs, reduction = "pca", dims = 1:x)
InjuredSCs <- FindClusters(InjuredSCs, resolution = 0.25)
InjuredSCs <- RunUMAP(InjuredSCs, reduction = "pca", dims = 1:x, seed.use = 42)

################################################################################
# FIGURE S4C: UMAP BY CLUSTER
################################################################################

cat("\n=== Generating Figure S4C ===\n")

DimPlot(InjuredSCs, reduction = "umap", group.by = c("seurat_clusters"))

################################################################################
# PREPARE METADATA FOR PLOTTING
################################################################################

cat("\n=== Preparing Metadata ===\n")

# Set label order
cat("Unique labels:\n")
print(unique(InjuredSCs@meta.data$label))

InjuredSCs$label <- factor(InjuredSCs$label, 
                           levels = c("2hpi-WT", "1dpi-WT", "1dpi-Sarm1-KO"))

cat("Unique Hashtags:\n")
print(unique(InjuredSCs$Hashtags))

InjuredSCs$Hashtags <- factor(InjuredSCs$Hashtags, 
                              levels = c("2hpi-WT-Female", "2hpi-WT-Male",
                                         "1dpi-WT-Female", "1dpi-WT-Male",
                                         "1dpi-Sarm1-KO-Female", "1dpi-Sarm1-KO-Male"))

filtered_data <- InjuredSCs@meta.data

################################################################################
# FIGURE S4F-G: COMPOSITION PLOTS
################################################################################

cat("\n=== Generating Figures S4F-G ===\n")

# Figure S4F: Cluster composition by label
ggplot(filtered_data, aes(x = label, fill = seurat_clusters)) + 
  geom_bar(position = "fill") +
  theme(panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank())

# Optional: by Hashtags (commented out)
# ggplot(filtered_data, aes(x = Hashtags, fill = seurat_clusters)) + 
#   geom_bar(position = "fill") +
#   theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank())

# Figure S4G: Label composition by cluster
ggplot(filtered_data, aes(x = seurat_clusters, fill = label)) + 
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("#004D40", "#1E88E5", "#D81B60")) +
  theme(panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank())

# Optional: by Hashtags (commented out)
# ggplot(filtered_data, aes(x = seurat_clusters, fill = Hashtags)) + 
#   geom_bar(position = "fill") +
#   theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank())

################################################################################
# FIND INJURED SCHWANN CELL CLUSTER MARKERS
################################################################################

cat("\n=== Finding Injured Schwann Cell Cluster Markers ===\n")

DefaultAssay(InjuredSCs) <- "RNA"
InjuredSCs[["RNA"]] <- JoinLayers(InjuredSCs[["RNA"]])

Idents(InjuredSCs) <- "seurat_clusters"
SC_Injury.integrated.markers <- FindAllMarkers(InjuredSCs, 
                                               only.pos = TRUE, 
                                               min.pct = 0.25,
                                               logfc.threshold = 0.25)

cat("Total markers found:", nrow(SC_Injury.integrated.markers), "\n")

# Save marker genes table (Table 5)
write.csv(SC_Injury.integrated.markers, file = injured_markers_output_path, row.names = FALSE)
cat("Injured SC markers saved to:", injured_markers_output_path, "\n")

################################################################################
# FIGURE S4D: HEATMAP OF TOP MARKERS
################################################################################

cat("\n=== Generating Figure S4D ===\n")

# Get top 5 markers per cluster (filtered by significance)
top5_Injury_SCs <- SC_Injury.integrated.markers %>% 
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)

# Set RNA assay and get average expression
DefaultAssay(InjuredSCs) <- "RNA"
cluster.averages_injury <- AverageExpression(InjuredSCs, 
                                             group.by = "seurat_clusters",
                                             assays = "RNA", 
                                             features = top5_Injury_SCs$gene)

# Convert the average expression data to a matrix
expr_mat_injury_SCs <- as.matrix(cluster.averages_injury$RNA)

# Order genes by cluster
ordered_genes_injured_SCs <- top5_Injury_SCs %>%
  arrange(cluster) %>%
  pull(gene)

# Reorder the matrix rows according to the ordered genes
expr_mat_injury_SCs <- expr_mat_injury_SCs[ordered_genes_injured_SCs, ]

# Generate heatmap
pheatmap(expr_mat_injury_SCs, 
         border_color = "NA", 
         scale = "row", 
         color = cm.colors(256), 
         show_rownames = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Injury Cluster−Average Expression of Top Genes for Each Schwann Cell Cluster")

################################################################################
# FIGURE S4E: SPLIT BY LABEL
################################################################################

cat("\n=== Generating Figure S4E ===\n")

# Create dataframe with UMAP coordinates and condition/label
df <- data.frame(UMAP1 = InjuredSCs@reductions$umap@cell.embeddings[, 1],
                 UMAP2 = InjuredSCs@reductions$umap@cell.embeddings[, 2],
                 label = InjuredSCs@meta.data$label)

# Function to create plot for a specific label with red highlighting
create_label_plot <- function(data, label_name) {
  # Create a new color column where cells of interest are red, others are light grey
  plot_data <- data %>%
    mutate(highlight_color = ifelse(label == label_name, "red", "lightgrey"))
  
  ggplot(data = plot_data) +
    geom_point(aes(x = UMAP1, y = UMAP2, color = highlight_color), size = 0.2, alpha = 0.5) +
    labs(title = label_name) +
    scale_color_identity(guide = "none") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 8),
          panel.grid = element_blank(),
          axis.line = element_line(),
          axis.title = element_text(size = 8))
}

# Create a list of plots for all unique labels
label_plots <- lapply(unique(df$label), function(label) {
  create_label_plot(df, label)
})

# Display all plots in a grid arrangement
grid.arrange(grobs = label_plots, ncol = 1)

################################################################################
# FIGURE S4H: PASC MARKER VIOLIN PLOTS
################################################################################

cat("\n=== Generating Figure S4H ===\n")

# Figure S4H - Part 1: By cluster
VlnPlot(InjuredSCs, features = c("Plp1", "Fth1", "Ptgds"), group.by = "seurat_clusters")

# Figure S4H - Part 2: By label
VlnPlot(InjuredSCs, features = c("Plp1", "Fth1", "Ptgds"), group.by = "label")

################################################################################
# FIGURE S4I: PASC MARKER FEATURE PLOTS
################################################################################

cat("\n=== Generating Figure S4I ===\n")

FeaturePlot(InjuredSCs, features = c("Plp1", "Fth1", "Ptgds"), ncol = 3)

################################################################################
# MITOCHONDRIAL GENES ANALYSIS
################################################################################

cat("\n=== Mitochondrial Genes Analysis ===\n")
cat("NOTE: Run this after running 'Figure 3' file to define oxphos_genes\n")

# NOTE: These sections require oxphos_genes variable from "Figure 3" analysis
# Uncomment and run after oxphos_genes is defined

# DefaultAssay(InjuredSCs) <- "RNA"
# Idents(InjuredSCs) <- "label"
# unique_genes_combined <- oxphos_genes
# 
# # Order clusters
# cat("Unique clusters:\n")
# print(unique(InjuredSCs$seurat_clusters))
# 
# InjuredSCs$seurat_clusters <- factor(InjuredSCs$seurat_clusters, 
#                                      levels = c("0", "1", "2", "3", '4', '5', '6', '7', '8'))
# 
# # Figure S4K: By label
# DotPlot(InjuredSCs, features = unique_genes_combined, group.by = "label") +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
#     axis.text.y = element_text(size = 10)
#   ) +
#   scale_x_discrete(label = function(x) gsub(",", "", x)) +
#   theme(plot.margin = margin(b = 50))
# 
# # Figure S4J: By cluster
# DotPlot(InjuredSCs, features = unique_genes_combined, group.by = "seurat_clusters") +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
#     axis.text.y = element_text(size = 10)
#   ) +
#   scale_x_discrete(label = function(x) gsub(",", "", x)) +
#   theme(plot.margin = margin(b = 50))

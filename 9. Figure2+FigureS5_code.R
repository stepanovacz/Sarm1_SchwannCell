################################################################################
# Schwann Cell Analysis and Pseudotime Trajectory - Figures 2 and supplemental Figure 5
# 
# Description: Subset Schwann cells from integrated dataset, perform
# reclustering, marker identification, and pseudotime trajectory analysis
# using Monocle3
#
# Requirements: 
# - Integrated Seurat object (Cond.combineds_integrated.rds)
# - R packages: Seurat, dplyr, tidyverse, ggplot2, future, RColorBrewer,
#   pheatmap, SeuratObject, gridExtra, monocle3, SeuratWrappers, ggridges,
#   SingleCellExperiment
# - Note: monocle3 and SeuratWrappers may need custom library paths
################################################################################

# Load required libraries
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(future)
library(RColorBrewer)
library(pheatmap)
library(SeuratObject)
library(gridExtra)
library(monocle3)
library(SeuratWrappers)
library(ggridges)
library(SingleCellExperiment)

################################################################################
# USER-DEFINED PATHS - MODIFY THESE FOR YOUR SYSTEM
################################################################################

# Path to integrated Seurat object
integrated_object_path <- "Data_Files/Files/Cond.combineds_integrated.rds"

# Path to save Schwann Cell subset
schwann_cells_output_path <- "Data_Files/Files/SchwannCells.rds"

# Path to save marker genes table
markers_output_path <- "Data_Files/Files/Table2.csv"

################################################################################
# LOAD DATA
################################################################################

cat("=== Loading Integrated Seurat Object ===\n")

# Read in all the cell data if already saved as a seurat object
Cond.combineds_integrated <- readRDS(integrated_object_path)

################################################################################
# IDENTIFY SCHWANN CELL CLUSTERS
################################################################################

cat("\n=== Identifying Schwann Cell Clusters ===\n")

DefaultAssay(Cond.combineds_integrated) <- "RNA"
Idents(Cond.combineds_integrated) <- "seurat_clusters"

# First figure out which clusters are Schwann Cells
# DotPlot(Cond.combineds_integrated, features = c("Mbp", "Ncam1", "Cadm2"))
# Expected: clusters 2, 5, 6, 7, 8

################################################################################
# FIGURE 2A, S3A, S4A: HIGHLIGHT SCHWANN CELLS ON FULL UMAP
################################################################################

cat("\n=== Generating Figures 2A, S3A, S4A ===\n")

# Highlighting only clusters 2, 5, 6, 7, 8 - Schwann Cells
df_SCs <- data.frame(UMAP1 = Cond.combineds_integrated@reductions$umap@cell.embeddings[, 1],
                     UMAP2 = Cond.combineds_integrated@reductions$umap@cell.embeddings[, 2],
                     seurat_clusters = Cond.combineds_integrated@meta.data$seurat_clusters)

# Create a new column for colors based on seurat_clusters
SchwannCells_highlight <- df_SCs %>%
  mutate(color = case_when(
    seurat_clusters == "2" ~ "black",
    seurat_clusters == "5" ~ "black",
    seurat_clusters == "8" ~ "black",
    seurat_clusters == "6" ~ "black",
    seurat_clusters == "7" ~ "black",
    TRUE ~ "lightgrey"
  ))

# Generate plot
ggplot(data = SchwannCells_highlight) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = color), size = 1) +
  labs(title = "Schwann Cells") +
  scale_color_identity(guide = "none") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_line())

################################################################################
# SUBSET SCHWANN CELLS
################################################################################

cat("\n=== Subsetting Schwann Cells ===\n")

# Subsetting Schwann Cells that are in clusters 2, 5, 6, 7, 8
Idents(Cond.combineds_integrated) <- "seurat_clusters"

SchwannCells <- subset(Cond.combineds_integrated, 
                       cells = WhichCells(Cond.combineds_integrated, 
                                          idents = c('2', "5", "6", '7', '8')))

# Double check that subsetting was done correctly by looking at SC markers
DotPlot(SchwannCells, features = c("Mbp", "Mpz", "Ncam1", "Cadm2"))

cat("Total Schwann cells:", length(WhichCells(SchwannCells)), "\n")
cat("Sham cells:", length(WhichCells(SchwannCells, expression = condition == "Sham")), "\n")
cat("1dpi cells:", length(WhichCells(SchwannCells, expression = condition == "1dpi")), "\n")
cat("2hpi cells:", length(WhichCells(SchwannCells, expression = condition == "2hpi")), "\n")

# Further subset to only include Sham, 1dpi, and 2hpi
Idents(SchwannCells) <- "condition"
SchwannCells <- subset(SchwannCells, 
                       cells = WhichCells(SchwannCells, 
                                          idents = c("Sham", '1dpi', "2hpi")))

cat("Final Schwann cells after filtering:", length(WhichCells(SchwannCells)), "\n")
cat("Conditions included:\n")
print(unique(SchwannCells$condition))

################################################################################
# RECLUSTERING SCHWANN CELLS
################################################################################

cat("\n=== Reclustering Schwann Cells ===\n")

DefaultAssay(SchwannCells) <- "RNA"

# Check batches present
cat("Batches present:\n")
print(unique(SchwannCells@meta.data$batch.orig))

# Create a unified layer
SchwannCells <- JoinLayers(SchwannCells)

# Split layers by batch for processing
SchwannCells[["RNA"]] <- split(SchwannCells[["RNA"]], f = SchwannCells$batch.orig)
cat("Schwann Cells object structure:\n")
print(SchwannCells)

# Normalization, feature selection, scaling, and PCA
SchwannCells <- NormalizeData(SchwannCells)
SchwannCells <- FindVariableFeatures(SchwannCells)
SchwannCells <- ScaleData(SchwannCells)
SchwannCells <- RunPCA(SchwannCells)

# Visualize elbow plot
ElbowPlot(SchwannCells)

# Clustering and UMAP
x <- 6
SchwannCells <- FindNeighbors(SchwannCells, reduction = "pca", dims = 1:x)
SchwannCells <- FindClusters(SchwannCells, resolution = 0.3)
SchwannCells <- RunUMAP(SchwannCells, dims = 1:x, reduction = "pca")

################################################################################
# FIGURE 2A: SCHWANN CELL UMAP WITH CLUSTERS
################################################################################

cat("\n=== Generating Figure 2A ===\n")

Idents(SchwannCells) <- "seurat_clusters"
DimPlot(SchwannCells, reduction = "umap", label = TRUE, pt.size = 0.05)

# Save the seurat object so you don't have to rerun the code above
saveRDS(SchwannCells, schwann_cells_output_path)
cat("Schwann Cells object saved to:", schwann_cells_output_path, "\n")

################################################################################
# FIGURE 2B (FIGURE 3SB): SHAM VS INJURY SAMPLES
################################################################################

cat("\n=== Generating Figure 2B (Figure 3SB) ===\n")

# Create a data frame with the necessary information
df <- data.frame(UMAP1 = SchwannCells@reductions$umap@cell.embeddings[, 1],
                 UMAP2 = SchwannCells@reductions$umap@cell.embeddings[, 2],
                 label = SchwannCells@meta.data$label)

# Create a new column to identify Sham samples
df$is_sham <- grepl("Sham", df$label)

ggplot(data = df) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = is_sham),
             size = 0.5) +
  labs(title = "Sham vs Injury Samples") +
  scale_color_manual(values = c("FALSE" = "lightgrey", "TRUE" = "blue"),
                     labels = c("Injury", "Sham"),
                     name = "Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_line())

################################################################################
# FIGURE S5A-B: UMAP BY LABEL AND BATCH
################################################################################

cat("\n=== Generating Figures S5A-B ===\n")

# Figure S5A
DimPlot(SchwannCells, reduction = "umap", label = FALSE, 
        pt.size = 0.05, group.by = "label")

# Figure S5B
DimPlot(SchwannCells, reduction = "umap", label = FALSE, 
        pt.size = 0.05, group.by = "batch.orig", ncol = 1)

################################################################################
# FIGURE 2C (SUPPLEMENTAL FIGURE 2D): SPLIT BY LABEL
################################################################################

cat("\n=== Generating Figure 2C (Supplemental Figure 2D) ===\n")

# Create dataframe with UMAP coordinates and condition/label
df <- data.frame(UMAP1 = SchwannCells@reductions$umap@cell.embeddings[, 1],
                 UMAP2 = SchwannCells@reductions$umap@cell.embeddings[, 2],
                 label = SchwannCells@meta.data$label)

# Define the order with the correct names
label_order <- c(
  "Sham-for-2hpi-WT",
  "2hpi-WT",
  "Sham-for-1dpi-WT",
  "1dpi-WT",
  "Sham-for-1dpi-Sarm1-KO",
  "1dpi-Sarm1-KO"
)

# Function to create plot for a specific label with conditional coloring
create_label_plot <- function(data, label_name) {
  # Create a new color column where "Naive" or "Sham" labels are blue, others are red
  plot_data <- data %>%
    mutate(highlight_color = ifelse(label == label_name, 
                                    ifelse(grepl("Naive|Sham", label_name), "blue", "red"), 
                                    "lightgrey"))
  
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

# Filter the order to only include labels that exist in the data
valid_labels <- label_order[label_order %in% unique(df$label)]

# Create plots in the specified order
label_plots <- lapply(valid_labels, function(label) {
  create_label_plot(df, label)
})

# Figure 2C
grid.arrange(grobs = label_plots, ncol = 2)

################################################################################
# FIND SCHWANN CELL CLUSTER MARKERS
################################################################################

cat("\n=== Finding Schwann Cell Cluster Markers ===\n")

DefaultAssay(SchwannCells) <- "RNA"
Idents(SchwannCells) <- "seurat_clusters"
SchwannCells[["RNA"]] <- JoinLayers(SchwannCells[["RNA"]])

SchwannCell.markers <- FindAllMarkers(SchwannCells, 
                                      only.pos = TRUE, 
                                      min.pct = 0.25,
                                      logfc.threshold = 0.25)

cat("Total markers found:", nrow(SchwannCell.markers), "\n")

# Save marker genes table (Table 2)
write.csv(SchwannCell.markers, file = markers_output_path, row.names = FALSE)
cat("Schwann Cell markers saved to:", markers_output_path, "\n")

################################################################################
# FIGURE 2D: HEATMAP OF TOP MARKERS
################################################################################

cat("\n=== Generating Figure 2D ===\n")

# Get top 10 markers per cluster
top10_SCs <- SchwannCell.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

# Set RNA assay and get average expression
DefaultAssay(SchwannCells) <- "RNA"
cluster.averages <- AverageExpression(SchwannCells, 
                                      group.by = c("seurat_clusters"),
                                      assays = "RNA", 
                                      features = top10_SCs$gene)

# Convert the average expression data to a matrix
expr_mat_SCs <- as.matrix(cluster.averages$RNA)

# Order genes by cluster
ordered_genes_SCs <- top10_SCs %>%
  arrange(cluster) %>%
  pull(gene)

# Reorder the matrix rows according to the ordered genes
expr_mat_SCs <- expr_mat_SCs[ordered_genes_SCs, ]

# Generate heatmap
pheatmap(expr_mat_SCs, 
         border_color = "NA", 
         scale = "row", 
         color = cm.colors(256), 
         show_rownames = TRUE,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         main = "Cluster−Average Expression of Top Genes for Each Schwann Cell Cluster",
         fontsize = 7)

################################################################################
# FIGURE 2G: FEATURE AND VIOLIN PLOTS
################################################################################

cat("\n=== Generating Figure 2G ===\n")

# Figure 2G - part 1
FeaturePlot(SchwannCells, 
            features = c("Camk1d", "Snhg11", "Mpz", "Ptgds", 
                         "Lsamp", "Stard13", "Ncam1", "Mbp"),
            combine = TRUE, pt.size = 0.1, ncol = 4)

# Figure 2G - part 2
VlnPlot(SchwannCells, 
        features = c("Camk1d", "Snhg11", "Mpz", "Ptgds", 
                     "Lsamp", "Stard13", "Ncam1", "Mbp"), 
        ncol = 4)

################################################################################
# FIGURE 2E-F: STACKED BAR PLOTS
################################################################################

cat("\n=== Generating Figures 2E-F ===\n")

# Filter the data to include only "Sham" in Hashtags
filtered_data_shams <- SchwannCells@meta.data %>%
  filter(grepl("Sham", Hashtags))

# Filter the data to include only injury samples
filtered_data_injury <- SchwannCells@meta.data %>%
  filter(!grepl("Sham", Hashtags))

# Create color palettes
color_palette_paired_sham <- brewer.pal(nlevels(factor(filtered_data_shams$Hashtags)), 
                                        "Paired")
color_palette_paired_injury <- brewer.pal(nlevels(factor(filtered_data_injury$Hashtags)), 
                                          "Paired")

# Define label order
label_order <- c(
  "Sham-for-2hpi-WT",
  "Sham-for-1dpi-WT",
  "Sham-for-1dpi-Sarm1-KO",
  "2hpi-WT",
  "1dpi-WT",
  "1dpi-Sarm1-KO"
)

# Convert label to factor with specific order
filtered_data_shams$label <- factor(filtered_data_shams$label, levels = label_order)
filtered_data_injury$label <- factor(filtered_data_injury$label, levels = label_order)

# Figure 2E
ggplot(filtered_data_shams, aes(x = label, fill = seurat_clusters)) + 
  geom_bar(position = "fill") +
  theme(panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Figure 2F
ggplot(filtered_data_injury, aes(x = label, fill = seurat_clusters)) + 
  geom_bar(position = "fill") +
  theme(panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

################################################################################
# PSEUDOTIME ANALYSIS - WILD TYPE
################################################################################

cat("\n=== Performing Pseudotime Analysis (WT) ===\n")

Idents(SchwannCells) <- "seurat_clusters"
SchwannCells <- JoinLayers(SchwannCells)

# Create logical vector for wild type samples
wild_type_logical <- SchwannCells$label %in% c("Sham-for-2hpi-WT", "Sham-for-1dpi-WT", 
                                               "1dpi-WT", "2hpi-WT")

# Subset the Seurat object to only wild type cells
SchwannCells_WT <- SchwannCells[, wild_type_logical]

# Convert the wild type subset to monocle3 cell_data_set
cds <- as.cell_data_set(SchwannCells_WT)

# Prepare gene metadata
fData(cds)$gene_short_name <- rownames(fData(cds))

# Assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

# Assign cluster information
list.cluster <- SchwannCells_WT@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

# Assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- 
  SchwannCells_WT@reductions$umap@cell.embeddings

# Plot before trajectory
cluster.before.traj <- plot_cells(cds, color_cells_by = "cluster", 
                                  label_groups_by_cluster = FALSE, 
                                  group_label_size = 5) + 
  theme(legend.position = "right")
print(cluster.before.traj)

# Learn Trajectory
cds <- learn_graph(cds, use_partition = FALSE)

# Define earliest principal node (starting point are sham cells)
get_earliest_principal_node <- function(cds, time_bin = c("Sham-for-2hpi-WT",
                                                          "Sham-for-1dpi-WT")) {
  cell_ids <- which(colData(cds)[, "label"] == time_bin)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(
    which.max(table(closest_vertex[cell_ids, ]))))]
  
  root_pr_nodes
}

cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds))

################################################################################
# FIGURE 2H: PSEUDOTIME TRAJECTORY (WT)
################################################################################

cat("\n=== Generating Figure 2H ===\n")

# Fully visible dots
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           trajectory_graph_color = "black",
           label_leaves = FALSE,
           label_branch_points = TRUE,
           graph_label_size = 5,
           trajectory_graph_segment_size = 1,
           group_label_size = 1,
           cell_size = 0.5,
           alpha = 1)

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

# Figure 2H - boxplot
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), 
                        fill = seurat_clusters)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  theme_ridges() +
  labs(x = 'Pseudotime', y = "") +
  theme(axis.text.x = element_text(hjust = 0.5), 
        axis.title.x = element_text(margin = margin(b = 5)))

################################################################################
# PSEUDOTIME ANALYSIS - SARM1-KO
################################################################################

cat("\n=== Performing Pseudotime Analysis (Sarm1-KO) ===\n")

SchwannCells <- JoinLayers(SchwannCells)

# Create logical vector for Sarm1-KO samples
Sarm1KO_logical <- SchwannCells$label %in% c("Sham-for-1dpi-Sarm1-KO", "1dpi-Sarm1-KO")

# Subset the Seurat object to only Sarm1-KO cells
SchwannCells_Sarm1KO <- SchwannCells[, Sarm1KO_logical]

# Convert to monocle3 cell_data_set
cds <- as.cell_data_set(SchwannCells_Sarm1KO)

# Prepare gene metadata
fData(cds)$gene_short_name <- rownames(fData(cds))

# Assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

# Assign cluster information
list.cluster <- SchwannCells_Sarm1KO@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

# Assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- 
  SchwannCells_Sarm1KO@reductions$umap@cell.embeddings

# Plot before trajectory
cluster.before.traj <- plot_cells(cds, color_cells_by = "cluster", 
                                  label_groups_by_cluster = FALSE, 
                                  group_label_size = 5) + 
  theme(legend.position = "right")
print(cluster.before.traj)

# Learn Trajectory
cds <- learn_graph(cds, use_partition = FALSE)

# Setting up Sham cells as starting timepoint
get_earliest_principal_node <- function(cds, time_bin = c("Sham-for-1dpi-Sarm1-KO")) {
  cell_ids <- which(colData(cds)[, "label"] == time_bin)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(
    which.max(table(closest_vertex[cell_ids, ]))))]
  
  root_pr_nodes
}

cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds))

################################################################################
# FIGURE 2I: PSEUDOTIME TRAJECTORY (SARM1-KO)
################################################################################

cat("\n=== Generating Figure 2I ===\n")

# Figure 2I - part 1
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           trajectory_graph_color = "black",
           label_leaves = FALSE,
           label_branch_points = TRUE,
           graph_label_size = 5,
           trajectory_graph_segment_size = 1,
           group_label_size = 1,
           cell_size = 0.5,
           alpha = 1)

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

# Figure 2I - part 2
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), 
                        fill = seurat_clusters)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  theme_ridges() +
  labs(x = 'Pseudotime', y = "") +
  theme(axis.text.x = element_text(hjust = 0.5), 
        axis.title.x = element_text(margin = margin(b = 5)))

################################################################################
# MITOCHONDRIAL GENES ANALYSIS
################################################################################

cat("\n=== Mitochondrial Genes Analysis ===\n")
cat("NOTE: Run this after running 'Figure 3' file to define oxphos_genes\n")

# Supplementary Figure 3D, 4F, 3E
# NOTE: These sections require oxphos_genes variable from "Figure 3" analysis
# Uncomment and run after oxphos_genes is defined

# DefaultAssay(SchwannCells) <- "RNA"
# Idents(SchwannCells) <- "label"
# unique_genes_combined <- oxphos_genes
# 
# # Order clusters
# SchwannCells$seurat_clusters <- factor(SchwannCells$seurat_clusters, 
#                                        levels = c("0", "1", "2", "3", '4', '5', '6', '7'))
# 
# # Supplemental Figure 3D
# DotPlot(SchwannCells, features = unique_genes_combined, group.by = "label") +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
#     axis.text.y = element_text(size = 10)
#   ) +
#   scale_x_discrete(label = function(x) gsub(",", "", x)) +
#   theme(plot.margin = margin(b = 50))
# 
# # Supplemental Figure 4F
# DotPlot(SchwannCells, features = unique_genes_combined, group.by = "seurat_clusters") +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
#     axis.text.y = element_text(size = 10)
#   ) +
#   scale_x_discrete(label = function(x) gsub(",", "", x)) +
#   theme(plot.margin = margin(b = 50))
# 
# # Supplemental Figure 3E
# # Filter to only include 1dpi-Sarm1-KO and 1dpi-WT
# SchwannCells_subset <- subset(SchwannCells, label %in% c("1dpi-Sarm1-KO", "1dpi-WT"))
# 
# # Create the combined variable for the subset
# SchwannCells_subset$cluster_label <- paste0("Cluster_", 
#                                             SchwannCells_subset$seurat_clusters, 
#                                             "_", 
#                                             SchwannCells_subset$label)
# 
# # Create the dot plot
# DotPlot(SchwannCells_subset, 
#         features = unique_genes_combined, 
#         group.by = "cluster_label") +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
#     axis.text.y = element_text(size = 10)
#   ) +
#   scale_y_discrete(limits = rev) +
#   theme(plot.margin = margin(b = 50))

cat("\n=== Schwann Cell Analysis Complete ===\n")
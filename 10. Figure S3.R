################################################################################
# Sham Schwann Cell Analysis and Pseudobulk Differential Expression
# 
# Description: Analysis of sham Schwann cells including reclustering,
# marker identification, and pseudobulk differential expression analysis
# comparing WT vs Sarm1-KO using Libra
#
# Requirements: 
# - Schwann Cells Seurat object (SchwannCells.rds)
# - R packages: Seurat, dplyr, tidyverse, ggplot2, pheatmap, RColorBrewer,
#   SeuratObject, edgeR, limma, Libra, plyr
# - Note: Libra may need custom library path
################################################################################

# Load required libraries
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(SeuratObject)
library(edgeR)
library(limma)
library(Libra, lib.loc="~/R/x86_64-pc-linux-gnu-library/4.4")

################################################################################
# USER-DEFINED PATHS - MODIFY THESE FOR YOUR SYSTEM
################################################################################

# Path to Schwann Cells Seurat object
schwann_cells_path <- "Data_Files/Files/SchwannCells.rds"

# Path to save marker genes table
sham_markers_output_path <- "Data_Files/Files/Table3.csv"

# Path to save pseudobulk results
pseudobulk_output_path <- "Data_Files/Files/Table4.csv"

################################################################################
# LOAD DATA
################################################################################

cat("=== Loading Schwann Cells Object ===\n")

# SchwannCells <- readRDS(schwann_cells_path)

################################################################################
# FIGURE S3B: IDENTIFY SHAM SCHWANN CELLS
################################################################################

cat("\n=== Generating Figure S3B ===\n")

# Create a data frame with the necessary information to identify injured SCs
df <- data.frame(UMAP1 = SchwannCells@reductions$umap@cell.embeddings[, 1],
                 UMAP2 = SchwannCells@reductions$umap@cell.embeddings[, 2],
                 label = SchwannCells@meta.data$label)

# Create a new column to identify Sham samples
cat("Unique labels:\n")
print(unique(df$label))

df$is_injury <- grepl("Sham", df$label)

# Figure S3B
ggplot(data = df) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = is_injury),
             size = 0.5) +
  labs(title = "Sham Schwann Cells") +
  scale_color_manual(values = c("FALSE" = "lightgrey", "TRUE" = "blue"),
                     labels = c("Injury", "Sham"),
                     name = "Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_line())

################################################################################
# SUBSET SHAM SCHWANN CELLS
################################################################################

cat("\n=== Subsetting Sham Schwann Cells ===\n")

DefaultAssay(SchwannCells) <- "RNA"
Idents(SchwannCells) <- "label"

cat("Total Schwann cells:", length(WhichCells(SchwannCells)), "\n")
cat("Unique labels:\n")
print(unique(SchwannCells@meta.data$label))

# Check if layers are split by batch
cat("Current layers:\n")
print(Layers(SchwannCells))

# Run only if the layers are not split by batch already
SchwannCells[["RNA"]] <- split(SchwannCells[["RNA"]], f = SchwannCells$batch.orig)

cat("Layers after splitting:\n")
print(Layers(SchwannCells))

# Subset only sham samples
ShamSCs <- subset(SchwannCells, 
                  cells = WhichCells(SchwannCells, 
                                     idents = c("Sham-for-1dpi-WT", 
                                                "Sham-for-1dpi-Sarm1-KO", 
                                                "Sham-for-2hpi-WT")))

cat("Sham Schwann cells:", length(WhichCells(ShamSCs)), "\n")
cat("Unique batches:\n")
print(unique(ShamSCs@meta.data$batch.orig))
cat("Unique labels:\n")
print(unique(ShamSCs@meta.data$label))

################################################################################
# RECLUSTER SHAM SCHWANN CELLS
################################################################################

cat("\n=== Reclustering Sham Schwann Cells ===\n")

DefaultAssay(ShamSCs) <- "RNA"

ShamSCs <- NormalizeData(ShamSCs)
ShamSCs <- FindVariableFeatures(ShamSCs)
ShamSCs <- ScaleData(ShamSCs)
ShamSCs <- RunPCA(ShamSCs)

# Visualize elbow plot
ElbowPlot(ShamSCs)

# Join layers for clustering
ShamSCs[["RNA"]] <- JoinLayers(ShamSCs[["RNA"]])

# Clustering
x <- 5
ShamSCs <- FindNeighbors(ShamSCs, reduction = "pca", dims = 1:x)
ShamSCs <- FindClusters(ShamSCs, resolution = 0.2)
ShamSCs <- RunUMAP(ShamSCs, dims = 1:x, reduction = "pca")

################################################################################
# FIGURE S3C: UMAP BY CLUSTER
################################################################################

cat("\n=== Generating Figure S3C ===\n")

DimPlot(ShamSCs, reduction = "umap", group.by = c("seurat_clusters"))

################################################################################
# FIGURE S3F: UMAP BY LABEL
################################################################################

cat("\n=== Generating Figure S3E ===\n")

DimPlot(ShamSCs, reduction = "umap", group.by = c("label"), 
        cols = c("#D81B60", "#1E88E5", "#004D40"))

################################################################################
# FIGURE S3F-G: COMPOSITION PLOTS
################################################################################

cat("\n=== Generating Figures S3F-G ===\n")

# Set label order
cat("Unique labels in ShamSCs:\n")
print(unique(ShamSCs@meta.data$label))

ShamSCs$label <- factor(ShamSCs$label, 
                        levels = c("Sham-for-2hpi-WT", "Sham-for-1dpi-WT",
                                   "Sham-for-1dpi-Sarm1-KO", "Sham-for-3dpi-WT"))

filtered_data <- ShamSCs@meta.data

# Figure S3F: Cluster composition by label
ggplot(filtered_data, aes(x = seurat_clusters, fill = label)) + 
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("#D81B60", "#1E88E5", "#004D40")) +
  theme(panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank())

# Figure S3G: Label composition by cluster
ggplot(filtered_data, aes(x = label, fill = seurat_clusters)) + 
  geom_bar(position = "fill") +
  theme(panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank())

################################################################################
# FIND SHAM SCHWANN CELL CLUSTER MARKERS
################################################################################

cat("\n=== Finding Sham Schwann Cell Cluster Markers ===\n")

DefaultAssay(ShamSCs) <- "RNA"
ShamSCs[["RNA"]] <- JoinLayers(ShamSCs[["RNA"]])

Idents(ShamSCs) <- "seurat_clusters"

SC_Sham.markers <- FindAllMarkers(ShamSCs, 
                                  only.pos = TRUE, 
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25)

cat("Total markers found:", nrow(SC_Sham.markers), "\n")

# Save marker genes table (Table 3)
write.csv(SC_Sham.markers, file = sham_markers_output_path, row.names = FALSE)
cat("Sham SC markers saved to:", sham_markers_output_path, "\n")

################################################################################
# FIGURE 3SD: VIOLIN PLOTS
################################################################################

cat("\n=== Generating Figure 3SD ===\n")

VlnPlot(ShamSCs, features = c("Mbp", "Mpz", "Prx", "Ncam1"), ncol = 2)

################################################################################
# FIGURE 3SE: HEATMAP OF TOP MARKERS
################################################################################

cat("\n=== Generating Figure 3SE ===\n")

# Get top 5 markers per cluster
top5_Sham_SCs <- SC_Sham.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)

# Set RNA assay and get average expression
DefaultAssay(ShamSCs) <- "RNA"
cluster.averages <- AverageExpression(ShamSCs, 
                                      group.by = "seurat_clusters",
                                      assays = "RNA", 
                                      features = top5_Sham_SCs$gene)

# Convert the average expression data to a matrix
expr_mat_sham_SCs <- as.matrix(cluster.averages$RNA)

# Order genes by cluster
ordered_genes_sham_SCs <- top5_Sham_SCs %>%
  arrange(cluster) %>%
  pull(gene)

# Reorder the matrix rows according to the ordered genes
expr_mat_sham_SCs <- expr_mat_sham_SCs[ordered_genes_sham_SCs, ]

# Generate heatmap
pheatmap(expr_mat_sham_SCs, 
         border_color = "NA", 
         scale = "row", 
         color = cm.colors(256), 
         show_rownames = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Sham Cluster−Average Expression of Top Genes for Each Schwann Cell Cluster")

################################################################################
# PSEUDOBULK DIFFERENTIAL EXPRESSION ANALYSIS
################################################################################

cat("\n=== Preparing for Pseudobulk Analysis ===\n")

DefaultAssay(SchwannCells) <- "RNA"

cat("Unique replicates:\n")
print(unique(SchwannCells@meta.data$replicate))

# Leave only 1dpi data
Idents(SchwannCells) <- "orig.ident"
cat("Unique orig.ident:\n")
print(unique(SchwannCells@meta.data$orig.ident))

SC.integrated_Pseudobulk <- subset(SchwannCells, 
                                   cells = WhichCells(SchwannCells, 
                                                      idents = c("Sham.singletB1", "Sham.singletB2",
                                                                 "Sham.singletB3", "Sham.singletB4")))

cat("Pseudobulk subset orig.ident:\n")
print(unique(SC.integrated_Pseudobulk@meta.data$orig.ident))

################################################################################
# ADD CELL TYPE METADATA
################################################################################

cat("\n=== Adding Cell Type Metadata ===\n")

# Add a new column "cell_type" that assigns the name of the cell type
Idents(SC.integrated_Pseudobulk) <- "seurat_clusters"
cat("Unique clusters:\n")
print(unique(SC.integrated_Pseudobulk@meta.data$seurat_clusters))

current.cluster.idsB <- c('0', '1', '2', '3', '4', "5", '6', '7')
new.cluster.idsB <- c('SC0', 'SC0', "SC0", "SC0", "SC0", 'SC0', 'SC0', 'SC0')
cell_type <- plyr::mapvalues(x = SC.integrated_Pseudobulk@meta.data$seurat_clusters, 
                             from = current.cluster.idsB, 
                             to = new.cluster.idsB)

names(cell_type) <- colnames(x = SC.integrated_Pseudobulk)
SC.integrated_Pseudobulk <- AddMetaData(object = SC.integrated_Pseudobulk,
                                        metadata = cell_type,
                                        col.name = 'cell_type')

# Check that all columns needed to run Libra are present
cat("Unique replicates:\n")
print(unique(SC.integrated_Pseudobulk@meta.data$replicate))
cat("Unique labels:\n")
print(unique(SC.integrated_Pseudobulk@meta.data$label))
cat("Unique cell types:\n")
print(unique(SC.integrated_Pseudobulk@meta.data$cell_type))

################################################################################
# RUN PSEUDOBULK ANALYSIS WITH LIBRA
################################################################################

cat("\n=== Running Pseudobulk Differential Expression Analysis ===\n")

SC.integrated_Pseudobulk <- JoinLayers(SC.integrated_Pseudobulk)

cat("Unique labels in pseudobulk object:\n")
print(unique(SC.integrated_Pseudobulk$label))

# Create the comparisons dataframe
comparisons <- data.frame(
  group1 = c("Sham-for-1dpi-WT"),
  group2 = c("Sham-for-1dpi-Sarm1-KO")
)

cat("Comparisons to perform:\n")
print(comparisons)

# Run differential expression analysis
results <- list()
for (i in 1:nrow(comparisons)) {
  Idents(SC.integrated_Pseudobulk) = SC.integrated_Pseudobulk$label
  group1 = comparisons[i, ]$group1
  group2 = comparisons[i, ]$group2
  
  cat("Running comparison:", group1, "vs", group2, "\n")
  
  sc0 = subset(SC.integrated_Pseudobulk, idents = c(group1, group2))
  
  de = run_de(sc0, n_threads = 8) %>%
    mutate(group1 = group1,
           group2 = group2)
  
  results[[i]] = de  # Store the result in the list
}

# Combine all results into a single table
SchwannCellGenes_Shams = do.call(bind_rows, results)

# Save the SchwannCellGenes in the global environment
assign("SchwannCellGenes_Shams", SchwannCellGenes_Shams, envir = .GlobalEnv)

cat("Total genes analyzed:", nrow(SchwannCellGenes_Shams), "\n")

# Save pseudobulk results (Table 4)
write.csv(SchwannCellGenes_Shams, file = pseudobulk_output_path, row.names = FALSE)
cat("Pseudobulk results saved to:", pseudobulk_output_path, "\n")

################################################################################
# FIGURE 3SI: SIGNIFICANT VS NON-SIGNIFICANT GENES
################################################################################

cat("\n=== Generating Figure 3SI ===\n")

# Create gene_summary using base R
sig_count <- sum(SchwannCellGenes_Shams$p_val_adj < 0.05, na.rm = TRUE)
nonsig_count <- sum(SchwannCellGenes_Shams$p_val_adj >= 0.05, na.rm = TRUE)

gene_summary <- data.frame(
  significance = c("Significant", "Non-significant"),
  n = c(sig_count, nonsig_count)
)

cat("Significant genes (padj < 0.05):", sig_count, "\n")
cat("Non-significant genes:", nonsig_count, "\n")

# Barplot to show significant vs non-significant genes
gene_summary %>%
  ggplot(aes(x = significance, y = n, fill = significance)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = n), vjust = -0.5, size = 5, fontface = "bold") +
  scale_x_discrete(limits = c("Significant", "Non-significant")) +
  scale_fill_manual(values = c("Significant" = "black", "Non-significant" = "#FF9999")) +
  labs(x = "Gene Significance", 
       y = "Number of Genes", 
       title = "WT vs Sarm1KO Sham: Significant vs Non-Significant Genes",
       subtitle = "Based on adjusted p-value < 0.05") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"), 
        axis.line = element_line(color = "black"), 
        plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, vjust = 0.5, size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  ylim(0, max(gene_summary$n) * 1.1)

cat("\n=== Sham Schwann Cell and Pseudobulk Analysis Complete ===\n")
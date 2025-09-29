#load the libraries
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(SeuratObject)
library(gridExtra)

#Read this in if have already SchwannCells.rds saved somewhere
#SchwannCells <- readRDS("Data_Files/Files/SchwannCells.rds")
Layers(SchwannCells)
# Create a data frame with the necessary information to identify injured SCs
df <- data.frame(UMAP1 = SchwannCells@reductions$umap@cell.embeddings[, 1],
                 UMAP2 = SchwannCells@reductions$umap@cell.embeddings[, 2],
                 label = SchwannCells@meta.data$label)

# Create a new column to identify Injured samples
unique(df$label)
df$is_injury <- !grepl("Sham|Naive", df$label)

#Figure S4B
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




#################Subset only Schwann Cell Clusters
#SchwannCells <- readRDS("Data_Files/Files/SchwannCells.rds")

DefaultAssay(SchwannCells) <- "RNA"
Idents(SchwannCells) <- "label"
length(WhichCells(SchwannCells))
unique(SchwannCells@meta.data$label)

Layers(SchwannCells) #check if layers are split by batch
SchwannCells[["RNA"]] <- split(SchwannCells[["RNA"]], f = SchwannCells$batch.orig) # Run only if the layers are not split by batch already
Layers(SchwannCells) #check if layers are split by batch


InjuredSCs<-subset(SchwannCells,cells=WhichCells(SchwannCells,idents=c("1dpi-Sarm1-KO","1dpi-WT","2hpi-WT")))
length(WhichCells(InjuredSCs))#7353
unique(InjuredSCs@meta.data$batch.orig)
unique(InjuredSCs@meta.data$label)

# Use your existing InjuredSCs object
DefaultAssay(InjuredSCs) <- "RNA"
InjuredSCs <- NormalizeData(InjuredSCs)
InjuredSCs <- FindVariableFeatures(InjuredSCs)
InjuredSCs <- ScaleData(InjuredSCs)
InjuredSCs <- RunPCA(InjuredSCs)



# re-join layers after integration
#InjuredSCs[["RNA"]] <- JoinLayers(InjuredSCs[["RNA"]])
ElbowPlot(InjuredSCs)
#number of principle components
x<-8

# Continue with downstream analysis
InjuredSCs <- FindNeighbors(InjuredSCs, reduction = "pca", dims = 1:x)
InjuredSCs <- FindClusters(InjuredSCs, resolution = 0.25)
InjuredSCs <- RunUMAP(InjuredSCs, reduction = "pca", dims = 1:x, seed.use = 42)
# Figure S4C
DimPlot(InjuredSCs, reduction = "umap", group.by = c("seurat_clusters"))


# Create the plot with the filtered data
unique(InjuredSCs@meta.data$label)
InjuredSCs$label <- factor(InjuredSCs$label, 
                           levels = c("2hpi-WT", "1dpi-WT",
                                      "1dpi-Sarm1-KO"))
unique(InjuredSCs$Hashtags)
InjuredSCs$Hashtags <- factor(InjuredSCs$Hashtags, 
                              levels = c("2hpi-WT-Female", "2hpi-WT-Male",
                                         "1dpi-WT-Female","1dpi-WT-Male",
                                         "1dpi-Sarm1-KO-Female","1dpi-Sarm1-KO-Male"))
filtered_data <- InjuredSCs@meta.data
#Figure S4F
ggplot(filtered_data, aes(x = label, fill = seurat_clusters)) + 
  geom_bar(position = "fill") +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank())

#ggplot(filtered_data, aes(x = Hashtags, fill = seurat_clusters)) + 
#  geom_bar(position = "fill") +
#  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank())

#Figure S4G
ggplot(filtered_data, aes(x = seurat_clusters, fill = label)) + 
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("#004D40", "#1E88E5", "#D81B60")) +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank())

#ggplot(filtered_data, aes(x = seurat_clusters, fill = Hashtags)) + 
#  geom_bar(position = "fill") +
#  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank())


#Markers and feature plots   
DefaultAssay(InjuredSCs) <- "RNA"
InjuredSCs[["RNA"]] <- JoinLayers(InjuredSCs[["RNA"]])

Idents(InjuredSCs) <- "seurat_clusters"
SC_Injury.integrated.markers <- FindAllMarkers(InjuredSCs, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25)
#Table 5
#write.csv(SC_Injury.markers, file="Data_Files/Files/Table5")


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
expr_mat_injury_SCs <- expr_mat_injury_SCs[ordered_genes_injured_SCs,]


#Figure S4D
pheatmap(expr_mat_injury_SCs, 
         border_color = "NA", 
         scale = "row", 
         color = cm.colors(256), 
         show_rownames = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Injury Cluster−Average Expression of Top Genes for Each Schwann Cell Cluster")


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
#Figure S4E
grid.arrange(grobs = label_plots, ncol = 1)  # Adjust ncol as needed

#Figure S4H - Part 1
VlnPlot(InjuredSCs, features=c("Plp1","Fth1","Ptgds"), group.by="seurat_clusters")
#Figure S4H - Part 2
VlnPlot(InjuredSCs, features=c("Plp1","Fth1","Ptgds"), group.by="label")

#Figure S4I
FeaturePlot(InjuredSCs, features=c("Plp1","Fth1","Ptgds"),ncol=3)





#####Run this after running file Figure 3


DefaultAssay(InjuredSCs) <- "RNA"
Idents(InjuredSCs) <- "label"
unique_genes_combined <- oxphos_genes

#order clusters
unique(InjuredSCs$seurat_clusters)
InjuredSCs$seurat_clusters <- factor(InjuredSCs$seurat_clusters, 
                                     levels = c("0","1","2","3",'4','5','6','7','8'))

#Figure S4H

DotPlot(InjuredSCs, features=unique_genes_combined, group.by="label")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),  # Adjust text size here
    axis.text.y = element_text(size = 10)
  ) +
  scale_x_discrete(label = function(x) gsub(",", "", x)) +
  theme(plot.margin = margin(b = 50))
#Figure S4I

DotPlot(InjuredSCs, features=unique_genes_combined_significant, group.by="seurat_clusters")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),  # Adjust text size here
    axis.text.y = element_text(size = 10)
  ) +
  scale_x_discrete(label = function(x) gsub(",", "", x)) +
  theme(plot.margin = margin(b = 50))

# Supplemental figure 4G
DotPlot(InjuredSCs, 
        features = unique_genes_combined_significant, 
        split.by = "label",
        group.by = "seurat_clusters",  # Add this parameter
        cols = c("blue", "blue", "blue", "blue", "blue")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  scale_x_discrete(label = function(x) gsub(",", "", x)) +
  # Add this to control gene ordering on y-axis
  scale_y_discrete(limits = rev) +
  theme(plot.margin = margin(b = 50))



# Supplemental figure 4F

PASC_genes<-c("Plp1","Fth1","Ptgds")
DotPlot(InjuredSCs, features=PASC_genes, group.by="seurat_clusters")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),  # Adjust text size here
    axis.text.y = element_text(size = 10)
  ) +
  scale_x_discrete(label = function(x) gsub(",", "", x)) +
  theme(plot.margin = margin(b = 50))

VlnPlot(InjuredSCs, features=PASC_genes, group.by = "seurat_clusters",ncol=1)
FeaturePlot(InjuredSCs, features=PASC_genes, ncol=3)


# Supplemental figure 4G
DotPlot(InjuredSCs, features=PASC_genes, group.by="label")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),  # Adjust text size here
    axis.text.y = element_text(size = 10)
  ) +
  scale_x_discrete(label = function(x) gsub(",", "", x)) +
  theme(plot.margin = margin(b = 50))

VlnPlot(InjuredSCs, features=PASC_genes, split.by = "label",ncol=1)
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(future)
library(RColorBrewer)
library(pheatmap)
library(SeuratObject)
library(gridExtra)
library(monocle3, lib.loc = "/sfs/gpfs/tardis/home/es3jk/R/library")
library(SeuratWrappers, lib.loc = "/sfs/gpfs/tardis/home/es3jk/R/library")
library(ggridges)
library(SingleCellExperiment)

#Read in all the cell data if already saved as a seurat object
#Cond.combineds_integrated <- readRDS("Data_Files/Files/Cond.combineds_integrated.rds")

#################Subsetting only Schwann Cell Clusters
DefaultAssay(Cond.combineds_integrated) <- "RNA"
Idents(Cond.combineds_integrated) <- "seurat_clusters"

#First figure out which clusters are SchwannCells
#DotPlot(Cond.combineds_integrated, features=c("Mbp","Ncam1", "Cadm2"))


####highlighting only clusters 2, 5,6,7,8 - Schwann Cells
df_SCs <- data.frame(UMAP1 = Cond.combineds_integrated@reductions$umap@cell.embeddings[, 1],
                     UMAP2 = Cond.combineds_integrated@reductions$umap@cell.embeddings[, 2],
                     seurat_clusters = Cond.combineds_integrated@meta.data$seurat_clusters)

# Create a new column for colors based on seurat_clusters
SchwannCells <- df_SCs %>%
  mutate(color = case_when(
    seurat_clusters == "2" ~ "black",
    seurat_clusters == "5" ~ "black",
    seurat_clusters == "8" ~ "black",
    seurat_clusters == "6" ~ "black",
    seurat_clusters == "7" ~ "black",
    TRUE ~ "lightgrey"
  ))



# Figure 2A, S3A,S4A
ggplot(data = SchwannCells) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = color), size = 1) +
  labs(title = "Schwann Cells") +
  scale_color_identity(guide = "none") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_line())




#Subsetting Schwann Cells that are in clusters 2, 5,6,7,8
Idents(Cond.combineds_integrated) <- "seurat_clusters"


SchwannCells<-subset(Cond.combineds_integrated,cells=WhichCells(Cond.combineds_integrated,idents=c('2',"5","6",'7','8')))
DotPlot(SchwannCells, features=c("Mbp","Mpz","Ncam1", "Cadm2")) #double check that subsetting was done correctly by looking at SC markers
length(WhichCells(SchwannCells))#14879
length(WhichCells(SchwannCells, expression = condition == "Sham"))#4724
length(WhichCells(SchwannCells, expression = condition == "1dpi"))#7049
length(WhichCells(SchwannCells, expression = condition == "2hpi"))#304
Idents(SchwannCells) <- "condition"
SchwannCells<-subset(SchwannCells,cells=WhichCells(SchwannCells,idents=c("Sham",'1dpi',"2hpi")))
unique(SchwannCells$condition)
length(WhichCells(SchwannCells))#12077


DefaultAssay(SchwannCells) <- "RNA"

##See how many batches 
unique(SchwannCells@meta.data$batch.orig)

# Create a unified layer
DefaultAssay(SchwannCells) <- "RNA"
SchwannCells <- JoinLayers(SchwannCells)

# Verify the names of each 10x lane
unique(SchwannCells@meta.data$batch.orig)

######
DefaultAssay(SchwannCells) <- "RNA"

SchwannCells[["RNA"]] <- split(SchwannCells[["RNA"]], f = SchwannCells$batch.orig)
SchwannCells
SchwannCells <- NormalizeData(SchwannCells)
SchwannCells <- FindVariableFeatures(SchwannCells)
SchwannCells <- ScaleData(SchwannCells)
SchwannCells <- RunPCA(SchwannCells)


ElbowPlot(SchwannCells)

x<-6
SchwannCells <- FindNeighbors(SchwannCells, reduction = "pca", dims = 1:x)
SchwannCells <- FindClusters(SchwannCells, resolution = 0.3) #0.2
SchwannCells <- RunUMAP(SchwannCells, dims = 1:x, reduction = "pca")

Idents(SchwannCells)<- "seurat_clusters"
#Figure 2A
DimPlot(SchwannCells, reduction = "umap", label=TRUE, pt.size = 0.05)
#Save the seurat object so you don't have to rerun the code above
saveRDS(SchwannCells, "Data_Files/Files/SchwannCells.rds")


##Figure 2B
# Create a data frame with the necessary information
df <- data.frame(UMAP1 = SchwannCells@reductions$umap@cell.embeddings[, 1],
                 UMAP2 = SchwannCells@reductions$umap@cell.embeddings[, 2],
                 label = SchwannCells@meta.data$label)
length(WhichCells(SchwannCells))

# Create a new column to identify Sham samples
df$is_sham <- grepl("Sham", df$label)

#Figure 3SB
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

#Figure S5A
DimPlot(SchwannCells, reduction = "umap", label=FALSE, pt.size = 0.05, group.by = "label")

#Figure S5B
DimPlot(SchwannCells, reduction = "umap", label=FALSE, pt.size = 0.05, group.by = "batch.orig", ncol = 1)

#Supplemental Figure 2D
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

# Create a list of plots for all labels in the specified order
# Filter the order to only include labels that exist in the data
valid_labels <- label_order[label_order %in% unique(df$label)]

# Create plots in the specified order
label_plots <- lapply(valid_labels, function(label) {
  create_label_plot(df, label)
})

#Figure 2C
library(gridExtra)
grid.arrange(grobs = label_plots, ncol = 2)



#####
#Markers and feature plots   
DefaultAssay(SchwannCells) <- "RNA"
Idents(SchwannCells) <- "seurat_clusters"
SchwannCells[["RNA"]] <- JoinLayers(SchwannCells[["RNA"]])

SchwannCell.markers <- FindAllMarkers(SchwannCells, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25)
View(SchwannCell.markers)
#Table 2
#write.csv(SchwannCell.markers, file="Data_Files/Files/Table2")


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
expr_mat_SCs <- expr_mat_SCs[ordered_genes_SCs,]


#Figure 2D
pheatmap(expr_mat_SCs, 
         border_color = "NA", 
         scale = "row", 
         color = cm.colors(256), 
         show_rownames = TRUE,
         cluster_rows = FALSE,
         cluster_cols = TRUE,  # Enable column clustering
         main = "Cluster−Average Expression of Top Genes for Each Schwann Cell Cluster",
         fontsize=7)


#Figure 2G - part 1
FeaturePlot(SchwannCells, 
            features = c("Camk1d","Snhg11","Mpz","Ptgds","Lsamp","Stard13","Ncam1","Mbp"),
            combine = TRUE, pt.size = 0.1,ncol=4)
#Figure 2G - part 2
VlnPlot(SchwannCells, 
        features = c("Camk1d","Snhg11","Mpz","Ptgds","Lsamp","Stard13","Ncam1","Mbp"), ncol=4)

###Stacked barplot to visualize composition of each cluster
##Plot Sham and Injury categories
unique(SchwannCells@meta.data$Hashtags)
# Filter the data to include only "Sham" in genotype
filtered_data_shams <- SchwannCells@meta.data %>%
  filter(grepl("Sham", Hashtags))
filtered_data_injury <- SchwannCells@meta.data %>%
  filter(!grepl("Sham", Hashtags))
# Create the color palette friendly to color blindness
color_palette_paired_sham <- brewer.pal(nlevels(factor(filtered_data_shams_and_naive$Hashtags)), "Paired")
color_palette_paired_injury <- brewer.pal(nlevels(factor(filtered_data_injury$Hashtags)), "Paired")


label_order <- c(
  "Sham-for-2hpi-WT",
  "Sham-for-1dpi-WT",
  "Sham-for-1dpi-Sarm1-KO",
  "2hpi-WT",
  "1dpi-WT",
  "1dpi-Sarm1-KO"
)


# Convert label to factor with specific order in both datasets
unique(filtered_data_shams_and_naive$label)
unique(filtered_data_injury$label)
filtered_data_shams_and_naive$label <- factor(filtered_data_shams_and_naive$label, 
                                              levels = label_order)
filtered_data_injury$label <- factor(filtered_data_injury$label, 
                                     levels = label_order)
# Figure2E
ggplot(filtered_data_shams_and_naive, aes(x = label, fill = seurat_clusters)) + 
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



#### Pseudotime analysis
Idents(SchwannCells) <- "seurat_clusters"
SchwannCells <- JoinLayers(SchwannCells)

# Create logical vector for wild type samples
wild_type_logical <- SchwannCells$label %in% c("Sham-for-2hpi-WT", "Sham-for-1dpi-WT", "1dpi-WT","2hpi-WT")

# Subset the Seurat object to only wild type cells
SchwannCells_WT <- SchwannCells[, wild_type_logical]

# Convert the wild type subset to monocle3 cell_data_set
cds <- as.cell_data_set(SchwannCells_WT)

fData(cds)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))
head(fData(cds))
head(counts(cds))

#Assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

#Assign cluster information
list.cluster <- SchwannCells_WT@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster


#Assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- SchwannCells_WT@reductions$umap@cell.embeddings
#Plot
cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
                                 group_label_size = 5) + theme(legend.position = "right")
cluster.before.traj

#Learn Trajectory
cds <- learn_graph(cds, use_partition = F)




#now lets find out what are our early and later time points
#starting point are sham cells
unique(SchwannCells_WT@meta.data$label)
get_earliest_principal_node <- function(cds, time_bin=c("Sham-for-2hpi-WT",
                                                        "Sham-for-2hpi-WT",
                                                        "Sham-for-1dpi-WT")){
  cell_ids <- which(colData(cds)[, "label"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

#fully visible dots
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           trajectory_graph_color = "black",
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=5,
           trajectory_graph_segment_size = 1,
           group_label_size = 1,
           cell_size = 0.5,
           alpha = 1)

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

#Figure 2H
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  theme_ridges() +
  labs(x = 'Pseudotime', y = "") +
  theme(axis.text.x = element_text(hjust = 0.5), 
        axis.title.x = element_text(margin = margin(b = 5)))



#Pseuotime on Sarm1KO Schwann Cells

####Idents(SchwannCells) <- "seurat_clusters"
SchwannCells <- JoinLayers(SchwannCells)
unique(SchwannCells$label)
# Create logical vector for wild type samples
Sarm1KO_logical <- SchwannCells$label %in% c("Sham-for-1dpi-Sarm1-KO", "1dpi-Sarm1-KO")

# Subset the Seurat object to only wild type cells
SchwannCells_Sarm1KO <- SchwannCells[, Sarm1KO_logical]

# Convert the wild type subset to monocle3 cell_data_set
cds <- as.cell_data_set(SchwannCells_Sarm1KO)

fData(cds)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))
head(fData(cds))
head(counts(cds))

#Assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

#Assign cluster information
list.cluster <- SchwannCells_Sarm1KO@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster


#Assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- SchwannCells_Sarm1KO@reductions$umap@cell.embeddings
#Plot
cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
                                 group_label_size = 5) + theme(legend.position = "right")
cluster.before.traj

#Learn Trajectory
cds <- learn_graph(cds, use_partition = F)


#Setting up Sham cells and starting timepoint
unique(SchwannCells_Sarm1KO@meta.data$label)
get_earliest_principal_node <- function(cds, time_bin=c("Sham-for-1dpi-Sarm1-KO")){
  cell_ids <- which(colData(cds)[, "label"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

#Figure 2I - part 1
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           trajectory_graph_color = "black",
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=5,
           trajectory_graph_segment_size = 1,
           group_label_size = 1,
           cell_size = 0.5,
           alpha = 1)

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

#Figure 2I - part 2
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  theme_ridges() +
  labs(x = 'Pseudotime', y = "") +
  theme(axis.text.x = element_text(hjust = 0.5), 
        axis.title.x = element_text(margin = margin(b = 5)))


##########Run this after running file "Figure 3"


### mitochondrial genes
###Supplementary Figure 3D
DefaultAssay(SchwannCells) <- "RNA"
Idents(SchwannCells) <- "label"
unique_genes_combined <- oxphos_genes

#order clusters
unique(SchwannCells$seurat_clusters)
SchwannCells$seurat_clusters <- factor(SchwannCells$seurat_clusters, 
                                       levels = c("0","1","2","3",'4','5','6','7'))

#
DotPlot(SchwannCells, features=unique_genes_combined, group.by="label")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),  # Adjust text size here
    axis.text.y = element_text(size = 10)
  ) +
  scale_x_discrete(label = function(x) gsub(",", "", x)) +
  theme(plot.margin = margin(b = 50))
# Supplemental figure 4F

DotPlot(SchwannCells, features=unique_genes_combined, group.by="seurat_clusters")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),  # Adjust text size here
    axis.text.y = element_text(size = 10)
  ) +
  scale_x_discrete(label = function(x) gsub(",", "", x)) +
  theme(plot.margin = margin(b = 50))

# Supplemental figure 3E

# Filter the Seurat object to only include the labels you want
SchwannCells_subset <- subset(SchwannCells, label %in% c("1dpi-Sarm1-KO", "1dpi-WT"))

# Create the combined variable for the subset
SchwannCells_subset$cluster_label <- paste0("Cluster_", SchwannCells_subset$seurat_clusters, "_", SchwannCells_subset$label)

# Create the dot plot
DotPlot(SchwannCells_subset, 
        features = unique_genes_combined, 
        group.by = "cluster_label") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  scale_y_discrete(limits = rev) +
  theme(plot.margin = margin(b = 50))


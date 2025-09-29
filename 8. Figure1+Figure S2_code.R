#load the libraries
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(SeuratObject)

#read in already saved object (can do this after have already generated the Cond.combined_integrated.rds seurat object)
#Cond.combineds_integrated <- readRDS("Data_Files/Files/Cond.combineds_integrated.rds")


#### Figure 1G - part 1#####
Idents(Cond.combineds_integrated)<- "seurat_clusters"
DimPlot(Cond.combineds_integrated, reduction = "umap", label=TRUE, pt.size = 0.01)


#validate that clustering is not happening on amount of genes detected in each cell
#FeaturePlot(Cond.combineds_integrated, features=c('nFeature_RNA'), pt.size=0.5)
#VlnPlot(Cond.combineds_integrated, features = c("nFeature_RNA"))


# Create a new column for plotting groups, preserving the original condition names
Cond.combineds_integrated$plot_groups <- as.character(Cond.combineds_integrated$condition)

# Reassign both Sham and Naive to the same value to combine them on one plot
Cond.combineds_integrated$plot_groups[Cond.combineds_integrated$plot_groups %in% c("Sham", "Naive")] <- "Sham_and_Naive"

# Convert to factor with specific level order to control plot arrangement
Cond.combineds_integrated$plot_groups <- factor(
  Cond.combineds_integrated$plot_groups,
  levels = c("Sham_and_Naive", "2hpi", "1dpi")
)

#### Figure 1G - part two#####
DimPlot(Cond.combineds_integrated, 
        reduction = "umap", 
        label = FALSE, 
        split.by = "plot_groups",
        ncol = 2)

# Find top significant genes enriched in each cluster
Idents(Cond.combineds_integrated) <- "seurat_clusters"
combined.integrated.markers <- FindAllMarkers(Cond.combineds_integrated, only.pos = TRUE, min.pct = 0.25, 
                                              logfc.threshold = 0.25)
#View(combined.integrated.markers)
#Table 1
#write.csv(combined.integrated.markers, file="Data_Files/Files/Table1")

# Get top 5 markers per cluster
top5 <- combined.integrated.markers %>% 
  filter(p_val_adj < 0.05) %>%          # Filter for significant p-values first
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)         # Then select top 5 by fold change

####Figure 1H
DefaultAssay(Cond.combineds_integrated) <- "RNA"
VlnPlot(Cond.combineds_integrated, 
        features=top5$gene,
        add.noise = FALSE,
        stack=TRUE, 
        flip=TRUE, 
        fill.by="ident") +
  theme(legend.position = "none",
        strip.text.x = element_text(angle = 90, 
                                    hjust = 0, 
                                    vjust = 0.5,
                                    face = "plain",
                                    size = 8))


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
expr_mat <- expr_mat[ordered_genes,]
#Figure S2B
pheatmap(expr_mat, 
         border_color = "NA", 
         scale = "row", 
         color = cm.colors(256), 
         show_rownames = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         angle_col = "45",
         main = "Cluster average expressoin of top 5 expressed markers for each cluster")


###Stacked barplot to visualize composition of each cluster
##Plot Sham and Injury categories
unique(Cond.combineds_integrated@meta.data$Hashtags)
# Filter the data to include only "Sham" in genotype
filtered_data_shams_and_naive <- Cond.combineds_integrated@meta.data %>%
  filter(grepl("Sham|Naive", Hashtags))
filtered_data_injury <- Cond.combineds_integrated@meta.data %>%
  filter(!grepl("Sham|Naive", Hashtags) & !grepl("Naive", Hashtags))
# Create the color palette friendly to color blindness
color_palette_paired_sham <- brewer.pal(nlevels(factor(filtered_data_shams_and_naive$Hashtags)), "Paired")
color_palette_paired_injury <- brewer.pal(nlevels(factor(filtered_data_injury$Hashtags)), "Paired")


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
unique(filtered_data_shams_and_naive$label)
unique(filtered_data_injury$label)
filtered_data_shams_and_naive$label <- factor(filtered_data_shams_and_naive$label, 
                                              levels = label_order)
filtered_data_injury$label <- factor(filtered_data_injury$label, 
                                     levels = label_order)
# Figure 1I
ggplot(filtered_data_shams_and_naive, aes(x = label, fill = seurat_clusters)) + 
  geom_bar(position = "fill") +
  theme(panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Figure 1J
ggplot(filtered_data_injury, aes(x = label, fill = seurat_clusters)) + 
  geom_bar(position = "fill") +
  theme(panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate total count from original dataset
total_count <- nrow(filtered_data_shams_and_naive)

#Figure S2J - part 1
#clusters 0-6 for sham cells
cluster0_6 <- filtered_data_shams_and_naive[filtered_data_shams_and_naive$seurat_clusters %in% 0:6, ]
ggplot(cluster0_6, aes(x = seurat_clusters, fill = label)) + 
  geom_bar(aes(y = after_stat(count)/total_count), position = "stack") +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
  scale_fill_manual(values = color_palette_paired_injury) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

#Figure S2J - part 2
#clusters 7-17 for sham cells
cluster7_17 <- filtered_data_shams_and_naive[filtered_data_shams_and_naive$seurat_clusters %in% 7:17, ]
ggplot(cluster7_17, aes(x = seurat_clusters, fill = label)) + 
  geom_bar(aes(y = after_stat(count)/total_count), position = "stack") +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
  scale_fill_manual(values = color_palette_paired_injury) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
# Calculate total count from original dataset
total_count <- nrow(filtered_data_injury)
#Figure S2K - part 1
#clusters 0-7 for injured SCs
cluster0_7 <- filtered_data_injury[filtered_data_injury$seurat_clusters %in% 0:7, ]

ggplot(cluster0_7, aes(x = seurat_clusters, fill = label)) + 
  geom_bar(aes(y = after_stat(count)/total_count), position = "stack") +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
  scale_fill_manual(values = color_palette_paired_injury) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

#Figure S2K - part 2
#clusters 8-17 for injured SCs
cluster8_17 <- filtered_data_injury[filtered_data_injury$seurat_clusters %in% 8:17, ]
ggplot(cluster8_17, aes(x = seurat_clusters, fill = label)) + 
  geom_bar(aes(y = after_stat(count)/total_count), position = "stack") +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
  scale_fill_manual(values = color_palette_paired_injury) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

#Which clusters are Schwann Cells based on the commonly used myelinating and non-myelinating Schwann Cell markers
#Figure S2I
DotPlot(Cond.combineds_integrated, features=c("Mbp", "Ncam1","Cadm2")) # clusters 2,5,6,7,8 


#Figure S2C
VlnPlot(Cond.combineds_integrated, features="nCount_RNA", group.by="seurat_clusters", pt.size=0)
#Figure S2D
VlnPlot(Cond.combineds_integrated, features="nFeature_RNA", group.by="seurat_clusters", pt.size=0)
#Figure S2E
VlnPlot(Cond.combineds_integrated, features="percent.mt", group.by="seurat_clusters", pt.size=0)

#Figure S2F
VlnPlot(Cond.combineds_integrated, features="nCount_RNA", group.by="batch.orig", pt.size=0)
#Figure S2G
VlnPlot(Cond.combineds_integrated, features="nFeature_RNA", group.by="batch.orig", pt.size=0)
#Figure S2H
VlnPlot(Cond.combineds_integrated, features="percent.mt", group.by="batch.orig", pt.size=0)


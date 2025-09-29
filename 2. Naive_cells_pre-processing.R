library(Seurat)
library(dplyr)
library(tidyverse)
library(BiocManager)
library(BUSpaRse)
library(ggplot2)
library(patchwork)
library(pheatmap)
#Set working directory
#cDNA
Naive.umis <- Read10X_h5('Data_Files/cDNA/KS220914_Force_Naives_filtered_feature_bc_matrix.h5')
dim(Naive.umis)
#HTO
Naive.htos <- BUSpaRse::read_count_output("Data_Files/HTOs/KS220914_sham_htos" , name = "cells_x_features", tcc = FALSE)
head(Naive.htos,n=6)
rownames(Naive.htos)
dim(Naive.htos)

# Join the RNA and HTO
#here is where we have to rename UMI columns to match HTO
head(colnames(Naive.htos))
head(colnames(Naive.umis))
colnames(Naive.umis) <- gsub("-1", "", colnames(Naive.umis))
tail(colnames(Naive.umis))

#all colnames had 1 but now we removed the -1
joint.bcs_Naive <- dplyr::intersect(colnames(Naive.htos), colnames(Naive.umis))

length(joint.bcs_Naive)
# Subset RNA and HTO counts by joint cell barcodes
Naive.umis <- Naive.umis[, joint.bcs_Naive]
Naive.htos<- as.matrix(Naive.htos[, joint.bcs_Naive])
head(colnames(Naive.htos))


dim(Naive.htos)
# Confirm that the HTO have the correct names
rownames(Naive.htos)

# Setup Seurat object
Naive.hashtag <- CreateSeuratObject(counts = Naive.umis, min.cells = 3)
# Normalize RNA data with log normalization
Naive.hashtag <- NormalizeData(Naive.hashtag)
# Find and scale variable features
Naive.hashtag <- FindVariableFeatures(Naive.hashtag, selection.method = "mean.var.plot")
Naive.hashtag <- ScaleData(Naive.hashtag, features = VariableFeatures(Naive.hashtag))
# Add HTO data as a new assay independent from RNA
Naive.hashtag[["HTO"]] <- CreateAssayObject(counts = Naive.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
Naive.hashtag <- NormalizeData(Naive.hashtag, assay = "HTO", normalization.method = "CLR", margin=2)
#Demultiplex the hashtags
Naive.hashtag <- MULTIseqDemux(Naive.hashtag, assay = "HTO", autoThresh = TRUE)

table(Naive.hashtag$MULTI_ID)

# Group cells based on the max HTO signal
Idents(Naive.hashtag) <- "MULTI_ID"
library(patchwork)
VlnPlot(Naive.hashtag, assay = "HTO", features = rownames(Naive.hashtag[["HTO"]])[1:4], ncol = 2)

RidgePlot(Naive.hashtag, assay = "HTO", features = rownames(Naive.hashtag[["HTO"]])[1:4],ncol = 2)


Idents(Naive.hashtag) <- "MULTI_ID"
VlnPlot(Naive.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)



# 1. Get the normalized HTO data
hto_data <- GetAssayData(Naive.hashtag, assay = "HTO", slot = "data")

# 2. Sample cells (to avoid plotting too many)
set.seed(42) # For reproducibility
cells_to_plot <- sample(colnames(Naive.hashtag), min(5000, ncol(Naive.hashtag)))

# 3. Order cells by their classification
cell_classifications <- Naive.hashtag$MULTI_ID[cells_to_plot]
cells_ordered <- cells_to_plot[order(cell_classifications)]

# 4. Create the heatmap
library(pheatmap)
hto_data_subset <- hto_data[, cells_ordered]
pheatmap(hto_data_subset, cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE, annotation_col = data.frame(ID = Naive.hashtag$MULTI_ID[cells_ordered],
                                                            row.names = cells_ordered))



# Extract the singlets
Naive.singlet <- subset(Naive.hashtag, idents = c("Doublet","Negative"), invert=TRUE)
#Name the project
current.project.id <- 'SeuratProject'
new.project.id <- 'Naive'
Naive.singlet@meta.data$orig.ident <- plyr::mapvalues(x = Naive.singlet@meta.data$orig.ident, from = current.project.id, to = new.project.id)

Naive.singlet[["percent.mt"]] <- PercentageFeatureSet(Naive.singlet, pattern = "^mt-")
VlnPlot(Naive.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
length(WhichCells(Naive.singlet))


#########
# Calculate means and SDs
feature_mean <- mean(Naive.singlet$nFeature_RNA)
count_mean <- mean(Naive.singlet$nCount_RNA)
feature_sd <- sd(Naive.singlet$nFeature_RNA) 
count_sd <- sd(Naive.singlet$nCount_RNA)

# Subset using 2 SD cutoffs
Naive.singlet <- subset(Naive.singlet, 
                        nFeature_RNA > feature_mean - 2*feature_sd &
                          nCount_RNA > count_mean - 2*count_sd &
                          percent.mt < 5)
#Naive.singlet <- subset(Naive.singlet,nFeature_RNA > 100 &nFeature_RNA < 1000 &nCount_RNA > 100 &nCount_RNA < 2000&percent.mt < 5)
# Visualize results
VlnPlot(Naive.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by="orig.ident")
length(WhichCells(Naive.singlet))
table(Naive.singlet$MULTI_ID)

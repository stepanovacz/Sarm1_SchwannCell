library(Seurat)
library(dplyr)
library(tidyverse)
library(BiocManager)
library(BUSpaRse)
library(ggplot2)
library(patchwork)
library(pheatmap)

####Process sham samples
#cDNA
ShamB1.umis <- Read10X_h5('Data_Files/cDNA/KS220601_ForceCells_Shams_filtered_feature_bc_matrix.h5')
dim(ShamB1.umis)
#HTOsSo 
ShamB1.htos <- BUSpaRse::read_count_output('Data_Files/HTOs/KS220601_Sham_counts_unfiltered', name = "cells_x_features", tcc = FALSE)
head(ShamB1.htos,n=6)
rownames(ShamB1.htos)
dim(ShamB1.htos)

# Join the RNA and HTO
#here is where we have to rename UMI columns to match HTO
head(colnames(ShamB1.htos))
head(colnames(ShamB1.umis))
colnames(ShamB1.umis) <- gsub("-1", "", colnames(ShamB1.umis))
tail(colnames(ShamB1.umis))

#all colnames had 1 but now we removed the -1
joint.bcsB1_Sham <- dplyr::intersect(colnames(ShamB1.htos), colnames(ShamB1.umis))
length(joint.bcsB1_Sham)

# Subset RNA and HTO counts by joint cell barcodes
ShamB1.umis <- ShamB1.umis[, joint.bcsB1_Sham]
ShamB1.htos<- as.matrix(ShamB1.htos[, joint.bcsB1_Sham])
head(colnames(ShamB1.htos))
dim(ShamB1.htos)
# Confirm that the HTO have the correct names
rownames(ShamB1.htos)

# Setup Seurat object
ShamB1.hashtag <- CreateSeuratObject(counts = ShamB1.umis, min.cells = 3)
# Normalize RNA data with log normalization
ShamB1.hashtag <- NormalizeData(ShamB1.hashtag)
# Find and scale variable features
ShamB1.hashtag <- FindVariableFeatures(ShamB1.hashtag, selection.method = "mean.var.plot")
ShamB1.hashtag <- ScaleData(ShamB1.hashtag, features = VariableFeatures(ShamB1.hashtag))
# Add HTO data as a new assay independent from RNA
ShamB1.hashtag[["HTO"]] <- CreateAssayObject(counts = ShamB1.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
ShamB1.hashtag <- NormalizeData(ShamB1.hashtag, assay = "HTO", normalization.method = "CLR", margin=2)
#Hashtag demultiplexing
ShamB1.hashtag <- MULTIseqDemux(ShamB1.hashtag, assay = "HTO", autoThresh = TRUE)
table(ShamB1.hashtag$MULTI_ID)

# Group cells based on the max HTO signal
Idents(ShamB1.hashtag) <- "MULTI_ID"
library(patchwork)
VlnPlot(ShamB1.hashtag, assay = "HTO", features = rownames(ShamB1.hashtag[["HTO"]])[1:4], ncol = 2)
RidgePlot(ShamB1.hashtag, assay = "HTO", features = rownames(ShamB1.hashtag[["HTO"]])[1:4],ncol = 2)
Idents(ShamB1.hashtag) <- "MULTI_ID"
VlnPlot(ShamB1.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)



# 1. Get the normalized HTO data
hto_data_ShamB1.hashtag <- GetAssayData(ShamB1.hashtag, assay = "HTO", slot = "data")

# 2. Sample cells (to avoid plotting too many)
set.seed(42) # For reproducibility
cells_to_plot_ShamB1.hashtag <- sample(colnames(ShamB1.hashtag), min(5000, ncol(ShamB1.hashtag)))

# 3. Order cells by their classification
cell_classifications_ShamB1.hashtag <- ShamB1.hashtag$MULTI_ID[cells_to_plot_ShamB1.hashtag]
cells_ordered_ShamB1.hashtag <- cells_to_plot_ShamB1.hashtag[order(cell_classifications_ShamB1.hashtag)]

# 4. Create the heatmap
hto_data_subset_ShamB1.hashtag <- hto_data_ShamB1.hashtag[, cells_ordered_ShamB1.hashtag]
pheatmap(hto_data_subset_ShamB1.hashtag, cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE, annotation_col = data.frame(ID = ShamB1.hashtag$MULTI_ID[cells_ordered_ShamB1.hashtag],
                                                            row.names = cells_ordered_ShamB1.hashtag))


#######Processing injured samples
#cDNA
InjuryB1.umis <- Read10X_h5('Data_Files/cDNA/KS220601_Forced_injury_filtered_feature_bc_matrix.h5')
dim(InjuryB1.umis)
#HTOs
InjuryB1.htos <- BUSpaRse::read_count_output("Data_Files/HTOs/KS220601_Undetermined_counts_unfiltered" , name = "cells_x_features", tcc = FALSE)
head(InjuryB1.htos, n = 4)
rownames(InjuryB1.htos)
dim(InjuryB1.htos)

# Join the RNA and HTO
#here is where we have to rename UMI columns to match HTO
head(colnames(InjuryB1.htos))
head(colnames(InjuryB1.umis))
colnames(InjuryB1.umis) <- gsub("-1", "", colnames(InjuryB1.umis))
tail(colnames(InjuryB1.umis))

#all colnames had 1 but now we removed the -1
joint.bcs_injuryB1 <- dplyr::intersect(colnames(InjuryB1.htos), colnames(InjuryB1.umis))
head(joint.bcs_injuryB1)
length(joint.bcs_injuryB1)
# Subset RNA and HTO counts by joint cell barcodes
InjuryB1.umis <- InjuryB1.umis[, joint.bcs_injuryB1]
InjuryB1.htos<- as.matrix(InjuryB1.htos[, joint.bcs_injuryB1])
head(colnames(InjuryB1.htos))


dim(InjuryB1.htos)
# Confirm that the HTO have the correct names
rownames(InjuryB1.htos)

# Setup Seurat object
InjuryB1_hashtag <- CreateSeuratObject(counts = InjuryB1.umis, min.cells = 3)
# Normalize RNA data with log normalization
InjuryB1_hashtag <- NormalizeData(InjuryB1_hashtag)
# Find and scale variable features
InjuryB1_hashtag <- FindVariableFeatures(InjuryB1_hashtag, selection.method = "mean.var.plot")
InjuryB1_hashtag <- ScaleData(InjuryB1_hashtag, features = VariableFeatures(InjuryB1_hashtag))
# Add HTO data as a new assay independent from RNA
InjuryB1_hashtag[["HTO"]] <- CreateAssayObject(counts = InjuryB1.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
InjuryB1_hashtag <- NormalizeData(InjuryB1_hashtag, assay = "HTO", normalization.method = "CLR", margin=2)
#Hashtag demultiplexing
InjuryB1_hashtag <- MULTIseqDemux(InjuryB1_hashtag, assay = "HTO", autoThresh = TRUE)
table(InjuryB1_hashtag$MULTI_ID)

# Group cells based on the max HTO signal
Idents(InjuryB1_hashtag) <- "MULTI_ID"
VlnPlot(InjuryB1_hashtag, assay = "HTO", features = rownames(InjuryB1_hashtag[["HTO"]])[1:4], ncol = 2)
RidgePlot(InjuryB1_hashtag, assay = "HTO", features = rownames(InjuryB1_hashtag[["HTO"]])[1:4],ncol = 2)
Idents(InjuryB1_hashtag) <- "MULTI_ID"
VlnPlot(InjuryB1_hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)



# 1. Get the normalized HTO data
hto_data_InjuryB1_hashtag <- GetAssayData(InjuryB1_hashtag, assay = "HTO", slot = "data")

# 2. Sample cells (to avoid plotting too many)
set.seed(42) # For reproducibility
cells_to_plot_InjuryB1_hashtag <- sample(colnames(InjuryB1_hashtag), min(5000, ncol(InjuryB1_hashtag)))

# 3. Order cells by their classification
cell_classifications_InjuryB1_hashtag <- InjuryB1_hashtag$MULTI_ID[cells_to_plot_InjuryB1_hashtag]
cells_ordered_InjuryB1_hashtag <- cells_to_plot_InjuryB1_hashtag[order(cell_classifications_InjuryB1_hashtag)]

# 4. Create the heatmap
library(pheatmap)
hto_data_subset_InjuryB1_hashtag <- hto_data_InjuryB1_hashtag[, cells_ordered_InjuryB1_hashtag]
pheatmap(hto_data_subset_InjuryB1_hashtag, cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE, annotation_col = data.frame(ID = InjuryB1_hashtag$MULTI_ID[cells_ordered_InjuryB1_hashtag],
                                                            row.names = cells_ordered_InjuryB1_hashtag))

# Cleaning up the data, making sure that the overlapping cell barcodes are removed

# First, identify overlapping barcodes for batch 1
sham_barcodes <- colnames(ShamB1.hashtag)
injury_barcodes <- colnames(InjuryB1_hashtag)
overlapping_barcodes <- intersect(sham_barcodes, injury_barcodes)

cat("Found", length(overlapping_barcodes), "overlapping barcodes\n")


# This ensures completely independent datasets
cat("\n=== OPTION 1: Remove from Injured dataset ===\n")
Sham.singletB1.clean <- ShamB1.hashtag
Injury.singletB1.clean <- InjuryB1_hashtag[, !colnames(InjuryB1_hashtag) %in% overlapping_barcodes]

# make sure that all overlapping cell barcodes are removed
sham_barcodes <- colnames(Sham.singletB1.clean)
injury_barcodes <- colnames(Injury.singletB1.clean)
overlapping_barcodes <- intersect(sham_barcodes, injury_barcodes)
cat("Found", length(overlapping_barcodes), "overlapping barcodes\n")


# Extract the singlets sham cells
Sham.singletB1 <- subset(Sham.singletB1.clean, idents = c("Doublet","Negative"), invert=TRUE)
#Name the project
current.project.id <- 'SeuratProject'
new.project.id <- 'Sham.singletB1'
Sham.singletB1@meta.data$orig.ident <- plyr::mapvalues(x = Sham.singletB1@meta.data$orig.ident, from = current.project.id, to = new.project.id)

Sham.singletB1[["percent.mt"]] <- PercentageFeatureSet(Sham.singletB1, pattern = "^mt-")
VlnPlot(Sham.singletB1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
length(WhichCells(Sham.singletB1))


#########
# Calculate means and SDs
feature_mean <- mean(Sham.singletB1$nFeature_RNA)
count_mean <- mean(Sham.singletB1$nCount_RNA)
feature_sd <- sd(Sham.singletB1$nFeature_RNA) 
count_sd <- sd(Sham.singletB1$nCount_RNA)

# Subset using 2 SD cutoffs
Sham.singletB1 <- subset(Sham.singletB1, 
                         nFeature_RNA > feature_mean - 2*feature_sd &
                           nCount_RNA > count_mean - 2*count_sd &
                           percent.mt < 5)
# Visualize results
VlnPlot(Sham.singletB1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by="orig.ident")
length(WhichCells(Sham.singletB1))
table(Sham.singletB1$MULTI_ID)

# Extract the singlets for injured samples
Injury.singletB1 <- subset(Injury.singletB1.clean, idents = c("Doublet","Negative"), invert=TRUE)
#Name the project
current.project.id <- 'SeuratProject'
new.project.id <- 'Injury.singletB1'
Injury.singletB1@meta.data$orig.ident <- plyr::mapvalues(x = Injury.singletB1@meta.data$orig.ident, from = current.project.id, to = new.project.id)

Injury.singletB1[["percent.mt"]] <- PercentageFeatureSet(Injury.singletB1, pattern = "^mt-")
VlnPlot(Injury.singletB1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
length(WhichCells(Injury.singletB1))
VlnPlot(Injury.singletB1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by="orig.ident")


#########
# Calculate means and SDs
feature_mean <- mean(Injury.singletB1$nFeature_RNA)
count_mean <- mean(Injury.singletB1$nCount_RNA)
feature_sd <- sd(Injury.singletB1$nFeature_RNA) 
count_sd <- sd(Injury.singletB1$nCount_RNA)

# Subset using 2 SD cutoffs
Injury.singletB1 <- subset(Injury.singletB1, 
                           nFeature_RNA > feature_mean - 2*feature_sd &
                             nCount_RNA > count_mean - 2*count_sd &
                             percent.mt < 5)
# Visualize results
VlnPlot(Injury.singletB1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by="orig.ident")
VlnPlot(Injury.singletB1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
length(WhichCells(Injury.singletB1))
table(Injury.singletB1$MULTI_ID)


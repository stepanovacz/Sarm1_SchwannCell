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
ShamB4.umis <- Read10X_h5('Data_Files/cDNA/KS230616_Force_Shams_filtered_feature_bc_matrix.h5')
#HTO
ShamB4.htos <- BUSpaRse::read_count_output('Data_Files/HTOs/KS230616_sham_counts_unfiltered', name = "cells_x_features", tcc = FALSE)
head(ShamB4.htos,n=6)
rownames(ShamB4.htos)
dim(ShamB4.htos)
ShamB4.htos
# Join the RNA and HTO
#here is where we have to rename UMI columns to match HTO
head(colnames(ShamB4.htos))
head(colnames(ShamB4.umis))
colnames(ShamB4.umis) <- gsub("-1", "", colnames(ShamB4.umis))
tail(colnames(ShamB4.umis))
#all colnames had 1 but now we removed the -1

#join cDNA and HTO data
joint.bcsB4_Sham <- dplyr::intersect(colnames(ShamB4.htos), colnames(ShamB4.umis))

length(joint.bcsB4_Sham)
# Subset RNA and HTO counts by joint cell barcodes
ShamB4.umis <- ShamB4.umis[, joint.bcsB4_Sham]
ShamB4.htos<- as.matrix(ShamB4.htos[, joint.bcsB4_Sham])
head(colnames(ShamB4.htos))


dim(ShamB4.htos)
# Confirm that the HTO have the correct names
rownames(ShamB4.htos)

# Setup Seurat object
ShamB4.hashtag <- CreateSeuratObject(counts = ShamB4.umis, min.cells = 3)
# Normalize RNA data with log normalization
ShamB4.hashtag <- NormalizeData(ShamB4.hashtag)
# Find and scale variable features
ShamB4.hashtag <- FindVariableFeatures(ShamB4.hashtag, selection.method = "mean.var.plot")
ShamB4.hashtag <- ScaleData(ShamB4.hashtag, features = VariableFeatures(ShamB4.hashtag))
# Add HTO data as a new assay independent from RNA
ShamB4.hashtag[["HTO"]] <- CreateAssayObject(counts = ShamB4.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
ShamB4.hashtag <- NormalizeData(ShamB4.hashtag, assay = "HTO", normalization.method = "CLR", margin=2)

#demultiplex hashtags
ShamB4.hashtag <- MULTIseqDemux(ShamB4.hashtag, assay = "HTO", autoThresh = TRUE)
table(ShamB4.hashtag$MULTI_ID)

# Group cells based on the max HTO signal
Idents(ShamB4.hashtag) <- "MULTI_ID"
VlnPlot(ShamB4.hashtag, assay = "HTO", features = rownames(ShamB4.hashtag[["HTO"]])[1:4], ncol = 2)
RidgePlot(ShamB4.hashtag, assay = "HTO", features = rownames(ShamB4.hashtag[["HTO"]])[1:4],ncol = 2)
Idents(ShamB4.hashtag) <- "MULTI_ID"
VlnPlot(ShamB4.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)



# 1. Get the normalized HTO data
hto_data_ShamB4.hashtag <- GetAssayData(ShamB4.hashtag, assay = "HTO", slot = "data")

# 2. Sample cells (to avoid plotting too many)
set.seed(42) # For reproducibility
cells_to_plot__ShamB4.hashtag <- sample(colnames(ShamB4.hashtag), min(5000, ncol(ShamB4.hashtag)))

# 3. Order cells by their classification
cell_classifications__ShamB4.hashtag <- ShamB4.hashtag$MULTI_ID[cells_to_plot__ShamB4.hashtag]
cells_ordered__ShamB4.hashtag <- cells_to_plot__ShamB4.hashtag[order(cell_classifications__ShamB4.hashtag)]

# 4. Create the heatmap
hto_data_subset__ShamB4.hashtag <- hto_data_ShamB4.hashtag[, cells_ordered__ShamB4.hashtag]
pheatmap(hto_data_subset__ShamB4.hashtag, cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE, annotation_col = data.frame(ID = ShamB4.hashtag$MULTI_ID[cells_ordered__ShamB4.hashtag],
                                                            row.names = cells_ordered__ShamB4.hashtag))


########processing the injured data
#read in the cDNA
InjuryB4.umis <- Read10X_h5('Data_Files/cDNA/KS230616_Force_Injury_filtered_feature_bc_matrix.h5')
dim(InjuryB4.umis)
#read in the hashtags
InjuryB4.htos <- BUSpaRse::read_count_output('Data_Files/HTOs/KS230616_Undetermined_counts_unfiltered', name = "cells_x_features", tcc = FALSE)
head(InjuryB4.htos, n = 4)
rownames(InjuryB4.htos)
dim(InjuryB4.htos)

# Join the RNA and HTO
#here is where we have to rename UMI columns to match HTO
head(colnames(InjuryB4.htos))
head(colnames(InjuryB4.umis))
colnames(InjuryB4.umis) <- gsub("-1", "", colnames(InjuryB4.umis))
tail(colnames(InjuryB4.umis))

#all colnames had 1 but now we removed the -1
joint.bcs_injuryB4 <- dplyr::intersect(colnames(InjuryB4.htos), colnames(InjuryB4.umis))
head(joint.bcs_injuryB4)
length(joint.bcs_injuryB4)
# Subset RNA and HTO counts by joint cell barcodes
InjuryB4.umis <- InjuryB4.umis[, joint.bcs_injuryB4]
InjuryB4.htos<- as.matrix(InjuryB4.htos[, joint.bcs_injuryB4])
head(colnames(InjuryB4.htos))
dim(InjuryB4.htos)
# Confirm that the HTO have the correct names
rownames(InjuryB4.htos)

# Setup Seurat object
InjuryB4_hashtag <- CreateSeuratObject(counts = InjuryB4.umis, min.cells = 3)
# Normalize RNA data with log normalization
InjuryB4_hashtag <- NormalizeData(InjuryB4_hashtag)
# Find and scale variable features
InjuryB4_hashtag <- FindVariableFeatures(InjuryB4_hashtag, selection.method = "mean.var.plot")
InjuryB4_hashtag <- ScaleData(InjuryB4_hashtag, features = VariableFeatures(InjuryB4_hashtag))
# Add HTO data as a new assay independent from RNA
InjuryB4_hashtag[["HTO"]] <- CreateAssayObject(counts = InjuryB4.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
InjuryB4_hashtag <- NormalizeData(InjuryB4_hashtag, assay = "HTO", normalization.method = "CLR", margin=2)

#Demultiplex the data
InjuryB4_hashtag <- MULTIseqDemux(InjuryB4_hashtag, assay = "HTO", autoThresh = TRUE)
table(InjuryB4_hashtag$MULTI_ID)

# Group cells based on the max HTO signal
Idents(InjuryB4_hashtag) <- "MULTI_ID"
VlnPlot(InjuryB4_hashtag, assay = "HTO", features = rownames(InjuryB4_hashtag[["HTO"]])[1:4], ncol = 2)

RidgePlot(InjuryB4_hashtag, assay = "HTO", features = rownames(InjuryB4_hashtag[["HTO"]])[1:4],ncol = 2)


Idents(InjuryB4_hashtag) <- "MULTI_ID"
VlnPlot(InjuryB4_hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


# 1. Get the normalized HTO data
hto_data_InjuryB4_hashtag <- GetAssayData(InjuryB4_hashtag, assay = "HTO", slot = "data")

# 2. Sample cells (to avoid plotting too many)
set.seed(42) # For reproducibility
cells_to_plot__InjuryB4_hashtag <- sample(colnames(InjuryB4_hashtag), min(5000, ncol(InjuryB4_hashtag)))

# 3. Order cells by their classification
cell_classifications__InjuryB4_hashtag <- InjuryB4_hashtag$MULTI_ID[cells_to_plot__InjuryB4_hashtag]
cells_ordered__InjuryB4_hashtag <- cells_to_plot__InjuryB4_hashtag[order(cell_classifications__InjuryB4_hashtag)]

# 4. Create the heatmap
hto_data_subset__InjuryB4_hashtag <- hto_data_InjuryB4_hashtag[, cells_ordered__InjuryB4_hashtag]
pheatmap(hto_data_subset__InjuryB4_hashtag, cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE, annotation_col = data.frame(ID = InjuryB4_hashtag$MULTI_ID[cells_ordered__InjuryB4_hashtag],
                                                            row.names = cells_ordered__InjuryB4_hashtag))


# First, identify overlapping barcodes
sham_barcodes <- colnames(ShamB4.hashtag)
injury_barcodes <- colnames(InjuryB4_hashtag)
overlapping_barcodes <- intersect(sham_barcodes, injury_barcodes)

cat("Found", length(overlapping_barcodes), "overlapping barcodes\n")

# This ensures completely independent dataset
cat("\ Remove from injured dataset\n")
ShamB4.hashtag.clean <- ShamB4.hashtag
InjuryB4.hashtag.clean <- InjuryB4_hashtag[, !colnames(InjuryB4_hashtag) %in% overlapping_barcodes]

# Check that overlapping celll barcodes have been removed
sham_barcodes <- colnames(ShamB4.hashtag.clean)
injury_barcodes <- colnames(InjuryB4.hashtag.clean)
overlapping_barcodes <- intersect(sham_barcodes, injury_barcodes)
cat("Found", length(overlapping_barcodes), "overlapping barcodes\n")


# Extract the singlets for shams
Sham.singletB4 <- subset(ShamB4.hashtag.clean, idents = c("Doublet","Negative"), invert=TRUE)
#Name the project
current.project.id <- 'SeuratProject'
new.project.id <- 'Sham.singletB4'
Sham.singletB4@meta.data$orig.ident <- plyr::mapvalues(x = Sham.singletB4@meta.data$orig.ident, from = current.project.id, to = new.project.id)

Sham.singletB4[["percent.mt"]] <- PercentageFeatureSet(Sham.singletB4, pattern = "^mt-")
VlnPlot(Sham.singletB4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
length(WhichCells(Sham.singletB4))


#########
# Calculate means and SDs
feature_mean <- mean(Sham.singletB4$nFeature_RNA)
count_mean <- mean(Sham.singletB4$nCount_RNA)
feature_sd <- sd(Sham.singletB4$nFeature_RNA) 
count_sd <- sd(Sham.singletB4$nCount_RNA)

# Subset using 2 SD cutoffs
Sham.singletB4 <- subset(Sham.singletB4, 
                         nFeature_RNA > feature_mean - 2*feature_sd &
                           nCount_RNA > count_mean - 2*count_sd &
                           percent.mt < 5)
# Visualize results
VlnPlot(Sham.singletB4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by="orig.ident")
length(WhichCells(Sham.singletB4))
table(Sham.singletB4$MULTI_ID)


# Extract the singlets for injury
InjuryB4.singlet <- subset(InjuryB4.hashtag.clean, idents = c("Doublet","Negative"), invert=TRUE)
#Name the project
current.project.id <- 'SeuratProject'
new.project.id <- 'Injury.singlet.B4'
InjuryB4.singlet@meta.data$orig.ident <- plyr::mapvalues(x = InjuryB4.singlet@meta.data$orig.ident, from = current.project.id, to = new.project.id)

InjuryB4.singlet[["percent.mt"]] <- PercentageFeatureSet(InjuryB4.singlet, pattern = "^mt-")
VlnPlot(InjuryB4.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
length(WhichCells(InjuryB4.singlet))
VlnPlot(InjuryB4.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by="orig.ident")


#########
# Calculate means and SDs
feature_mean <- mean(InjuryB4.singlet$nFeature_RNA)
count_mean <- mean(InjuryB4.singlet$nCount_RNA)
feature_sd <- sd(InjuryB4.singlet$nFeature_RNA) 
count_sd <- sd(InjuryB4.singlet$nCount_RNA)

# Subset using 2 SD cutoffs
Injury.singletB4 <- subset(InjuryB4.singlet, 
                           nFeature_RNA > feature_mean - 2*feature_sd &
                             nCount_RNA > count_mean - 2*count_sd &
                             percent.mt < 5)
# Visualize results
VlnPlot(Injury.singletB4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by="orig.ident")
VlnPlot(Injury.singletB4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(Injury.singletB4$MULTI_ID)
length(WhichCells(Injury.singletB4))
dim(Injury.singletB4)

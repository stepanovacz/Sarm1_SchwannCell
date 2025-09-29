library(Seurat)
library(dplyr)
library(tidyverse)
library(BiocManager)
library(BUSpaRse)
library(ggplot2)
library(patchwork)
library(pheatmap)

#read in the data
#cDNA
ShamB3.umis <- Read10X_h5('Data_Files/cDNA/KS230314_Force_Shams_filtered_feature_bc_matrix.h5')
dim(ShamB3.umis)
#Hashtag data
ShamB3.htos <- BUSpaRse::read_count_output('Data_Files/HTOs/KS230314_sham_htos', name = "cells_x_features", tcc = FALSE)
head(ShamB3.htos,n=6)
rownames(ShamB3.htos)
dim(ShamB3.htos)

# Join the RNA and HTO
#here is where we have to rename UMI columns to match HTO
head(colnames(ShamB3.htos))
head(colnames(ShamB3.umis))
colnames(ShamB3.umis) <- gsub("-1", "", colnames(ShamB3.umis))
tail(colnames(ShamB3.umis))
#all colnames had 1 but now we removed the -1

#join the hashtag and cDNA data
joint.bcsB3_Sham <- dplyr::intersect(colnames(ShamB3.htos), colnames(ShamB3.umis))

length(joint.bcsB3_Sham)
# Subset RNA and HTO counts by joint cell barcodes
ShamB3.umis <- ShamB3.umis[, joint.bcsB3_Sham]
ShamB3.htos<- as.matrix(ShamB3.htos[, joint.bcsB3_Sham])
head(colnames(ShamB3.htos))
dim(ShamB3.htos)
# Confirm that the HTO have the correct names
rownames(ShamB3.htos)

# Setup Seurat object
ShamB3.hashtag <- CreateSeuratObject(counts = ShamB3.umis, min.cells = 3)
# Normalize RNA data with log normalization
ShamB3.hashtag <- NormalizeData(ShamB3.hashtag)
# Find and scale variable features
ShamB3.hashtag <- FindVariableFeatures(ShamB3.hashtag, selection.method = "mean.var.plot")
ShamB3.hashtag <- ScaleData(ShamB3.hashtag, features = VariableFeatures(ShamB3.hashtag))
# Add HTO data as a new assay independent from RNA
ShamB3.hashtag[["HTO"]] <- CreateAssayObject(counts = ShamB3.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
ShamB3.hashtag <- NormalizeData(ShamB3.hashtag, assay = "HTO", normalization.method = "CLR", margin=2)
#Demultiplex the hashtags
ShamB3.hashtag <- MULTIseqDemux(ShamB3.hashtag, assay = "HTO", autoThresh = TRUE)

table(ShamB3.hashtag$MULTI_ID)

# Group cells based on the max HTO signal
Idents(ShamB3.hashtag) <- "MULTI_ID"
VlnPlot(ShamB3.hashtag, assay = "HTO", features = rownames(ShamB3.hashtag[["HTO"]])[1:4], ncol = 2)
RidgePlot(ShamB3.hashtag, assay = "HTO", features = rownames(ShamB3.hashtag[["HTO"]])[1:4],ncol = 2)


Idents(ShamB3.hashtag) <- "MULTI_ID"
VlnPlot(ShamB3.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


# 1. Get the normalized HTO data
hto_data_ShamB3.hashtag <- GetAssayData(ShamB3.hashtag, assay = "HTO", slot = "data")

# 2. Sample cells (to avoid plotting too many)
set.seed(42) # For reproducibility
cells_to_plot_ShamB3.hashtag <- sample(colnames(ShamB3.hashtag), min(5000, ncol(ShamB3.hashtag)))

# 3. Order cells by their classification
cell_classifications_ShamB3.hashtag <- ShamB3.hashtag$MULTI_ID[cells_to_plot_ShamB3.hashtag]
cells_ordered_ShamB3.hashtag <- cells_to_plot_ShamB3.hashtag[order(cell_classifications_ShamB3.hashtag)]

# 4. Create the heatmap
hto_data_subset_ShamB3.hashtag <- hto_data_ShamB3.hashtag[, cells_ordered_ShamB3.hashtag]
pheatmap(hto_data_subset_ShamB3.hashtag, cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE, annotation_col = data.frame(ID = ShamB3.hashtag$MULTI_ID[cells_ordered_ShamB3.hashtag],
                                                            row.names = cells_ordered_ShamB3.hashtag))

##Now process the injured data
#Read in cDNA data
InjuryB3.umis <- Read10X_h5('Data_Files/cDNA/KS230314_Expect_Injured_filtered_feature_bc_matrix.h5')
dim(InjuryB3.umis)
#Read in HTO data
InjuryB3.htos <- BUSpaRse::read_count_output("Data_Files/HTOs/KS230314_injury_htos" , name = "cells_x_features", tcc = FALSE)
head(InjuryB3.htos, n = 4)
rownames(InjuryB3.htos)

dim(InjuryB3.htos)
InjuryB3.htos
# Join the RNA and HTO
#here is where we have to rename UMI columns to match HTO
head(colnames(InjuryB3.htos))
head(colnames(InjuryB3.umis))
colnames(InjuryB3.umis) <- gsub("-1", "", colnames(InjuryB3.umis))
tail(colnames(InjuryB3.umis))
#all colnames had 1 but now we removed the -1

#overlap cDNA and HTO data

joint.bcs_injuryB3 <- dplyr::intersect(colnames(InjuryB3.htos), colnames(InjuryB3.umis))
head(joint.bcs_injuryB3)
length(joint.bcs_injuryB3)
# Subset RNA and HTO counts by joint cell barcodes
InjuryB3.umis <- InjuryB3.umis[, joint.bcs_injuryB3]
InjuryB3.htos<- as.matrix(InjuryB3.htos[, joint.bcs_injuryB3])
head(colnames(InjuryB3.htos))


dim(InjuryB3.htos)
# Confirm that the HTO have the correct names
rownames(InjuryB3.htos)

# Setup Seurat object
InjuryB3_hashtag <- CreateSeuratObject(counts = InjuryB3.umis, min.cells = 3)
# Normalize RNA data with log normalization
InjuryB3_hashtag <- NormalizeData(InjuryB3_hashtag)
# Find and scale variable features
InjuryB3_hashtag <- FindVariableFeatures(InjuryB3_hashtag, selection.method = "mean.var.plot")
InjuryB3_hashtag <- ScaleData(InjuryB3_hashtag, features = VariableFeatures(InjuryB3_hashtag))
# Add HTO data as a new assay independent from RNA
InjuryB3_hashtag[["HTO"]] <- CreateAssayObject(counts = InjuryB3.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
InjuryB3_hashtag <- NormalizeData(InjuryB3_hashtag, assay = "HTO", normalization.method = "CLR", margin=2)

#Demultiplex hashtags
InjuryB3_hashtag <- MULTIseqDemux(InjuryB3_hashtag, assay = "HTO", autoThresh = TRUE)

table(InjuryB3_hashtag$MULTI_ID)

# Group cells based on the max HTO signal
Idents(InjuryB3_hashtag) <- "MULTI_ID"
VlnPlot(InjuryB3_hashtag, assay = "HTO", features = rownames(InjuryB3_hashtag[["HTO"]])[1:4], ncol = 2)
RidgePlot(InjuryB3_hashtag, assay = "HTO", features = rownames(InjuryB3_hashtag[["HTO"]])[1:4],ncol = 2)
Idents(InjuryB3_hashtag) <- "MULTI_ID"
VlnPlot(InjuryB3_hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


# 1. Get the normalized HTO data
hto_data_InjuryB3_hashtag <- GetAssayData(InjuryB3_hashtag, assay = "HTO", slot = "data")

# 2. Sample cells (to avoid plotting too many)
set.seed(42) # For reproducibility
cells_to_plot_InjuryB3_hashtag <- sample(colnames(InjuryB3_hashtag), min(5000, ncol(InjuryB3_hashtag)))

# 3. Order cells by their classification
cell_classifications_InjuryB3_hashtag <- InjuryB3_hashtag$MULTI_ID[cells_to_plot_InjuryB3_hashtag]
cells_ordered_InjuryB3_hashtag <- cells_to_plot_InjuryB3_hashtag[order(cell_classifications_InjuryB3_hashtag)]

# 4. Create the heatmap
hto_data_subset_InjuryB3_hashtag <- hto_data_InjuryB3_hashtag[, cells_ordered_InjuryB3_hashtag]
pheatmap(hto_data_subset_InjuryB3_hashtag, cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE, annotation_col = data.frame(ID = InjuryB3_hashtag$MULTI_ID[cells_ordered_InjuryB3_hashtag],
                                                            row.names = cells_ordered_InjuryB3_hashtag))

# Cleaning up the data, making sure that the overlapping cell barcodes are removed

sham_barcodes <- colnames(ShamB3.hashtag)
injury_barcodes <- colnames(InjuryB3_hashtag)
overlapping_barcodes <- intersect(sham_barcodes, injury_barcodes)

cat("Found", length(overlapping_barcodes), "overlapping barcodes\n")

#Remove overlapping cells from injured dataset
# This ensures completely independent datasetscat("\nRemove from injured dataset\n")
Sham.B3.clean <- ShamB3.hashtag
Injury.B3.clean <- InjuryB3_hashtag[, !colnames(InjuryB3_hashtag) %in% overlapping_barcodes]

#Check that the barcodes were removed
sham_barcodes <- colnames(Sham.B3.clean)
injury_barcodes <- colnames(Injury.B3.clean)
overlapping_barcodes <- intersect(sham_barcodes, injury_barcodes)

cat("Found", length(overlapping_barcodes), "overlapping barcodes\n")



# Extract the singlets
Sham.singletB3.clean <- subset(Sham.B3.clean, idents = c("Doublet","Negative"), invert=TRUE)
#Name the project
current.project.id <- 'SeuratProject'
new.project.id <- 'Sham.singletB3'
Sham.singletB3.clean@meta.data$orig.ident <- plyr::mapvalues(x = Sham.singletB3.clean@meta.data$orig.ident, from = current.project.id, to = new.project.id)

Sham.singletB3.clean[["percent.mt"]] <- PercentageFeatureSet(Sham.singletB3.clean, pattern = "^mt-")
VlnPlot(Sham.singletB3.clean, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
length(WhichCells(Sham.singletB3.clean))


#########
# Calculate means and SDs
feature_mean <- mean(Sham.singletB3.clean$nFeature_RNA)
count_mean <- mean(Sham.singletB3.clean$nCount_RNA)
feature_sd <- sd(Sham.singletB3.clean$nFeature_RNA) 
count_sd <- sd(Sham.singletB3.clean$nCount_RNA)

# Subset using 2 SD cutoffs
Sham.singletB3 <- subset(Sham.singletB3.clean,
                         nFeature_RNA > feature_mean - 2*feature_sd &
                           nCount_RNA > count_mean - 2*count_sd &
                           percent.mt < 5)
# Visualize results
VlnPlot(Sham.singletB3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by="orig.ident")
length(WhichCells(Sham.singletB3))
table(Sham.singletB3$MULTI_ID)


# Extract the singlets from the injured sample
Injury.singletB3.clean <- subset(Injury.B3.clean, idents = c("Doublet","Negative"), invert=TRUE)

#Name the project
current.project.id <- 'SeuratProject'
new.project.id <- 'Injury.singletB3'
Injury.singletB3.clean@meta.data$orig.ident <- plyr::mapvalues(x = Injury.singletB3.clean@meta.data$orig.ident, from = current.project.id, to = new.project.id)

Injury.singletB3.clean[["percent.mt"]] <- PercentageFeatureSet(Injury.singletB3.clean, pattern = "^mt-")
VlnPlot(Injury.singletB3.clean, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
length(WhichCells(Injury.singletB3.clean))
VlnPlot(Injury.singletB3.clean, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by="orig.ident")

#########
# Calculate means and SDs
feature_mean <- mean(Injury.singletB3.clean$nFeature_RNA)
count_mean <- mean(Injury.singletB3.clean$nCount_RNA)
feature_sd <- sd(Injury.singletB3.clean$nFeature_RNA) 
count_sd <- sd(Injury.singletB3.clean$nCount_RNA)

# Subset using 2 SD cutoffs
Injury.singletB3 <- subset(Injury.singletB3.clean, 
                           nFeature_RNA > feature_mean - 2*feature_sd &
                             nCount_RNA > count_mean - 2*count_sd &
                             percent.mt < 5)
# Visualize results
VlnPlot(Injury.singletB3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by="orig.ident")
VlnPlot(Injury.singletB3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(Injury.singletB3$MULTI_ID)
length(WhichCells(Injury.singletB3))


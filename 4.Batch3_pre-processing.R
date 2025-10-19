################################################################################
# Single-Cell RNA-seq Analysis with Cell Hashing - Batch 3
# 
# Description: Processing and demultiplexing of single-cell RNA-seq data using
# hashtag oligonucleotides (HTOs) for sham and injured samples (Batch 3)
#
# Requirements: 
# - 10X Genomics filtered feature barcode matrix (H5 format)
# - Kallisto/BUStools count output for HTOs
# - R packages: Seurat, dplyr, tidyverse, BiocManager, BUSpaRse, 
#   ggplot2, patchwork, pheatmap, plyr
################################################################################

# Load required libraries
library(Seurat)
library(dplyr)
library(tidyverse)
library(BiocManager)
library(BUSpaRse)
library(ggplot2)
library(patchwork)
library(pheatmap)

################################################################################
# USER-DEFINED PATHS - MODIFY THESE FOR YOUR SYSTEM
################################################################################

# Path to cDNA data (10X H5 files)
shamB3_cdna_path <- "Data_Files/cDNA/KS230314_Sham_filtered_feature_bc_matrix.h5"
injuryB3_cdna_path <- "Data_Files/cDNA/KS230314_Injury_filtered_feature_bc_matrix.h5"

# Path to HTO count data (kallisto/BUStools output)
shamB3_hto_path <- "Data_Files/HTOs/KS230314_sham_htos"
injuryB3_hto_path <- "Data_Files/HTOs/KS230314_injury_htos"

################################################################################
# PROCESS SHAM SAMPLES - BATCH 3
################################################################################

cat("=== Processing Sham Samples (Batch 3) ===\n")

# Read in cDNA data
ShamB3.umis <- Read10X_h5(shamB3_cdna_path)
cat("Sham B3 cDNA dimensions:", dim(ShamB3.umis), "\n")

# Read in hashtag data
ShamB3.htos <- BUSpaRse::read_count_output(shamB3_hto_path, 
                                           name = "cells_x_features", 
                                           tcc = FALSE)
cat("Sham B3 HTO dimensions:", dim(ShamB3.htos), "\n")
cat("HTO names:\n")
print(head(ShamB3.htos, n = 6))
print(rownames(ShamB3.htos))

# Join the RNA and HTO
# Rename UMI columns to match HTO (remove -1 suffix)
cat("\nBefore modification:\n")
cat("cDNA barcodes:", head(colnames(ShamB3.umis)), "\n")
cat("HTO barcodes:", head(colnames(ShamB3.htos)), "\n")

colnames(ShamB3.umis) <- gsub("-1", "", colnames(ShamB3.umis))

cat("After modification:\n")
cat("cDNA barcodes:", tail(colnames(ShamB3.umis)), "\n")

# Join the hashtag and cDNA data
joint.bcsB3_Sham <- dplyr::intersect(colnames(ShamB3.htos), colnames(ShamB3.umis))
cat("Joint barcodes found:", length(joint.bcsB3_Sham), "\n")

# Subset RNA and HTO counts by joint cell barcodes
ShamB3.umis <- ShamB3.umis[, joint.bcsB3_Sham]
ShamB3.htos <- as.matrix(ShamB3.htos[, joint.bcsB3_Sham])
cat("Final Sham B3 HTO dimensions:", dim(ShamB3.htos), "\n")
cat("HTO row names confirmed:\n")
print(rownames(ShamB3.htos))

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
ShamB3.hashtag <- NormalizeData(ShamB3.hashtag, assay = "HTO", 
                                normalization.method = "CLR", margin = 2)

# Demultiplex the hashtags
ShamB3.hashtag <- MULTIseqDemux(ShamB3.hashtag, assay = "HTO", autoThresh = TRUE)
cat("\nSham B3 demultiplexing results:\n")
print(table(ShamB3.hashtag$MULTI_ID))

# Group cells based on the max HTO signal
Idents(ShamB3.hashtag) <- "MULTI_ID"

# Visualizations
VlnPlot(ShamB3.hashtag, assay = "HTO", 
        features = rownames(ShamB3.hashtag[["HTO"]])[1:4], ncol = 2)
RidgePlot(ShamB3.hashtag, assay = "HTO", 
          features = rownames(ShamB3.hashtag[["HTO"]])[1:4], ncol = 2)
VlnPlot(ShamB3.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# Create HTO heatmap
# 1. Get the normalized HTO data
hto_data_ShamB3.hashtag <- GetAssayData(ShamB3.hashtag, assay = "HTO", slot = "data")

# 2. Sample cells (to avoid plotting too many)
set.seed(42) # For reproducibility
cells_to_plot_ShamB3.hashtag <- sample(colnames(ShamB3.hashtag), 
                                       min(5000, ncol(ShamB3.hashtag)))

# 3. Order cells by their classification
cell_classifications_ShamB3.hashtag <- ShamB3.hashtag$MULTI_ID[cells_to_plot_ShamB3.hashtag]
cells_ordered_ShamB3.hashtag <- cells_to_plot_ShamB3.hashtag[order(cell_classifications_ShamB3.hashtag)]

# 4. Create the heatmap
hto_data_subset_ShamB3.hashtag <- hto_data_ShamB3.hashtag[, cells_ordered_ShamB3.hashtag]
pheatmap(hto_data_subset_ShamB3.hashtag, cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE, 
         annotation_col = data.frame(ID = ShamB3.hashtag$MULTI_ID[cells_ordered_ShamB3.hashtag],
                                     row.names = cells_ordered_ShamB3.hashtag))

################################################################################
# PROCESS INJURED SAMPLES - BATCH 3
################################################################################

cat("\n=== Processing Injured Samples (Batch 3) ===\n")

# Read in cDNA data
InjuryB3.umis <- Read10X_h5(injuryB3_cdna_path)
cat("Injury B3 cDNA dimensions:", dim(InjuryB3.umis), "\n")

# Read in HTO data
InjuryB3.htos <- BUSpaRse::read_count_output(injuryB3_hto_path, 
                                             name = "cells_x_features", 
                                             tcc = FALSE)
cat("Injury B3 HTO dimensions:", dim(InjuryB3.htos), "\n")
cat("HTO names:\n")
print(head(InjuryB3.htos, n = 4))
print(rownames(InjuryB3.htos))

# Join the RNA and HTO
# Rename UMI columns to match HTO
cat("\nBefore modification:\n")
cat("cDNA barcodes:", head(colnames(InjuryB3.umis)), "\n")
cat("HTO barcodes:", head(colnames(InjuryB3.htos)), "\n")

colnames(InjuryB3.umis) <- gsub("-1", "", colnames(InjuryB3.umis))

cat("After modification:\n")
cat("cDNA barcodes:", tail(colnames(InjuryB3.umis)), "\n")

# Overlap cDNA and HTO data
joint.bcs_injuryB3 <- dplyr::intersect(colnames(InjuryB3.htos), colnames(InjuryB3.umis))
cat("Joint barcodes found:", length(joint.bcs_injuryB3), "\n")

# Subset RNA and HTO counts by joint cell barcodes
InjuryB3.umis <- InjuryB3.umis[, joint.bcs_injuryB3]
InjuryB3.htos <- as.matrix(InjuryB3.htos[, joint.bcs_injuryB3])
cat("Final Injury B3 HTO dimensions:", dim(InjuryB3.htos), "\n")
cat("HTO row names confirmed:\n")
print(rownames(InjuryB3.htos))

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
InjuryB3_hashtag <- NormalizeData(InjuryB3_hashtag, assay = "HTO", 
                                  normalization.method = "CLR", margin = 2)

# Demultiplex hashtags
InjuryB3_hashtag <- MULTIseqDemux(InjuryB3_hashtag, assay = "HTO", autoThresh = TRUE)
cat("\nInjury B3 demultiplexing results:\n")
print(table(InjuryB3_hashtag$MULTI_ID))

# Group cells based on the max HTO signal
Idents(InjuryB3_hashtag) <- "MULTI_ID"

# Visualizations
VlnPlot(InjuryB3_hashtag, assay = "HTO", 
        features = rownames(InjuryB3_hashtag[["HTO"]])[1:4], ncol = 2)
RidgePlot(InjuryB3_hashtag, assay = "HTO", 
          features = rownames(InjuryB3_hashtag[["HTO"]])[1:4], ncol = 2)
VlnPlot(InjuryB3_hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# Create HTO heatmap
# 1. Get the normalized HTO data
hto_data_InjuryB3_hashtag <- GetAssayData(InjuryB3_hashtag, assay = "HTO", slot = "data")

# 2. Sample cells (to avoid plotting too many)
set.seed(42) # For reproducibility
cells_to_plot_InjuryB3_hashtag <- sample(colnames(InjuryB3_hashtag), 
                                         min(5000, ncol(InjuryB3_hashtag)))

# 3. Order cells by their classification
cell_classifications_InjuryB3_hashtag <- InjuryB3_hashtag$MULTI_ID[cells_to_plot_InjuryB3_hashtag]
cells_ordered_InjuryB3_hashtag <- cells_to_plot_InjuryB3_hashtag[order(cell_classifications_InjuryB3_hashtag)]

# 4. Create the heatmap
hto_data_subset_InjuryB3_hashtag <- hto_data_InjuryB3_hashtag[, cells_ordered_InjuryB3_hashtag]
pheatmap(hto_data_subset_InjuryB3_hashtag, cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE, 
         annotation_col = data.frame(ID = InjuryB3_hashtag$MULTI_ID[cells_ordered_InjuryB3_hashtag],
                                     row.names = cells_ordered_InjuryB3_hashtag))

################################################################################
# REMOVE OVERLAPPING BARCODES
################################################################################

cat("\n=== Cleaning Data: Removing Overlapping Barcodes ===\n")

# First, identify overlapping barcodes
sham_barcodes <- colnames(ShamB3.hashtag)
injury_barcodes <- colnames(InjuryB3_hashtag)
overlapping_barcodes <- intersect(sham_barcodes, injury_barcodes)

cat("Found", length(overlapping_barcodes), "overlapping barcodes\n")

# Remove overlapping cells from injured dataset
# This ensures completely independent datasets
cat("\nRemove from injured dataset\n")
Sham.B3.clean <- ShamB3.hashtag
Injury.B3.clean <- InjuryB3_hashtag[, !colnames(InjuryB3_hashtag) %in% overlapping_barcodes]

# Check that the barcodes were removed
sham_barcodes <- colnames(Sham.B3.clean)
injury_barcodes <- colnames(Injury.B3.clean)
overlapping_barcodes <- intersect(sham_barcodes, injury_barcodes)

cat("After cleanup, found", length(overlapping_barcodes), "overlapping barcodes\n")

################################################################################
# EXTRACT SINGLETS AND PERFORM QC - SHAM SAMPLES
################################################################################

cat("\n=== Processing Sham B3 Singlets ===\n")

# Extract the singlets (remove doublets and negatives)
Sham.singletB3.clean <- subset(Sham.B3.clean, idents = c("Doublet", "Negative"), invert = TRUE)

# Rename the project
current.project.id <- 'SeuratProject'
new.project.id <- 'Sham.singletB3'
Sham.singletB3.clean@meta.data$orig.ident <- plyr::mapvalues(x = Sham.singletB3.clean@meta.data$orig.ident, 
                                                             from = current.project.id, 
                                                             to = new.project.id)

# Calculate mitochondrial percentage
Sham.singletB3.clean[["percent.mt"]] <- PercentageFeatureSet(Sham.singletB3.clean, pattern = "^mt-")
VlnPlot(Sham.singletB3.clean, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
cat("Cells before QC filtering:", length(WhichCells(Sham.singletB3.clean)), "\n")

# Calculate means and SDs for QC filtering
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
VlnPlot(Sham.singletB3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, group.by = "orig.ident")
cat("Cells after QC filtering:", length(WhichCells(Sham.singletB3)), "\n")
cat("Sham B3 samples by HTO:\n")
print(table(Sham.singletB3$MULTI_ID))

################################################################################
# EXTRACT SINGLETS AND PERFORM QC - INJURED SAMPLES
################################################################################

cat("\n=== Processing Injured B3 Singlets ===\n")

# Extract the singlets from the injured sample (remove doublets and negatives)
Injury.singletB3.clean <- subset(Injury.B3.clean, idents = c("Doublet", "Negative"), invert = TRUE)

# Rename the project
current.project.id <- 'SeuratProject'
new.project.id <- 'Injury.singletB3'
Injury.singletB3.clean@meta.data$orig.ident <- plyr::mapvalues(x = Injury.singletB3.clean@meta.data$orig.ident, 
                                                               from = current.project.id, 
                                                               to = new.project.id)

# Calculate mitochondrial percentage
Injury.singletB3.clean[["percent.mt"]] <- PercentageFeatureSet(Injury.singletB3.clean, pattern = "^mt-")
VlnPlot(Injury.singletB3.clean, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
cat("Cells before QC filtering:", length(WhichCells(Injury.singletB3.clean)), "\n")
VlnPlot(Injury.singletB3.clean, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, group.by = "orig.ident")

# Calculate means and SDs for QC filtering
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
VlnPlot(Injury.singletB3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, group.by = "orig.ident")
VlnPlot(Injury.singletB3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
cat("Cells after QC filtering:", length(WhichCells(Injury.singletB3)), "\n")
cat("Injured B3 samples by HTO:\n")
print(table(Injury.singletB3$MULTI_ID))

cat("\n=== Processing Complete ===\n")
cat("Final objects created:\n")
cat("  - Sham.singletB3: QC-filtered sham B3 singlets\n")
cat("  - Injury.singletB3: QC-filtered injured B3 singlets\n")
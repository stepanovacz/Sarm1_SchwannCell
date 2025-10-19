################################################################################
# Single-Cell RNA-seq Analysis with Cell Hashing - Batch 2
# 
# Description: Processing and demultiplexing of single-cell RNA-seq data using
# hashtag oligonucleotides (HTOs) for sham and injured samples (Batch 2)
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
shamB2_cdna_path <- "Data_Files/cDNA/KS221027_Sham_filtered_feature_bc_matrix.h5"
injuryB2_cdna_path <- "Data_Files/cDNA/KS221027_Injury_filtered_feature_bc_matrix.h5"

# Path to HTO count data (kallisto/BUStools output)
shamB2_hto_path <- "Data_Files/HTOs/KS231003_sham_counts_unfiltered"
injuryB2_hto_path <- "Data_Files/HTOs/KS221003_Undetermined_counts_unfiltered"

################################################################################
# PROCESS SHAM SAMPLES - BATCH 2
################################################################################

cat("=== Processing Sham Samples (Batch 2) ===\n")

# Load cDNA data
ShamB2.umis <- Read10X_h5(shamB2_cdna_path)
cat("Sham B2 cDNA dimensions:", dim(ShamB2.umis), "\n")

# Load HTOs
ShamB2.htos <- BUSpaRse::read_count_output(shamB2_hto_path, 
                                           name = "cells_x_features", 
                                           tcc = FALSE)
cat("Sham B2 HTO dimensions:", dim(ShamB2.htos), "\n")
cat("HTO names:\n")
print(head(ShamB2.htos, n = 6))
print(rownames(ShamB2.htos))

# Join the RNA and HTO
# Rename UMI columns to match HTO (remove -1 suffix)
cat("\nBefore modification:\n")
cat("cDNA barcodes:", head(colnames(ShamB2.umis)), "\n")
cat("HTO barcodes:", head(colnames(ShamB2.htos)), "\n")

colnames(ShamB2.umis) <- gsub("-1", "", colnames(ShamB2.umis))

cat("After modification:\n")
cat("cDNA barcodes:", tail(colnames(ShamB2.umis)), "\n")

# Join the hashtag and cDNA data
joint.bcsB2_Sham <- dplyr::intersect(colnames(ShamB2.htos), colnames(ShamB2.umis))
cat("Joint barcodes found:", length(joint.bcsB2_Sham), "\n")

# Subset RNA and HTO counts by joint cell barcodes
ShamB2.umis <- ShamB2.umis[, joint.bcsB2_Sham]
ShamB2.htos <- as.matrix(ShamB2.htos[, joint.bcsB2_Sham])
cat("Final Sham B2 HTO dimensions:", dim(ShamB2.htos), "\n")
cat("HTO row names confirmed:\n")
print(rownames(ShamB2.htos))

# Setup Seurat object
ShamB2.hashtag <- CreateSeuratObject(counts = ShamB2.umis, min.cells = 3)

# Normalize RNA data with log normalization
ShamB2.hashtag <- NormalizeData(ShamB2.hashtag)

# Find and scale variable features
ShamB2.hashtag <- FindVariableFeatures(ShamB2.hashtag, selection.method = "mean.var.plot")
ShamB2.hashtag <- ScaleData(ShamB2.hashtag, features = VariableFeatures(ShamB2.hashtag))

# Add HTO data as a new assay independent from RNA
ShamB2.hashtag[["HTO"]] <- CreateAssayObject(counts = ShamB2.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
ShamB2.hashtag <- NormalizeData(ShamB2.hashtag, assay = "HTO", 
                                normalization.method = "CLR", margin = 2)

# Hashtags demultiplexing
ShamB2.hashtag <- MULTIseqDemux(ShamB2.hashtag, assay = "HTO", autoThresh = TRUE)
cat("\nSham B2 demultiplexing results:\n")
print(table(ShamB2.hashtag$MULTI_ID))

# Group cells based on the max HTO signal
Idents(ShamB2.hashtag) <- "MULTI_ID"

# Visualizations
VlnPlot(ShamB2.hashtag, assay = "HTO", 
        features = rownames(ShamB2.hashtag[["HTO"]])[1:4], ncol = 2)
RidgePlot(ShamB2.hashtag, assay = "HTO", 
          features = rownames(ShamB2.hashtag[["HTO"]])[1:4], ncol = 2)
VlnPlot(ShamB2.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# Create HTO heatmap
# 1. Get the normalized HTO data
hto_data_ShamB2.hashtag <- GetAssayData(ShamB2.hashtag, assay = "HTO", slot = "data")

# 2. Sample cells (to avoid plotting too many)
set.seed(42) # For reproducibility
cells_to_plot_ShamB2.hashtag <- sample(colnames(ShamB2.hashtag), 
                                       min(5000, ncol(ShamB2.hashtag)))

# 3. Order cells by their classification
cell_classifications_ShamB2.hashtag <- ShamB2.hashtag$MULTI_ID[cells_to_plot_ShamB2.hashtag]
cells_ordered_ShamB2.hashtag <- cells_to_plot_ShamB2.hashtag[order(cell_classifications_ShamB2.hashtag)]

# 4. Create the heatmap
hto_data_subset_ShamB2.hashtag <- hto_data_ShamB2.hashtag[, cells_ordered_ShamB2.hashtag]
pheatmap(hto_data_subset_ShamB2.hashtag, cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE, 
         annotation_col = data.frame(ID = ShamB2.hashtag$MULTI_ID[cells_ordered_ShamB2.hashtag],
                                     row.names = cells_ordered_ShamB2.hashtag))

################################################################################
# PROCESS INJURED SAMPLES - BATCH 2
################################################################################

cat("\n=== Processing Injured Samples (Batch 2) ===\n")

# Load cDNA data
InjuryB2.umis <- Read10X_h5(injuryB2_cdna_path)
cat("Injury B2 cDNA dimensions:", dim(InjuryB2.umis), "\n")

# Load HTOs
InjuryB2.htos <- BUSpaRse::read_count_output(injuryB2_hto_path, 
                                             name = "cells_x_features", 
                                             tcc = FALSE)
cat("Injury B2 HTO dimensions:", dim(InjuryB2.htos), "\n")
cat("HTO names:\n")
print(head(InjuryB2.htos, n = 4))
print(rownames(InjuryB2.htos))

# Join the RNA and HTO
# Rename UMI columns to match HTO
cat("\nBefore modification:\n")
cat("cDNA barcodes:", head(colnames(InjuryB2.umis)), "\n")
cat("HTO barcodes:", head(colnames(InjuryB2.htos)), "\n")

colnames(InjuryB2.umis) <- gsub("-1", "", colnames(InjuryB2.umis))

cat("After modification:\n")
cat("cDNA barcodes:", tail(colnames(InjuryB2.umis)), "\n")

# Join the hashtag and cDNA data
joint.bcs_injuryB2 <- dplyr::intersect(colnames(InjuryB2.htos), colnames(InjuryB2.umis))
cat("Joint barcodes found:", length(joint.bcs_injuryB2), "\n")

# Subset RNA and HTO counts by joint cell barcodes
InjuryB2.umis <- InjuryB2.umis[, joint.bcs_injuryB2]
InjuryB2.htos <- as.matrix(InjuryB2.htos[, joint.bcs_injuryB2])
cat("Final Injury B2 HTO dimensions:", dim(InjuryB2.htos), "\n")
cat("HTO row names confirmed:\n")
print(rownames(InjuryB2.htos))

# Setup Seurat object
InjuryB2_hashtag <- CreateSeuratObject(counts = InjuryB2.umis, min.cells = 3)

# Normalize RNA data with log normalization
InjuryB2_hashtag <- NormalizeData(InjuryB2_hashtag)

# Find and scale variable features
InjuryB2_hashtag <- FindVariableFeatures(InjuryB2_hashtag, selection.method = "mean.var.plot")
InjuryB2_hashtag <- ScaleData(InjuryB2_hashtag, features = VariableFeatures(InjuryB2_hashtag))

# Add HTO data as a new assay independent from RNA
InjuryB2_hashtag[["HTO"]] <- CreateAssayObject(counts = InjuryB2.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
InjuryB2_hashtag <- NormalizeData(InjuryB2_hashtag, assay = "HTO", 
                                  normalization.method = "CLR", margin = 2)

# Hashtag demultiplexing
InjuryB2_hashtag <- MULTIseqDemux(InjuryB2_hashtag, assay = "HTO", autoThresh = TRUE)
cat("\nInjury B2 demultiplexing results:\n")
print(table(InjuryB2_hashtag$MULTI_ID))

# Group cells based on the max HTO signal
Idents(InjuryB2_hashtag) <- "MULTI_ID"

# Visualizations
VlnPlot(InjuryB2_hashtag, assay = "HTO", 
        features = rownames(InjuryB2_hashtag[["HTO"]])[1:4], ncol = 2)
RidgePlot(InjuryB2_hashtag, assay = "HTO", 
          features = rownames(InjuryB2_hashtag[["HTO"]])[1:4], ncol = 2)
VlnPlot(InjuryB2_hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# Create HTO heatmap
# 1. Get the normalized HTO data
hto_data_InjuryB2_hashtag <- GetAssayData(InjuryB2_hashtag, assay = "HTO", slot = "data")

# 2. Sample cells (to avoid plotting too many)
set.seed(42) # For reproducibility
cells_to_plot_InjuryB2_hashtag <- sample(colnames(InjuryB2_hashtag), 
                                         min(5000, ncol(InjuryB2_hashtag)))

# 3. Order cells by their classification
cell_classifications_InjuryB2_hashtag <- InjuryB2_hashtag$MULTI_ID[cells_to_plot_InjuryB2_hashtag]
cells_ordered_InjuryB2_hashtag <- cells_to_plot_InjuryB2_hashtag[order(cell_classifications_InjuryB2_hashtag)]

# 4. Create the heatmap
hto_data_subset_InjuryB2_hashtag <- hto_data_InjuryB2_hashtag[, cells_ordered_InjuryB2_hashtag]
pheatmap(hto_data_subset_InjuryB2_hashtag, cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE, 
         annotation_col = data.frame(ID = InjuryB2_hashtag$MULTI_ID[cells_ordered_InjuryB2_hashtag],
                                     row.names = cells_ordered_InjuryB2_hashtag))

################################################################################
# REMOVE OVERLAPPING BARCODES
################################################################################

cat("\n=== Cleaning Data: Removing Overlapping Barcodes ===\n")

# First, identify overlapping barcodes
sham_barcodes <- colnames(ShamB2.hashtag)
injury_barcodes <- colnames(InjuryB2_hashtag)
overlapping_barcodes <- intersect(sham_barcodes, injury_barcodes)

cat("Found", length(overlapping_barcodes), "overlapping barcodes\n")

# Remove overlapping cells from injured dataset
# This ensures completely independent datasets
cat("\n=== Remove from Injured datasets ===\n")
Sham.singletB2.clean <- ShamB2.hashtag
Injury.singletB2.clean <- InjuryB2_hashtag[, !colnames(InjuryB2_hashtag) %in% overlapping_barcodes]

# Verify that all overlapping cell barcodes are removed
sham_barcodes <- colnames(Sham.singletB2.clean)
injury_barcodes <- colnames(Injury.singletB2.clean)
overlapping_barcodes <- intersect(sham_barcodes, injury_barcodes)

cat("After cleanup, found", length(overlapping_barcodes), "overlapping barcodes\n")

################################################################################
# EXTRACT SINGLETS AND PERFORM QC - SHAM SAMPLES
################################################################################

cat("\n=== Processing Sham B2 Singlets ===\n")

# Extract the singlets (remove doublets and negatives)
Sham.singletB2 <- subset(Sham.singletB2.clean, idents = c("Doublet", "Negative"), invert = TRUE)

# Rename the project
current.project.id <- 'SeuratProject'
new.project.id <- 'Sham.singletB2'
Sham.singletB2@meta.data$orig.ident <- plyr::mapvalues(x = Sham.singletB2@meta.data$orig.ident, 
                                                       from = current.project.id, 
                                                       to = new.project.id)

# Calculate mitochondrial percentage
Sham.singletB2[["percent.mt"]] <- PercentageFeatureSet(Sham.singletB2, pattern = "^mt-")
VlnPlot(Sham.singletB2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
cat("Cells before QC filtering:", length(WhichCells(Sham.singletB2)), "\n")

# Calculate means and SDs for QC filtering
feature_mean <- mean(Sham.singletB2$nFeature_RNA)
count_mean <- mean(Sham.singletB2$nCount_RNA)
feature_sd <- sd(Sham.singletB2$nFeature_RNA) 
count_sd <- sd(Sham.singletB2$nCount_RNA)

# Subset using 2 SD cutoffs
Sham.singletB2 <- subset(Sham.singletB2, 
                         nFeature_RNA > feature_mean - 2*feature_sd &
                           nCount_RNA > count_mean - 2*count_sd &
                           percent.mt < 5)

# Visualize results
VlnPlot(Sham.singletB2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, group.by = "orig.ident")
cat("Cells after QC filtering:", length(WhichCells(Sham.singletB2)), "\n")
cat("Sham B2 samples by HTO:\n")
print(table(Sham.singletB2$MULTI_ID))

################################################################################
# EXTRACT SINGLETS AND PERFORM QC - INJURED SAMPLES
################################################################################

cat("\n=== Processing Injured B2 Singlets ===\n")

# Extract the singlets (remove doublets and negatives)
InjuryB2_hashtag_singlet <- subset(Injury.singletB2.clean, idents = c("Doublet", "Negative"), invert = TRUE)

# Rename the project
current.project.id <- 'SeuratProject'
new.project.id <- 'Injury.singletB2'
InjuryB2_hashtag_singlet@meta.data$orig.ident <- plyr::mapvalues(x = InjuryB2_hashtag_singlet@meta.data$orig.ident, 
                                                                 from = current.project.id, 
                                                                 to = new.project.id)

# Calculate mitochondrial percentage
InjuryB2_hashtag_singlet[["percent.mt"]] <- PercentageFeatureSet(InjuryB2_hashtag_singlet, pattern = "^mt-")
VlnPlot(InjuryB2_hashtag_singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
cat("Cells before QC filtering:", length(WhichCells(InjuryB2_hashtag_singlet)), "\n")
VlnPlot(InjuryB2_hashtag_singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, group.by = "orig.ident")

# Calculate means and SDs for QC filtering
feature_mean <- mean(InjuryB2_hashtag_singlet$nFeature_RNA)
count_mean <- mean(InjuryB2_hashtag_singlet$nCount_RNA)
feature_sd <- sd(InjuryB2_hashtag_singlet$nFeature_RNA) 
count_sd <- sd(InjuryB2_hashtag_singlet$nCount_RNA)

# Subset using 2 SD cutoffs
Injury.singletB2 <- subset(InjuryB2_hashtag_singlet, 
                           nFeature_RNA > feature_mean - 2*feature_sd &
                             nCount_RNA > count_mean - 2*count_sd &
                             percent.mt < 5)

# Visualize results
VlnPlot(Injury.singletB2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, group.by = "orig.ident")
VlnPlot(Injury.singletB2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
cat("Cells after QC filtering:", length(WhichCells(Injury.singletB2)), "\n")
cat("Injured B2 samples by HTO:\n")
print(table(Injury.singletB2$MULTI_ID))

cat("\n=== Processing Complete ===\n")
cat("Final objects created:\n")
cat("  - Sham.singletB2: QC-filtered sham B2 singlets\n")
cat("  - Injury.singletB2: QC-filtered injured B2 singlets\n")
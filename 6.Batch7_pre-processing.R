################################################################################
# Single-Cell RNA-seq Analysis with Cell Hashing - 2 Hours Post Injury - Batch 7
# 
# Description: Processing and demultiplexing of single-cell RNA-seq data using
# hashtag oligonucleotides (HTOs) for sham and 2hpi (2 hours post injury) samples
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
sham_for_2hpi_cdna_path <- "Data_Files/cDNA/ex028_Sham_filtered_feature_bc_matrix.h5"
two_hpi_cdna_path <- "Data_Files/cDNA/ex028_2hpis_filtered_feature_bc_matrix.h5"

# Path to HTO count data (kallisto/BUStools output)
sham_for_2hpi_hto_path <- "Data_Files/HTOs/ex028_sham_counts_unfiltered"
two_hpi_hto_path <- "Data_Files/HTOs/ex028_Undetermined_counts_unfiltered"

################################################################################
# PROCESS SHAM SAMPLES (FOR 2HPI COMPARISON)
################################################################################

cat("=== Processing Sham Samples (for 2hpi comparison) ===\n")

# Load cDNA data
Sham_for_TwoHpi_umis <- Read10X_h5(sham_for_2hpi_cdna_path)
cat("Sham for 2hpi cDNA dimensions:", dim(Sham_for_TwoHpi_umis), "\n")

# Load HTO data
Sham_for_TwoHpi_htos <- BUSpaRse::read_count_output(sham_for_2hpi_hto_path, 
                                                    name = "cells_x_features", 
                                                    tcc = FALSE)
cat("Sham for 2hpi HTO dimensions:", dim(Sham_for_TwoHpi_htos), "\n")
cat("HTO names:\n")
print(head(Sham_for_TwoHpi_htos, n = 6))
print(rownames(Sham_for_TwoHpi_htos))

# Join the RNA and HTO
# Rename UMI columns to match HTO (remove -1 suffix)
cat("\nBefore modification:\n")
cat("cDNA barcodes:", head(colnames(Sham_for_TwoHpi_umis)), "\n")
cat("HTO barcodes:", head(colnames(Sham_for_TwoHpi_htos)), "\n")

colnames(Sham_for_TwoHpi_umis) <- gsub("-1", "", colnames(Sham_for_TwoHpi_umis))

cat("After modification:\n")
cat("cDNA barcodes:", tail(colnames(Sham_for_TwoHpi_umis)), "\n")

# Overlap cDNA and HTO data
joint.bcs_sham_for_TwoHpi <- dplyr::intersect(colnames(Sham_for_TwoHpi_htos), 
                                              colnames(Sham_for_TwoHpi_umis))
cat("Joint barcodes found:", length(joint.bcs_sham_for_TwoHpi), "\n")

# Subset RNA and HTO counts by joint cell barcodes
Sham_for_TwoHpi_umis <- Sham_for_TwoHpi_umis[, joint.bcs_sham_for_TwoHpi]
Sham_for_TwoHpi_htos <- as.matrix(Sham_for_TwoHpi_htos[, joint.bcs_sham_for_TwoHpi])
cat("Final Sham for 2hpi HTO dimensions:", dim(Sham_for_TwoHpi_htos), "\n")
cat("HTO row names confirmed:\n")
print(rownames(Sham_for_TwoHpi_htos))

# Setup Seurat object
Sham_for_two_hpi <- CreateSeuratObject(counts = Sham_for_TwoHpi_umis, min.cells = 3)

# Normalize RNA data with log normalization
Sham_for_two_hpi <- NormalizeData(Sham_for_two_hpi)

# Find and scale variable features
Sham_for_two_hpi <- FindVariableFeatures(Sham_for_two_hpi, selection.method = "mean.var.plot")
Sham_for_two_hpi <- ScaleData(Sham_for_two_hpi, features = VariableFeatures(Sham_for_two_hpi))

# Add HTO data as a new assay independent from RNA
Sham_for_two_hpi[["HTO"]] <- CreateAssayObject(counts = Sham_for_TwoHpi_htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
Sham_for_two_hpi <- NormalizeData(Sham_for_two_hpi, assay = "HTO", 
                                  normalization.method = "CLR", margin = 2)

# Demultiplex the data
Sham_for_two_hpi <- MULTIseqDemux(Sham_for_two_hpi, assay = "HTO", autoThresh = TRUE)
cat("\nSham for 2hpi demultiplexing results:\n")
print(table(Sham_for_two_hpi$MULTI_ID))

# Group cells based on the max HTO signal
Idents(Sham_for_two_hpi) <- "MULTI_ID"

# Visualizations
VlnPlot(Sham_for_two_hpi, assay = "HTO", 
        features = rownames(Sham_for_two_hpi[["HTO"]])[1:4], ncol = 2)
RidgePlot(Sham_for_two_hpi, assay = "HTO", 
          features = rownames(Sham_for_two_hpi[["HTO"]])[1:4], ncol = 2)
VlnPlot(Sham_for_two_hpi, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# Create HTO heatmap to check hashtag specificity
# 1. Get the normalized HTO data
hto_data_Sham_for_two_hpi <- GetAssayData(Sham_for_two_hpi, assay = "HTO", slot = "data")

# 2. Sample cells (to avoid plotting too many)
set.seed(42) # For reproducibility
cells_to_plot_Sham_for_two_hpi <- sample(colnames(Sham_for_two_hpi), 
                                         min(5000, ncol(Sham_for_two_hpi)))

# 3. Order cells by their classification
cell_classifications_Sham_for_two_hpi <- Sham_for_two_hpi$MULTI_ID[cells_to_plot_Sham_for_two_hpi]
cells_ordered_Sham_for_two_hpi <- cells_to_plot_Sham_for_two_hpi[order(cell_classifications_Sham_for_two_hpi)]

# 4. Create the heatmap
hto_data_subset_Sham_for_two_hpi <- hto_data_Sham_for_two_hpi[, cells_ordered_Sham_for_two_hpi]
pheatmap(hto_data_subset_Sham_for_two_hpi, cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE, 
         annotation_col = data.frame(ID = Sham_for_two_hpi$MULTI_ID[cells_ordered_Sham_for_two_hpi],
                                     row.names = cells_ordered_Sham_for_two_hpi))

################################################################################
# PROCESS 2HPI (2 HOURS POST INJURY) SAMPLES
################################################################################

cat("\n=== Processing 2hpi (2 Hours Post Injury) Samples ===\n")

# Load cDNA data
TwoHpi_umis <- Read10X_h5(two_hpi_cdna_path)
cat("2hpi cDNA dimensions:", dim(TwoHpi_umis), "\n")

# Load hashtags
TwoHpi_htos <- BUSpaRse::read_count_output(two_hpi_hto_path, 
                                           name = "cells_x_features", 
                                           tcc = FALSE)
cat("2hpi HTO dimensions:", dim(TwoHpi_htos), "\n")
cat("HTO names:\n")
print(head(TwoHpi_htos, n = 4))
print(rownames(TwoHpi_htos))

# Join the RNA and HTO
# Rename UMI columns to match HTO
cat("\nBefore modification:\n")
cat("cDNA barcodes:", head(colnames(TwoHpi_umis)), "\n")
cat("HTO barcodes:", head(colnames(TwoHpi_htos)), "\n")

colnames(TwoHpi_umis) <- gsub("-1", "", colnames(TwoHpi_umis))

cat("After modification:\n")
cat("cDNA barcodes:", tail(colnames(TwoHpi_umis)), "\n")

# Overlap cDNA and HTO data
joint.bcs_TwoHpi <- dplyr::intersect(colnames(TwoHpi_htos), colnames(TwoHpi_umis))
cat("Joint barcodes found:", length(joint.bcs_TwoHpi), "\n")

# Subset RNA and HTO counts by joint cell barcodes
TwoHpi_umis <- TwoHpi_umis[, joint.bcs_TwoHpi]
TwoHpi_htos <- as.matrix(TwoHpi_htos[, joint.bcs_TwoHpi])
cat("Final 2hpi HTO dimensions:", dim(TwoHpi_htos), "\n")
cat("HTO row names confirmed:\n")
print(rownames(TwoHpi_htos))

# Setup Seurat object
Injury_Two_hpi_hashtag <- CreateSeuratObject(counts = TwoHpi_umis, min.cells = 3)

# Normalize RNA data with log normalization
Injury_Two_hpi_hashtag <- NormalizeData(Injury_Two_hpi_hashtag)

# Find and scale variable features
Injury_Two_hpi_hashtag <- FindVariableFeatures(Injury_Two_hpi_hashtag, 
                                               selection.method = "mean.var.plot")
Injury_Two_hpi_hashtag <- ScaleData(Injury_Two_hpi_hashtag, 
                                    features = VariableFeatures(Injury_Two_hpi_hashtag))

# Add HTO data as a new assay independent from RNA
Injury_Two_hpi_hashtag[["HTO"]] <- CreateAssayObject(counts = TwoHpi_htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
Injury_Two_hpi_hashtag <- NormalizeData(Injury_Two_hpi_hashtag, assay = "HTO", 
                                        normalization.method = "CLR", margin = 2)

# Demultiplex the hashtags
Injury_Two_hpi_hashtag <- MULTIseqDemux(Injury_Two_hpi_hashtag, assay = "HTO", 
                                        autoThresh = TRUE)
cat("\n2hpi demultiplexing results:\n")
print(table(Injury_Two_hpi_hashtag$MULTI_ID))

# Group cells based on the max HTO signal
Idents(Injury_Two_hpi_hashtag) <- "MULTI_ID"

# Visualizations
VlnPlot(Injury_Two_hpi_hashtag, assay = "HTO", 
        features = rownames(Injury_Two_hpi_hashtag[["HTO"]])[1:4], ncol = 2)
RidgePlot(Injury_Two_hpi_hashtag, assay = "HTO", 
          features = rownames(Injury_Two_hpi_hashtag[["HTO"]])[1:4], ncol = 2)
VlnPlot(Injury_Two_hpi_hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# Create HTO heatmap
# 1. Get the normalized HTO data
hto_data_Injury_Two_hpi_hashtag <- GetAssayData(Injury_Two_hpi_hashtag, 
                                                assay = "HTO", slot = "data")

# 2. Sample cells (to avoid plotting too many)
set.seed(42) # For reproducibility
cells_to_plot_Injury_Two_hpi_hashtag <- sample(colnames(Injury_Two_hpi_hashtag), 
                                               min(5000, ncol(Injury_Two_hpi_hashtag)))

# 3. Order cells by their classification
cell_classifications_Injury_Two_hpi_hashtag <- Injury_Two_hpi_hashtag$MULTI_ID[cells_to_plot_Injury_Two_hpi_hashtag]
cells_ordered_Injury_Two_hpi_hashtag <- cells_to_plot_Injury_Two_hpi_hashtag[order(cell_classifications_Injury_Two_hpi_hashtag)]

# 4. Create the heatmap
hto_data_subset_Injury_Two_hpi_hashtag <- hto_data_Injury_Two_hpi_hashtag[, cells_ordered_Injury_Two_hpi_hashtag]
pheatmap(hto_data_subset_Injury_Two_hpi_hashtag, cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE, 
         annotation_col = data.frame(ID = Injury_Two_hpi_hashtag$MULTI_ID[cells_ordered_Injury_Two_hpi_hashtag],
                                     row.names = cells_ordered_Injury_Two_hpi_hashtag))

################################################################################
# REMOVE OVERLAPPING BARCODES
################################################################################

cat("\n=== Cleaning Data: Removing Overlapping Barcodes ===\n")

# First, identify overlapping barcodes between both datasets
sham_barcodes <- colnames(Sham_for_two_hpi)
injury_barcodes <- colnames(Injury_Two_hpi_hashtag)
overlapping_barcodes <- intersect(sham_barcodes, injury_barcodes)

cat("Found", length(overlapping_barcodes), "overlapping barcodes\n")

# Remove overlapping cells from injured dataset
Sham_for2hpi_hashtag.clean <- Sham_for_two_hpi
Two_hpi_hashtag.clean <- Injury_Two_hpi_hashtag[, !colnames(Injury_Two_hpi_hashtag) %in% overlapping_barcodes]

# Check that the duplicate cell barcodes have been removed
sham_barcodes <- colnames(Sham_for2hpi_hashtag.clean)
injury_barcodes <- colnames(Two_hpi_hashtag.clean)
overlapping_barcodes <- intersect(sham_barcodes, injury_barcodes)
cat("After cleanup, found", length(overlapping_barcodes), "overlapping barcodes\n")

################################################################################
# EXTRACT SINGLETS AND PERFORM QC - SHAM SAMPLES
################################################################################

cat("\n=== Processing Sham for 2hpi Singlets ===\n")

# Extract the singlets (remove doublets and negatives)
Sham_for_TwoHpi_singlet.clean <- subset(Sham_for2hpi_hashtag.clean, 
                                        idents = c("Doublet", "Negative"), invert = TRUE)

# Rename the project
current.project.id <- 'SeuratProject'
new.project.id <- 'Sham_for2hpi.singlet'
Sham_for_TwoHpi_singlet.clean@meta.data$orig.ident <- plyr::mapvalues(
  x = Sham_for_TwoHpi_singlet.clean@meta.data$orig.ident, 
  from = current.project.id, 
  to = new.project.id
)

# Calculate mitochondrial percentage
Sham_for_TwoHpi_singlet.clean[["percent.mt"]] <- PercentageFeatureSet(
  Sham_for_TwoHpi_singlet.clean, pattern = "^mt-"
)
VlnPlot(Sham_for_TwoHpi_singlet.clean, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
cat("Cells before QC filtering:", length(WhichCells(Sham_for_TwoHpi_singlet.clean)), "\n")

# Calculate means and SDs for QC filtering
feature_mean <- mean(Sham_for_TwoHpi_singlet.clean$nFeature_RNA)
count_mean <- mean(Sham_for_TwoHpi_singlet.clean$nCount_RNA)
feature_sd <- sd(Sham_for_TwoHpi_singlet.clean$nFeature_RNA) 
count_sd <- sd(Sham_for_TwoHpi_singlet.clean$nCount_RNA)

# Subset using 2 SD cutoffs
Sham_for2hpi_hashtag_singlets <- subset(Sham_for_TwoHpi_singlet.clean, 
                                        nFeature_RNA > feature_mean - 2*feature_sd &
                                          nCount_RNA > count_mean - 2*count_sd &
                                          percent.mt < 5)

# Visualize results
VlnPlot(Sham_for2hpi_hashtag_singlets, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, group.by = "orig.ident")
cat("Cells after QC filtering:", length(WhichCells(Sham_for2hpi_hashtag_singlets)), "\n")
cat("Sham for 2hpi samples by HTO:\n")
print(table(Sham_for2hpi_hashtag_singlets$MULTI_ID))

################################################################################
# EXTRACT SINGLETS AND PERFORM QC - 2HPI SAMPLES
################################################################################

cat("\n=== Processing 2hpi Singlets ===\n")

# Extract the singlets (remove doublets and negatives)
Two_hpi_hashtag_singlets.clean <- subset(Two_hpi_hashtag.clean, 
                                         idents = c("Doublet", "Negative"), invert = TRUE)

# Rename the project
current.project.id <- 'SeuratProject'
new.project.id <- '2hpi.singlet'
Two_hpi_hashtag_singlets.clean@meta.data$orig.ident <- plyr::mapvalues(
  x = Two_hpi_hashtag_singlets.clean@meta.data$orig.ident, 
  from = current.project.id, 
  to = new.project.id
)

# Calculate mitochondrial percentage
Two_hpi_hashtag_singlets.clean[["percent.mt"]] <- PercentageFeatureSet(
  Two_hpi_hashtag_singlets.clean, pattern = "^mt-"
)
VlnPlot(Two_hpi_hashtag_singlets.clean, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
cat("Cells before QC filtering:", length(WhichCells(Two_hpi_hashtag_singlets.clean)), "\n")

# Calculate means and SDs for QC filtering
feature_mean <- mean(Two_hpi_hashtag_singlets.clean$nFeature_RNA)
count_mean <- mean(Two_hpi_hashtag_singlets.clean$nCount_RNA)
feature_sd <- sd(Two_hpi_hashtag_singlets.clean$nFeature_RNA) 
count_sd <- sd(Two_hpi_hashtag_singlets.clean$nCount_RNA)

# Subset using 2 SD cutoffs
Two_hpi_hashtag_singlets.clean <- subset(Two_hpi_hashtag_singlets.clean, 
                                         nFeature_RNA > feature_mean - 2*feature_sd &
                                           nCount_RNA > count_mean - 2*count_sd &
                                           percent.mt < 5)

# Visualize results
VlnPlot(Two_hpi_hashtag_singlets.clean, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, group.by = "orig.ident")
cat("Cells after QC filtering:", length(WhichCells(Two_hpi_hashtag_singlets.clean)), "\n")
cat("2hpi samples by HTO:\n")
print(table(Two_hpi_hashtag_singlets.clean$MULTI_ID))

cat("\n=== Processing Complete ===\n")
cat("Final objects created:\n")
cat("  - Sham_for2hpi_hashtag_singlets: QC-filtered sham singlets for 2hpi comparison\n")
cat("  - Two_hpi_hashtag_singlets.clean: QC-filtered 2hpi singlets\n")
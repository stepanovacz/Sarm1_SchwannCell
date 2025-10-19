################################################################################
# Single-Cell RNA-seq Analysis with Cell Hashing - Naive Samples
# 
# Description: Processing and demultiplexing of naive single-cell RNA-seq data 
# using hashtag oligonucleotides (HTOs)
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

# Path to cDNA data (10X H5 file)
naive_cdna_path <- "Data_Files/cDNA/KS220914_Naives_filtered_feature_bc_matrix.h5"

# Path to HTO count data (kallisto/BUStools output)
naive_hto_path <- "Data_Files/HTOs/KS220914_sham_htos"

################################################################################
# LOAD DATA
################################################################################

cat("=== Loading Naive Sample Data ===\n")

# Load cDNA data
Naive.umis <- Read10X_h5(naive_cdna_path)
cat("Naive cDNA dimensions:", dim(Naive.umis), "\n")

# Load HTO data
Naive.htos <- BUSpaRse::read_count_output(naive_hto_path, 
                                          name = "cells_x_features", 
                                          tcc = FALSE)
cat("Naive HTO dimensions:", dim(Naive.htos), "\n")
cat("HTO names:\n")
print(head(Naive.htos, n = 6))
print(rownames(Naive.htos))

################################################################################
# JOIN RNA AND HTO DATA
################################################################################

cat("\n=== Joining RNA and HTO Data ===\n")

# Rename UMI columns to match HTO (remove -1 suffix)
cat("Before modification:\n")
cat("cDNA barcodes:", head(colnames(Naive.umis)), "\n")
cat("HTO barcodes:", head(colnames(Naive.htos)), "\n")

colnames(Naive.umis) <- gsub("-1", "", colnames(Naive.umis))

cat("\nAfter modification:\n")
cat("cDNA barcodes:", tail(colnames(Naive.umis)), "\n")

# Find joint barcodes
joint.bcs_Naive <- dplyr::intersect(colnames(Naive.htos), colnames(Naive.umis))
cat("Joint barcodes found:", length(joint.bcs_Naive), "\n")

# Subset RNA and HTO counts by joint cell barcodes
Naive.umis <- Naive.umis[, joint.bcs_Naive]
Naive.htos <- as.matrix(Naive.htos[, joint.bcs_Naive])
cat("Final HTO dimensions:", dim(Naive.htos), "\n")
cat("HTO row names confirmed:\n")
print(rownames(Naive.htos))

################################################################################
# CREATE AND PROCESS SEURAT OBJECT
################################################################################

cat("\n=== Creating Seurat Object ===\n")

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
Naive.hashtag <- NormalizeData(Naive.hashtag, assay = "HTO", 
                               normalization.method = "CLR", margin = 2)

################################################################################
# HASHTAG DEMULTIPLEXING
################################################################################

cat("\n=== Performing Hashtag Demultiplexing ===\n")

# Demultiplex the hashtags
Naive.hashtag <- MULTIseqDemux(Naive.hashtag, assay = "HTO", autoThresh = TRUE)

cat("Demultiplexing results:\n")
print(table(Naive.hashtag$MULTI_ID))

################################################################################
# VISUALIZE DEMULTIPLEXING RESULTS
################################################################################

# Group cells based on the max HTO signal
Idents(Naive.hashtag) <- "MULTI_ID"

# Violin plots
VlnPlot(Naive.hashtag, assay = "HTO", 
        features = rownames(Naive.hashtag[["HTO"]])[1:4], ncol = 2)

# Ridge plots
RidgePlot(Naive.hashtag, assay = "HTO", 
          features = rownames(Naive.hashtag[["HTO"]])[1:4], ncol = 2)

# Count plot
VlnPlot(Naive.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# Create HTO heatmap
# 1. Get the normalized HTO data
hto_data <- GetAssayData(Naive.hashtag, assay = "HTO", slot = "data")

# 2. Sample cells (to avoid plotting too many)
set.seed(42) # For reproducibility
cells_to_plot <- sample(colnames(Naive.hashtag), min(5000, ncol(Naive.hashtag)))

# 3. Order cells by their classification
cell_classifications <- Naive.hashtag$MULTI_ID[cells_to_plot]
cells_ordered <- cells_to_plot[order(cell_classifications)]

# 4. Create the heatmap
hto_data_subset <- hto_data[, cells_ordered]
pheatmap(hto_data_subset, cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE, 
         annotation_col = data.frame(ID = Naive.hashtag$MULTI_ID[cells_ordered],
                                     row.names = cells_ordered))

################################################################################
# EXTRACT SINGLETS AND PERFORM QC
################################################################################

cat("\n=== Extracting Singlets and Performing QC ===\n")

# Extract the singlets (remove doublets and negatives)
Naive.singlet <- subset(Naive.hashtag, idents = c("Doublet", "Negative"), invert = TRUE)

# Rename the project
current.project.id <- 'SeuratProject'
new.project.id <- 'Naive'
Naive.singlet@meta.data$orig.ident <- plyr::mapvalues(x = Naive.singlet@meta.data$orig.ident, 
                                                      from = current.project.id, 
                                                      to = new.project.id)

# Calculate mitochondrial percentage
Naive.singlet[["percent.mt"]] <- PercentageFeatureSet(Naive.singlet, pattern = "^mt-")
VlnPlot(Naive.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
cat("Cells before QC filtering:", length(WhichCells(Naive.singlet)), "\n")

# Calculate means and SDs for QC filtering
feature_mean <- mean(Naive.singlet$nFeature_RNA)
count_mean <- mean(Naive.singlet$nCount_RNA)
feature_sd <- sd(Naive.singlet$nFeature_RNA) 
count_sd <- sd(Naive.singlet$nCount_RNA)

# Subset using 2 SD cutoffs
Naive.singlet <- subset(Naive.singlet, 
                        nFeature_RNA > feature_mean - 2*feature_sd &
                          nCount_RNA > count_mean - 2*count_sd &
                          percent.mt < 5)

# Visualize results
VlnPlot(Naive.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, group.by = "orig.ident")
cat("Cells after QC filtering:", length(WhichCells(Naive.singlet)), "\n")
cat("Naive samples by HTO:\n")
print(table(Naive.singlet$MULTI_ID))

cat("\n=== Processing Complete ===\n")
cat("Final object created: Naive.singlet (QC-filtered naive singlets)\n")

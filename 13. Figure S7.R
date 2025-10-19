################################################################################
# Bulk RNA-seq Analysis: Neuregulin Treatment in WT and SARM1-KO Samples
# 
# Description: Differential expression analysis comparing neuregulin treatment
# effects in wild-type and SARM1-knockout samples at 24 hours
#
# Requirements: 
# - Salmon quantification files in subdirectories (named *_quant/)
# - tx2gene mapping file (transcript ID to gene ID)
# - Pseudobulk results file (for comparison analysis)
# - R packages: tximport, DESeq2, AnnotationDbi, org.Mm.eg.db, 
#   pheatmap, viridis, reshape2, dplyr, ggplot2, gridExtra
################################################################################

# Load required libraries
library(tximport)
library(DESeq2)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(pheatmap)
library(viridis)
library(reshape2)
library(dplyr)
library(ggplot2)
library(gridExtra)

################################################################################
# USER-DEFINED PATHS - MODIFY THESE FOR YOUR SYSTEM
################################################################################

# Path to tx2gene file (transcript to gene ID mapping)
tx2gene_path <- "path/to/tx2gene_corrected.csv"
# Output directory for results
output_dir <- "./results"

# Path to pseudobulk results (for comparison analysis) from snRNAseq
pseudobulk_path <- "path/to/pseudobulk_results.csv"
# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Set working directory to your quants folder (or modify path below)
setwd("path/to/quants/folder")
# List all directories in the quants folder
sample_dirs <- list.dirs(path = ".", recursive = FALSE)

# Exclude the results directory
sample_dirs <- sample_dirs[!grepl("results$", sample_dirs)]

print("Found these sample directories:")
print(sample_dirs)

# Get paths to all quant.sf files
quant_files <- file.path(sample_dirs, "quant.sf")
print("\nConstructed these paths to quant.sf files:")
print(quant_files)

# Read tx2gene file
tx2gene <- read.csv(tx2gene_path, 
                    header = TRUE, 
                    col.names = c("transcript_id", "gene_id"))
# Look at the first few rows
print("First few rows of tx2gene:")
print(head(tx2gene))

# Add .1 to transcript IDs (since they need versions to match quant.sf)
tx2gene$transcript_id <- paste0(tx2gene$transcript_id, ".1")

# Check the modification
print("First few rows of modified tx2gene:")
print(head(tx2gene))

# Let's first look at just the clean sample names
sample_names <- gsub("_quant$", "", basename(sample_dirs))
print("Sample names:")
print(sample_names)

# Get clean sample names
sample_names <- gsub("_quant$", "", basename(sample_dirs))

# Create metadata with correct pattern matching
sample_metadata <- data.frame(
  sample = sample_names,
  # Match "Sarm1-KO" instead of "SARM1-KO"
  genotype = factor(ifelse(startsWith(sample_names, "Sarm1"), "SARM1-KO", "WT")),
  # Extract treatment
  treatment = factor(ifelse(grepl("[Pp]lus", sample_names), "plus_NRG", "minus_NRG")),
  # Extract batch number
  batch = factor(gsub(".*Batch-([0-9]+).*", "\\1", sample_names))
)

rownames(sample_metadata) <- sample_metadata$sample

print("\nSample metadata:")
print(sample_metadata)

print("\nNumber of samples per condition:")
print(table(sample_metadata$genotype, sample_metadata$treatment))

# First, confirm our tx2gene and quant files are ready
print("Number of transcript-gene mappings:")
print(nrow(tx2gene))

print("\nChecking first quant file exists:")
print(file.exists(quant_files[1]))

# Remove version numbers from tx2gene
tx2gene$transcript_id <- gsub("\\.\\d+$", "", tx2gene$transcript_id)

# Now import with tximport, ignoring transcript versions
txi <- tximport(quant_files, 
                type = "salmon", 
                tx2gene = tx2gene, 
                ignoreTxVersion = TRUE,  # Now ignore versions completely
                ignoreAfterBar = TRUE)

# Check the dimensions
print("\nDimensions of imported count data:")
dim(txi$counts)

# Look at the first few genes and their counts
print("\nFirst few genes and their counts:")
head(txi$counts)

##
# First, build our DESeq object with an interaction term
dds <- DESeqDataSetFromTximport(txi,
                                colData = sample_metadata,
                                design = ~ batch + genotype + treatment + genotype:treatment)

# Filter and run DESeq
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

# Get results for NRG effect in WT
res_nrg_wt <- results(dds, contrast=list(c("treatment_plus_NRG_vs_minus_NRG")))
head(res_nrg_wt)

# Get results for NRG effect in SARM1-KO
res_nrg_sarm <- results(dds, contrast=list(
  c("treatment_plus_NRG_vs_minus_NRG",
    "genotypeWT.treatmentplus_NRG")))
head(res_nrg_sarm)

# Order by adjusted p-value
res_nrg_wt_ordered <- res_nrg_wt[order(res_nrg_wt$padj),]
res_nrg_sarm_ordered <- res_nrg_sarm[order(res_nrg_sarm$padj),]
length(res_nrg_wt_ordered)
length(res_nrg_sarm_ordered)
# Print summaries
print("NRG effect in WT summary:")
summary(res_nrg_wt)

print("\nNRG effect in SARM1-KO summary:")
summary(res_nrg_sarm)

# Look at top genes
print("\nTop 10 genes affected by NRG in WT:")
head(res_nrg_wt_ordered, 10)

print("\nTop 10 genes affected by NRG in SARM1-KO:")
head(res_nrg_sarm_ordered, 10)

##Differential expression analysis
# Correct the subsetting for each genotype
wt_samples <- sample_metadata$genotype == "WT"
ko_samples <- sample_metadata$genotype == "SARM1-KO"  # Note the hyphen

# Create txi list for WT samples
txi_wt <- list(
  abundance = txi$abundance[, wt_samples],
  counts = txi$counts[, wt_samples],
  length = txi$length[, wt_samples],
  countsFromAbundance = txi$countsFromAbundance
)

# Create txi list for KO samples
txi_ko <- list(
  abundance = txi$abundance[, ko_samples],
  counts = txi$counts[, ko_samples],
  length = txi$length[, ko_samples],
  countsFromAbundance = txi$countsFromAbundance
)
# Create DESeq datasets for WT and KO samples
#WT data
dds_wt <- DESeqDataSetFromTximport(txi_wt,
                                   colData = sample_metadata[wt_samples,],
                                   design = ~ treatment)
#filter out low counts for WT
#smallestGroupSize <- 3
keep <- rowSums(counts(dds_wt)) >= 3
dds_wt <- dds_wt[keep,]
head(dds_wt)
#KO data

dds_ko <- DESeqDataSetFromTximport(txi_ko,
                                   colData = sample_metadata[ko_samples,],
                                   design = ~ treatment)
levels(dds_ko$treatment)
#filter out low counts for KO
#smallestGroupSize <- 3
keep <- rowSums(counts(dds_ko)) >= 3
dds_ko <- dds_ko[keep,]
head(dds_ko)


#Differential expression analysys for WT
WT_dds <- DESeq(dds_wt)
KO_dds <- DESeq(dds_ko)

#results. here comparign treatment : plus/minus NRG since
res_wt <- results(WT_dds, alpha=0.1)
res_wt #check that the right stuff was compared
res_ko <- results(KO_dds, alpha=0.1)
res_ko #check that the right stuff was compared



res_wt <- results(WT_dds, contrast=c("treatment", "plus_NRG", "minus_NRG"))
res_ko <- results(KO_dds, contrast=c("treatment", "plus_NRG", "minus_NRG"))
head(res_ko)
head(res_wt)
resultsNames(WT_dds)
resultsNames(KO_dds)

######
res_treatment <- results(dds, contrast=c("treatment", "plus_NRG", "minus_NRG"))
res_treatment_ordered <- res_treatment[order(res_treatment$padj),]
head(res_treatment_ordered)
# Map genes to symbols
gene_symbols <- mapIds(org.Mm.eg.db,
                       keys = rownames(res_treatment_ordered),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Check what we got
print("First few gene symbols:")
head(gene_symbols)

# Prepare data for heatmap
vsd <- vst(dds, blind=FALSE)
top_genes <- head(order(res_treatment$padj), 50)
mat <- assay(vsd)[top_genes, ]

# Create more informative row names (combining symbol and ID)
row_labels <- ifelse(!is.na(gene_symbols[rownames(mat)]),
                     paste0(gene_symbols[rownames(mat)], " (", rownames(mat), ")"),
                     rownames(mat))
rownames(mat) <- row_labels

# Create annotation dataframe
annotation_col <- data.frame(
  Genotype = dds$genotype,
  Treatment = dds$treatment,
  row.names = colnames(mat)
)

# Create heatmap
# Define your genes of interest


# Define genes of interest and sort alphabetically
#Same genes as detected in KEGG pseudobulk analysis in single-nucleus RNA sequencing in Sarm1KO Schwann Cells 1 day after entry. 
genes_of_interest <- c("Cycs","Cox7b2","Cox5b","Ndufc1","Ndufa4l2","Cox7a1","Uqcr11","Ndufb3","Uqcr10","Ndufb11","Ndufs8","Ndufb5","Atp6v0a4",
                       "Ndufb8", "Uqcrq", "Ndufb4", "Atp6v1g1" ,"Ndufs6",  "Ndufs2", "Cox6b2"  , "Ndufv1"  , "Cox6b1"  , "Cox7a2" ,  "Atp6v1f" , "Ndufa12" , "Cox7b"  , 
                       "Cox4i1" , "Cox8a"  ,  "Ndufa13" , "Ndufa4" ,  "Uqcrfs1" , "Ppa1" ,    "Atp6v0c" , "Cox6a1" ,  "Ndufa1" ,  "Ndufb10" , "Ndufb2" ,  "Cox6c"  ,  "Cox11"  , 
                       "Ndufb9" ,  "Ndufa11" , "Ndufs7" ,  "Ndufa9" ,  "Ndufb6" ,  "Ndufa2" ,  "Ndufa5"  , "Uqcrc2"  , "Ndufb7" ,  "Ndufa8"  , "Cox5a" ,   "Ndufa7" ,  "Cyc1"   , 
                       "Atp6v1d",  "Ndufa10" , "Ndufab1" , "Cox17"  ,  "Uqcrb"  ,  "Ndufs5" ,  "Ndufc2",   "Atp6v0b" , "Atp6v0e2" ,"Atp6v1e1" ,"Cox10"   , "Ndufa6" ,  "Ndufs1" , 
                       "Ndufs3" ,  "Atp6v1c1", "Lhpp" ,    "Uqcrh"   , "Cox7a2l" , "Atp6v1g2")


# Sort genes alphabetically
genes_of_interest_sorted <- sort(genes_of_interest)

# Get ENSEMBL IDs for these genes (maintain the sorted order)
ensembl_ids_sorted <- names(gene_symbols)[match(genes_of_interest_sorted, gene_symbols)]

# Remove any NA values (genes not found)
valid_indices <- !is.na(ensembl_ids_sorted)
ensembl_ids_sorted <- ensembl_ids_sorted[valid_indices]
genes_of_interest_sorted <- genes_of_interest_sorted[valid_indices]

# Create matrix of fold changes for selected genes
fc_matrix <- matrix(0, nrow = length(ensembl_ids_sorted), ncol = 2)
colnames(fc_matrix) <- c("WT", "SARM1-KO")
rownames(fc_matrix) <- ensembl_ids_sorted

# Fill in fold changes
fc_matrix[, "WT"] <- res_wt$log2FoldChange[match(ensembl_ids_sorted, rownames(res_wt))]
fc_matrix[, "SARM1-KO"] <- res_ko$log2FoldChange[match(ensembl_ids_sorted, rownames(res_ko))]

# Replace row names with gene symbols (now in alphabetical order)
rownames(fc_matrix) <- genes_of_interest_sorted

# Create annotation
annotation_col_mean <- data.frame(
  Condition = c("WT", "SARM1-KO")
)
rownames(annotation_col_mean) <- colnames(fc_matrix)

# Verify the order
print(head(rownames(fc_matrix), 10))

# Convert fc_matrix to long format for ggplot
fc_long <- melt(fc_matrix, varnames = c("Gene", "Sample"), value.name = "LogFC")

#Figure S7D
ggplot(fc_long, aes(x = Sample, y = Gene, fill = LogFC)) +
  geom_tile(color = "white", size = 0.6) +
  scale_fill_gradient2(
    low = "#FED976",     
    mid = "white",       
    high = "#08519C",     
    midpoint = 0,
    name = "Log2FC"
  ) +
  labs(x = "Sample", 
       y = "Gene",
       title = "Log2 Fold Changes after 24hr Neuregulin treatment in Mitochondrial Respiration Genes") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  ) +
  # Keep the original order (no clustering)
  scale_x_discrete(limits = colnames(fc_matrix)) +
  scale_y_discrete(limits = rev(rownames(fc_matrix)))


#Prepare data for Supplemental Figure 7E

# First, get gene symbols for all genes in both result sets
gene_symbols_wt <- mapIds(org.Mm.eg.db,
                          keys = rownames(res_wt),
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")

gene_symbols_ko <- mapIds(org.Mm.eg.db,
                          keys = rownames(res_ko),
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")

# Create comprehensive dataframes for WT results
wt_results_df <- data.frame(
  ENSEMBL_ID = rownames(res_wt),
  GeneSymbol = gene_symbols_wt[rownames(res_wt)],
  log2FoldChange = res_wt$log2FoldChange,
  padj = res_wt$padj,
  pvalue = res_wt$pvalue,
  baseMean = res_wt$baseMean,
  lfcSE = res_wt$lfcSE,
  stat = res_wt$stat,
  stringsAsFactors = FALSE
)

# Create comprehensive dataframes for SARM1-KO results
ko_results_df <- data.frame(
  ENSEMBL_ID = rownames(res_ko),
  GeneSymbol = gene_symbols_ko[rownames(res_ko)],
  log2FoldChange = res_ko$log2FoldChange,
  padj = res_ko$padj,
  pvalue = res_ko$pvalue,
  baseMean = res_ko$baseMean,
  lfcSE = res_ko$lfcSE,
  stat = res_ko$stat,
  stringsAsFactors = FALSE
)

# Handle missing gene symbols (replace NA with ENSEMBL ID)
wt_results_df$GeneSymbol[is.na(wt_results_df$GeneSymbol)] <- wt_results_df$ENSEMBL_ID[is.na(wt_results_df$GeneSymbol)]
ko_results_df$GeneSymbol[is.na(ko_results_df$GeneSymbol)] <- ko_results_df$ENSEMBL_ID[is.na(ko_results_df$GeneSymbol)]

# Order by adjusted p-value
wt_results_df <- wt_results_df[order(wt_results_df$padj, na.last = TRUE), ]
ko_results_df <- ko_results_df[order(ko_results_df$padj, na.last = TRUE), ]

# Print summaries
cat("WT Results Summary:\n")
cat("Total genes:", nrow(wt_results_df), "\n")
cat("Significant genes (padj < 0.05):", sum(wt_results_df$padj < 0.05, na.rm = TRUE), "\n")
cat("Significant genes (padj < 0.1):", sum(wt_results_df$padj < 0.1, na.rm = TRUE), "\n\n")

cat("SARM1-KO Results Summary:\n")
cat("Total genes:", nrow(ko_results_df), "\n")
cat("Significant genes (padj < 0.05):", sum(ko_results_df$padj < 0.05, na.rm = TRUE), "\n")
cat("Significant genes (padj < 0.1):", sum(ko_results_df$padj < 0.1, na.rm = TRUE), "\n\n")

# Show top 10 results
cat("Top 10 WT results:\n")
print(head(wt_results_df[, c("GeneSymbol", "log2FoldChange", "padj")], 10))

cat("\nTop 10 SARM1-KO results:\n")
print(head(ko_results_df[, c("GeneSymbol", "log2FoldChange", "padj")], 10))

# Write to CSV files
write.csv(wt_results_df, 
          file = file.path(output_dir, "WT_NRG_treatment_results_complete.csv"), 
          row.names = FALSE)

write.csv(ko_results_df, 
          file = file.path(output_dir, "SARM1KO_NRG_treatment_results_complete.csv"), 
          row.names = FALSE)

# Also create simplified versions with just the three columns you requested
wt_simple <- wt_results_df[, c("GeneSymbol", "log2FoldChange", "padj")]
ko_simple <- ko_results_df[, c("GeneSymbol", "log2FoldChange", "padj")]

write.csv(wt_simple, 
          file = file.path(output_dir, "WT_NRG_treatment_results_simple.csv"), 
          row.names = FALSE)

write.csv(ko_simple, 
          file = file.path(output_dir, "SARM1KO_NRG_treatment_results_simple.csv"), 
          row.names = FALSE)

# Read the data
pseudobulk <- read.csv(pseudobulk_path)
bulk_wt <- read.csv(file.path(output_dir, "WT_NRG_treatment_results_simple.csv"))
bulk_sarm <- read.csv(file.path(output_dir, "SARM1KO_NRG_treatment_results_simple.csv"))

# Process WT data
wt_pseudo <- pseudobulk %>%
  filter(group1 == "Sham-for-1dpi-WT",
         group2 == "1dpi-WT",
         gene %in% genes_of_interest) %>%
  dplyr::select(gene, avg_logFC, p_val_adj)

wt_bulk <- bulk_wt %>%
  filter(GeneSymbol %in% genes_of_interest) %>%
  dplyr::select(GeneSymbol, log2FoldChange, padj)

wt_combined <- wt_pseudo %>%
  inner_join(wt_bulk, by = c("gene" = "GeneSymbol")) %>%
  mutate(
    quadrant = case_when(
      avg_logFC >= 0 & log2FoldChange >= 0 ~ "Q1",
      avg_logFC < 0 & log2FoldChange >= 0 ~ "Q2",
      avg_logFC < 0 & log2FoldChange < 0 ~ "Q3",
      avg_logFC >= 0 & log2FoldChange < 0 ~ "Q4"
    ),
    dataset = "WT"
  )
head(wt_combined)

# Process Sarm1KO data
sarm_pseudo <- pseudobulk %>%
  filter(group1 == "Sham-for-1dpi-Sarm1-KO",
         group2 == "1dpi-Sarm1-KO",
         gene %in% genes_of_interest) %>%
  dplyr::select(gene, avg_logFC, p_val_adj)

sarm_bulk <- bulk_sarm %>%
  filter(GeneSymbol %in% genes_of_interest) %>%
  dplyr::select(GeneSymbol, log2FoldChange, padj)

sarm_combined <- sarm_pseudo %>%
  inner_join(sarm_bulk, by = c("gene" = "GeneSymbol")) %>%
  mutate(
    quadrant = case_when(
      avg_logFC >= 0 & log2FoldChange >= 0 ~ "Q1",
      avg_logFC < 0 & log2FoldChange >= 0 ~ "Q2",
      avg_logFC < 0 & log2FoldChange < 0 ~ "Q3",
      avg_logFC >= 0 & log2FoldChange < 0 ~ "Q4"
    ),
    dataset = "Sarm1KO"
  )


# Process WT - all genes
wt_pseudo_all <- pseudobulk %>%
  filter(group1 == "Sham-for-1dpi-WT",
         group2 == "1dpi-WT") %>%
  dplyr::select(gene, avg_logFC, p_val_adj)

wt_bulk_all <- bulk_wt %>%
  dplyr::select(GeneSymbol, log2FoldChange, padj)

# Inner join keeps only genes present in both datasets
wt_all_combined <- wt_pseudo_all %>%
  inner_join(wt_bulk_all, by = c("gene" = "GeneSymbol")) %>%
  filter(!is.na(avg_logFC) & !is.na(log2FoldChange)) %>%
  mutate(
    quadrant = case_when(
      avg_logFC >= 0 & log2FoldChange >= 0 ~ "Q1",
      avg_logFC < 0 & log2FoldChange >= 0 ~ "Q2",
      avg_logFC < 0 & log2FoldChange < 0 ~ "Q3",
      avg_logFC >= 0 & log2FoldChange < 0 ~ "Q4"
    ),
    dataset = "WT"
  )

# Process Sarm1KO - all genes
sarm_pseudo_all <- pseudobulk %>%
  filter(group1 == "Sham-for-1dpi-Sarm1-KO",
         group2 == "1dpi-Sarm1-KO") %>%
  dplyr::select(gene, avg_logFC, p_val_adj)

sarm_bulk_all <- bulk_sarm %>%
  dplyr::select(GeneSymbol, log2FoldChange, padj)

sarm_all_combined <- sarm_pseudo_all %>%
  inner_join(sarm_bulk_all, by = c("gene" = "GeneSymbol")) %>%
  filter(!is.na(avg_logFC) & !is.na(log2FoldChange)) %>%
  mutate(
    quadrant = case_when(
      avg_logFC >= 0 & log2FoldChange >= 0 ~ "Q1",
      avg_logFC < 0 & log2FoldChange >= 0 ~ "Q2",
      avg_logFC < 0 & log2FoldChange < 0 ~ "Q3",
      avg_logFC >= 0 & log2FoldChange < 0 ~ "Q4"
    ),
    dataset = "Sarm1KO"
  )

# ============================================================================
#Concordance in top N genes
# Shows overlap in top genes between methods
# ============================================================================

analyze_top_gene_overlap <- function(data, thresholds = c(100, 500, 1000, 5000, 10000)) {
  results <- data.frame()
  
  for (n in thresholds) {
    pseudo_top <- data %>%
      arrange(desc(abs(avg_logFC))) %>%
      head(n) %>%
      pull(gene)
    
    bulk_top <- data %>%
      arrange(desc(abs(log2FoldChange))) %>%
      head(n) %>%
      pull(gene)
    
    overlap <- length(intersect(pseudo_top, bulk_top))
    jaccard <- overlap / length(union(pseudo_top, bulk_top))
    percent_overlap <- (overlap / n) * 100
    
    results <- rbind(results, data.frame(
      n_genes = n,
      overlap = overlap,
      jaccard_index = jaccard,
      percent_overlap = percent_overlap
    ))
  }
  
  return(results)
}

overlap_wt <- analyze_top_gene_overlap(wt_all_combined)
overlap_sarm <- analyze_top_gene_overlap(sarm_all_combined)

cat("\n=== TOP GENE OVERLAP ANALYSIS ===\n")
cat("\nWT:\n")
print(overlap_wt)
cat("\nSarm1KO:\n")
print(overlap_sarm)

# Visualize overlap
plot_overlap <- function(overlap_data, title) {
  p <- ggplot(overlap_data, aes(x = n_genes, y = percent_overlap)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "blue", size = 3) +
    geom_text(aes(label = sprintf("%.1f%%", percent_overlap)), 
              vjust = -0.5, size = 3) +
    scale_x_continuous(trans = "log10", breaks = overlap_data$n_genes) +
    labs(title = title,
         x = "Number of Top Genes",
         y = "% Overlap Between Methods") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

p_overlap_wt <- plot_overlap(overlap_wt, "WT - Top Gene Overlap")
p_overlap_sarm <- plot_overlap(overlap_sarm, "Sarm1KO - Top Gene Overlap")

#Supplemental Figure 7E
grid.arrange(p_overlap_wt, p_overlap_sarm, ncol = 2)

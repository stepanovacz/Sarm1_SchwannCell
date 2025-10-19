################################################################################
# Pseudobulk Differential Expression with GO and KEGG Enrichment Analysis
# 
# Description: Comprehensive pseudobulk analysis of Schwann cells comparing
# WT vs Sarm1-KO at different time points, followed by GO and KEGG pathway
# enrichment analysis, and sex-specific comparisons
#
# Requirements: 
# - Schwann Cells Seurat object (SchwannCells.rds)
# - R packages: edgeR, limma, Seurat, tidyverse, dplyr, clusterProfiler,
#   org.Mm.eg.db, ggplot2, ggrepel, Libra, plyr
# - Note: Libra may need custom library path
################################################################################

# Load required libraries
library(edgeR)
library(limma)
library(Seurat)
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(ggrepel)
library(Libra)

################################################################################
# USER-DEFINED PATHS - MODIFY THESE FOR YOUR SYSTEM
################################################################################

# Path to Schwann Cells Seurat object
schwann_cells_path <- "Data_Files/Files/SchwannCells.rds"

# Paths to save output tables
table6_path <- "Data_Files/Files/Table6.csv"   # Main pseudobulk results
table7_path <- "Data_Files/Files/Table7.csv"   # GO BP WT
table8_path <- "Data_Files/Files/Table8.csv"   # GO CC WT
table9_path <- "Data_Files/Files/Table9.csv"   # GO BP Sarm1-KO
table10_path <- "Data_Files/Files/Table10.csv" # GO CC Sarm1-KO
table11_path <- "Data_Files/Files/Table11.csv" # KEGG WT
table12_path <- "Data_Files/Files/Table12.csv" # KEGG Sarm1-KO
table13_path <- "Data_Files/Files/Table13.csv" # Sex-split results
table14_path <- "Data_Files/Files/Table14.csv" # All time points

################################################################################
# LOAD DATA
################################################################################

cat("=== Loading Schwann Cells Object ===\n")

# SchwannCells <- readRDS(schwann_cells_path)

DefaultAssay(SchwannCells) <- "RNA"

cat("Unique replicates:\n")
print(unique(SchwannCells@meta.data$replicate))
cat("Unique MULTI_IDs:\n")
print(unique(SchwannCells@meta.data$MULTI_ID))

################################################################################
# SUBSET TO 1DPI DATA
################################################################################

cat("\n=== Subsetting to 1dpi Data ===\n")

Idents(SchwannCells) <- "orig.ident"
cat("Unique orig.ident:\n")
print(unique(SchwannCells@meta.data$orig.ident))

SchwannCells_oneDPI <- subset(SchwannCells, 
                              cells = WhichCells(SchwannCells, 
                                                 idents = c("Sham.singletB1", "Injury.singletB1",
                                                            "Sham.singletB2", "Injury.singletB2",
                                                            "Sham.singletB3", "Injury.singletB3",
                                                            "Sham.singletB4", "Injury.singlet.B4")))

cat("1dpi subset orig.ident:\n")
print(unique(SchwannCells_oneDPI@meta.data$orig.ident))

################################################################################
# ADD CELL TYPE METADATA
################################################################################

cat("\n=== Adding Cell Type Metadata ===\n")

Idents(SchwannCells_oneDPI) <- "seurat_clusters"
cat("Unique clusters:\n")
print(unique(SchwannCells_oneDPI@meta.data$seurat_clusters))

current.cluster.idsB <- c('0', '1', '2', '3', '4', '5', '6', '7')
new.cluster.idsB <- c('SC0', 'SC0', "SC0", "SC0", "SC0", "SC0", "SC0", "SC0")
cell_type <- plyr::mapvalues(x = SchwannCells_oneDPI@meta.data$seurat_clusters, 
                             from = current.cluster.idsB, 
                             to = new.cluster.idsB)

names(cell_type) <- colnames(x = SchwannCells_oneDPI)
SchwannCells_oneDPI <- AddMetaData(object = SchwannCells_oneDPI,
                                   metadata = cell_type,
                                   col.name = 'cell_type')

# Check that all columns needed to run Libra are present
cat("Unique replicates:\n")
print(unique(SchwannCells_oneDPI@meta.data$replicate))
cat("Unique labels:\n")
print(unique(SchwannCells_oneDPI@meta.data$label))
cat("Unique cell types:\n")
print(unique(SchwannCells_oneDPI@meta.data$cell_type))

################################################################################
# RUN PSEUDOBULK ANALYSIS - 1DPI
################################################################################

cat("\n=== Running Pseudobulk Differential Expression (1dpi) ===\n")

SC.integrated_Pseudobulk <- JoinLayers(SchwannCells_oneDPI)

cat("Unique labels:\n")
print(unique(SC.integrated_Pseudobulk$label))

# Create the comparisons dataframe
comparisons <- data.frame(
  group1 = c("Sham-for-1dpi-WT", "Sham-for-1dpi-Sarm1-KO"),
  group2 = c("1dpi-WT", "1dpi-Sarm1-KO")
)

cat("Comparisons to perform:\n")
print(comparisons)

# Run differential expression analysis
results <- list()
for (i in 1:nrow(comparisons)) {
  Idents(SC.integrated_Pseudobulk) = SC.integrated_Pseudobulk$label
  group1 = comparisons[i, ]$group1
  group2 = comparisons[i, ]$group2
  
  cat("Running comparison:", group1, "vs", group2, "\n")
  
  sc0 = subset(SC.integrated_Pseudobulk, idents = c(group1, group2))
  
  de = run_de(sc0, n_threads = 8) %>%
    mutate(group1 = group1,
           group2 = group2)
  
  results[[i]] = de
}

# Combine all results into a single table
SchwannCellGenes = do.call(bind_rows, results)
assign("SchwannCellGenes", SchwannCellGenes, envir = .GlobalEnv)

cat("Total genes analyzed:", nrow(SchwannCellGenes), "\n")

# Save results (Table 6)
write.csv(SchwannCellGenes, file = table6_path, row.names = FALSE)
cat("Pseudobulk results saved to:", table6_path, "\n")

################################################################################
# FIGURE 3A: VOLCANO PLOT - WT
################################################################################

cat("\n=== Generating Figure 3A (WT Volcano Plot) ===\n")

comparison_volcano <- SchwannCellGenes %>% 
  filter(group1 == "Sham-for-1dpi-WT", group2 == "1dpi-WT")

cat("Significant genes (padj < 0.05):", 
    nrow(comparison_volcano %>% filter(p_val_adj < 0.05)), "\n")
cat("Non-significant genes:", 
    nrow(comparison_volcano %>% filter(p_val_adj > 0.05)), "\n")

# Get ALL significant upregulated and downregulated genes
all_pos_genes <- comparison_volcano %>%
  filter(p_val_adj < 0.05, avg_logFC > 0)

all_neg_genes <- comparison_volcano %>%
  filter(p_val_adj < 0.05, avg_logFC < 0)

# Get top 10 for labeling
top_pos_genes <- all_pos_genes %>%
  arrange(p_val_adj) %>%
  head(10)

top_neg_genes <- all_neg_genes %>%
  arrange(p_val_adj) %>%
  head(10)

# Figure 3A
comparison_volcano %>%
  ggplot(aes(x = avg_logFC, y = -log10(p_val_adj))) + 
  geom_point(color = ifelse(comparison_volcano$p_val_adj < 0.05, "black", "#FF9999")) +
  geom_text_repel(data = top_pos_genes, aes(label = gene), color = "black", size = 4, 
                  direction = "both", box.padding = 0.3, point.padding = 0.1,
                  max.overlaps = Inf, force = 2, force_pull = 1) +
  geom_text_repel(data = top_neg_genes, aes(label = gene), color = "black", size = 4, 
                  direction = "both", box.padding = 0.3, point.padding = 0.1,
                  max.overlaps = Inf, force = 2, force_pull = 1) +
  labs(x = "avg_logFC", y = "log10(p_val_adj)", 
       title = "WT Sham vs Injury differentially expressed genes") +
  theme(panel.background = element_rect(fill = "white"), 
        axis.line = element_line(color = "black"), 
        plot.title = element_text(hjust = 0.5, vjust = 0.5),
        panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = -log10(0.05), color = "darkred")

################################################################################
# FIGURE 3B: VOLCANO PLOT - SARM1-KO
################################################################################

cat("\n=== Generating Figure 3B (Sarm1-KO Volcano Plot) ===\n")

comparison_volcano <- SchwannCellGenes %>% 
  filter(group1 == "Sham-for-1dpi-Sarm1-KO", group2 == "1dpi-Sarm1-KO")

cat("Significant genes (padj < 0.05):", 
    nrow(comparison_volcano %>% filter(p_val_adj < 0.05)), "\n")

# Get top genes for labeling
all_pos_genes <- comparison_volcano %>%
  filter(p_val_adj < 0.05, avg_logFC > 0)

all_neg_genes <- comparison_volcano %>%
  filter(p_val_adj < 0.05, avg_logFC < 0)

top_pos_genes <- all_pos_genes %>%
  arrange(p_val_adj) %>%
  head(10)

top_neg_genes <- all_neg_genes %>%
  arrange(p_val_adj) %>%
  head(10)

# Figure 3B
comparison_volcano %>%
  ggplot(aes(x = avg_logFC, y = -log10(p_val_adj))) + 
  geom_point(color = ifelse(comparison_volcano$p_val_adj < 0.05, "black", "#FF9999")) +
  geom_text_repel(data = top_pos_genes, aes(label = gene), color = "black", size = 4, 
                  direction = "both", box.padding = 0.3, point.padding = 0.1,
                  max.overlaps = Inf, force = 2, force_pull = 1) +
  geom_text_repel(data = top_neg_genes, aes(label = gene), color = "black", size = 4, 
                  direction = "both", box.padding = 0.3, point.padding = 0.1,
                  max.overlaps = Inf, force = 2, force_pull = 1) +
  labs(x = "avg_logFC", y = "log10(p_val_adj)", 
       title = "Sarm1KO Sham vs Injury differentially expressed genes") +
  theme(panel.background = element_rect(fill = "white"), 
        axis.line = element_line(color = "black"), 
        plot.title = element_text(hjust = 0.5, vjust = 0.5),
        panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = -log10(0.05), color = "darkred")

################################################################################
# FIGURE 3C: WT VS SARM1-KO COMPARISON
################################################################################

cat("\n=== Generating Figure 3C (WT vs Sarm1-KO Comparison) ===\n")

# Prepare data for comparison
wt_data <- SchwannCellGenes %>% 
  filter(group2 == "1dpi-WT") %>%
  dplyr::select(gene, wt_logFC = avg_logFC, wt_p_val_adj = p_val_adj)

ko_data <- SchwannCellGenes %>% 
  filter(group2 == "1dpi-Sarm1-KO") %>%
  dplyr::select(gene, ko_logFC = avg_logFC, ko_p_val_adj = p_val_adj)

# Full join and filter
plot_data <- full_join(wt_data, ko_data, by = "gene") %>%
  filter(!is.na(wt_logFC) & !is.na(ko_logFC))

# Filter for significant genes in AT LEAST one condition
sig_genes <- plot_data %>%
  filter(wt_p_val_adj < 0.05 | ko_p_val_adj < 0.05) %>%
  mutate(
    sig_in_wt = wt_p_val_adj < 0.05,
    sig_in_ko = ko_p_val_adj < 0.05,
    sig_status = case_when(
      sig_in_wt & sig_in_ko ~ "Both",
      sig_in_wt & !sig_in_ko ~ "WT only",
      !sig_in_wt & sig_in_ko ~ "KO only",
      TRUE ~ "Neither"
    ),
    regulation = case_when(
      (wt_logFC > 0 & ko_logFC > 0) ~ "Same (Up)",
      (wt_logFC < 0 & ko_logFC < 0) ~ "Same (Down)",
      (wt_logFC > 0 & ko_logFC < 0) ~ "Opposite (WT Up, KO Down)",
      (wt_logFC < 0 & ko_logFC > 0) ~ "Opposite (WT Down, KO Up)",
      TRUE ~ "Other"
    ),
    regulation_type = ifelse(grepl("Same", regulation), "Same direction", "Opposite direction"),
    discrepancy = abs(wt_logFC - ko_logFC)
  )

# Get top 20 genes for each opposite regulation quadrant
top_left_genes <- sig_genes %>%
  filter(regulation_type == "Opposite direction", 
         wt_logFC < 0, ko_logFC > 0) %>%
  arrange(desc(abs(wt_logFC - ko_logFC))) %>%
  head(20)

bottom_right_genes <- sig_genes %>%
  filter(regulation_type == "Opposite direction",
         wt_logFC > 0, ko_logFC < 0) %>%
  arrange(desc(abs(wt_logFC - ko_logFC))) %>%
  head(20)

highlight_genes <- rbind(top_left_genes, bottom_right_genes)

# Figure 3C
ggplot(sig_genes, aes(x = wt_logFC, y = ko_logFC, 
                      color = regulation_type, 
                      shape = sig_status)) +
  geom_point(alpha = 0.7, size = 0.25) +
  geom_point(data = highlight_genes, 
             aes(x = wt_logFC, y = ko_logFC, color = regulation_type),
             size = 0.5, alpha = 0.9) +
  geom_text_repel(
    data = top_left_genes,
    aes(label = gene),
    color = "black",
    size = 2,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "black",
    min.segment.length = 0,
    segment.size = 0.1,
    force = 10,
    force_pull = 1,
    direction = "both"
  ) +
  geom_text_repel(
    data = bottom_right_genes,
    aes(label = gene),
    color = "black",
    size = 2,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "black",
    min.segment.length = 0,
    segment.size = 0.1,
    force = 10,
    force_pull = 1,
    direction = "both"
  ) +
  scale_color_manual(
    values = c("Same direction" = "purple", "Opposite direction" = "orange"),
    name = "Regulation Pattern"
  ) +
  scale_shape_manual(
    values = c("Both" = 19, "WT only" = 0, "KO only" = 25, "Neither" = 8),
    name = "Significance"
  ) +
  labs(
    title = "Gene Expression Changes in WT vs Sarm1-KO Schwann Cells",
    subtitle = "Color: regulation pattern | Shape: significance status",
    x = "Log Fold Change in WT (1dpi vs Sham)",
    y = "Log Fold Change in Sarm1-KO (1dpi vs Sham)"
  ) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))

################################################################################
# PREPARE OPPOSITE REGULATION GENES FOR GO/KEGG
################################################################################

cat("\n=== Preparing Opposite Regulation Genes ===\n")

# Get all genes with opposite regulation
left_genes_all <- sig_genes %>%
  filter(wt_logFC < 0 & ko_logFC > 0)

right_genes_all <- sig_genes %>%
  filter(wt_logFC > 0 & ko_logFC < 0)

all_genes_opposite <- c(left_genes_all$gene, right_genes_all$gene)

cat("Total opposite regulation genes:", length(all_genes_opposite), "\n")

# Create dataset with only opposite regulation genes
all_genes_opposite_dataset <- SchwannCellGenes[SchwannCellGenes$gene %in% all_genes_opposite, ]

################################################################################
# GO ENRICHMENT ANALYSIS - WT
################################################################################

cat("\n=== Running GO Enrichment Analysis (WT) ===\n")

# Filter for WT opposite regulation genes
opposites_WTs <- all_genes_opposite_dataset %>%
  filter(group1 == "Sham-for-1dpi-WT", group2 == "1dpi-WT") %>%
  filter(!avg_logFC %in% 0) %>%
  filter(p_val_adj < 0.05)

# Convert gene symbols to ENTREZ IDs
ids <- bitr(opposites_WTs$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]), ]
sig_genes <- dedup_ids$ENTREZID

cat("Mapped", length(sig_genes), "genes for GO analysis\n")

# Run GO enrichment analysis for Biological Process
go_bp <- enrichGO(
  gene = sig_genes,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# Figure 3D
dotplot(go_bp, showCategory = 15, title = "GO Biological Process Enrichment WT 1dpi")

# Save results (Table 7)
simple_results_go_BP_WT <- data.frame(
  GO_ID = go_bp@result$ID,
  pValue = go_bp@result$pvalue
)
write.csv(simple_results_go_BP_WT, file = table7_path, row.names = FALSE)
cat("GO BP WT results saved to:", table7_path, "\n")

# Run GO enrichment analysis for Cellular Component
go_cc_WT <- enrichGO(
  gene = sig_genes,
  OrgDb = org.Mm.eg.db,
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# Figure S6C
dotplot(go_cc_WT, showCategory = 15, title = "GO Cellular Component ontology WT 1dpi")

# Save results (Table 8)
simple_results_cc_WT <- data.frame(
  GO_ID = go_cc_WT@result$ID,
  pValue = go_cc_WT@result$pvalue
)
write.csv(simple_results_cc_WT, file = table8_path, row.names = FALSE)
cat("GO CC WT results saved to:", table8_path, "\n")

################################################################################
# GO ENRICHMENT ANALYSIS - SARM1-KO
################################################################################

cat("\n=== Running GO Enrichment Analysis (Sarm1-KO) ===\n")

# Filter for Sarm1-KO opposite regulation genes
opposites_Sarm1KO <- all_genes_opposite_dataset %>%
  filter(group1 == "Sham-for-1dpi-Sarm1-KO", group2 == "1dpi-Sarm1-KO") %>%
  filter(!avg_logFC %in% 0) %>%
  filter(p_val_adj < 0.05)

# Convert gene symbols to ENTREZ IDs
ids <- bitr(opposites_Sarm1KO$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]), ]
sig_genes <- dedup_ids$ENTREZID

cat("Mapped", length(sig_genes), "genes for GO analysis\n")

# Run GO enrichment analysis for Biological Process
go_bp_Sarm1KO <- enrichGO(
  gene = sig_genes,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# Figure 3E
dotplot(go_bp_Sarm1KO, showCategory = 15, title = "GO Biological Processes ontology Sarm1KO 1dpi")

# Save results (Table 9)
simple_results_go_BP_Sarm1KO <- data.frame(
  GO_ID = go_bp_Sarm1KO@result$ID,
  pValue = go_bp_Sarm1KO@result$pvalue
)
write.csv(simple_results_go_BP_Sarm1KO, file = table9_path, row.names = FALSE)
cat("GO BP Sarm1-KO results saved to:", table9_path, "\n")

# Run GO enrichment analysis for Cellular Component
go_cc_Sarm1KO <- enrichGO(
  gene = sig_genes,
  OrgDb = org.Mm.eg.db,
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# Figure S6D
dotplot(go_cc_Sarm1KO, showCategory = 15, title = "GO Cellular Component ontology Sarm1KO 1dpi")

# Save results (Table 10)
simple_results_go_CC_Sarm1KO <- data.frame(
  GO_ID = go_cc_Sarm1KO@result$ID,
  pValue = go_cc_Sarm1KO@result$pvalue
)
write.csv(simple_results_go_CC_Sarm1KO, file = table10_path, row.names = FALSE)
cat("GO CC Sarm1-KO results saved to:", table10_path, "\n")

################################################################################
# KEGG PATHWAY ANALYSIS - WT
################################################################################

cat("\n=== Running KEGG Pathway Analysis (WT) ===\n")

# Filter for WT 1dpi comparison
WT_1dpi <- SchwannCellGenes %>%
  filter(group1 == "Sham-for-1dpi-WT", group2 == "1dpi-WT") %>%
  filter(!avg_logFC %in% 0)

# Create ranked gene list
ranked_genes <- WT_1dpi$avg_logFC
names(ranked_genes) <- WT_1dpi$gene
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Prepare significant gene list
sig_genes <- WT_1dpi %>%
  arrange(desc(abs(avg_logFC))) %>%
  filter(p_val_adj < 0.05) %>%
  pull(gene)

# Convert to ENTREZ IDs
entrez_ids <- bitr(sig_genes, 
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)

cat("Successfully mapped", nrow(entrez_ids), "out of", length(sig_genes), "genes\n")

# Convert ranked list to ENTREZ IDs
ranked_entrez <- bitr(names(ranked_genes),
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Mm.eg.db)

ranked_genes_entrez <- ranked_genes[ranked_entrez$SYMBOL]
names(ranked_genes_entrez) <- ranked_entrez$ENTREZID

# Run GSEA
gc()
set.seed(123)
kk_WT <- gseKEGG(
  geneList = ranked_genes_entrez,
  organism = 'mmu',
  minGSSize = 3,
  pvalueCutoff = 0.05
)

# Figure 3F
dotplot(kk_WT, showCategory = 5) +
  ggtitle("KEGG Wild Type Schwann Cells 1 day post injury")

# Save results (Table 11)
write.csv(kk_WT, file = table11_path, row.names = FALSE)
cat("KEGG WT results saved to:", table11_path, "\n")

################################################################################
# KEGG PATHWAY ANALYSIS - SARM1-KO
################################################################################

cat("\n=== Running KEGG Pathway Analysis (Sarm1-KO) ===\n")

# Filter for Sarm1KO 1dpi comparison
Sarm1KO_1dpi <- SchwannCellGenes %>%
  filter(group1 == "Sham-for-1dpi-Sarm1-KO", group2 == "1dpi-Sarm1-KO") %>%
  filter(!avg_logFC %in% 0)

# Create ranked gene list
ranked_genes <- Sarm1KO_1dpi$avg_logFC
names(ranked_genes) <- Sarm1KO_1dpi$gene
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Prepare significant gene list
sig_genes <- Sarm1KO_1dpi %>%
  arrange(desc(abs(avg_logFC))) %>%
  filter(p_val_adj < 0.05) %>%
  pull(gene)

# Convert to ENTREZ IDs
entrez_ids <- bitr(sig_genes, 
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)

cat("Successfully mapped", nrow(entrez_ids), "out of", length(sig_genes), "genes\n")

# Convert ranked list to ENTREZ IDs
ranked_entrez <- bitr(names(ranked_genes),
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Mm.eg.db)

ranked_genes_entrez <- ranked_genes[ranked_entrez$SYMBOL]
names(ranked_genes_entrez) <- ranked_entrez$ENTREZID

# Run GSEA
gc()
set.seed(123)
kk_Sarm1KO <- gseKEGG(
  geneList = ranked_genes_entrez,
  organism = 'mmu',
  minGSSize = 3,
  pvalueCutoff = 0.05
)

# Figure 3E
dotplot(kk_Sarm1KO, showCategory = 5) +
  ggtitle("KEGG Sarm1 Knockout Schwann Cells 1 day post injury")

# Save results (Table 12)
write.csv(kk_Sarm1KO, file = table12_path, row.names = FALSE)
cat("KEGG Sarm1-KO results saved to:", table12_path, "\n")

################################################################################
# EXTRACT OXIDATIVE PHOSPHORYLATION GENES
################################################################################

cat("\n=== Extracting Oxidative Phosphorylation Genes ===\n")

# Function to convert ENTREZ IDs to symbols
convert_entrez_to_symbols <- function(entrez_ids) {
  entrez_ids <- gsub(" ", "", entrez_ids)
  entrez_ids <- strsplit(entrez_ids, ",")[[1]]
  
  genes_df <- bitr(entrez_ids, 
                   fromType = "ENTREZID",
                   toType = "SYMBOL",
                   OrgDb = org.Mm.eg.db)
  
  return(genes_df)
}

# Find the specific pathway "Oxidative phosphorylation"
oxphos_index <- grep("Oxidative phosphorylation", kk_Sarm1KO@result$Description)
if (length(oxphos_index) > 0) {
  pathway_name <- kk_Sarm1KO@result$Description[oxphos_index]
  genes <- kk_Sarm1KO@result$core_enrichment[oxphos_index]
  
  # Check if these are ENTREZ IDs (numeric) or symbols
  if (grepl("^\\d+", strsplit(genes, "/")[[1]][1])) {
    genes_df <- convert_entrez_to_symbols(gsub("/", ", ", genes))
    gene_symbols <- genes_df$SYMBOL
  } else {
    gene_symbols <- strsplit(genes, "/")[[1]]
  }
  
  # Save into a vector
  oxphos_genes <- gene_symbols
  
  cat("\n===", pathway_name, "===\n")
  cat(paste0('"', gene_symbols, '"', collapse = ", "), "\n")
  cat("\nVector saved as 'oxphos_genes'\n")
  
} else {
  cat("Oxidative phosphorylation pathway not found in the results.\n")
}

################################################################################
# FIGURE 3H: HEATMAP OF OXPHOS GENES
################################################################################

cat("\n=== Generating Figure 3H (OxPhos Heatmap) ===\n")

# Sort the genes alphabetically
all_oxphos_genes <- sort(oxphos_genes)

SchwannCellGenes_NotSplit_By_Sex <- SchwannCellGenes %>%
  mutate(gene_order = factor(gene, levels = all_oxphos_genes))

# Filter out NA rows from the 'gene_order' column
heatmap_genes <- SchwannCellGenes_NotSplit_By_Sex %>%
  filter(!is.na(gene_order))

selected_combinations <- heatmap_genes %>%
  filter(
    (group1 == "Sham-for-2hpi-WT" & group2 == "2hpi-WT") |
      (group1 == "Sham-for-1dpi-Sarm1-KO" & group2 == "1dpi-Sarm1-KO") |
      (group1 == "Sham-for-1dpi-WT" & group2 == "1dpi-WT") |
      (group1 == "Sham-for-3dpi-WT" & group2 == "3dpi-Sarm1-KO") |
      (group1 == "Sham-for-3dpi-WT" & group2 == "3dpi-WT")
  ) %>%
  mutate(combined_group = interaction(group1, group2, sep = "_")) %>%
  group_by(group1, group2) %>%
  ungroup()

# Filter based on selected combinations
filtered_SchwannCellGenes <- heatmap_genes %>%
  semi_join(selected_combinations, by = c("group1", "group2"))

filtered_SchwannCellGenes$group2 <- factor(filtered_SchwannCellGenes$group2, 
                                           levels = c("2hpi-WT", "1dpi-WT", "1dpi-Sarm1-KO",
                                                      "3dpi-WT", "3dpi-Sarm1-KO"))

# Figure 3H
ggplot(filtered_SchwannCellGenes,
       mapping = aes(x = gene_order,
                     y = group2,
                     fill = avg_logFC)) +
  geom_tile(color = "white", size = 0.6) +
  scale_fill_gradient2(
    low = "#FED976",
    mid = "white",
    high = "#08519C",
    midpoint = 0,
    name = "avg_logFC"
  ) +
  scale_y_discrete(limits = rev) +
  labs(x = "", y = "Condition") +
  ggtitle(label = "Top Regulated Genes after Injury\nfor Each Condition") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.line = element_blank(),
        axis.ticks = element_blank()) +
  geom_text(data = filtered_SchwannCellGenes,
            aes(label = case_when(
              p_val_adj < 0.0001 ~ "***",
              p_val_adj < 0.001 ~ "**", 
              p_val_adj < 0.05 ~ "*",
              TRUE ~ ""
            )),
            size = 4,
            hjust = 0.5, vjust = 0.5,
            color = "black",
            fontface = "bold")

################################################################################
# FIGURE S6B: HEATMAP ALTERNATIVE VIEW
################################################################################

cat("\n=== Generating Figure S6D (OxPhos Heatmap Alternative) ===\n")

ggplot(filtered_SchwannCellGenes,
       mapping = aes(x = group2,
                     y = gene_order,
                     fill = avg_logFC)) +
  geom_tile(color = "white", size = 0.6) +
  scale_fill_gradient2(
    low = "#FED976",
    mid = "white",
    high = "#08519C",
    midpoint = 0,
    name = "avg_logFC"
  ) +
  scale_x_discrete(limits = function(x) c("1dpi-WT", setdiff(x, "1dpi-WT"))) +
  scale_y_discrete(limits = function(y) {
    genes_with_A <- y[substr(y, 1, 1) == "A"]
    genes_without_A <- y[substr(y, 1, 1) != "A"]
    rev(c(sort(genes_with_A), sort(genes_without_A)))
  }) +
  labs(x = "Condition", y = "") +
  ggtitle(label = "Top Regulated Genes after Injury\nfor Each Condition") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.line = element_blank(),
        axis.ticks = element_blank()) +
  geom_text(data = filtered_SchwannCellGenes,
            aes(label = case_when(
              p_val_adj < 0.0001 ~ "***",
              p_val_adj < 0.001 ~ "**", 
              p_val_adj < 0.05 ~ "*",
              TRUE ~ ""
            )),
            size = 4,
            hjust = 0.5, vjust = 0.5,
            color = "black",
            fontface = "bold")

################################################################################
# SEX-SPECIFIC ANALYSIS
################################################################################

cat("\n=== Running Sex-Specific Pseudobulk Analysis ===\n")

SC.integrated_Pseudobulk_split_sex <- subset(
  SchwannCells,
  subset = orig.ident %in% c("Sham.singletB1", "Injury.singletB1",
                             "Sham.singletB2", "Injury.singletB2",
                             "Sham.singletB3", "Injury.singletB3",
                             "Sham.singletB4", "Injury.singlet.B4")
)

# Add cell type metadata
Idents(SC.integrated_Pseudobulk_split_sex) <- "seurat_clusters"
current.cluster.idsB <- c('0', '1', '2', '3', '4', '5', '6', '7')
new.cluster.idsB <- c('SC0', 'SC0', "SC0", "SC0", "SC0", "SC0", "SC0", "SC0")
cell_type <- plyr::mapvalues(x = SC.integrated_Pseudobulk_split_sex@meta.data$seurat_clusters, 
                             from = current.cluster.idsB, 
                             to = new.cluster.idsB)

names(cell_type) <- colnames(x = SC.integrated_Pseudobulk_split_sex)
SC.integrated_Pseudobulk_split_sex <- AddMetaData(
  object = SC.integrated_Pseudobulk_split_sex,
  metadata = cell_type,
  col.name = 'cell_type'
)

# Use Hashtags as label for sex-specific analysis
SC.integrated_Pseudobulk_split_sex@meta.data$label <- SC.integrated_Pseudobulk_split_sex@meta.data$Hashtags

SC.integrated_Pseudobulk_split_sex <- JoinLayers(SC.integrated_Pseudobulk_split_sex)

# Create sex-specific comparisons
comparisons <- data.frame(
  group1 = c("Sham-WT-Male", "Sham-WT-Female", "Sham-Sarm1-KO-Male", "Sham-Sarm1-KO-Female"),
  group2 = c("1dpi-WT-Male", "1dpi-WT-Female", "1dpi-Sarm1-KO-Male", "1dpi-Sarm1-KO-Female")
)

cat("Sex-specific comparisons:\n")
print(comparisons)

# Run differential expression analysis
results <- list()
for (i in 1:nrow(comparisons)) {
  Idents(SC.integrated_Pseudobulk_split_sex) = SC.integrated_Pseudobulk_split_sex$label
  group1 = comparisons[i, ]$group1
  group2 = comparisons[i, ]$group2
  
  cat("Running comparison:", group1, "vs", group2, "\n")
  
  sc0 = subset(SC.integrated_Pseudobulk_split_sex, idents = c(group1, group2))
  
  de = run_de(sc0, n_threads = 8) %>%
    mutate(group1 = group1,
           group2 = group2)
  
  results[[i]] = de
}

# Combine all results
SchwannCellGenes_sex_split = do.call(bind_rows, results)
assign("SchwannCellGenes_sex_split", SchwannCellGenes_sex_split, envir = .GlobalEnv)

# Save results (Table 13)
write.csv(SchwannCellGenes_sex_split, file = table13_path, row.names = FALSE)
cat("Sex-split results saved to:", table13_path, "\n")

################################################################################
# FIGURE S6A: SEX-SPECIFIC HEATMAP
################################################################################

cat("\n=== Generating Figure S6A (Sex-Specific Heatmap) ===\n")

# Prepare data for heatmap
SchwannCellGenes_NotSplit_By_Sex_GenesofInterest_new <- SchwannCellGenes_sex_split %>%
  mutate(gene_order = factor(gene, levels = all_oxphos_genes))

heatmap_genes <- SchwannCellGenes_NotSplit_By_Sex_GenesofInterest_new %>%
  filter(!is.na(gene_order))

selected_combinations <- heatmap_genes %>%
  filter(
    (group1 == "Sham-WT-Male" & group2 == "1dpi-WT-Male") |
      (group1 == "Sham-WT-Female" & group2 == "1dpi-WT-Female") |
      (group1 == "Sham-Sarm1-KO-Male" & group2 == "1dpi-Sarm1-KO-Male") |
      (group1 == "Sham-Sarm1-KO-Female" & group2 == "1dpi-Sarm1-KO-Female")
  ) %>%
  mutate(combined_group = interaction(group1, group2, sep = "_")) %>%
  group_by(group1, group2) %>%
  ungroup()

filtered_SchwannCellGenes <- heatmap_genes %>%
  semi_join(selected_combinations, by = c("group1", "group2"))

filtered_SchwannCellGenes$group2 <- factor(filtered_SchwannCellGenes$group2, 
                                           levels = c("1dpi-WT-Male", 
                                                      "1dpi-WT-Female", 
                                                      "1dpi-Sarm1-KO-Male", 
                                                      "1dpi-Sarm1-KO-Female"))

# Figure S6A
ggplot(filtered_SchwannCellGenes,
       mapping = aes(x = group2,
                     y = gene_order,
                     fill = avg_logFC)) +
  geom_tile(color = "white", size = 0.6) +
  scale_fill_gradient2(
    low = "#FED976",
    mid = "white",
    high = "#08519C",
    midpoint = 0,
    name = "avg_logFC"
  ) +
  labs(x = "Condition", y = "") +
  ggtitle(label = "Top Regulated Genes after Injury\nfor Each Condition") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.line = element_blank(),
        axis.ticks = element_blank()) +
  geom_text(data = filtered_SchwannCellGenes,
            aes(label = case_when(
              p_val_adj < 0.0001 ~ "***",
              p_val_adj < 0.001 ~ "**", 
              p_val_adj < 0.05 ~ "*",
              TRUE ~ ""
            )),
            size = 2,
            hjust = 0.5, vjust = 0.5,
            color = "black",
            fontface = "bold")

################################################################################
# FIGURES S6E-F: SEX COMPARISON PLOTS
################################################################################

cat("\n=== Generating Figures S6E-F (Sex Comparison) ===\n")

# Figure S6E: WT Male vs Female
wt_male_data <- SchwannCellGenes_sex_split %>% 
  filter(group2 == "1dpi-WT-Male") %>%
  dplyr::select(gene, male_logFC = avg_logFC, male_p_val_adj = p_val_adj)

wt_female_data <- SchwannCellGenes_sex_split %>% 
  filter(group2 == "1dpi-WT-Female") %>%
  dplyr::select(gene, female_logFC = avg_logFC, female_p_val_adj = p_val_adj)

wt_plot_data <- full_join(wt_male_data, wt_female_data, by = "gene") %>%
  filter(!is.na(male_logFC) & !is.na(female_logFC))

wt_sig_genes <- wt_plot_data %>%
  filter(male_p_val_adj < 0.05 & female_p_val_adj < 0.05) %>%
  mutate(
    sig_status = "Both",
    regulation = case_when(
      (male_logFC > 0 & female_logFC > 0) ~ "Same (Up)",
      (male_logFC < 0 & female_logFC < 0) ~ "Same (Down)",
      (male_logFC > 0 & female_logFC < 0) ~ "Opposite (Male Up, Female Down)",
      (male_logFC < 0 & female_logFC > 0) ~ "Opposite (Male Down, Female Up)",
      TRUE ~ "Other"
    ),
    regulation_type = ifelse(grepl("Same", regulation), "Same direction", "Opposite direction"),
    discrepancy = abs(male_logFC - female_logFC)
  )

wt_top_left_genes <- wt_sig_genes %>%
  filter(regulation_type == "Opposite direction", 
         male_logFC < 0, female_logFC > 0) %>%
  arrange(desc(abs(male_logFC - female_logFC))) %>%
  head(10)

wt_bottom_right_genes <- wt_sig_genes %>%
  filter(regulation_type == "Opposite direction",
         male_logFC > 0, female_logFC < 0) %>%
  arrange(desc(abs(male_logFC - female_logFC))) %>%
  head(10)

wt_highlight_genes <- rbind(wt_top_left_genes, wt_bottom_right_genes)

# Figure S6E
ggplot(wt_sig_genes, aes(x = male_logFC, y = female_logFC, 
                         color = regulation_type)) +
  geom_point(alpha = 0.7, size = 1) +
  geom_point(data = wt_highlight_genes, 
             aes(x = male_logFC, y = female_logFC, color = regulation_type),
             size = 0.5, alpha = 0.9) +
  geom_text_repel(
    data = wt_top_left_genes,
    aes(label = gene),
    color = "black",
    size = 2,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "black",
    min.segment.length = 0,
    segment.size = 0.1,
    force = 10,
    force_pull = 1,
    direction = "both"
  ) +
  geom_text_repel(
    data = wt_bottom_right_genes,
    aes(label = gene),
    color = "black",
    size = 2,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "black",
    min.segment.length = 0,
    segment.size = 0.1,
    force = 10,
    force_pull = 1,
    direction = "both"
  ) +
  scale_color_manual(
    values = c("Same direction" = "purple", "Opposite direction" = "orange"),
    name = "Regulation Pattern"
  ) +
  labs(
    title = "Gene Expression Changes: WT Male vs WT Female Schwann Cells",
    subtitle = "Only genes significant in BOTH sexes (p < 0.05)",
    x = "Log Fold Change in WT Male (1dpi vs Sham)",
    y = "Log Fold Change in WT Female (1dpi vs Sham)"
  ) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))

# Figure S6F: Sarm1-KO Male vs Female (similar code structure)
ko_male_data <- SchwannCellGenes_sex_split %>% 
  filter(group2 == "1dpi-Sarm1-KO-Male") %>%
  dplyr::select(gene, male_logFC = avg_logFC, male_p_val_adj = p_val_adj)

ko_female_data <- SchwannCellGenes_sex_split %>% 
  filter(group2 == "1dpi-Sarm1-KO-Female") %>%
  dplyr::select(gene, female_logFC = avg_logFC, female_p_val_adj = p_val_adj)

ko_plot_data <- full_join(ko_male_data, ko_female_data, by = "gene") %>%
  filter(!is.na(male_logFC) & !is.na(female_logFC))

ko_sig_genes <- ko_plot_data %>%
  filter(male_p_val_adj < 0.05 & female_p_val_adj < 0.05) %>%
  mutate(
    sig_status = "Both",
    regulation = case_when(
      (male_logFC > 0 & female_logFC > 0) ~ "Same (Up)",
      (male_logFC < 0 & female_logFC < 0) ~ "Same (Down)",
      (male_logFC > 0 & female_logFC < 0) ~ "Opposite (Male Up, Female Down)",
      (male_logFC < 0 & female_logFC > 0) ~ "Opposite (Male Down, Female Up)",
      TRUE ~ "Other"
    ),
    regulation_type = ifelse(grepl("Same", regulation), "Same direction", "Opposite direction"),
    discrepancy = abs(male_logFC - female_logFC)
  )

ko_top_left_genes <- ko_sig_genes %>%
  filter(regulation_type == "Opposite direction", 
         male_logFC < 0, female_logFC > 0) %>%
  arrange(desc(abs(male_logFC - female_logFC))) %>%
  head(10)

ko_bottom_right_genes <- ko_sig_genes %>%
  filter(regulation_type == "Opposite direction",
         male_logFC > 0, female_logFC < 0) %>%
  arrange(desc(abs(male_logFC - female_logFC))) %>%
  head(10)

ko_highlight_genes <- rbind(ko_top_left_genes, ko_bottom_right_genes)

# Figure S6F
ggplot(ko_sig_genes, aes(x = male_logFC, y = female_logFC, 
                         color = regulation_type)) +
  geom_point(alpha = 0.7, size = 1) +
  geom_point(data = ko_highlight_genes, 
             aes(x = male_logFC, y = female_logFC, color = regulation_type),
             size = 0.5, alpha = 0.9) +
  geom_text_repel(
    data = ko_top_left_genes,
    aes(label = gene),
    color = "black",
    size = 2,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "black",
    min.segment.length = 0,
    segment.size = 0.1,
    force = 10,
    force_pull = 1,
    direction = "both"
  ) +
  geom_text_repel(
    data = ko_bottom_right_genes,
    aes(label = gene),
    color = "black",
    size = 2,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "black",
    min.segment.length = 0,
    segment.size = 0.1,
    force = 10,
    force_pull = 1,
    direction = "both"
  ) +
  scale_color_manual(
    values = c("Same direction" = "purple", "Opposite direction" = "orange"),
    name = "Regulation Pattern"
  ) +
  labs(
    title = "Gene Expression Changes: Sarm1-KO Male vs Sarm1-KO Female Schwann Cells",
    subtitle = "Only genes significant in BOTH sexes (p < 0.05)",
    x = "Log Fold Change in Sarm1-KO Male (1dpi vs Sham)",
    y = "Log Fold Change in Sarm1-KO Female (1dpi vs Sham)"
  ) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))

################################################################################
# ALL TIME POINTS ANALYSIS
################################################################################

cat("\n=== Running All Time Points Pseudobulk Analysis ===\n")

SC.integrated_Pseudobulk <- SchwannCells

# Add cell type metadata
Idents(SC.integrated_Pseudobulk) <- "seurat_clusters"
current.cluster.idsB <- c('0', '1', '2', '3', '4', '5', '6', '7')
new.cluster.idsB <- c('SC0', 'SC0', "SC0", "SC0", "SC0", "SC0", "SC0", "SC0")
cell_type <- plyr::mapvalues(x = SC.integrated_Pseudobulk@meta.data$seurat_clusters, 
                             from = current.cluster.idsB, 
                             to = new.cluster.idsB)

names(cell_type) <- colnames(x = SC.integrated_Pseudobulk)
SC.integrated_Pseudobulk <- AddMetaData(
  object = SC.integrated_Pseudobulk,
  metadata = cell_type,
  col.name = 'cell_type'
)

SC.integrated_Pseudobulk <- JoinLayers(SC.integrated_Pseudobulk)

# Create all time points comparisons
comparisons <- data.frame(
  group1 = c("Sham-for-1dpi-WT", "Sham-for-1dpi-Sarm1-KO", "Sham-for-2hpi-WT"),
  group2 = c("1dpi-WT", "1dpi-Sarm1-KO", "2hpi-WT")
)

cat("All time points comparisons:\n")
print(comparisons)

# Run differential expression analysis
results <- list()
for (i in 1:nrow(comparisons)) {
  Idents(SC.integrated_Pseudobulk) = SC.integrated_Pseudobulk$label
  group1 = comparisons[i, ]$group1
  group2 = comparisons[i, ]$group2
  
  cat("Running comparison:", group1, "vs", group2, "\n")
  
  sc0 = subset(SC.integrated_Pseudobulk, idents = c(group1, group2))
  
  de = run_de(sc0, n_threads = 8) %>%
    mutate(group1 = group1,
           group2 = group2)
  
  results[[i]] = de
}

# Combine all results
SchwannCellGenes_all_time_points = do.call(bind_rows, results)
assign("SchwannCellGenes_all_time_points", SchwannCellGenes_all_time_points, envir = .GlobalEnv)

# Save results (Table 14)
write.csv(SchwannCellGenes_all_time_points, file = table14_path, row.names = FALSE)
cat("All time points results saved to:", table14_path, "\n")

################################################################################
# FIGURE S6B: ALL TIME POINTS HEATMAP
################################################################################

cat("\n=== Generating Figure S6B (All Time Points Heatmap) ===\n")

# NOTE: Code structure similar to Figure 3H but with all_time_points data
# Prepare the heatmap with all time points...

cat("\n=== Analysis Complete ===\n")
cat("All tables and figures have been generated\n")

positiveWTs_2hpi <- SchwannCellGenes_all_time_points %>%
  filter(group1 == "Sham-for-2hpi-WT", group2 == "2hpi-WT")

negativeWTs_2hpi <- SchwannCellGenes_all_time_points %>%
  filter(group1 == "Sham-for-2hpi-WT", group2 == "2hpi-WT")

positiveWTs_1dpi <- SchwannCellGenes_all_time_points %>%
  filter(group1 == "Sham-for-1dpi-WT", group2 == "1dpi-WT")

negativeWTs_1dpi <- SchwannCellGenes_all_time_points %>%
  filter(group1 == "Sham-for-1dpi-WT", group2 == "1dpi-WT")

positiveSarm1s_1dpi <- SchwannCellGenes_all_time_points %>%
  filter(group1 == "Sham-for-1dpi-Sarm1-KO", group2 == "1dpi-Sarm1-KO")

negativeSarm1s_1dpi <- SchwannCellGenes_all_time_points %>%
  filter(group1 == "Sham-for-1dpi-Sarm1-KO", group2 == "1dpi-Sarm1-KO")


combined <- rbind(
  positiveWTs_2hpi,
  negativeWTs_2hpi,
  positiveWTs_1dpi,
  negativeWTs_1dpi,
  positiveSarm1s_1dpi,
  negativeSarm1s_1dpi
)



# Sort the genes alphabetically
all_oxphos_genes <- sort(oxphos_genes) #need to run 


# Create a new variable 'gene_order' to indicate the order the genes were pulled from 'combined'. 
#this allows only plotting on the y axis of ggplot only the genes that pop up in the combined dataset, while still 
#retaining information about this gene for other groups
SchwannCellGenes_NotSplit_By_Sex_GenesofInterest_new <- SchwannCellGenes_all_time_points %>%
  mutate(gene_order = factor(gene, levels = all_oxphos_genes))

# Filter out NA rows from the 'gene_order' column
heatmap_genes<- SchwannCellGenes_NotSplit_By_Sex_GenesofInterest_new %>%
  filter(!is.na(gene_order))

selected_combinations <- heatmap_genes %>%
  filter(
    (group1 == "Sham-for-2hpi-WT" & group2 == "2hpi-WT")|
      (group1 == "Sham-for-1dpi-Sarm1-KO" & group2 == "1dpi-Sarm1-KO")|
      (group1 == "Sham-for-1dpi-WT" & group2 == "1dpi-WT")
  ) %>%
  mutate(combined_group = interaction(group1, group2, sep = "_")) %>%
  group_by(group1, group2) %>%
  ungroup()
tail(selected_combinations)
tail(heatmap_genes)
# Filter SchwannCellGenes based on specific group1 and group2 combinations from selected_combinations
filtered_SchwannCellGenes <- heatmap_genes %>%
  semi_join(selected_combinations, by = c("group1", "group2"))
tail(filtered_SchwannCellGenes)


filtered_SchwannCellGenes$group2 <- factor(filtered_SchwannCellGenes$group2, 
                                           levels = c("2hpi-WT", "1dpi-WT","1dpi-Sarm1-KO"))

#View(filtered_SchwannCellGenes)

#Figure S6B
ggplot(filtered_SchwannCellGenes,
       mapping = aes(x = group2,
                     y = gene_order,
                     fill = avg_logFC)) +
  geom_tile(color = "white", size = 0.6) +
  scale_fill_gradient2(
    low = "#FED976",     # dark blue
    mid = "white",       # white
    high = "#08519C",     # Very dark teal
    midpoint = 0,
    name = "avg_logFC"
  ) +
  labs(x = "Condition", y = "") +
  ggtitle(label = "Top Regulated Genes after Injury\nfor Each Condition") +
  theme_classic() +  # Different theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.line = element_blank(),
        axis.ticks = element_blank()) +
  geom_text(data = filtered_SchwannCellGenes,
            aes(label = case_when(
              p_val_adj < 0.0001 ~ "***",
              p_val_adj < 0.001 ~ "**", 
              p_val_adj < 0.05 ~ "*",
              TRUE ~ ""
            )),
            size = 4,
            hjust = 0.5, vjust = 0.5,
            color = "black",
            fontface = "bold")

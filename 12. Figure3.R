library(edgeR)
library(limma)
library(Seurat)
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(ggrepel)


DefaultAssay(SchwannCells) <- "RNA"

unique(SchwannCells@meta.data$replicate)
unique(SchwannCells@meta.data$MULTI_ID)


#lets leave only 1dpi data
Idents(SchwannCells) <- "orig.ident"
unique(SchwannCells@meta.data$orig.ident)
SchwannCells<-subset(SchwannCells,cells=WhichCells(SchwannCells,idents=c("Sham.singletB1","Injury.singletB1","Sham.singletB2","Injury.singletB2",
                                                                                     "Sham.singletB3","Injury.singletB3","Sham.singletB4","Injury.singlet.B4")))
unique(SchwannCells@meta.data$orig.ident)


##add a new column "cell_type" that assigns the name of the cell type
Idents(SchwannCells) <- "seurat_clusters"
unique(SchwannCells@meta.data$seurat_clusters)

current.cluster.idsB <- c('0', '1', '2','3', '4','5','6','7')
new.cluster.idsB <- c('SC0','SC0',"SC0","SC0", "SC0","SC0","SC0","SC0")
cell_type <- plyr::mapvalues(x = SchwannCells@meta.data$seurat_clusters, from = current.cluster.idsB, to = new.cluster.idsB)
tail(x = SchwannCells[[]])


names(cell_type) <- colnames(x = SchwannCells)
SchwannCells <- AddMetaData(
  object = SchwannCells,
  metadata = cell_type,
  col.name = 'cell_type'
)
tail(x = SchwannCells[[]])

#now check that all columns need to run Libra for analysis are present
unique(SchwannCells@meta.data$replicate)
unique(SchwannCells@meta.data$label)
unique(SchwannCells@meta.data$cell_type)

##################
#devtools::install_github("neurorestore/Libra")
#sometimes linear algebra package "Libra" is being installed. Make sure the correct libra is being loaded. Incorrect loading can be avoided
#by loading edgeR and limma libraries in advance
library(Libra, lib.loc="~/R/x86_64-pc-linux-gnu-library/4.4")
SC.integrated_Pseudobulk <- JoinLayers(SC.integrated_Pseudobulk)

unique(SC.integrated_Pseudobulk$label)
# Create the comparisons dataframe
comparisons <- data.frame(
  group1 = c("Sham-for-1dpi-WT", "Sham-for-1dpi-Sarm1-KO"),
  group2 = c("1dpi-WT", "1dpi-Sarm1-KO")
)


# View the valid comparisons
print(comparisons)
results <- list()
for (i in 1:nrow(comparisons)) {
  Idents(SC.integrated_Pseudobulk) = SC.integrated_Pseudobulk$label
  group1 = comparisons[i,]$group1
  group2 = comparisons[i,]$group2
  
  sc0 = subset(SC.integrated_Pseudobulk, idents = c(group1, group2))
  
  de = run_de(sc0, n_threads = 8) %>%
    mutate(group1 = group1,
           group2 = group2)
  
  results[[i]] = de  # Store the result in the list
}


# Combine all results into a single table
SchwannCellGenes = do.call(bind_rows, results)

# Save the SchwannCellGenes in the global environment
assign("SchwannCellGenes", SchwannCellGenes, envir = .GlobalEnv)

#View(SchwannCellGenes)
#Table 6
write.csv(SchwannCellGenes, file="Data_Files/Files/Table6")


comparison_volcano <- SchwannCellGenes %>% filter(group1=="Sham-for-1dpi-WT", group2=="1dpi-WT")

#number of significantly different genes
comparison_volcano %>%
  filter(p_val_adj < 0.05) %>%
  nrow()

#number of non-significant genes
comparison_volcano %>%
  filter(p_val_adj > 0.05) %>%
  nrow()

# Get ALL significant upregulated and downregulated genes for plotting
all_pos_genes <- comparison_volcano %>%
  filter(p_val_adj < 0.05, avg_logFC > 0)

all_neg_genes <- comparison_volcano %>%
  filter(p_val_adj < 0.05, avg_logFC < 0)

# Get top 10 upregulated genes (positive logFC) for LABELING - most significant first
top_pos_genes <- all_pos_genes %>%
  arrange(p_val_adj) %>%  # Sort by most significant (smallest p-value)
  head(10)  # Take top 10 for labeling

# Get top 10 downregulated genes (negative logFC) for LABELING - most significant first
top_neg_genes <- all_neg_genes %>%
  arrange(p_val_adj) %>%  # Sort by most significant (smallest p-value)
  head(10)  # Take top 10 for labeling

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
  labs(x = "avg_logFC", y = "log10(p_val_adj)", title = "WT Sham vs Injury differentially expressed genes") +
  theme(panel.background = element_rect(fill = "white"), 
        axis.line = element_line(color = "black"), 
        plot.title = element_text(hjust = 0.5, vjust = 0.5),
        panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = -log10(0.05), color = "darkred")

#Figure S6A

comparison_volcano <- SchwannCellGenes %>% filter(group1=="Sham-for-1dpi-Sarm1-KO", group2=="1dpi-Sarm1-KO")

#number of significantly different genes
comparison_volcano %>%
  filter(p_val_adj < 0.05) %>%
  nrow()

#number of non-significant genes
comparison_volcano %>%
  filter(p_val_adj > 0.05) %>%
  nrow()

# Get ALL significant upregulated and downregulated genes for plotting
all_pos_genes <- comparison_volcano %>%
  filter(p_val_adj < 0.05, avg_logFC > 0)

all_neg_genes <- comparison_volcano %>%
  filter(p_val_adj < 0.05, avg_logFC < 0)

# Get top 10 upregulated genes (positive logFC) for LABELING - most significant first
top_pos_genes <- all_pos_genes %>%
  arrange(p_val_adj) %>%  # Sort by most significant (smallest p-value)
  head(10)  # Take top 10 for labeling

# Get top 10 downregulated genes (negative logFC) for LABELING - most significant first
top_neg_genes <- all_neg_genes %>%
  arrange(p_val_adj) %>%  # Sort by most significant (smallest p-value)
  head(10)  # Take top 10 for labeling

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
  labs(x = "avg_logFC", y = "log10(p_val_adj)", title = "Sarm1KO Sham vs Injury differentially expressed genes") +
  theme(panel.background = element_rect(fill = "white"), 
        axis.line = element_line(color = "black"), 
        plot.title = element_text(hjust = 0.5, vjust = 0.5),
        panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = -log10(0.05), color = "darkred")



wt_data <- SchwannCellGenes %>% 
  filter(group2 == "1dpi-WT") %>%
  dplyr::select(gene, wt_logFC = avg_logFC, wt_p_val_adj = p_val_adj)

ko_data <- SchwannCellGenes %>% 
  filter(group2 == "1dpi-Sarm1-KO") %>%
  dplyr::select(gene, ko_logFC = avg_logFC, ko_p_val_adj = p_val_adj)

# Do a full join to ensure we capture all genes
plot_data <- full_join(wt_data, ko_data, by = "gene")

# Handle potential NA values - exclude genes missing either dataset
plot_data <- plot_data %>%
  filter(!is.na(wt_logFC) & !is.na(ko_logFC))

# Filter for significant genes in AT LEAST one condition
sig_genes <- plot_data %>%
  filter(wt_p_val_adj < 0.05 | ko_p_val_adj < 0.05)

# Create information about which genes are significant in which conditions
sig_genes <- sig_genes %>%
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
    # Calculate discrepancy score (higher = bigger difference between conditions)
    discrepancy = abs(wt_logFC - ko_logFC)
  )

# Get top 25 genes for top-left quadrant (WT down, KO up) - OPPOSITE DIRECTION
top_left_genes <- sig_genes %>%
  filter(wt_logFC < 0 & ko_logFC > 0) %>%
  arrange(desc(discrepancy)) %>%
  head(25) %>%
  mutate(highlight_group = "Top Left")

# Get top 25 genes for bottom-right quadrant (WT up, KO down) - OPPOSITE DIRECTION
bottom_right_genes <- sig_genes %>%
  filter(wt_logFC > 0 & ko_logFC < 0) %>%
  arrange(desc(discrepancy)) %>%
  head(25) %>%
  mutate(highlight_group = "Bottom Right")

# Combine the highlighted genes
highlight_genes <- bind_rows(top_left_genes, bottom_right_genes)
############
# Get genes for top-left quadrant (WT down, KO up) - OPPOSITE DIRECTION
left_genes_all <- sig_genes %>%
  filter(wt_logFC < 0 & ko_logFC > 0) %>%
  mutate(highlight_group = "Top Left")

# Get genes for bottom-right quadrant (WT up, KO down) - OPPOSITE DIRECTION
right_genes_all <- sig_genes %>%
  filter(wt_logFC > 0 & ko_logFC < 0) %>%
  mutate(highlight_group = "Bottom Right")

# Combine the highlighted genes
all_genes_opposite <- bind_rows(left_genes_all, right_genes_all)
all_genes_opposite <- all_genes_opposite$gene

####Create dataset with only opposite regulation genes
all_genes_opposite_dataset <- SchwannCellGenes[SchwannCellGenes$gene %in% all_genes_opposite, ]



# Select exactly top 20 genes for each opposite regulation quadrant
# You'll need to define how you're selecting these

# Top 20 genes in TOP LEFT quadrant (high KO, low WT - opposite regulation)
top_left_genes <- sig_genes %>%
  filter(regulation_type == "Opposite direction", 
         wt_logFC < 0, ko_logFC > 0) %>%  # Adjust filter criteria as needed
  arrange(desc(abs(wt_logFC - ko_logFC))) %>%  # Sort by magnitude of difference
  head(20)  # Take exactly top 20

# Top 20 genes in BOTTOM RIGHT quadrant (low KO, high WT - opposite regulation)  
bottom_right_genes <- sig_genes %>%
  filter(regulation_type == "Opposite direction",
         wt_logFC > 0, ko_logFC < 0) %>%  # Adjust filter criteria as needed
  arrange(desc(abs(wt_logFC - ko_logFC))) %>%  # Sort by magnitude of difference
  head(20)  # Take exactly top 20

# Combine for highlighting
highlight_genes <- rbind(top_left_genes, bottom_right_genes)

#Figure 3C
ggplot(sig_genes, aes(x = wt_logFC, y = ko_logFC, 
                      color = regulation_type, 
                      shape = sig_status)) +
  # Single geom_point layer with both color and shape
  geom_point(alpha = 0.7, size = 0.25) +
  
  # Highlight the top discrepancy genes
  geom_point(data = highlight_genes, 
             aes(x = wt_logFC, y = ko_logFC, color = regulation_type),
             size = 0.5, alpha = 0.9) +
  
  # Add labels for TOP LEFT quadrant - FORCE to show all 20
  geom_text_repel(
    data = top_left_genes,
    aes(label = gene),
    color = "black",
    size = 2,
    max.overlaps = Inf,  # Force all labels to show
    box.padding = 0.3,   # Reduced padding
    point.padding = 0.2,  # Reduced padding
    segment.color = "black",
    min.segment.length = 0,
    segment.size = 0.1,
    force = 10,          # Increased force
    force_pull = 1,      # Added force_pull
    direction = "both"
  ) +
  
  # Add labels for BOTTOM RIGHT quadrant - FORCE to show all 20
  geom_text_repel(
    data = bottom_right_genes,
    aes(label = gene),
    color = "black",
    size = 2,
    max.overlaps = Inf,  # Force all labels to show
    box.padding = 0.3,   # Reduced padding
    point.padding = 0.2,  # Reduced padding
    segment.color = "black",
    min.segment.length = 0,
    segment.size = 0.1,
    force = 10,          # Increased force
    force_pull = 1,      # Added force_pull
    direction = "both"
  ) +
  
  # Set colors for regulation pattern
  scale_color_manual(
    values = c("Same direction" = "purple", "Opposite direction" = "orange"),
    name = "Regulation Pattern"
  ) +
  
  # Set shapes for significance status
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
  # Expand the plot area to make room for labels
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))



head(all_genes_opposite)
length(all_genes_opposite)
#Filter out only genes that are regulated in the opposite manner 1 day after injury between WT and Sarm1KO Schwann Cells
all_genes_opposite_dataset <- SchwannCellGenes[SchwannCellGenes$gene %in% all_genes_opposite, ]
opposites_WTs <- all_genes_opposite_dataset %>%
  filter(group1 == "Sham-for-1dpi-WT", group2 == "1dpi-WT") %>%
  filter(!avg_logFC %in% 0) %>%
  filter(p_val_adj < 0.05)


# Convert gene symbols to ENTREZ IDs
ids <- bitr(opposites_WTs$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]),]

# Get significant gene ENTREZ IDs
sig_genes <- dedup_ids$ENTREZID


# Run GO enrichment analysis for Biological Process
go_bp <- enrichGO(
  gene = sig_genes,  # Your list of significant gene ENTREZ IDs
  OrgDb = org.Mm.eg.db,
  ont = "BP",  # Biological Process ontology
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# Figure 3D
dotplot(go_bp, showCategory = 15, title = "GO Biological Process Enrichment WT 1dpi")
# Extract only GO_ID and pValue
simple_results_go_BP_WT <- data.frame(
  GO_ID = go_bp@result$ID,
  pValue = go_bp@result$pvalue
)

# View the first few rows
head(simple_results_go_BP_WT)

#Table7
write.csv(simple_results_go_BP_WT, file="Data_Files/Files/Table7")

# Run GO enrichment analysis for Cellular Component

go_cc_WT <- enrichGO(
  gene = sig_genes,
  OrgDb = org.Mm.eg.db,
  ont = "CC",  # Cellular Component ontology
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
# Extract only GO_ID and pValue
simple_results_cc_WT <- data.frame(
  GO_ID = go_cc_WT@result$ID,
  pValue = go_cc_WT@result$pvalue
)

# View the first few rows
head(simple_results_cc_WT)

#Table8
write.csv(simple_results_cc_WT, file="Data_Files/Files/Table8")

# Figure S6C
dotplot(go_cc_WT, showCategory = 15,title="GO Cellular Component ontology WT 1dpi")


#Look at the Sarm1KO Schwann Cells genes opposite regulation
opposites_Sarm1KO <- all_genes_opposite_dataset %>%
  filter(group1 == "Sham-for-1dpi-Sarm1-KO", group2 == "1dpi-Sarm1-KO") %>%
  filter(!avg_logFC %in% 0) %>%
  filter(p_val_adj < 0.05)
opposites_Sarm1KO
#View(opposites_Sarm1KO)

# Create a named vector for gene list ranking
ranked_genes <- opposites_Sarm1KO$avg_logFC
names(ranked_genes) <- opposites_Sarm1KO$gene
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Prepare gene list for standard enrichment
sig_genes <- opposites_Sarm1KO %>%
  arrange(desc(abs(avg_logFC))) %>%
  pull(gene)

# Convert gene symbols to ENTREZ IDs
ids <- bitr(opposites_Sarm1KO$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]),]

# Get significant gene ENTREZ IDs
sig_genes <- dedup_ids$ENTREZID

# Run GO enrichment analysis for Biological Process
go_bp_Sarm1KO <- enrichGO(
  gene = sig_genes,  # Your list of significant gene ENTREZ IDs
  OrgDb = org.Mm.eg.db,
  ont = "BP",  # Biological Process ontology
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# Figure 3E
dotplot(go_bp_Sarm1KO, showCategory = 15, title="GO Biological Processes ontology Sarm1KO 1dpi")
# Extract only GO_ID and pValue
simple_results_go_BP_Sarm1KO <- data.frame(
  GO_ID = go_bp_Sarm1KO@result$ID,
  pValue = go_bp_Sarm1KO@result$pvalue
)

# View the first few rows
head(simple_results_go_BP_Sarm1KO)

#Table9
write.csv(simple_results_go_BP_Sarm1KO, file="Data_Files/Files/Table9")


go_cc_Sarm1KO <- enrichGO(
  gene = sig_genes,
  OrgDb = org.Mm.eg.db,
  ont = "CC",  # Cellular Component ontology
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
# Extract only GO_ID and pValue
simple_results_go_CC_Sarm1KO <- data.frame(
  GO_ID = go_cc_Sarm1KO@result$ID,
  pValue = go_cc_Sarm1KO@result$pvalue
)

# View the first few rows
head(simple_results_go_CC_Sarm1KO)


#Figure S6D
dotplot(go_cc_Sarm1KO, showCategory = 15, title="GO Cellular Component ontology Sarm1KO 1dpi")

#Table10
write.csv(simple_results_go_CC_Sarm1KO, file="Data_Files/Files/Table10")


###
# Filter for WT 1dpi comparison
WT_1dpi <- SchwannCellGenes %>%
  filter(group1 == "Sham-for-1dpi-WT", group2 == "1dpi-WT") %>%
  filter(!avg_logFC %in% 0)

# Create a named vector for gene list ranking
# This fixes the gseKEGG error by properly formatting the gene list
ranked_genes <- WT_1dpi$avg_logFC
names(ranked_genes) <- WT_1dpi$gene
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Prepare gene list for standard enrichment
sig_genes <- WT_1dpi %>%
  arrange(desc(abs(avg_logFC))) %>%
  filter(p_val_adj < 0.05) %>%
  pull(gene)

# Convert gene symbols to ENTREZ IDs
entrez_ids <- bitr(sig_genes, 
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)

# Print mapping statistics
print(paste("Successfully mapped", nrow(entrez_ids), "out of", length(sig_genes), "genes"))



# Run GSEA KEGG analysis with properly formatted gene list
# Convert gene names to ENTREZ IDs for the ranked list
ranked_entrez <- bitr(names(ranked_genes),
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Mm.eg.db)

# Create new ranked list with ENTREZ IDs
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

#Figure 3F

dotplot(kk_WT, showCategory=5) +
  ggtitle("KEGG Wild Type Schwann Cells 1 day post injury")


#Table11
write.csv(kk_WT, file="Data_Files/Files/Table11")


# Filter for Sarm1KO 1dpi comparison
Sarm1KO_1dpi <- SchwannCellGenes %>%
  filter(group1 == "Sham-for-1dpi-Sarm1-KO", group2 == "1dpi-Sarm1-KO") %>%
  filter(!avg_logFC %in% 0)

# Create a named vector for gene list ranking
# This fixes the gseKEGG error by properly formatting the gene list
ranked_genes <- Sarm1KO_1dpi$avg_logFC
names(ranked_genes) <- Sarm1KO_1dpi$gene
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Prepare gene list for standard enrichment
sig_genes <- Sarm1KO_1dpi %>%
  arrange(desc(abs(avg_logFC))) %>%
  filter(p_val_adj < 0.05) %>%
  pull(gene)

# Convert gene symbols to ENTREZ IDs
entrez_ids <- bitr(sig_genes, 
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)

# Print mapping statistics
print(paste("Successfully mapped", nrow(entrez_ids), "out of", length(sig_genes), "genes"))



# Run GSEA KEGG analysis with properly formatted gene list
# Convert gene names to ENTREZ IDs for the ranked list
ranked_entrez <- bitr(names(ranked_genes),
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Mm.eg.db)

# Create new ranked list with ENTREZ IDs
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

# Display results
#View(kk_Sarm1KO@result)

#Figure 4E
dotplot(kk_Sarm1KO, showCategory=5) +
  ggtitle("KEGG Sarm1 Knockout Schwann Cells 1 day post injury")


#Table12
write.csv(kk_Sarm1KO, file="Data_Files/Files/Table12")


###Create the function ot extract the genes
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
  # Get pathway name
  pathway_name <- kk_Sarm1KO@result$Description[oxphos_index]
  
  # Get genes
  genes <- kk_Sarm1KO@result$core_enrichment[oxphos_index]
  
  # Check if these are ENTREZ IDs (numeric) or symbols
  if (grepl("^\\d+", strsplit(genes, "/")[[1]][1])) {
    # Convert ENTREZ IDs to symbols using your function
    genes_df <- convert_entrez_to_symbols(gsub("/", ", ", genes))
    gene_symbols <- genes_df$SYMBOL
  } else {
    # Already symbols, just split by "/"
    gene_symbols <- strsplit(genes, "/")[[1]]
  }
  
  # Save into a vector
  oxphos_genes <- gene_symbols
  
  # Print the results
  cat("\n===", pathway_name, "===\n")
  cat(paste0('"', gene_symbols, '"', collapse = ", "), "\n")
  cat("\nVector saved as 'oxphos_genes'\n")
  
} else {
  cat("Oxidative phosphorylation pathway not found in the results.\n")
}

oxphos_genes

WT_2hpis <- SchwannCellGenes %>%
  filter(group1 == "Sham-for-2hpi-WT", group2 == "2hpi-WT")

WTs_1dpi <- SchwannCellGenes %>%
  filter(group1 == "Sham-for-1dpi-WT", group2 == "1dpi-WT")

Sarm1KOs_1dpi <- SchwannCellGenes %>%
  filter(group1 == "Sham-for-1dpi-Sarm1-KO", group2 == "1dpi-Sarm1-KO")


combined <- rbind(
  WT_2hpis,
  WTs_1dpi,
  Sarm1KOs_1dpi
)



# Sort the genes alphabetically.
all_oxphos_genes <- sort(oxphos_genes)

SchwannCellGenes_NotSplit_By_Sex <- SchwannCellGenes %>%
  mutate(gene_order = factor(gene, levels = all_oxphos_genes))

# Filter out NA rows from the 'gene_order' column
heatmap_genes<- SchwannCellGenes_NotSplit_By_Sex %>%
  filter(!is.na(gene_order))

selected_combinations <- heatmap_genes %>%
  filter(
    (group1 == "Sham-for-2hpi-WT" & group2 == "2hpi-WT")|
      (group1 == "Sham-for-1dpi-Sarm1-KO" & group2 == "1dpi-Sarm1-KO")|
      (group1 == "Sham-for-1dpi-WT" & group2 == "1dpi-WT")|
      (group1 == "Sham-for-3dpi-WT" & group2 == "3dpi-Sarm1-KO")|
      (group1 == "Sham-for-3dpi-WT" & group2 == "3dpi-WT")
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
                                           levels = c("2hpi-WT", "1dpi-WT","1dpi-Sarm1-KO",
                                                      "3dpi-WT", 
                                                      "3dpi-Sarm1-KO"))

#View(filtered_SchwannCellGenes)
#Figure 3H
ggplot(filtered_SchwannCellGenes,
       mapping = aes(x = gene_order,
                     y = group2,
                     fill = avg_logFC)) +
  geom_tile(color = "white", size = 0.6) +
  scale_fill_gradient2(
    low = "#FED976",     # dark blue
    mid = "white",       # white
    high = "#08519C",     # Very dark teal
    midpoint = 0,
    name = "avg_logFC"
  ) +
  scale_y_discrete(limits = rev) +  # This reverses the y-axis order to put WT on top
  labs(x = "", y = "Condition") +
  ggtitle(label = "Top Regulated Genes after Injury\nfor Each Condition") +
  theme_classic() +  # Different theme
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

##Figure S7D
ggplot(filtered_SchwannCellGenes,
       mapping = aes(x = group2,           # Conditions now on x-axis
                     y = gene_order,       # Genes now on y-axis
                     fill = avg_logFC)) +
  geom_tile(color = "white", size = 0.6) +
  scale_fill_gradient2(
    low = "#FED976",     # dark blue
    mid = "white",       # white
    high = "#08519C",     # Very dark teal
    midpoint = 0,
    name = "avg_logFC"
  ) +
  scale_x_discrete(limits = function(x) c("1dpi-WT", setdiff(x, "1dpi-WT"))) +  # Put 1dpi-WT first, keep all others
  scale_y_discrete(limits = function(y) {
    # Assuming gene_order contains gene names - put genes starting with A on top
    genes_with_A <- y[substr(y, 1, 1) == "A"]
    genes_without_A <- y[substr(y, 1, 1) != "A"]
    rev(c(sort(genes_with_A), sort(genes_without_A)))  # rev() to put A genes at top
  }) +
  labs(x = "Condition", y = "") +   # Swapped axis labels
  ggtitle(label = "Top Regulated Genes after Injury\nfor Each Condition") +
  theme_classic() +  # Different theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Conditions at angle
        axis.text.y = element_text(size = 10),                         # Gene names horizontal
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


###Now let's look at sex differences

SC.integrated_Pseudobulk_split_sex<-subset(SchwannCells,cells=WhichCells(SchwannCells,idents=c("Sham.singletB1","Injury.singletB1","Sham.singletB2","Injury.singletB2",
                                                                                               "Sham.singletB3","Injury.singletB3","Sham.singletB4","Injury.singlet.B4")))
unique(SC.integrated_Pseudobulk_split_sex@meta.data$orig.ident)



##add a new column "cell_type" that assigns the name of the cell type
Idents(SC.integrated_Pseudobulk_split_sex) <- "seurat_clusters"
unique(SC.integrated_Pseudobulk_split_sex@meta.data$seurat_clusters)

current.cluster.idsB <- c('0', '1', '2','3', '4','5','6','7')
new.cluster.idsB <- c('SC0','SC0',"SC0","SC0", "SC0","SC0", "SC0","SC0")
cell_type <- plyr::mapvalues(x = SC.integrated_Pseudobulk_split_sex@meta.data$seurat_clusters, from = current.cluster.idsB, to = new.cluster.idsB)
tail(x = SC.integrated_Pseudobulk_split_sex[[]])


names(cell_type) <- colnames(x = SC.integrated_Pseudobulk_split_sex)
SC.integrated_Pseudobulk_split_sex <- AddMetaData(
  object = SC.integrated_Pseudobulk_split_sex,
  metadata = cell_type,
  col.name = 'cell_type'
)
tail(x = SC.integrated_Pseudobulk_split_sex[[]])

#now check that all columns need to run Libra for analysis are present
unique(SC.integrated_Pseudobulk_split_sex@meta.data$replicate)
unique(SC.integrated_Pseudobulk_split_sex@meta.data$Hashtags)
SC.integrated_Pseudobulk_split_sex@meta.data$label<-SC.integrated_Pseudobulk_split_sex@meta.data$Hashtags
unique(SC.integrated_Pseudobulk_split_sex@meta.data$label)

unique(SC.integrated_Pseudobulk_split_sex@meta.data$cell_type)

SC.integrated_Pseudobulk_split_sex <- JoinLayers(SC.integrated_Pseudobulk_split_sex)

unique(SC.integrated_Pseudobulk_split_sex$label)
# Create the comparisons dataframe
comparisons <- data.frame(
  group1 = c("Sham-WT-Male","Sham-WT-Female","Sham-Sarm1-KO-Male", "Sham-Sarm1-KO-Female"),
  group2 = c("1dpi-WT-Male" ,"1dpi-WT-Female","1dpi-Sarm1-KO-Male", "1dpi-Sarm1-KO-Female")
)


# View the valid comparisons
print(comparisons)
results <- list()
for (i in 1:nrow(comparisons)) {
  Idents(SC.integrated_Pseudobulk_split_sex) = SC.integrated_Pseudobulk_split_sex$label
  group1 = comparisons[i,]$group1
  group2 = comparisons[i,]$group2
  
  sc0 = subset(SC.integrated_Pseudobulk_split_sex, idents = c(group1, group2))
  
  de = run_de(sc0, n_threads = 8) %>%
    mutate(group1 = group1,
           group2 = group2)
  
  results[[i]] = de  # Store the result in the list
}


# Combine all results into a single table
SchwannCellGenes_sex_split = do.call(bind_rows, results)

# Save the SchwannCellGenes_sex_split in the global environment
assign("SchwannCellGenes_sex_split", SchwannCellGenes_sex_split, envir = .GlobalEnv)

#View(SchwannCellGenes_sex_split)

#Table13
write.csv(SchwannCellGenes_sex_split, file="Data_Files/Files/Table13")


WT_male_1dpi <- SchwannCellGenes_sex_split %>%
  filter(group1 == "Sham-WT-Male", group2 == "1dpi-WT-Male")

WT_female_1dpi <- SchwannCellGenes_sex_split %>%
  filter(group1 == "Sham-WT-Female", group2 == "1dpi-WT-Female")

Sarm1KO_female_1dpi <- SchwannCellGenes_sex_split %>%
  filter(group1 == "Sham-Sarm1-KO-Female", group2 == "1dpi-Sarm1-KO-Female")

Sarm1KO_male_1dpi <- SchwannCellGenes_sex_split %>%
  filter(group1 == "Sham-Sarm1-KO-Male", group2 == "1dpi-Sarm1-KO-Male")




combined <- rbind(
  WT_male_1dpi,
  WT_female_1dpi,
  Sarm1KO_female_1dpi,
  Sarm1KO_male_1dpi)



# Sort the genes alphabetically.
all_oxphos_genes <- sort(oxphos_genes) #need to run 

# Create a new variable 'gene_order' to indicate the order the genes were pulled from 'combined'. 
#this allows only plotting on the y axis of ggplot only the genes that pop up in the combined dataset, while still 
#retaining information about this gene for other groups
SchwannCellGenes_NotSplit_By_Sex_GenesofInterest_new <- SchwannCellGenes_sex_split %>%
  mutate(gene_order = factor(gene, levels = all_oxphos_genes))

# Filter out NA rows from the 'gene_order' column
heatmap_genes<- SchwannCellGenes_NotSplit_By_Sex_GenesofInterest_new %>%
  filter(!is.na(gene_order))

selected_combinations <- heatmap_genes %>%
  filter(
    (group1 == "Sham-WT-Male" & group2 == "1dpi-WT-Male")|
      (group1 == "Sham-WT-Female" & group2 == "1dpi-WT-Female")|
      (group1 == "Sham-Sarm1-KO-Male" & group2 == "1dpi-Sarm1-KO-Male")|
      (group1 == "Sham-Sarm1-KO-Female" & group2 == "1dpi-Sarm1-KO-Female")
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

unique(filtered_SchwannCellGenes$group2)
filtered_SchwannCellGenes$group2 <- factor(filtered_SchwannCellGenes$group2, 
                                           levels = c("1dpi-WT-Male", 
                                                      "1dpi-WT-Female", 
                                                      "1dpi-Sarm1-KO-Male", 
                                                      "1dpi-Sarm1-KO-Female"))

#View(filtered_SchwannCellGenes)
#Figure S6A
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
            size = 2,
            hjust = 0.5, vjust = 0.5,
            color = "black",
            fontface = "bold")

######
SC.integrated_Pseudobulk<-SchwannCells
unique(SC.integrated_Pseudobulk@meta.data$orig.ident)
##add a new column "cell_type" that assigns the name of the cell type
Idents(SC.integrated_Pseudobulk) <- "seurat_clusters"
unique(SC.integrated_Pseudobulk@meta.data$seurat_clusters)

current.cluster.idsB <- c('0', '1', '2','3', '4','5','6','7')
new.cluster.idsB <- c('SC0','SC0',"SC0","SC0", "SC0","SC0", "SC0","SC0")
cell_type <- plyr::mapvalues(x = SC.integrated_Pseudobulk@meta.data$seurat_clusters, from = current.cluster.idsB, to = new.cluster.idsB)
tail(x = SC.integrated_Pseudobulk[[]])


names(cell_type) <- colnames(x = SC.integrated_Pseudobulk)
SC.integrated_Pseudobulk <- AddMetaData(
  object = SC.integrated_Pseudobulk,
  metadata = cell_type,
  col.name = 'cell_type'
)
tail(x = SC.integrated_Pseudobulk[[]])

#now check that all columns need to run Libra for analysis are present
unique(SC.integrated_Pseudobulk@meta.data$replicate)
unique(SC.integrated_Pseudobulk@meta.data$label)
unique(SC.integrated_Pseudobulk@meta.data$cell_type)

SC.integrated_Pseudobulk <- JoinLayers(SC.integrated_Pseudobulk)

unique(SC.integrated_Pseudobulk$label)
# Create the comparisons dataframe
comparisons <- data.frame(
  group1 = c("Sham-for-1dpi-WT", "Sham-for-1dpi-Sarm1-KO","Sham-for-2hpi-WT"),
  group2 = c("1dpi-WT", "1dpi-Sarm1-KO", "2hpi-WT")
)


# View the valid comparisons
print(comparisons)
results <- list()
for (i in 1:nrow(comparisons)) {
  Idents(SC.integrated_Pseudobulk) = SC.integrated_Pseudobulk$label
  group1 = comparisons[i,]$group1
  group2 = comparisons[i,]$group2
  
  sc0 = subset(SC.integrated_Pseudobulk, idents = c(group1, group2))
  
  de = run_de(sc0, n_threads = 8) %>%
    mutate(group1 = group1,
           group2 = group2)
  
  results[[i]] = de  # Store the result in the list
}


# Combine all results into a single table
SchwannCellGenes_all_time_points = do.call(bind_rows, results)

# Save the SchwannCellGenes_all_time_points in the global environment
assign("SchwannCellGenes_all_time_points", SchwannCellGenes_all_time_points, envir = .GlobalEnv)

#View(SchwannCellGenes_all_time_points)

#Table14
write.csv(SchwannCellGenes_all_time_points, file="Data_Files/Files/Table14")

#Plot this
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

positiveWTs_3dpi <- SchwannCellGenes_all_time_points %>%
  filter(group1 == "Sham-for-3dpi-WT", group2 == "3dpi-WT")

negativeWTs_3dpi <- SchwannCellGenes_all_time_points %>%
  filter(group1 == "Sham-for-3dpi-WT", group2 == "3dpi-WT")

positiveSarm1s_3dpi <- SchwannCellGenes_all_time_points %>%
  filter(group1 == "Sham-for-3dpi-Sarm1-KO", group2 == "3dpi-Sarm1-KO")

negativeSarm1s_3dpi <- SchwannCellGenes_all_time_points %>%
  filter(group1 == "Sham-for-3dpi-Sarm1-KO", group2 == "3dpi-Sarm1-KO")

combined <- rbind(
  positiveWTs_2hpi,
  negativeWTs_2hpi,
  positiveWTs_1dpi,
  negativeWTs_1dpi,
  positiveSarm1s_1dpi,
  negativeSarm1s_1dpi,
  positiveWTs_3dpi,
  negativeWTs_3dpi,
  positiveSarm1s_3dpi,
  negativeSarm1s_3dpi
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
      (group1 == "Sham-for-1dpi-WT" & group2 == "1dpi-WT")|
      (group1 == "Sham-for-3dpi-WT" & group2 == "3dpi-Sarm1-KO")|
      (group1 == "Sham-for-3dpi-WT" & group2 == "3dpi-WT")
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
                                           levels = c("2hpi-WT", "1dpi-WT","1dpi-Sarm1-KO",
                                                      "3dpi-WT", 
                                                      "3dpi-Sarm1-KO"))

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


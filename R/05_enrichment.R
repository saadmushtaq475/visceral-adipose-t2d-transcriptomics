# ============================================================
# PROJECT: Visceral Adipose Tissue Transcriptomic Analysis
# DATASET: GSE78721
# ANALYSIS: Functional Enrichment Analysis
#
# PURPOSE:
# Perform Gene Ontology (GO) and KEGG pathway enrichment
# analysis using exploratory candidate genes identified from
# the limma differential expression analysis.
#
# INPUT:
# - Candidate gene list from differential expression analysis
#
# OUTPUT:
# - GO Biological Process
# - GO Molecular Function
# - GO Cellular Component
# - KEGG pathway enrichment
# - Enrichment result tables
# - Enrichment plots
# ============================================================


# ============================================================
# 01. LOAD REQUIRED PACKAGES
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(writexl)


# ============================================================
# 02. CREATE OUTPUT DIRECTORIES
# ============================================================

dir.create("results", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)


# ============================================================
# 03. LOAD CANDIDATE GENES
# ============================================================

# Candidate genes were defined as genes with:
# |log2FC| >= 1
# P-value < 0.05
#
# IMPORTANT:
# These are exploratory candidate genes.
# They are NOT considered FDR-significant DEGs unless
# adjusted P-value < 0.05.

candidate_file <- "results/tables/05_Candidate_Genes_for_STRING.xlsx"

candidate_annotated <- readxl::read_xlsx(candidate_file)

# Inspect data
head(candidate_annotated)

# Extract unique gene symbols
gene_symbols <- unique(
  candidate_annotated$Gene_Symbol
)

# Remove missing/empty values
gene_symbols <- gene_symbols[
  !is.na(gene_symbols) &
    gene_symbols != ""
]

cat(
  "Number of candidate genes:",
  length(gene_symbols),
  "\n"
)


# ============================================================
# 04. CONVERT GENE SYMBOLS TO ENTREZ IDs
# ============================================================

entrez <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Remove duplicate mappings
entrez <- entrez %>%
  distinct(SYMBOL, .keep_all = TRUE)

cat(
  "Genes successfully mapped to Entrez IDs:",
  nrow(entrez),
  "\n"
)

# Check unmapped genes
unmapped_genes <- setdiff(
  gene_symbols,
  entrez$SYMBOL
)

cat(
  "Unmapped genes:",
  length(unmapped_genes),
  "\n"
)


# ============================================================
# 05. GO BIOLOGICAL PROCESS
# ============================================================

GO_BP <- enrichGO(
  gene = entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

GO_BP_results <- as.data.frame(GO_BP)

head(GO_BP_results)


# ============================================================
# 06. GO MOLECULAR FUNCTION
# ============================================================

GO_MF <- enrichGO(
  gene = entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

GO_MF_results <- as.data.frame(GO_MF)

head(GO_MF_results)


# ============================================================
# 07. GO CELLULAR COMPONENT
# ============================================================

GO_CC <- enrichGO(
  gene = entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

GO_CC_results <- as.data.frame(GO_CC)

head(GO_CC_results)


# ============================================================
# 08. KEGG PATHWAY ENRICHMENT
# ============================================================

KEGG <- enrichKEGG(
  gene = entrez$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

KEGG_results <- as.data.frame(KEGG)

head(KEGG_results)


# ============================================================
# 09. SAVE ENRICHMENT RESULTS
# ============================================================

write_xlsx(
  list(
    GO_Biological_Process = GO_BP_results,
    GO_Molecular_Function = GO_MF_results,
    GO_Cellular_Component = GO_CC_results,
    KEGG_Pathways = KEGG_results,
    Gene_ID_Mapping = entrez
  ),
  "results/tables/06_GO_KEGG_Enrichment.xlsx"
)


# ============================================================
# 10. SAVE UNMAPPED GENES
# ============================================================

write.table(
  unmapped_genes,
  "results/tables/06_Unmapped_Genes.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)


# ============================================================
# 11. GO BIOLOGICAL PROCESS PLOT
# ============================================================

if (nrow(GO_BP_results) > 0) {

  p_GO_BP <- barplot(
    GO_BP,
    showCategory = 15,
    title = "GO Biological Process"
  )

  ggsave(
    "results/figures/06_GO_Biological_Process.png",
    p_GO_BP,
    width = 10,
    height = 7,
    dpi = 300
  )
}


# ============================================================
# 12. GO MOLECULAR FUNCTION PLOT
# ============================================================

if (nrow(GO_MF_results) > 0) {

  p_GO_MF <- barplot(
    GO_MF,
    showCategory = 15,
    title = "GO Molecular Function"
  )

  ggsave(
    "results/figures/06_GO_Molecular_Function.png",
    p_GO_MF,
    width = 10,
    height = 7,
    dpi = 300
  )
}


# ============================================================
# 13. GO CELLULAR COMPONENT PLOT
# ============================================================

if (nrow(GO_CC_results) > 0) {

  p_GO_CC <- barplot(
    GO_CC,
    showCategory = 15,
    title = "GO Cellular Component"
  )

  ggsave(
    "results/figures/06_GO_Cellular_Component.png",
    p_GO_CC,
    width = 10,
    height = 7,
    dpi = 300
  )
}


# ============================================================
# 14. KEGG PLOT
# ============================================================

if (nrow(KEGG_results) > 0) {

  p_KEGG <- barplot(
    KEGG,
    showCategory = 15,
    title = "KEGG Pathway Enrichment"
  )

  ggsave(
    "results/figures/06_KEGG_Pathways.png",
    p_KEGG,
    width = 10,
    height = 7,
    dpi = 300
  )
}


# ============================================================
# 15. DOT PLOTS
# ============================================================

if (nrow(GO_BP_results) > 0) {

  p_GO_BP_dot <- dotplot(
    GO_BP,
    showCategory = 15
  )

  ggsave(
    "results/figures/06_GO_BP_Dotplot.png",
    p_GO_BP_dot,
    width = 10,
    height = 7,
    dpi = 300
  )
}


if (nrow(KEGG_results) > 0) {

  p_KEGG_dot <- dotplot(
    KEGG,
    showCategory = 15
  )

  ggsave(
    "results/figures/06_KEGG_Dotplot.png",
    p_KEGG_dot,
    width = 10,
    height = 7,
    dpi = 300
  )
}


# ============================================================
# 16. SUMMARY
# ============================================================

cat("\n")
cat("============================================\n")
cat("FUNCTIONAL ENRICHMENT ANALYSIS COMPLETED\n")
cat("============================================\n")

cat(
  "Candidate genes:",
  length(gene_symbols),
  "\n"
)

cat(
  "GO Biological Process terms:",
  nrow(GO_BP_results),
  "\n"
)

cat(
  "GO Molecular Function terms:",
  nrow(GO_MF_results),
  "\n"
)

cat(
  "GO Cellular Component terms:",
  nrow(GO_CC_results),
  "\n"
)

cat(
  "KEGG pathways:",
  nrow(KEGG_results),
  "\n"
)

cat("============================================\n")

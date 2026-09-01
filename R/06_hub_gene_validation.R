# ============================================================
# PROJECT: Visceral Adipose Tissue Transcriptomic Analysis
# DATASET: GSE78721
# ANALYSIS: Hub Gene Validation
#
# PURPOSE:
# Validate the expression-level evidence for the eight hub genes
# identified through STRING/Cytoscape network analysis.
#
# IMPORTANT:
# Hub status was determined from network analysis and does not
# automatically mean that a gene is a statistically significant
# DEG.
#
# Expression evidence is therefore classified as:
# 1. FDR-significant DEG
# 2. Nominally significant candidate
# 3. Network-derived hub only
# ============================================================


# ============================================================
# 01. LOAD REQUIRED PACKAGES
# ============================================================

library(dplyr)
library(readxl)
library(writexl)


# ============================================================
# 02. CREATE OUTPUT DIRECTORIES
# ============================================================

dir.create("results", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE)


# ============================================================
# 03. LOAD ANNOTATED LIMMA RESULTS
# ============================================================

# This file should contain the complete limma results
# annotated using GPL15207.

deg_annotated <- read_xlsx(
  "results/tables/04_All_DEGs_Annotated_GPL15207.xlsx"
)

head(deg_annotated)


# ============================================================
# 04. DEFINE THE EIGHT HUB GENES
# ============================================================

hub_genes <- c(
  "COL1A1",
  "VTN",
  "TIMP1",
  "VCAM1",
  "THY1",
  "KRT19",
  "KRT18",
  "KRT8"
)


# ============================================================
# 05. EXTRACT HUB-GENE STATISTICAL RESULTS
# ============================================================

hub_validation <- deg_annotated %>%
  filter(
    Gene_Symbol %in% hub_genes
  )

# Inspect results
hub_validation


# ============================================================
# 06. SAVE COMPLETE HUB-GENE RESULTS
# ============================================================

write_xlsx(
  hub_validation,
  "results/tables/07_Eight_Hub_Genes_Validation.xlsx"
)


# ============================================================
# 07. CREATE ONE SUMMARY PER GENE
# ============================================================

hub_summary <- hub_validation %>%
  group_by(Gene_Symbol) %>%
  summarise(
    Best_logFC = max(logFC, na.rm = TRUE),
    Best_P = min(P.Value, na.rm = TRUE),
    Best_FDR = min(adj.P.Val, na.rm = TRUE),
    Number_of_probes = n(),
    .groups = "drop"
  )


# ============================================================
# 08. CLASSIFY EXPRESSION EVIDENCE
# ============================================================

hub_summary <- hub_summary %>%
  mutate(
    Expression_Evidence = case_when(

      Best_FDR < 0.05 &
        abs(Best_logFC) >= 1 ~
        "FDR-significant DEG",

      Best_P < 0.05 &
        Best_FDR >= 0.05 &
        abs(Best_logFC) >= 1 ~
        "Nominally significant candidate",

      TRUE ~
        "Network-derived hub only"
    )
  )


# ============================================================
# 09. DETERMINE DIRECTION OF EXPRESSION
# ============================================================

hub_summary <- hub_summary %>%
  mutate(
    Direction = case_when(
      Best_logFC > 0 ~ "Upregulated",
      Best_logFC < 0 ~ "Downregulated",
      TRUE ~ "No change"
    )
  )


# ============================================================
# 10. SORT RESULTS
# ============================================================

hub_summary <- hub_summary %>%
  arrange(Best_P)


# Display final summary
hub_summary


# ============================================================
# 11. SAVE HUB-GENE SUMMARY
# ============================================================

write_xlsx(
  hub_summary,
  "results/tables/07_Hub_Gene_Summary.xlsx"
)


# ============================================================
# 12. CHECK FDR-SIGNIFICANT HUBS
# ============================================================

FDR_hubs <- hub_summary %>%
  filter(
    Best_FDR < 0.05
  )

cat(
  "Number of FDR-significant hub genes:",
  nrow(FDR_hubs),
  "\n"
)


# ============================================================
# 13. CHECK NOMINALLY SIGNIFICANT HUBS
# ============================================================

nominal_hubs <- hub_summary %>%
  filter(
    Best_P < 0.05,
    Best_FDR >= 0.05
  )

cat(
  "Number of nominally significant hub genes:",
  nrow(nominal_hubs),
  "\n"
)


# ============================================================
# 14. CHECK NETWORK-DERIVED HUBS
# ============================================================

network_only_hubs <- hub_summary %>%
  filter(
    Best_P >= 0.05
  )

cat(
  "Number of network-derived hubs without nominal
  expression significance:",
  nrow(network_only_hubs),
  "\n"
)


# ============================================================
# 15. FINAL SUMMARY
# ============================================================

cat("\n")
cat("============================================\n")
cat("HUB GENE VALIDATION COMPLETED\n")
cat("============================================\n")

print(hub_summary)

cat("============================================\n")

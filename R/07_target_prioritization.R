# ============================================================
# PROJECT: Visceral Adipose Tissue Transcriptomic Analysis
# DATASET: GSE78721
# ANALYSIS: Target Prioritization
#
# PURPOSE:
# Prioritize hub genes for downstream drug-repurposing analysis
# based on expression evidence, network relevance, biological
# interpretation, druggability, and structural availability.
#
# IMPORTANT:
# Network hub status does not by itself establish therapeutic
# relevance. Target prioritization is therefore treated as a
# multi-criteria exploratory step.
#
# DOWNSTREAM:
# DrugRep receptor-based virtual screening
# ============================================================


# ============================================================
# 01. LOAD REQUIRED PACKAGES
# ============================================================

library(dplyr)
library(readxl)
library(writexl)


# ============================================================
# 02. CREATE OUTPUT DIRECTORY
# ============================================================

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 03. LOAD HUB-GENE SUMMARY
# ============================================================

hub_summary <- read_xlsx(
  "results/tables/07_Hub_Gene_Summary.xlsx"
)

hub_summary


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
# 05. BIOLOGICAL INTERPRETATION
# ============================================================

target_annotation <- data.frame(

  Gene = hub_genes,

  Biological_Role = c(
    "Extracellular matrix collagen; associated with adipose tissue fibrosis and tissue remodeling",
    "Vitronectin; extracellular matrix adhesion glycoprotein involved in cell adhesion and integrin signaling",
    "Tissue inhibitor of metalloproteinases; regulates extracellular matrix turnover and metalloproteinase activity",
    "Vascular cell adhesion molecule; involved in endothelial adhesion and inflammatory cell recruitment",
    "THY1/CD90; cell-surface glycoprotein involved in cell adhesion and extracellular matrix signaling",
    "Keratin 19; epithelial/intermediate filament structural protein",
    "Keratin 18; intermediate filament structural protein",
    "Keratin 8; intermediate filament structural protein"
  ),

  T2D_Relevance = c(
    "ECM remodeling and adipose tissue fibrosis",
    "ECM remodeling, adhesion and integrin signaling",
    "ECM turnover and fibrosis",
    "Inflammation and immune-cell adhesion",
    "Cell adhesion and stromal signaling",
    "Cellular structural/remodeling changes",
    "Cellular structural/remodeling changes",
    "Cellular structural/remodeling changes"
  ),

  DrugRepurposing_Priority = c(
    "Moderate",
    "High",
    "High",
    "High",
    "Moderate",
    "Low",
    "Low",
    "Low"
  ),

  stringsAsFactors = FALSE
)


# ============================================================
# 06. STRUCTURAL / DRUGGABILITY ASSESSMENT
# ============================================================

target_annotation$Druggability <- c(
  "Low",
  "High",
  "Moderate-High",
  "High",
  "Moderate",
  "Low",
  "Low",
  "Low"
)

target_annotation$Structural_Consideration <- c(
  "Fibrous structural protein; lacks a conventional small-molecule pocket",
  "Defined extracellular interaction surfaces; suitable for receptor-based screening",
  "High-resolution unbound structure available; suitable for exploratory screening",
  "Defined extracellular interaction surface; structural data available",
  "Extracellular interaction surface; direct small-molecule precedent limited",
  "Structural protein; limited direct small-molecule targeting rationale",
  "Structural protein; limited direct small-molecule targeting rationale",
  "Structural protein; limited direct small-molecule targeting rationale"
)


# ============================================================
# 07. STRUCTURAL PRIORITY
# ============================================================

target_annotation$Docking_Priority <- c(
  "Exploratory",
  "Primary",
  "Primary",
  "Primary",
  "Secondary",
  "Not prioritized",
  "Not prioritized",
  "Not prioritized"
)


# ============================================================
# 08. ASSIGN TARGET RANK
# ============================================================

target_annotation$Priority_Rank <- c(
  5,
  2,
  3,
  1,
  4,
  NA,
  NA,
  NA
)


# ============================================================
# 09. COMBINE EXPRESSION AND TARGET INFORMATION
# ============================================================

target_prioritization <- target_annotation %>%
  left_join(
    hub_summary,
    by = c("Gene" = "Gene_Symbol")
  ) %>%
  arrange(
    is.na(Priority_Rank),
    Priority_Rank
  )


# ============================================================
# 10. VIEW FINAL TARGET PRIORITIZATION
# ============================================================

target_prioritization


# ============================================================
# 11. EXPORT COMPLETE TARGET PRIORITIZATION
# ============================================================

write_xlsx(
  target_prioritization,
  "results/tables/08_Target_Prioritization.xlsx"
)


# ============================================================
# 12. EXPORT PRIORITIZED TARGET LIST
# ============================================================

prioritized_targets <- target_prioritization %>%
  filter(
    Docking_Priority %in%
      c("Primary", "Secondary")
  ) %>%
  arrange(Priority_Rank)

prioritized_targets


# ============================================================
# 13. SAVE TARGET LIST FOR DOWNSTREAM ANALYSIS
# ============================================================

write.table(
  prioritized_targets$Gene,
  "results/tables/08_Prioritized_Targets.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)


# ============================================================
# 14. FINAL PRIMARY TARGETS
# ============================================================

primary_targets <- target_prioritization %>%
  filter(
    Docking_Priority == "Primary"
  ) %>%
  arrange(Priority_Rank)

cat("\n")
cat("============================================\n")
cat("PRIMARY TARGETS\n")
cat("============================================\n")

print(primary_targets$Gene)


# ============================================================
# 15. SECONDARY TARGETS
# ============================================================

secondary_targets <- target_prioritization %>%
  filter(
    Docking_Priority == "Secondary"
  ) %>%
  arrange(Priority_Rank)

cat("\n")
cat("============================================\n")
cat("SECONDARY TARGETS\n")
cat("============================================\n")

print(secondary_targets$Gene)


# ============================================================
# 16. FINAL SUMMARY
# ============================================================

cat("\n")
cat("============================================\n")
cat("TARGET PRIORITIZATION COMPLETED\n")
cat("============================================\n")

cat(
  "Primary targets:",
  paste(primary_targets$Gene, collapse = ", "),
  "\n"
)

cat(
  "Secondary targets:",
  paste(secondary_targets$Gene, collapse = ", "),
  "\n"
)

cat(
  "Exploratory structural target:",
  "COL1A1",
  "\n"
)

cat("============================================\n")


# ============================================================
# END OF R-BASED TARGET PRIORITIZATION
#
# NEXT STEP:
# DrugRep receptor-based virtual screening
#
# DrugRep results should be documented separately from the
# R transcriptomic workflow.
# ============================================================

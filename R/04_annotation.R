# ============================================================
# PROJECT: Visceral Adipose Tissue Transcriptomics in T2D
# DATASET: GSE78721
# PLATFORM: GPL15207
#
# STEP 04: PROBE-TO-GENE ANNOTATION
# ============================================================


# ============================================================
# 01. LOAD REQUIRED PACKAGES
# ============================================================

library(GEOquery)
library(dplyr)
library(writexl)


# ============================================================
# 02. PROJECT DIRECTORIES
# ============================================================

data_dir <- "data"
results_dir <- "results"
tables_dir <- file.path(results_dir, "tables")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 03. LOAD GEO EXPRESSION DATA
# ============================================================

gse <- getGEO(
  filename = file.path(
    data_dir,
    "GSE78721_series_matrix.txt.gz"
  ),
  getGPL = FALSE
)

meta <- pData(gse)
expr <- exprs(gse)


# ============================================================
# 04. FILTER FOR VISCERAL ADIPOSE TISSUE
# ============================================================

visceral_samples <- rownames(meta)[
  meta$`tissue:ch1` == "Visceral"
]

expr_visceral <- expr[, visceral_samples]


# ============================================================
# 05. LOAD AUTHORITATIVE GPL15207 ANNOTATION
# ============================================================

gpl <- getGEO(
  filename = file.path(
    data_dir,
    "GPL15207.txt"
  )
)

gpl_table <- Table(gpl)


# ============================================================
# 06. CREATE PROBE → GENE MAPPING
# ============================================================

probe_to_gene <- gpl_table %>%
  dplyr::select(
    ID,
    `Gene Symbol`
  ) %>%
  dplyr::rename(
    ProbeID = ID,
    Gene_Symbol = `Gene Symbol`
  )


# ============================================================
# 07. ANNOTATE PROBE-LEVEL LIMMA RESULTS
# ============================================================
#
# The limma results are loaded from the previous analysis.
# This script assumes that the DEG results have been saved
# as an R object or regenerated when the complete workflow
# is run.
#
# For a standalone GitHub workflow, we regenerate limma
# results here to maintain reproducibility.


library(limma)


# Create diagnosis variable
group <- factor(
  meta$`diagnosis:ch1`[visceral_samples],
  levels = c(
    "Normal",
    "Type-2 Diabetes"
  )
)


# Design matrix
design <- model.matrix(
  ~ group
)

colnames(design) <- c(
  "Normal",
  "T2D_vs_Normal"
)


# Run limma
fit <- lmFit(
  expr_visceral,
  design
)

fit <- eBayes(
  fit
)


# Extract complete results
deg_results <- topTable(
  fit,
  coef = "T2D_vs_Normal",
  number = Inf,
  sort.by = "P"
)


# Add probe IDs
deg_results$ProbeID <- rownames(deg_results)


# ============================================================
# 08. MERGE LIMMA RESULTS WITH GPL ANNOTATION
# ============================================================

deg_annotated <- deg_results %>%
  left_join(
    probe_to_gene,
    by = "ProbeID"
  )


# ============================================================
# 09. CLEAN GENE SYMBOLS
# ============================================================

deg_annotated <- deg_annotated %>%
  mutate(
    Gene_Symbol = trimws(Gene_Symbol)
  ) %>%
  filter(
    !is.na(Gene_Symbol),
    Gene_Symbol != ""
  )


# ============================================================
# 10. CREATE EXPLORATORY CANDIDATE GENE SET
# ============================================================
#
# Criteria:
# |log2FC| >= 1
# P-value < 0.05
#
# These are exploratory candidates and are NOT called
# FDR-significant DEGs.


candidate_annotated <- deg_annotated %>%
  filter(
    abs(logFC) >= 1,
    P.Value < 0.05
  )


# ============================================================
# 11. PREPARE GENE SYMBOL LIST FOR STRING
# ============================================================

string_genes <- candidate_annotated %>%
  distinct(Gene_Symbol) %>%
  pull(Gene_Symbol)


# ============================================================
# 12. REMOVE OBVIOUSLY INVALID ENTRIES
# ============================================================

string_genes <- string_genes[
  !is.na(string_genes) &
    string_genes != ""
]


# ============================================================
# 13. SAVE ANNOTATED RESULTS
# ============================================================

write_xlsx(
  list(
    Annotated_Limma = deg_annotated,
    Candidate_Genes = candidate_annotated
  ),
  file.path(
    tables_dir,
    "04_Annotated_DEG_Results.xlsx"
  )
)


# ============================================================
# 14. EXPORT STRING GENE LIST
# ============================================================

write.table(
  string_genes,
  file.path(
    tables_dir,
    "Genes_for_STRING.txt"
  ),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)


# ============================================================
# 15. SUMMARY
# ============================================================

cat(
  "Annotated probes:",
  nrow(deg_annotated),
  "\n"
)

cat(
  "Exploratory candidate genes:",
  nrow(candidate_annotated),
  "\n"
)

cat(
  "Unique genes exported for STRING:",
  length(string_genes),
  "\n"
)


# ============================================================
# END OF SCRIPT
# ============================================================

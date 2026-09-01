# ============================================================
# PROJECT: Visceral Adipose Tissue Transcriptomics in T2D
# DATASET: GSE78721
# PLATFORM: GPL15207
# ANALYSIS: Normal vs Type 2 Diabetes
#
# STEP 01: DATA LOADING AND VISCERAL TISSUE FILTERING
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

dir.create(data_dir, showWarnings = FALSE)
dir.create(results_dir, showWarnings = FALSE)
dir.create(tables_dir, showWarnings = FALSE)


# ============================================================
# 03. LOAD GEO DATASET
# ============================================================

# Downloaded GEO series matrix file:
# GSE78721_series_matrix.txt.gz

gse <- getGEO(
  filename = file.path(
    data_dir,
    "GSE78721_series_matrix.txt.gz"
  ),
  getGPL = FALSE
)


# ============================================================
# 04. EXTRACT SAMPLE METADATA
# ============================================================

meta <- pData(gse)

# Inspect available metadata
colnames(meta)

# Tissue distribution
table(meta$`tissue:ch1`)

# Diagnosis distribution
table(meta$`diagnosis:ch1`)

# Gender distribution
table(meta$`gender:ch1`)


# ============================================================
# 05. EXTRACT EXPRESSION MATRIX
# ============================================================

expr <- exprs(gse)

cat(
  "Expression matrix dimensions:",
  nrow(expr),
  "probes x",
  ncol(expr),
  "samples\n"
)


# ============================================================
# 06. FILTER FOR VISCERAL ADIPOSE TISSUE
# ============================================================

visceral_samples <- rownames(meta)[
  meta$`tissue:ch1` == "Visceral"
]

cat(
  "Number of visceral adipose samples:",
  length(visceral_samples),
  "\n"
)

# Filter metadata
meta_visceral <- meta[visceral_samples, ]

# Filter expression matrix
expr_visceral <- expr[, visceral_samples]


# ============================================================
# 07. CREATE CLEAN METADATA
# ============================================================

meta_visceral_clean <- meta_visceral %>%
  dplyr::select(
    geo_accession,
    title,
    `tissue:ch1`,
    `diagnosis:ch1`,
    `gender:ch1`
  ) %>%
  dplyr::rename(
    Sample_ID = geo_accession,
    Title = title,
    Tissue = `tissue:ch1`,
    Diagnosis = `diagnosis:ch1`,
    Gender = `gender:ch1`
  )


# ============================================================
# 08. CHECK FILTERED DATA
# ============================================================

cat("\nDiagnosis distribution:\n")
print(table(meta_visceral_clean$Diagnosis))

cat("\nGender distribution:\n")
print(table(meta_visceral_clean$Gender))

cat("\nGender x Diagnosis:\n")
print(
  table(
    meta_visceral_clean$Gender,
    meta_visceral_clean$Diagnosis
  )
)


# ============================================================
# 09. SAVE FILTERED DATA
# ============================================================

expr_visceral_df <- as.data.frame(expr_visceral)

expr_visceral_df <- cbind(
  Gene_ID = rownames(expr_visceral_df),
  expr_visceral_df
)

rownames(expr_visceral_df) <- NULL

write_xlsx(
  list(
    Metadata = meta_visceral_clean,
    Expression_Matrix = expr_visceral_df
  ),
  file.path(
    tables_dir,
    "01_GSE78721_Visceral_Data.xlsx"
  )
)


# ============================================================
# END OF SCRIPT
# ============================================================

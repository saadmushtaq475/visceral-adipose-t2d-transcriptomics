# ============================================================
# PROJECT: Visceral Adipose Tissue T2D Transcriptomics
# PURPOSE: Install required R packages
# ============================================================

# CRAN packages
install.packages(c(
  "dplyr",
  "ggplot2",
  "pheatmap",
  "writexl",
  "readxl"
))

# Install BiocManager if required
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Bioconductor packages
BiocManager::install(c(
  "GEOquery",
  "limma",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
))

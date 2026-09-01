# ============================================================
# PROJECT: Visceral Adipose Tissue Transcriptomics in T2D
# DATASET: GSE78721
# PLATFORM: GPL15207
#
# STEP 02: QUALITY CONTROL AND EXPLORATORY ANALYSIS
# ============================================================


# ============================================================
# 01. LOAD REQUIRED PACKAGES
# ============================================================

library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)
library(pheatmap)


# ============================================================
# 02. PROJECT DIRECTORIES
# ============================================================

data_dir <- "data"
results_dir <- "results"
figures_dir <- file.path(results_dir, "figures")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 03. LOAD GEO DATA
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

meta_visceral <- meta[visceral_samples, ]
expr_visceral <- expr[, visceral_samples]


# ============================================================
# 05. CLEAN METADATA
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
# 06. EXPRESSION DISTRIBUTION — BOXPLOT
# ============================================================

png(
  filename = file.path(
    figures_dir,
    "01_expression_boxplot.png"
  ),
  width = 2400,
  height = 1600,
  res = 300
)

boxplot(
  expr_visceral,
  las = 2,
  outline = FALSE,
  main = "Expression Distribution of Visceral Adipose Samples",
  xlab = "Samples",
  ylab = "Expression"
)

dev.off()


# ============================================================
# 07. EXPRESSION DENSITY
# ============================================================

png(
  filename = file.path(
    figures_dir,
    "02_expression_density.png"
  ),
  width = 2400,
  height = 1600,
  res = 300
)

plotDensities(
  expr_visceral,
  legend = FALSE,
  main = "Expression Density of Visceral Adipose Samples"
)

dev.off()


# ============================================================
# 08. PRINCIPAL COMPONENT ANALYSIS (PCA)
# ============================================================

pca <- prcomp(
  t(expr_visceral),
  scale. = TRUE
)

pca_df <- as.data.frame(pca$x)

pca_df$Diagnosis <- meta_visceral_clean$Diagnosis
pca_df$Gender <- meta_visceral_clean$Gender
pca_df$Sample_ID <- meta_visceral_clean$Sample_ID


# Calculate percentage variance
percent_variance <- 100 * (pca$sdev^2 / sum(pca$sdev^2))


pca_plot <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = Diagnosis
  )
) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCA of Visceral Adipose Tissue Samples",
    x = paste0(
      "PC1 (",
      round(percent_variance[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(percent_variance[2], 1),
      "%)"
    )
  )


ggsave(
  filename = file.path(
    figures_dir,
    "03_PCA.png"
  ),
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 09. SAMPLE CORRELATION
# ============================================================

cor_matrix <- cor(expr_visceral)


# Create sample annotation
annotation_col <- data.frame(
  Diagnosis = meta_visceral_clean$Diagnosis,
  Gender = meta_visceral_clean$Gender
)

rownames(annotation_col) <-
  meta_visceral_clean$Sample_ID


png(
  filename = file.path(
    figures_dir,
    "04_sample_correlation_heatmap.png"
  ),
  width = 2400,
  height = 2000,
  res = 300
)

pheatmap(
  cor_matrix,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  main = "Sample Correlation Heatmap"
)

dev.off()


# ============================================================
# 10. BATCH INFORMATION ASSESSMENT
# ============================================================

# Inspect metadata columns potentially related to batch,
# series, source, or platform.

batch_columns <- meta_visceral[
  ,
  grep(
    "batch|series|source|platform",
    colnames(meta_visceral),
    ignore.case = TRUE
  ),
  drop = FALSE
]

print(batch_columns)


# ============================================================
# 11. EXPORT PCA DATA
# ============================================================

write.csv(
  pca_df,
  file.path(
    results_dir,
    "PCA_sample_coordinates.csv"
  ),
  row.names = TRUE
)


# ============================================================
# END OF SCRIPT
# ============================================================

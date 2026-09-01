# ============================================================
# PROJECT: Visceral Adipose Tissue Transcriptomics in T2D
# DATASET: GSE78721
# PLATFORM: GPL15207
#
# STEP 03: DIFFERENTIAL EXPRESSION ANALYSIS
# METHOD: limma
# COMPARISON: Type 2 Diabetes vs Normal
# ============================================================


# ============================================================
# 01. LOAD REQUIRED PACKAGES
# ============================================================

library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)
library(writexl)


# ============================================================
# 02. PROJECT DIRECTORIES
# ============================================================

data_dir <- "data"
results_dir <- "results"
figures_dir <- file.path(results_dir, "figures")
tables_dir <- file.path(results_dir, "tables")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)


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
# 05. CREATE DIAGNOSIS VARIABLE
# ============================================================

group <- factor(
  meta_visceral$`diagnosis:ch1`,
  levels = c(
    "Normal",
    "Type-2 Diabetes"
  )
)

print(table(group))


# ============================================================
# 06. DESIGN MATRIX
# ============================================================

design <- model.matrix(
  ~ group
)

colnames(design) <- c(
  "Normal",
  "T2D_vs_Normal"
)

print(design)


# ============================================================
# 07. RUN LIMMA
# ============================================================

fit <- lmFit(
  expr_visceral,
  design
)

fit <- eBayes(
  fit
)


# ============================================================
# 08. EXTRACT DIFFERENTIAL EXPRESSION RESULTS
# ============================================================

deg_results <- topTable(
  fit,
  coef = "T2D_vs_Normal",
  number = Inf,
  sort.by = "P"
)


# ============================================================
# 09. BASIC STATISTICAL SUMMARY
# ============================================================

number_fdr_significant <- sum(
  deg_results$adj.P.Val < 0.05,
  na.rm = TRUE
)

minimum_fdr <- min(
  deg_results$adj.P.Val,
  na.rm = TRUE
)

number_fc_fdr <- sum(
  abs(deg_results$logFC) >= 1 &
    deg_results$adj.P.Val < 0.05,
  na.rm = TRUE
)

cat(
  "Number of FDR-significant genes:",
  number_fdr_significant,
  "\n"
)

cat(
  "Minimum adjusted P-value:",
  minimum_fdr,
  "\n"
)

cat(
  "Number with |log2FC| >= 1 and FDR < 0.05:",
  number_fc_fdr,
  "\n"
)


# ============================================================
# 10. STRICT FDR-SIGNIFICANT DEGs
# ============================================================

FDR_DEGs <- deg_results %>%
  filter(
    abs(logFC) >= 1,
    adj.P.Val < 0.05
  )


# ============================================================
# 11. EXPLORATORY CANDIDATE GENES
# ============================================================
#
# IMPORTANT:
# These genes are NOT called statistically significant DEGs.
#
# They meet:
#   |log2FC| >= 1
#   P-value < 0.05
#
# but do not necessarily meet FDR < 0.05.
#
# They are therefore used as exploratory candidates for
# downstream network analysis.


candidate_genes <- deg_results %>%
  filter(
    abs(logFC) >= 1,
    P.Value < 0.05
  )


cat(
  "FDR-significant DEGs:",
  nrow(FDR_DEGs),
  "\n"
)

cat(
  "Exploratory candidate genes:",
  nrow(candidate_genes),
  "\n"
)


# ============================================================
# 12. SAVE DIFFERENTIAL EXPRESSION RESULTS
# ============================================================

write_xlsx(
  list(
    Complete_Limma = deg_results,
    FDR_DEGs = FDR_DEGs,
    Exploratory_Candidates = candidate_genes
  ),
  file.path(
    tables_dir,
    "02_DEG_Results.xlsx"
  )
)


# ============================================================
# 13. PREPARE VOLCANO PLOT
# ============================================================

deg_results$Category <- "Not Significant"

deg_results$Category[
  deg_results$logFC >= 1 &
    deg_results$adj.P.Val < 0.05
] <- "Upregulated"

deg_results$Category[
  deg_results$logFC <= -1 &
    deg_results$adj.P.Val < 0.05
] <- "Downregulated"


# ============================================================
# 14. VOLCANO PLOT
# ============================================================

volcano_plot <- ggplot(
  deg_results,
  aes(
    x = logFC,
    y = -log10(P.Value),
    color = Category
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  theme_minimal() +
  labs(
    title = "Differential Expression: T2D vs Normal",
    subtitle = "Visceral Adipose Tissue — GSE78721",
    x = "log2 Fold Change",
    y = "-log10(P-value)",
    color = "Category"
  )


# ============================================================
# 15. SAVE VOLCANO PLOT
# ============================================================

ggsave(
  filename = file.path(
    figures_dir,
    "05_volcano_plot.png"
  ),
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 16. EXPORT TOP EXPLORATORY CANDIDATES
# ============================================================

top_candidates <- candidate_genes %>%
  arrange(P.Value)

write_xlsx(
  top_candidates,
  file.path(
    tables_dir,
    "03_exploratory_candidate_genes.xlsx"
  )
)


# ============================================================
# END OF SCRIPT
# ============================================================

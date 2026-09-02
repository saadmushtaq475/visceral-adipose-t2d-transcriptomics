# Results

This directory contains the main outputs generated from the transcriptomic,
functional enrichment, network analysis, target prioritization, and
drug-repurposing workflow.

## Directory Structure

```text
results/
├── README.md
├── figures/
└── tables/


### One important point

I deliberately used **"selected graphical outputs"** rather than saying that all figures are stored here. That's consistent with your decision to keep only the Volcano plot while leaving the plotting code in the R scripts.

Also, the statement about **no FDR-significant DEGs** is supported by your actual analysis code, which explicitly separates `FDR_DEGs` from exploratory `candidate_genes`.

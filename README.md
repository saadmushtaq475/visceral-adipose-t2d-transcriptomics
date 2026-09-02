# Visceral Adipose Tissue Transcriptomics & Drug Repurposing in Insulin Resistance

## Overview
This repository contains the full computational pipeline for identifying hub genes
associated with insulin resistance in visceral adipose tissue (Type 2 Diabetes vs.
Normal), and a structure-based drug repurposing screen against these targets using
FDA-approved compounds.

The study uses publicly available gene expression data (GEO accession **GSE78721**)
and combines differential expression analysis, functional enrichment, protein-protein
interaction network analysis, and molecular docking to prioritize both disease-relevant
genes and candidate repurposable drugs.

## Pipeline

GEO (GSE78721) → Visceral tissue filtering → Quality control → Differential expression
(limma) → Gene annotation (GPL15207) → GO/KEGG enrichment → STRING PPI network →
MCODE + cytoHubba hub gene analysis → Hub gene validation → Druggability assessment →
Drug-target interaction mining → Molecular docking (AutoDock Vina) → Candidate selection

## Key Results

| Stage | Result |
|---|---|
| Samples analyzed | 35 visceral adipose samples (16 Normal, 19 T2D) |
| Candidate genes | 139 (126 upregulated, 13 downregulated) |
| FDR-significant DEGs | None (adjusted P < 0.05) |
| Hub genes identified | COL1A1, VTN, TIMP1, VCAM1, THY1, KRT19, KRT18, KRT8 |
| Prioritized docking targets | VCAM1, VTN, TIMP1 |
| Final drug candidates | Lumacaftor (TIMP1), Flufenamic acid (VTN) |

The candidate genes were selected using nominal statistical significance
(P < 0.05) and fold-change criteria. Because no genes passed the predefined
FDR-adjusted significance threshold, downstream network and drug-repurposing
analyses should be considered exploratory.


## Repository Structure

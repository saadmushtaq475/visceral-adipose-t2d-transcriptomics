# Integrated Transcriptomic Analysis and Drug Repurposing of Visceral Adipose Tissue in Type 2 Diabetes

An end-to-end bioinformatics workflow for identifying transcriptional signatures, functional pathways, network hub genes, candidate therapeutic targets, and drug-repurposing candidates associated with Type 2 Diabetes (T2D) in visceral adipose tissue.

---

## Project Overview

Type 2 Diabetes (T2D) is a complex metabolic disorder involving insulin resistance, inflammation, extracellular matrix remodeling, and altered adipose tissue biology.

This project investigates gene-expression changes in **visceral adipose tissue** from individuals with normal glucose status and Type 2 Diabetes using publicly available transcriptomic data.

The workflow integrates:

**GEO → Quality Control → Differential Expression → Probe Annotation → GO/KEGG → STRING → Cytoscape → MCODE → CytoHubba → Hub-Gene Validation → Target Prioritization → Drug Repurposing → Molecular Docking**

The project was developed as a reproducible bioinformatics research and portfolio project using R, Bioconductor, Cytoscape, and DrugRep.

---

## Research Objective

The main objective was to identify biologically relevant candidate genes and therapeutic targets associated with T2D in visceral adipose tissue and subsequently investigate potential drug-repurposing opportunities through receptor-based virtual screening.

---

## Dataset

| Feature | Information |
|---|---|
| GEO accession | **GSE78721** |
| Platform | **GPL15207** |
| Organism | *Homo sapiens* |
| Tissue | Visceral adipose tissue |
| Comparison | Normal vs Type 2 Diabetes |
| Analysis type | Transcriptomic gene-expression analysis |

The original GEO dataset was filtered to retain samples corresponding specifically to visceral adipose tissue.

---

## Analysis Workflow

### 1. Data Acquisition and Preprocessing

- Downloaded and processed the GSE78721 series matrix.
- Extracted sample metadata.
- Identified visceral adipose tissue samples.
- Constructed a clean sample metadata table.
- Prepared the expression matrix for downstream analysis.

### 2. Quality Control

Quality-control analyses included:

- Expression distribution assessment
- Density plots
- Principal Component Analysis (PCA)
- Sample-to-sample correlation analysis
- Correlation heatmap
- Assessment of potential batch-related variables

### 3. Differential Expression Analysis

Differential expression was performed using the **limma** framework.

The primary comparison was:

**Type 2 Diabetes vs Normal**

Two categories were retained:

- **FDR-significant DEGs:** |log2FC| ≥ 1 and adjusted P-value < 0.05
- **Exploratory candidate genes:** |log2FC| ≥ 1 and nominal P-value < 0.05

Importantly, nominally significant genes were **not** treated as FDR-significant DEGs.

### 4. Probe-to-Gene Annotation

Probe identifiers were mapped to gene symbols using the authoritative **GPL15207** platform annotation.

### 5. Functional Enrichment

Candidate genes were investigated using:

- Gene Ontology (GO)
  - Biological Process
  - Molecular Function
  - Cellular Component
- KEGG pathway enrichment

The enrichment analysis highlighted biological themes involving extracellular matrix organization, cell adhesion, immune processes, and related adipose-tissue biology.

### 6. Protein-Protein Interaction Network

Candidate genes were submitted to **STRING** to construct a protein-protein interaction network.

The resulting network was subsequently analyzed using **Cytoscape**, including:

- MCODE
- CytoHubba
- Degree centrality
- Betweenness centrality
- Closeness centrality

Multiple network-based approaches were used to identify highly connected and central genes.

---

## Hub Gene Identification

The final hub-gene set consisted of:

- **COL1A1**
- **VTN**
- **TIMP1**
- **VCAM1**
- **THY1**
- **KRT19**
- **KRT18**
- **KRT8**

These genes represent network-derived candidate hubs and were subsequently evaluated using their transcriptomic evidence and biological relevance.

### Important Statistical Interpretation

The analysis did **not** identify genes meeting the conventional genome-wide FDR threshold of adjusted P-value < 0.05.

Therefore, the hub genes are not described as statistically significant DEGs.

Instead, they were interpreted as **network-derived candidate genes**, with nominal expression evidence considered separately from FDR significance.

This distinction was maintained throughout the analysis to avoid overinterpretation of the results.

---

## Target Prioritization

The eight hub genes were evaluated according to:

- Network centrality
- MCODE participation
- Expression direction and magnitude
- Nominal statistical evidence
- Functional enrichment
- Biological relevance to T2D
- Protein structure availability
- Druggability
- Suitability for structure-based analysis

Three targets were prioritized for receptor-based virtual screening:

| Target | PDB structure | Priority |
|---|---|---|
| **TIMP1** | **9SOS** | Primary |
| **VTN** | **1S4G** | Primary |
| **VCAM1** | **1VCA** | Primary screening target |

---

# Drug Repurposing

## DrugRep Receptor-Based Virtual Screening

The prioritized protein structures were subjected to receptor-based virtual screening using **DrugRep**.

The available drug library was screened independently against each selected protein structure.

Hits were ranked according to predicted Vina docking scores and subsequently evaluated using additional criteria including:

- Drug/therapeutic status
- Chemical properties
- Molecular size
- Lipinski-related characteristics
- ADME/toxicity considerations
- Binding-site interpretation
- Overall suitability for drug repurposing

---

## Top Screening Results

### TIMP1 — PDB 9SOS

| Rank | Drug | Vina score |
|---:|---|---:|
| 1 | Conivaptan | −9.8 |
| 2 | **Lumacaftor** | **−9.5** |
| 3 | Entrectinib | −9.5 |
| 4 | Nilotinib | −9.4 |
| 5 | Rimegepant | −9.3 |

Following candidate filtering and prioritization, **Lumacaftor** was retained as the primary computational candidate for TIMP1.

### Final candidate

**TIMP1 → Lumacaftor**

Predicted Vina score:

**−9.5 kcal/mol**

---

### VTN — PDB 1S4G

| Rank | Drug | Vina score |
|---:|---|---:|
| 1 | **Flufenamic acid** | **−9.6** |
| 2 | Triclocarban | −9.5 |
| 3 | Isocarboxazid | −9.3 |
| 4 | Ensulizole | −9.2 |
| 5 | Granisetron | −8.9 |

After filtering unsuitable compounds, **Flufenamic acid** was retained as the primary computational candidate for VTN.

### Final candidate

**VTN → Flufenamic acid**

Predicted Vina score:

**−9.6 kcal/mol**

---

### VCAM1 — PDB 1VCA

| Rank | Drug | Vina score |
|---:|---|---:|
| 1 | Cyproheptadine | −7.0 |
| 2 | Quinupramine | −6.9 |
| 3 | Dantron | −6.9 |
| 4 | Anthralin | −6.8 |
| 5 | Morphine | −6.8 |

No VCAM1 compound was retained as a confident final repurposing candidate.

The predicted compounds showed limited target-specific biological rationale, and the identified binding region may not represent the most biologically relevant VCAM1 interaction surface.

Therefore, VCAM1 was classified as **inconclusive** and is recommended for binding-site refinement and future re-docking.

---

## Final Computational Candidates

| Target | Drug | PDB | Vina score | Status |
|---|---|---|---:|---|
| **TIMP1** | **Lumacaftor** | 9SOS | −9.5 kcal/mol | Computationally prioritized |
| **VTN** | **Flufenamic acid** | 1S4G | −9.6 kcal/mol | Computationally prioritized |
| VCAM1 | None selected | 1VCA | −7.0 kcal/mol top hit | Inconclusive |

These candidates should be regarded as **computationally prioritized drug-repurposing hypotheses**, not experimentally validated drug-target interactions.

---

# Molecular Interaction / Pose Validation

The final drug-target pairs selected for structural evaluation were:

- **Lumacaftor–TIMP1**
- **Flufenamic acid–VTN**

The structural analysis was used to evaluate the predicted docking poses and molecular interactions within the respective protein structures.

Docking results provide hypotheses regarding potential molecular interactions and require experimental validation.

---

# Key Findings

1. Visceral adipose tissue samples from GSE78721 were successfully analyzed using an integrated transcriptomic workflow.
2. No genes passed the conventional FDR threshold of adjusted P-value < 0.05.
3. Nominally associated genes were therefore treated as exploratory candidates rather than statistically significant DEGs.
4. Network analysis identified eight hub genes:
   **COL1A1, VTN, TIMP1, VCAM1, THY1, KRT19, KRT18, and KRT8.**
5. Functional analysis highlighted processes related to extracellular matrix organization, cell adhesion, and immune/inflammatory biology.
6. TIMP1, VTN, and VCAM1 were selected for receptor-based virtual screening.
7. DrugRep screening identified **Lumacaftor** as the prioritized candidate for TIMP1.
8. DrugRep screening identified **Flufenamic acid** as the prioritized candidate for VTN.
9. VCAM1 screening was considered inconclusive because of weak target-specific interpretation of the predicted hits.
10. The final results provide computational hypotheses that require experimental validation.

---

# Limitations

- The study uses a relatively limited number of visceral adipose tissue samples.
- No genes reached the conventional FDR significance threshold.
- Candidate-gene selection therefore involves exploratory nominal statistical evidence and network-based prioritization.
- PPI networks are database-derived and may contain indirect or context-independent interactions.
- Docking scores are computational predictions and do not demonstrate experimental binding.
- DrugRep virtual screening depends on the selected protein structure and predicted binding site.
- The VCAM1 binding-site result requires further structural validation.
- The proposed drug-target relationships require biochemical, cellular, and ultimately clinical investigation before therapeutic conclusions can be made.

---

# Reproducibility

The analysis scripts are organized sequentially in the `R/` directory:

```text
00_setup.R
01_data_loading.R
02_quality_control.R
03_differential_expression.R
04_annotation.R
05_enrichment.R
06_hub_gene_validation.R
07_target_prioritization.R
08_drug_repurposing.R
09_molecular_interaction_validation.R

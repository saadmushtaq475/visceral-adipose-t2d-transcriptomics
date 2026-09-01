# ============================================================
# 08. DRUG REPURPOSING
# ============================================================
#
# Project:
# Integrated Bioinformatics Analysis and Drug Repurposing
# of Visceral Adipose Tissue Gene Expression in Type 2 Diabetes
#
# Method:
# Receptor-based virtual screening using DrugRep
#
# Prioritized targets:
# TIMP1, VTN, VCAM1
#
# ============================================================


# ============================================================
# 1. TARGET PROTEINS AND PDB STRUCTURES
# ============================================================

target_structures <- data.frame(
  Target = c("TIMP1", "VTN", "VCAM1"),
  PDB_ID = c("9SOS", "1S4G", "1VCA"),
  stringsAsFactors = FALSE
)

target_structures


# ============================================================
# 2. DRUGREP RECEPTOR-BASED VIRTUAL SCREENING
# ============================================================

# DrugRep was used for receptor-based virtual screening
# against the available drug library.
#
# Screening was performed independently for each selected
# protein target using the corresponding PDB structure.
#
# DrugRep output was ranked primarily according to the
# predicted Vina docking score.


# ============================================================
# 3. TOP DRUGREP SCREENING RESULTS
# ============================================================

# TIMP1 — PDB 9SOS

TIMP1_top_hits <- data.frame(
  Rank = 1:5,
  Drug = c(
    "Conivaptan",
    "Lumacaftor",
    "Entrectinib",
    "Nilotinib",
    "Rimegepant"
  ),
  Vina_Score = c(
    -9.8,
    -9.5,
    -9.5,
    -9.4,
    -9.3
  )
)

TIMP1_top_hits


# VTN — PDB 1S4G

VTN_top_hits <- data.frame(
  Rank = 1:5,
  Drug = c(
    "Flufenamic acid",
    "Triclocarban",
    "Isocarboxazid",
    "Ensulizole",
    "Granisetron"
  ),
  Vina_Score = c(
    -9.6,
    -9.5,
    -9.3,
    -9.2,
    -8.9
  )
)

VTN_top_hits


# VCAM1 — PDB 1VCA

VCAM1_top_hits <- data.frame(
  Rank = 1:5,
  Drug = c(
    "Cyproheptadine",
    "Quinupramine",
    "Dantron",
    "Anthralin",
    "Morphine"
  ),
  Vina_Score = c(
    -7.0,
    -6.9,
    -6.9,
    -6.8,
    -6.8
  )
)

VCAM1_top_hits


# ============================================================
# 4. FILTERING OF SCREENING HITS
# ============================================================

# Initial filtering considered:
#
# 1. Drug/therapeutic status
# 2. Suitability for drug repurposing
# 3. Molecular size and chemical properties
# 4. Lipinski-related characteristics
# 5. ADME/toxicity considerations
# 6. Biological relevance to the selected target
# 7. Quality and interpretability of the predicted binding site
#
# Compounds such as cosmetic compounds, UV filters,
# preservatives, and unsuitable biologics were not retained
# as repurposing candidates.


# ============================================================
# 5. PRIORITIZED DRUG-TARGET PAIRS
# ============================================================

final_candidates <- data.frame(
  Target = c("TIMP1", "VTN"),
  PDB_ID = c("9SOS", "1S4G"),
  Drug = c(
    "Lumacaftor",
    "Flufenamic acid"
  ),
  Vina_Score = c(
    -9.5,
    -9.6
  ),
  Priority = c(
    "Primary candidate",
    "Primary candidate"
  ),
  Confidence = c(
    "Moderate",
    "Moderate"
  ),
  stringsAsFactors = FALSE
)

final_candidates


# ============================================================
# 6. VCAM1 INTERPRETATION
# ============================================================

# VCAM1 produced relatively weaker top screening scores
# compared with TIMP1 and VTN.
#
# The identified compounds were predominantly tricyclic
# amines and opioid/antihistamine-like compounds, without
# convincing target-specific biological rationale.
#
# Therefore, no VCAM1 drug was selected as a final
# repurposing candidate from this screening run.
#
# VCAM1 should be considered inconclusive and may require
# binding-site refinement and/or re-docking in future work.


# ============================================================
# 7. FINAL DRUG REPURPOSING TABLE
# ============================================================

write.csv(
  final_candidates,
  "R/data/results/tables/08_final_drug_repurposing_candidates.csv",
  row.names = FALSE
)

write.csv(
  TIMP1_top_hits,
  "R/data/results/tables/08_TIMP1_top_hits.csv",
  row.names = FALSE
)

write.csv(
  VTN_top_hits,
  "R/data/results/tables/08_VTN_top_hits.csv",
  row.names = FALSE
)

write.csv(
  VCAM1_top_hits,
  "R/data/results/tables/08_VCAM1_top_hits.csv",
  row.names = FALSE
)


# ============================================================
# 8. FINAL INTERPRETATION
# ============================================================

# Final computationally prioritized candidates:
#
# TIMP1 — Lumacaftor
# Vina score = -9.5 kcal/mol
#
# VTN — Flufenamic acid
# Vina score = -9.6 kcal/mol
#
# These compounds represent computationally prioritized
# drug-repurposing candidates and require experimental
# validation.
#
# Docking score alone does not establish:
# - direct experimental binding
# - target inhibition
# - therapeutic efficacy
# - cellular activity
# - clinical effectiveness


# ============================================================
# END OF DRUG REPURPOSING ANALYSIS
# ============================================================

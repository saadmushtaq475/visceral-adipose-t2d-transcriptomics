# ============================================================
# 09. MOLECULAR INTERACTION / POSE VALIDATION
# ============================================================
#
# Project:
# Integrated Bioinformatics Analysis and Drug Repurposing
# of Visceral Adipose Tissue Gene Expression in Type 2 Diabetes
#
# Purpose:
# Document structural validation of the prioritized
# drug-target docking results obtained after DrugRep
# receptor-based virtual screening.
#
# Prioritized drug-target pairs:
#
# TIMP1  -> Lumacaftor
# VTN    -> Flufenamic acid
#
# VCAM1 was not retained as a final docking candidate because
# the receptor-based screening result did not provide a
# sufficiently convincing biologically relevant binding site.
#
# ============================================================


# ============================================================
# 1. PRIORITIZED TARGETS
# ============================================================

final_targets <- data.frame(
  Target = c("TIMP1", "VTN"),
  PDB_ID = c("9SOS", "1S4G"),
  Drug = c("Lumacaftor", "Flufenamic acid"),
  stringsAsFactors = FALSE
)

final_targets


# ============================================================
# 2. STRUCTURAL VALIDATION INFORMATION
# ============================================================

# TIMP1
# PDB: 9SOS
# Drug: Lumacaftor
#
# VTN
# PDB: 1S4G
# Drug: Flufenamic acid


# ============================================================
# 3. DRUGREP VIRTUAL SCREENING SCORES
# ============================================================

docking_results <- data.frame(
  Target = c("TIMP1", "VTN"),
  Drug = c("Lumacaftor", "Flufenamic acid"),
  PDB_ID = c("9SOS", "1S4G"),
  Vina_Score_kcal_mol = c(-9.5, -9.6),
  stringsAsFactors = FALSE
)

docking_results


# ============================================================
# 4. FINAL INTERPRETATION
# ============================================================

final_interpretation <- data.frame(
  Target = c("TIMP1", "VTN"),
  Drug = c("Lumacaftor", "Flufenamic acid"),
  Interpretation = c(
    "Strong predicted binding from receptor-based virtual screening; requires experimental validation.",
    "Strong predicted binding from receptor-based virtual screening; requires experimental validation."
  ),
  stringsAsFactors = FALSE
)

final_interpretation


# ============================================================
# 5. EXPORT FINAL DOCKING SUMMARY
# ============================================================

write.csv(
  docking_results,
  "R/data/results/tables/09_final_docking_results.csv",
  row.names = FALSE
)

write.csv(
  final_interpretation,
  "R/data/results/tables/09_final_docking_interpretation.csv",
  row.names = FALSE
)


# ============================================================
# 6. IMPORTANT SCIENTIFIC NOTE
# ============================================================

# Docking scores represent predicted binding affinity.
#
# A favorable docking score does NOT establish:
# - direct biological binding,
# - target inhibition,
# - therapeutic efficacy,
# - cellular activity,
# - clinical effectiveness.
#
# Therefore, Lumacaftor-TIMP1 and Flufenamic acid-VTN
# are reported as computationally prioritized repurposing
# candidates requiring experimental validation.


# ============================================================
# END OF ANALYSIS
# ============================================================

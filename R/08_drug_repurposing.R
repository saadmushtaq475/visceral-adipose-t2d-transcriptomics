# ============================================================
# 08. DRUG REPURPOSING
# Receptor-based virtual screening using DrugRep
# ============================================================

# Selected targets:
# 1. TIMP1
# 2. VTN
# 3. VCAM1

# Protein structures:
# TIMP1  -> PDB: 9SOS
# VTN    -> PDB: 1S4G
# VCAM1  -> PDB: 1VCA

# DrugRep was used for receptor-based virtual screening
# against the available drug library.
#
# Top-ranked compounds were subsequently filtered based on:
# - Drug status
# - Chemical suitability
# - Lipinski properties
# - ADME/Toxicity considerations
# - Binding score
#
# Final candidates:
# TIMP1 -> Lumacaftor
# VTN   -> Flufenamic acid
# VCAM1 -> Inconclusive; requires binding-site validation

# DrugRep output tables should be stored in:
# R/data/results/tables/

# DrugRep figures should be stored in:
# R/data/results/figures/

## Data

The analysis uses the publicly available GEO dataset GSE78721
and platform GPL15207.

The repository includes selected processed/input files required
for downstream analysis and molecular interaction validation.

### Included data

- `01_GSE78721_Visceral_Data.xlsx` — filtered visceral adipose tissue dataset.
- `1S4G.pdb` — VTN protein structure used for molecular docking.
- `1VCA.pdb` — VCAM1 protein structure used for molecular docking.
- `9SOS.pdb` — TIMP1 protein structure used for molecular docking.

### Data not included

The original GSE78721 expression matrix and GPL15207 platform
annotation files are not included because of their large file size.
These data are publicly available through NCBI GEO and can be
retrieved using the analysis workflow provided in the `R/` directory.

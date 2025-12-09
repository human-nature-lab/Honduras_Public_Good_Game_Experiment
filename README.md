# Honduras_Public_Good_Game_Experiment
Analysis and replication code for “Group Cooperation Diverges onto Durable Low versus High Paths: Public Goods Experiments in 134 Honduran Villages,” reproducing all results in the main text and Supplementary Information. De-identified data available on request from the Human Nature Lab.

# Replication Code — Honduras Public Goods Game (PGG)

This repository contains the analysis scripts used for the paper:

Papamichalis, M., Christakis, N. A., & Fu, F. (2025).
"Group Cooperation Diverges onto Durable Low versus High Paths:
Public Goods Experiments in 134 Honduran Villages."

## Repository Structure
.
├── Main_Text/
│   ├── Plot_Figure_3.R
│   ├── Section_3_3_Initial_Analysis.R
│   ├── Section_5_1_Moran_Process.R
│   ├── Section_5_2_Clustering.R
│   ├── Section_5_2_Heterogeneous_Analysis.R
│   └── Section_5_3_IV_Analysis.R
└── Supplementary Information/
    ├── SI_B.R
    ├── SI_B1.R
    ├── SI_C1_F.R
    ├── SI_C2.R
    ├── SI_D.R
    ├── SI_E2.R
    ├── SI_G_Robustness_Sensitivity.R
    ├── SI_G4_Alpha_Robustness.R
    ├── SI_H_Plots.R
    ├── SI_I.R
    ├── SI_J_Calibration.R
    └── SI_L.R

All scripts expect a single CSV file named:

- `data_set.csv`

Each script loads the dataset autonomously using:

```r
pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)



# Replication Code — Honduras Public Goods Game (PGG)

Important: Because scripts use the relative path "data_set.csv", you must set your
working directory to the repository root (the folder containing this CSV).
Each .R file is written to run autonomously, without depending on any other script.

## How to Run

### Option A: Command Line

From the repository root:

Rscript Main_Text/Section_3_3_Initial_Analysis.R
Rscript Main_Text/Section_5_1_Moran_Process.R
Rscript Main_Text/Section_5_2_Clustering.R
Rscript Main_Text/Section_5_2_Heterogeneous_Analysis.R
Rscript Main_Text/Section_5_3_IV_Analysis.R
Rscript Main_Text/Plot_Figure_3.R

Rscript "Supplementary Information/SI_B.R"
# ...repeat for other SI_*.R files as needed.

### Option B: RStudio

1. Open the repository as an RStudio Project
   (or set your working directory to the repository root).
2. Open any .R file and click "Source" (or run it line-by-line).

## Autonomous Execution

- Each .R file runs standalone.
- No script depends on output from any other script.
- All scripts load data_set.csv internally.

## Script Map

### Main_Text/

- Section_3_3_Initial_Analysis.R — Initial empirical analyses.
- Section_5_1_Moran_Process.R — Evolutionary modeling & Moran calibration.
- Section_5_2_Clustering.R — Clustering analyses for contribution-path heterogeneity.
- Section_5_2_Heterogeneous_Analysis.R — Heterogeneous-effects models.
- Section_5_3_IV_Analysis.R — Instrumental-variables analyses for peer effects.
- Plot_Figure_3.R — Code reproducing Figure 3.

### Supplementary Information/

Scripts SI_*.R reproduce all Supplementary Information figures, tables, and robustness
checks. Filenames correspond to appendix sections (e.g., SI_J_Calibration.R = Appendix J).

## Dependencies

Each script loads the R packages it needs via library(...).
If a package is missing, install it in R with:

install.packages("PACKAGE_NAME")

## Outputs

Figures, tables, and intermediate results are saved to disk according to output paths
defined inside each script. Create any required folders if they do not already exist.

## Data Access & Contact

**The dataset is not included in this repository.**
- Data are available upon request to the Human Nature Lab (Yale University).
- Author contact emails appear in the manuscript front matter.
- The manuscript notes that a fully de-identified replication package (data + scripts)
  will be made public in an AEA-compliant repository prior to publication.

## Citation

If you use this code, please cite:

Papamichalis, M., Christakis, N. A., & Fu, F. (2025).
"Group Cooperation Diverges onto Durable Low versus High Paths:
Public Goods Experiments in 134 Honduran Villages."

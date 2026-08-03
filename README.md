
# SwitchClass
<img src="inst/logo.png" align="right" width="110" height="110" />

`SwitchClass` is an R package for quantifying baseline-aligned and perturbation-aligned molecular features across biological perturbations using a label-switch classification framework.
It provides a unified workflow for identifying *baseline-aligned* versus *perturbation-aligned* molecular features in longitudinal or comparative omics datasets.

The package was developed as part of the study **Dissecting Molecular Feature Alignment via a Label-Switch Classification Framework**, which systematically maps molecular trajectories across biological or therapeutic conditions.

---

## Overview

<img src="inst/framework.png" align="right" width="440" height="220" />

`SwitchClass` implements the following core components:

1. **Label-switch classification**  
   A random-forest-based approach that compares importance profiles between two label schemes to compute a per-feature *directional importance score* (`delta = importance_baseline - importance_perturbation`).

2. **Visualization utilities**  
   Functions for visualizing molecular states via UMAP embeddings, quadrant-based scatterplots, feature-level boxplots, and annotated heatmaps.

3. **Downstream interpretation**  
   Tools for pathway enrichment (Reactome) and quadrant-based biological annotation.

4. **Example datasets and vignettes**  
   Includes demonstration data from colorectal cancer (CRC).

---

## Installation

Install the development version from GitHub:

```r
# install dependencies if not already installed
# BiocManager::install("PhosR")
# BiocManager::install("reactome.db")
# BiocManager::install("org.Hs.eg.db")

devtools::install_github("PYangLab/SwitchClass",
                         build_vignettes = TRUE,
                         dependencies = TRUE)

```

## Vignette 

Please find our vignette: 
```r
browseVignettes("SwitchClass")
```

## Contact us

If you have any enquiries about SwitchClass, please contact d.xiao@sydney.edu.au. We are also happy to receive any suggestions and comments.

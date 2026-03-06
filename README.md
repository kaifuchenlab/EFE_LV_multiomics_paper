# EFE_multiomics

Single-cell multiomics data analysis for the **Endocardial Fibroelastosis (EFE)** multiomics project.

---

## Overview

This repository supports the analysis of the **Endocardial Fibroelastosis (EFE)** multiomics project, including:

- Single nucleus multiome and snRNA-seq data integration
- Single nucleus multiomics (sn-multiomics and snRNA-seq) analysis
- Spatial transcriptomics analysis (MERFISH)

---

## Repository Structure

---

## Supplemental Materials

| File | Description |
|-----|-----|
| `EFE_vs_CTRL_props.txt` | Summary statistics comparing cell proportions between **EFE (Children)** and **CTRL (Children)** groups |
| `EFE_vs_Fetal_props.txt` | Summary statistics comparing cell proportions between **EFE (Children)** and **CTRL (Fetal)** groups |
| `EFE_vs_Young_props.txt` | Summary statistics comparing cell proportions between **EFE (Children)** and **CTRL (Young)** groups |
| `4_Reviewer4_Q6_FigureE` | Boxplot comparing cell proportions in EFE to age-matched CTRL samples |
| `5_ReviewMaterials.xls` | Detailed gene lists including: (1) MERFISH gene probes, (2) endothelial cell marker genes for Figure G, and (3) fibrosis-related genes for Figure G |
| `6_Methods_Online_030426.docx` | Comprehensive Materials and Methods section for the study |

---


## Analysis Scripts

| Script | Description |
|------|------|
| `1_Sample_Integrations_030126.R` | Integration of single-nucleus multiomics data with public datasets |
| `1_SingleCell_Analysis_030126.R` | Core single-cell analysis pipeline |
| `1_CellChat_GO_030126.R` | Cell–cell communication analysis and Gene Ontology (GO) enrichment analysis |
| `1_MERFISH_analysis_030126.ipynb` | Spatial transcriptomics analysis (MERFISH), including cell type annotation and neighborhood enrichment |

---



## Citation

If you use this repository or associated analysis, please cite the corresponding publication (to be added upon publication).

---

## Contact

For questions regarding the analysis or data, please open a GitHub issue or contact the repository maintainer.

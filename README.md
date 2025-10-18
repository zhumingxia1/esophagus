# Radiomics-Based Prognostic Modeling Workflow

## Overview
This repository provides a complete R-based workflow for radiomics feature preprocessing, survival modeling, and biological interpretation.  
It includes step-by-step scripts for reproducible analyses covering normalization, feature selection, survival modeling, validation, and immune/functional profiling.

---

## 📂 Scripts Overview

| No. | Script | Description |
|:--|:--|:--|
| [1](#1-z-score--icc) | **z-score & ICC.R** | Performs z-score normalization of radiomics features and computes interclass correlation coefficients (ICC) to assess reproducibility. |
| [2](#2-unicox-score-clir) | **uniCox score cli.R** | Runs univariate Cox proportional hazards regression for each variable (risk and clinical) to identify prognostic factors. |
| [3](#3-feature-forestr) | **feature forest.R** | Draws a forest plot displaying hazard ratios (HRs), 95% confidence intervals, and p-values for each variable. |
| [4](#4-km-curvesr) | **KM curves.R** | Generates Kaplan–Meier survival curves for training and testing cohorts using optimal cutoff values. |
| [5](#5-rocr) | **ROC.R** | Produces time-dependent ROC curves (1-, 2-, and 3-year AUCs) and compares AUCs among models and clinical variables. |
| [6](#6-c-indexr) | **c-index.R** | Calculates time-dependent concordance indices (C-index) via bootstrap cross-validation to evaluate model discrimination. |
| [7](#7-nomor) | **Nomo.R** | Builds a Cox-based nomogram and visualizes 1-, 2-, and 3-year survival probabilities. |
| [8](#8-dcar) | **DCA.R** | Conducts Decision Curve Analysis (DCA) to evaluate clinical net benefit across multiple predictive models. |
| [9](#9-tcga-gsea--ssgsea--cibersortr) | **TCGA (GSEA & ssGSEA & cibersort).R** | Performs transcriptomic functional enrichment (GSEA), single-sample GSEA, and immune deconvolution (CIBERSORT). |

---

## ⚙️ Requirements

- **R version:** ≥ 4.1  
- **Operating system:** Windows / macOS / Linux  

### Main R Packages
```
survival, survminer, timeROC, pec, regplot, ggDCA,
limma, clusterProfiler, enrichplot, GSVA, GSEABase, ggplot2
```

Install missing dependencies:
```r
install.packages(c("survival", "survminer", "timeROC", "pec", 
                   "regplot", "ggDCA", "limma", "clusterProfiler", 
                   "enrichplot", "GSVA", "GSEABase", "ggplot2"))
```

---

## 🚀 Usage

1. Place all input data files (e.g., `trainorig.csv`, `score.txt`, `mRNA.txt`) in the working directory.
2. Run scripts sequentially to reproduce the analysis:
   ```bash
   Rscript "z-score & ICC.R"
   Rscript "uniCox score cli.R"
   Rscript "feature forest.R"
   Rscript "KM curves.R"
   Rscript "ROC.R"
   Rscript "c-index.R"
   Rscript "Nomo.R"
   Rscript "DCA.R"
   Rscript "TCGA (GSEA & ssGSEA & cibersort).R"
   ```
3. All outputs (tables and plots) will be automatically saved as `.txt` or `.pdf` files in the same directory.

---

## 📊 Output Examples

| Script | Output | Description |
|:--|:--|:--|
| z-score & ICC.R | `train.csv`, `test.csv`, `icc_inter.csv` | Normalized datasets and ICC values |
| uniCox score cli.R | `uniCox score cli.txt`, `uni scorecli.txt` | Univariate Cox results and significant feature list |
| feature forest.R | `feature forest.pdf` | HR forest plot |
| KM curves.R | `KM curves (train/test).pdf` | Kaplan–Meier survival plots |
| ROC.R | `ROC ici score.pdf`, `cliROC ici 2year.pdf` | Time-dependent ROC analyses |
| c-index.R | `crt time-dependent Cindex.pdf` | Time-dependent C-index curves |
| Nomo.R | `nomogram.pdf` | Cox-based nomogram |
| DCA.R | `DCA crt 2year.pdf`, `DCA ici 2year.pdf` | Decision Curve Analysis |
| TCGA (GSEA & ssGSEA & cibersort).R | `GO.txt`, `GO.pdf`, `Hallmark score.txt`, `CIBERSORT_results.txt` | Functional and immune enrichment results |

---

## 🧠 Notes
- Ensure input tables follow the expected format (headers, row names, and data types).
- Adjust working directories or file paths before running each script if necessary.
- For reproducibility, set random seeds (`set.seed()`) where applicable.
- All scripts are modular — they can be executed independently if inputs are available.

---


## 🔗 Contact
**Corresponding author:** Xin Zhou , M.D., Ph.D.  
Department of Oncology, First Affiliated Hospital of Nanjing Medical University  
Email: zhouxin5523@jsph.org.cn

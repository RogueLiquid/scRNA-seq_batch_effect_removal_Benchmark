# 🫁 Lung Atlas Integration Benchmark (Scenario 1)

This repository contains code and data for benchmarking atlas-level data integration methods on single-cell RNA-seq data using the Lung Atlas dataset. We compare Seurat v5, Harmony, scVI, and fastMNN under a controlled subset setting (Scenario 1: Batches 1–3 only).

---

## 🔁 Reproducibility Guide

### 1. 🔽 Download the Dataset

You can download the h5ad dataset with get_data.py. Data should be automatically put into folder scenario1 and scenario2.

### 2. 🧼 Quality Control (Python)

Filter cells based on gene count and mitochondrial content.

```bash
python "quality control.py"
```

This generates:
- `Lung_atlas_QCfiltered.h5ad`

### 3. ⚙️ Run scVI (Python)

Before evaluating scVI results in R, latent embedding must be generated:

```bash
python "run_scvi.py"
```

This generates:
- `lung_scvi_latent.h5ad`

---

## 🧪 Evaluation & Metric Calculation

Each integration method is evaluated independently using the following R scripts:

| Method      | Script                   | Notes |
|-------------|--------------------------|-------|
| Seurat v5   | `Seurat_evaluation.R`    | Also evaluates raw and calculates metrics |
| Harmony     | `Harmony_evaluation.R`   | UMAP and metric outputs |
| fastMNN     | `fastMNN_evaluation.R`   | Integration and evaluation |
| MNN (failed)| `MNN_evaluation.R`       | Optional (very long runtime) |
| scVI        | `scvi_evaluation.R`      | Uses output from `run_scvi.py` |

All metric results are collected in `benchmark_results.csv`.

---

## 📊 UMAP Visualization

To generate consistent UMAP plots across all methods:

```R
source("plot_UMAP.R")
```

This generates:
- `UMAP_all_celltype.png`
- `UMAP_all_batch.png`

---

## 📈 Evaluation Metrics

Metrics used in this study:
- **iLISI / cLISI** (Integration & Cell Type Local Inverse Simpson's Index)
- **ASW (batch & cell type)** (Average Silhouette Width)
- **kNN Graph Connectivity**
- **Isolated Label Score (ILS)**
- **Adjusted Rand Index (ARI)**
- **PCR R² (Principal Component Regression)**

---

## 🛠 Dependencies

### R
- `zellkonverter`
- `Seurat`
- `SeuratDisk`
- `SingleCellExperiment`
- `harmony`
- `batchelor`
- `aricode`, `cluster`, `Polychrome`
- `colorspace`, `dplyr`, `ggplot2`, `RColorBrewer`
- `kBET`, `lisi`, `mclust`, `FNN`
- `patchwork`, `readr`, `scales`
- `scater`, `scran`, `Matrix`

### Python
- `scanpy`
- `scvi-tools`
- `anndata`, `numpy`, `pandas`

---

## 📜 Citation

Original dataset:  
**Theis Lab (2021). Benchmarking atlas-level data integration in single-cell genomics.**  
🔗 [DOI: 10.6084/m9.figshare.12420968](https://doi.org/10.6084/m9.figshare.12420968)
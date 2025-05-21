import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns                  # nicer histograms

# --------------------------------------------------------------------------
# 1.  Read the full AnnData
# --------------------------------------------------------------------------
adata = sc.read_h5ad("Lung_atlas_public.h5ad")   # adjust path if needed

# The raw UMI counts are in layers["counts"]; copy to .X so qc metrics use it
adata.X = adata.layers["counts"]

# --------------------------------------------------------------------------
# 2.  Add standard QC metrics
# --------------------------------------------------------------------------
adata.var['mt'] = adata.var_names.str.startswith('MT-')  #  human mt-genes
sc.pp.calculate_qc_metrics(
    adata, 
    qc_vars=['mt'],       # computes pct_counts_mt
    percent_top=None,
    inplace=True
)

# `adata.obs` now contains:
#   - total_counts   (nUMI)
#   - n_genes_by_counts  (nGene)
#   - pct_counts_mt

# --------------------------------------------------------------------------
# 3.  Visualise distributions
# --------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

sns.histplot(adata.obs['total_counts'], bins=50, ax=axes[0])
axes[0].set_title('UMI counts per cell')
axes[0].set_xlabel('total_counts')

sns.histplot(adata.obs['n_genes_by_counts'], bins=50, ax=axes[1])
axes[1].set_title('Genes detected per cell')
axes[1].set_xlabel('n_genes_by_counts')

sns.histplot(adata.obs['pct_counts_mt'], bins=50, ax=axes[2])
axes[2].set_title('% mitochondrial counts per cell')
axes[2].set_xlabel('pct_counts_mt')

plt.tight_layout()
plt.show()

# --------------------------------------------------------------------------
# 4.  Choose thresholds (example values – adjust after seeing plots)
# --------------------------------------------------------------------------
min_genes   = 200
max_genes   = 6000
max_mt      = 15          # percent

initial_n   = adata.n_obs

adata = adata[(
    (adata.obs['n_genes_by_counts'] >= min_genes) &
    (adata.obs['n_genes_by_counts'] <= max_genes) &
    (adata.obs['pct_counts_mt']      <= max_mt)
)].copy()

print(f"Kept {adata.n_obs} / {initial_n} cells "
      f"({adata.n_obs/initial_n:.1%}) after QC filtering.")

# --------------------------------------------------------------------------
# 5.  Save the clean subset (optional)
# --------------------------------------------------------------------------
adata.write_h5ad("Lung_atlas_QCfiltered.h5ad")

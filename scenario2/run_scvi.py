import scanpy as sc
import scvi
import os

scvi.settings.seed = 0
adata = sc.read("Lung_atlas_QCfiltered.h5ad")

sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=3000,
    layer="counts",
    batch_key="batch",
    subset=True,
)

scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
model_path = "lung_model"
if os.path.exists(model_path):
    # 已有保存模型，直接加载
    model = scvi.model.SCVI.load(model_path, adata=adata)
    print("📂 Loaded existing scVI model.")
else:
    # 无模型，创建并训练
    model = scvi.model.SCVI(
        adata,
        n_layers=2,
        n_latent=30,
        gene_likelihood="nb"
    )
    model.train()
    model.save(model_path, overwrite=True)
    print("✅ Trained and saved new scVI model.")

SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(adata)

# finally: save everything (counts, obs, latent, UMAP) for R
sc.write("lung_scvi_latent.h5ad", adata, compression="gzip")
print("✅ Saved lung_scvi_latent.h5ad")
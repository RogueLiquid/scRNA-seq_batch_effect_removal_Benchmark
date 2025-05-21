import scanpy as sc
import os
import requests

# === Step 1: Download the file from Figshare ===
url = "https://figshare.com/ndownloader/files/24539942"
output_dir = "scenario2"
output_file = os.path.join(output_dir, "Lung_atlas_public.h5ad")

os.makedirs(output_dir, exist_ok=True)

if not os.path.exists(output_file):
    print("📥 Downloading Lung Atlas H5AD file...")
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(output_file, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    print("✅ Download complete.")
else:
    print("✅ File already exists. Skipping download.")

# === Step 2: Read the downloaded AnnData file ===
print("📖 Reading H5AD...")
adata = sc.read_h5ad(output_file)

# === Step 3: Subset batches 1, 2, and 3 ===
batches_to_keep = ["1", "2", "3"]
adata_small = adata[adata.obs["batch"].isin(batches_to_keep)].copy()

print("Original shape :", adata.shape)
print("Subset  shape  :", adata_small.shape)
print("Batches kept   :", adata_small.obs['batch'].unique().tolist())

# === Step 4: Save the subset to scenario1 ===
os.makedirs("scenario1", exist_ok=True)
subset_output = os.path.join("scenario1", "Lung_atlas_batch123.h5ad")
adata_small.write_h5ad(subset_output, compression="gzip")
print(f"💾 Subset saved to: {subset_output}")

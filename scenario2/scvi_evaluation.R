setwd("D:/college/year_3/CMML/ICA2/scenario2")

library(zellkonverter)
library(Seurat)
library(dplyr)
library(patchwork)
# source your benchmark_metrics() function here -----------------------------

# ---------- read the new .h5ad --------------------------------------------
adata_scvi <- readH5AD("lung_scvi_latent.h5ad")

# Convert to Seurat; keep raw counts
seurat_scvi <- as.Seurat(adata_scvi,
                         counts = "counts",   # raw counts layer
                         data   = NULL)       # we’ll log-normalise fresh if needed
Reductions(seurat_scvi)
embedding_scvi <- Embeddings(seurat_scvi, "X_scVI")

# Rename assay so code downstream can reuse nFeature_RNA, etc.
seurat_scvi <- RenameAssays(seurat_scvi, originalexp = "RNA")
DefaultAssay(seurat_scvi) <- "RNA"

# The latent embedding is imported as a reduction called "X_scVI"
# (zellkonverter converts every obsm slot to a DimReduc)
Embeddings(seurat_scvi, "X_scVI")[1:5, 1:5]  # quick sanity check

# Leiden clusters created in Python arrive as a column:
table(head(seurat_scvi@meta.data$leiden))

# -------------- build the meta + embedding list ---------------------------
meta_scvi <- seurat_scvi@meta.data |>
  mutate(
    cell_id   = colnames(seurat_scvi),
    dataset   = batch,                 # or whatever column marks the batches
    nGene     = nFeature_RNA,          # now it exists because we renamed assay
    cell_type = cell_type,             # already present in the .h5ad
    cluster   = leiden            # use scvi-based clustering
  )

cell_lines_scvi <- list(
  meta_data = meta_scvi,
  scaled_pcs = Embeddings(seurat_scvi, "X_scVI")  # 30-D scVI space
)

benchmark_metrics <- function(cell_lines, label = "Method", output_csv = "benchmark_results.csv") {
  cat("\n====", label, "====\n")
  embedding <- cell_lines$scaled_pcs
  meta <- cell_lines$meta_data[match(rownames(embedding), cell_lines$meta_data$cell_id), ]
  k <- 30
  
  lisi <- compute_lisi(embedding, meta, c("dataset", "cell_type"))
  mean_ilisi <- round(mean(lisi$dataset), 3)
  mean_clisi <- round(mean(lisi$cell_type), 3)
  
  sil_ct <- silhouette(as.numeric(factor(meta$cell_type)), dist(embedding))
  sil_batch <- silhouette(as.numeric(factor(meta$dataset)), dist(embedding))
  asw_ct <- round(mean(sil_ct[, "sil_width"]), 3)
  asw_batch <- round(mean(sil_batch[, "sil_width"]), 3)
  
  knn <- get.knn(embedding, k = k)
  knn_mix <- round(mean(sapply(seq_len(nrow(embedding)), function(i) {
    mean(meta$dataset[knn$nn.index[i, ]] != meta$dataset[i])
  })), 3)
  
  r2_pcr <- round(mean(sapply(seq_len(ncol(embedding)), function(j) {
    summary(lm(embedding[, j] ~ meta$dataset))$r.squared
  })), 3)
  
  isolated_label <- round(mean(sapply(seq_len(nrow(embedding)), function(i) {
    mean(meta$cell_type[knn$nn.index[i, ]] == meta$cell_type[i])
  })), 3)
  
  nmi <- round(NMI(meta$cell_type, meta$cluster, variant = "sqrt"), 3)
  ari <- round(adjustedRandIndex(meta$cell_type, meta$cluster), 3)
  
  # Print results
  cat("🧪 mean iLISI:", mean_ilisi, "\n")
  cat("🧬 mean cLISI:", mean_clisi, "\n")
  cat("📐 ASW (cell type):", asw_ct, "\n")
  cat("🧯 ASW (batch):", asw_batch, "\n")
  cat("🔗 kNN Graph Connectivity:", knn_mix, "\n")
  cat("🧪 Mean R² from PCR:", r2_pcr, "\n")
  cat("🧬 Isolated Label Score:", isolated_label, "\n")
  cat("📊 NMI (cell type vs. cluster):", nmi, "\n")
  cat("📊 ARI (cell type vs. cluster):", ari, "\n")
  
  # Store in data frame
  res <- data.frame(
    Method = label,
    mean_iLISI = mean_ilisi,
    mean_cLISI = mean_clisi,
    ASW_celltype = asw_ct,
    ASW_batch = asw_batch,
    kNN_batch_mix = knn_mix,
    PCR_R2 = r2_pcr,
    Isolated_Label = isolated_label,
    NMI = nmi,
    ARI = ari,
    stringsAsFactors = FALSE
  )
  
  # Write (append if exists)
  if (!file.exists(output_csv)) {
    write.csv(res, output_csv, row.names = FALSE)
  } else {
    write.table(res, output_csv, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}

# -------------- run your metrics ------------------------------------------
benchmark_metrics(cell_lines_scvi, label = "scVI Integrated")

# ------------ 1. Does a UMAP already exist? ---------------------------------
if (! any(grepl("umap", Reductions(seurat_scvi), ignore.case = TRUE))) {
  
  message("No UMAP found – computing UMAP on the scVI latent space …")
  
  # scVI latent was imported as "XscVI" (key auto-renamed by Seurat)
  seurat_scvi <- RunUMAP(
    seurat_scvi,
    reduction      = "X_scVI",
    dims           = 1:30,
    reduction.name = "umap_scvi",   # give it a clear name
    reduction.key  = "umap_scvi_"
  )
} else {
  message("UMAP reduction already present: ",
          paste(Reductions(seurat_scvi), collapse = ", "))
}

# identify the UMAP DimReduc we’ll use
umap_name <- grep("umap", Reductions(seurat_scvi), value = TRUE, ignore.case = TRUE)[1]

# ------------ 2. Save the coordinates for later plotting --------------------
save_umap <- function(obj, file, reduction, extra.cols = c("batch", "cell_type")) {
  key <- Key(obj[[reduction]])
  df  <- FetchData(obj,
                   vars = c(paste0(key, "1"),
                            paste0(key, "2"),
                            extra.cols))
  colnames(df)[1:2] <- c("UMAP_1", "UMAP_2")
  df$method <- basename(tools::file_path_sans_ext(file))
  saveRDS(df, file = file, compress = "xz")
  invisible(df)
}

save_umap(
  seurat_scvi,
  file       = "umap_scvi.rds",
  reduction  = umap_name,
  extra.cols = c(if ("batch" %in% colnames(seurat_scvi@meta.data))
    "batch" else "X_scvi_batch",
    "cell_type")
)

# ------------ 3. Quick visual check ----------------------------------------
p_batch <- DimPlot(seurat_scvi, reduction = umap_name, group.by = if ("batch" %in% colnames(seurat_scvi@meta.data)) "batch" else "X_scvi_batch") +
  ggtitle("scVI – batch") + theme_void()
p_type  <- DimPlot(seurat_scvi, reduction = umap_name, group.by = "cell_type") +
  ggtitle("scVI – cell type") + theme_void()




##############################################################################
## Harmony batch-correction pipeline (Seurat v5)
##  – dataset: Lung_atlas_public.h5ad
##  – integration key: "batch"
##############################################################################
setwd("D:/college/year_3/CMML/ICA2/scenario1")

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)          # install.packages("harmony") if missing
  library(zellkonverter)    # readH5AD
  library(dplyr)
  library(Matrix)
  library(patchwork)
  ## metrics dependencies
  library(FNN);  library(lisi);  library(kBET)
  library(mclust); library(aricode)
  library(SeuratDisk)
})

# ---------- 1. read h5ad and convert to Seurat ------------------------------
adata <- readH5AD("Lung_atlas_QCfiltered.h5ad")
seurat <- as.Seurat(adata, counts = "counts", data = NULL)
DefaultAssay(seurat) <- "originalexp"          # use imported assay
rm(adata); gc()

# ---------- 2. minimal preprocessing ----------------------------------------
seurat <- FindVariableFeatures(seurat, nfeatures = 3000)
seurat <- NormalizeData(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, npcs = 30)

# choose batch column
batch_col <- if ("batch" %in% colnames(seurat@meta.data)) "batch" else "dataset"

# ---------- 3. Harmony integration ------------------------------------------
seurat <- RunHarmony(
  object          = seurat,
  group.by.vars   = batch_col,
  reduction.save  = "harmony",    # name of the new DimReduc
  plot_convergence = FALSE
)

# ---------- 4. UMAP / clustering on Harmony space ---------------------------
seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:30)
seurat <- FindClusters(seurat, resolution = 0.5)

# ---------- 5. build cell_lines list for benchmark_metrics ------------------
meta_harm <- seurat@meta.data |>
  mutate(
    cell_id   = colnames(seurat),
    dataset   = .data[[batch_col]],
    nGene     = nFeature_originalexp,
    cell_type = cell_type,
    cluster   = seurat_clusters
  )

cell_lines_harmony <- list(
  meta_data  = meta_harm,
  scaled_pcs = Embeddings(seurat, "harmony")
)

# ---------- 6. run metrics ---------------------------------------------------
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

benchmark_metrics(cell_lines_harmony, label = "Harmony Integrated")

# ---------- 7. save tiny UMAP table -----------------------------------------
save_umap <- function(obj, file, reduction, extra.cols = c("batch","cell_type")) {
  key <- Key(obj[[reduction]])
  df  <- FetchData(obj, vars = c(paste0(key, "1"),
                                 paste0(key, "2"),
                                 extra.cols))
  colnames(df)[1:2] <- c("UMAP_1", "UMAP_2")
  df$method <- basename(tools::file_path_sans_ext(file))
  saveRDS(df, file = file, compress = "xz")
}

save_umap(
  seurat,
  file       = "umap_harmony.rds",
  reduction  = "umap_harmony",
  extra.cols = c(batch_col, "cell_type")
)

# ---------- 8. slim & save Seurat object (optional) --------------------------
harm_slim <- DietSeurat(
  object  = seurat,
  assays  = c("originalexp", "RNA"),       # RNA created by NormalizeData()
  layers  = list(RNA = "data",             # log-norm
                 originalexp = "counts"),  # raw counts
  dimreducs = c("harmony", "umap_harmony")
)

SaveH5Seurat(harm_slim, "processed/lung_harmony_slim.h5seurat", overwrite = TRUE)

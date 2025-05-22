setwd("D:/college/year_3/CMML/ICA2/scenario2")

suppressPackageStartupMessages({
  library(zellkonverter)         # readH5AD
  library(SingleCellExperiment)
  library(batchelor)             # fastMNN
  library(scater)                # QC & logNorm
  library(scran)                 # HVG
  library(Seurat)                # RunUMAP
  library(SeuratDisk)           # Conversion
  library(dplyr)
  library(FNN)
  library(lisi)
  library(kBET)
  library(mclust)
  library(aricode)
  library(patchwork)
  library(cluster)
  ##  metrics deps: FNN, lisi, kBET, mclust, aricode already loaded earlier
})

# ---------- 1. read AnnData --------------------------------------------------
sce <- readH5AD("Lung_atlas_QCfiltered.h5ad")
batch_col <- if ("batch" %in% colnames(colData(sce))) "batch" else "dataset"

# ---------- 2. select 3 000 highly-variable genes ---------------------------
dec     <- modelGeneVar(sce, assay.type = "counts")
top_hvgs<- getTopHVGs(dec, n = 3000)

# ---------- 3. fastMNN integration ------------------------------------------
set.seed(0)

sce <- logNormCounts(sce)

sce_mnn <- fastMNN(
  sce,
  batch      = sce[[batch_col]],
  subset.row = top_hvgs,
  d          = 50,     # 50 corrected PCs by default
  k          = 20,
  assay.type = "logcounts"
)

# ---------- 4. UMAP on corrected PCs using Seurat --------------------------
assay(sce_mnn, "reconstructed") <- as(assay(sce_mnn, "reconstructed"), "dgCMatrix")
seurat_mnn <- as.Seurat(
  sce_mnn,
  counts = NULL,
  data   = "reconstructed"   # batch-corrected expression
)

cell_meta <- cell_lines_mnn$meta_data
seurat_mnn@meta.data$cell_type <- cell_meta$cell_type[match(rownames(seurat_mnn@meta.data), cell_meta$cell_id)]
seurat_mnn@meta.data$batch     <- cell_meta$dataset[match(rownames(seurat_mnn@meta.data), cell_meta$cell_id)]

# Add corrected PCs
seurat_mnn[["pca_corrected"]] <- CreateDimReducObject(
  embeddings = reducedDim(sce_mnn, "corrected"),
  key        = "mnnPC_",
  assay      = "RNA"
)

# Use Seurat's RunUMAP instead of scater's runUMAP
seurat_mnn <- RunUMAP(
  seurat_mnn,
  reduction      = "pca_corrected",
  dims           = 1:30,
  reduction.name = "umap_mnn",
  reduction.key  = "umap_mnn_"
)

# ---------- 5. Louvain clustering on corrected graph --------------------------
set.seed(0)
g <- buildSNNGraph(sce_mnn, use.dimred = "corrected", k = 20)
colLabels(sce_mnn, "cluster_fastmnn") <- igraph::cluster_louvain(g)$membership
colData(sce_mnn)$cluster_fastmnn <- colLabels(sce_mnn, "cluster_fastmnn")

# ---------- 6. Metadata prep --------------------------------------------------
if (! "nGene" %in% colnames(colData(sce_mnn))) {
  if ("counts" %in% assayNames(sce)) {
    colData(sce_mnn)$nGene <- Matrix::colSums(assay(sce, "counts") > 0)
  } else if ("counts" %in% assayNames(sce_mnn)) {
    colData(sce_mnn)$nGene <- Matrix::colSums(assay(sce_mnn, "counts") > 0)
  } else {
    colData(sce_mnn)$nGene <- NA_real_
    warning("nGene not found; set to NA")
  }
}

if (! "cell_type" %in% colnames(colData(sce_mnn))) {
  colData(sce_mnn)$cell_type <-
    colData(sce)$cell_type[ match(colnames(sce_mnn), colnames(sce)) ]
}

meta_mnn <- as.data.frame(colData(sce_mnn)) |>
  mutate(
    cell_id   = colnames(sce_mnn),
    dataset   = .data[[batch_col]],
    nGene     = nGene,                 # already in colData
    cell_type = cell_type,
    cluster   = cluster_fastmnn
  )

cell_lines_mnn <- list(
  meta_data  = meta_mnn,
  scaled_pcs = reducedDim(sce_mnn, "corrected")
)

# ---------- 7. save UMAP result -----------------------------------------
umap_coords <- Embeddings(seurat_mnn, reduction = "umap_mnn")

umap_df <- as.data.frame(umap_coords)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")

# Add metadata
umap_df$batch      <- seurat_mnn@meta.data$batch
umap_df$cell_type  <- seurat_mnn@meta.data$cell_type
umap_df$method     <- "fastMNN"

saveRDS(umap_df, file = "umap_fastmnn.rds", compress = "xz")

# ---------- 8. Run benchmarking metrics -----------------------------------
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

benchmark_metrics(cell_lines_mnn, label = "fastMNN Integrated")

##############################################################################
##  Classic MNN (mnnCorrect) integration pipeline
##############################################################################
setwd("D:/college/year_3/CMML/ICA2/scenario2")

suppressPackageStartupMessages({
  library(zellkonverter)     # readH5AD
  library(SingleCellExperiment)
  library(batchelor)         # mnnCorrect
  library(scater)            # runPCA, runUMAP
  library(scran)             # modelGeneVar
  library(Matrix)
  ## metrics deps
  library(dplyr); library(FNN); library(lisi)
  library(kBET); library(mclust); library(aricode)
})

# ---------- 1. load AnnData --------------------------------------------------
sce <- readH5AD("Lung_atlas_QCfiltered.h5ad")

batch_col <- if ("batch" %in% colnames(colData(sce))) "batch" else "dataset"

# ---------- 2. HVG selection & log-normalisation ----------------------------
sce <- logNormCounts(sce)                                   # adds 'logcounts'
dec  <- modelGeneVar(sce, assay.type = "logcounts")
hvg  <- getTopHVGs(dec, n = 3000)

# ---------- 3. classic MNN correction ---------------------------------------
set.seed(0)
sce_mnn <- mnnCorrect(
  sce,
  batch        = sce[[batch_col]],
  subset.row   = hvg,
  k            = 20,
  assay.type   = "logcounts"
)

# mnnCorrect returns a SingleCellExperiment with assay "corrected"
# ---------- 4. PCA & UMAP on corrected expression ---------------------------
sce_mnn <- runPCA(sce_mnn, ncomponents = 50,
                  exprs_values = "corrected", name = "mnn_pca")
sce_mnn <- runUMAP(sce_mnn, dimred = "mnn_pca", n_neighbors = 30)

# ---------- 5. graph-based clustering (Louvain) -----------------------------
library(igraph); library(scran)
g <- buildSNNGraph(sce_mnn, use.dimred = "mnn_pca", k = 20)
colLabels(sce_mnn, "cluster_mnn") <- igraph::cluster_louvain(g)$membership
colData(sce_mnn)$cluster_mnn <- colLabels(sce_mnn, "cluster_mnn")

# ---------- 6. ensure meta columns exist ------------------------------------
if (! "nGene" %in% colnames(colData(sce_mnn))) {
  colData(sce_mnn)$nGene <- Matrix::colSums(assay(sce, "counts") > 0)
}
if (! "cell_type" %in% colnames(colData(sce_mnn))) {
  colData(sce_mnn)$cell_type <- colData(sce)$cell_type[
    match(colnames(sce_mnn), colnames(sce))]
}

# ---------- 7. build cell_lines & run metrics -------------------------------
meta_mnn <- as.data.frame(colData(sce_mnn)) |>
  mutate(cell_id = colnames(sce_mnn),
         dataset = .data[[batch_col]],
         cluster = cluster_mnn)

cell_lines_mnn <- list(
  meta_data  = meta_mnn,
  scaled_pcs = reducedDim(sce_mnn, "mnn_pca")
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

benchmark_metrics(cell_lines_mnn, label = "MNN Integrated")

# ---------- 8. save tiny UMAP table -----------------------------------------
save_umap_sce <- function(sce, file,
                          umap_name = "UMAP",
                          extra.cols = c("dataset","cell_type")) {
  umap <- reducedDim(sce, umap_name)[, 1:2]
  df   <- as.data.frame(umap)
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df   <- cbind(df, as.data.frame(colData(sce)[, extra.cols, drop = FALSE]))
  df$method <- basename(tools::file_path_sans_ext(file))
  saveRDS(df, file = file, compress = "xz")
}

save_umap_sce(
  sce_mnn,
  file       = "umap_mnn.rds",
  extra.cols = c(batch_col, "cell_type")
)

# ---------- 9. (optional) convert to Seurat for DimPlot ---------------------
# Comment out if you don’t need Seurat visualisation
if (requireNamespace("Seurat", quietly = TRUE)) {
  library(Seurat)
  
  # Convert corrected assay to sparse matrix for Seurat
  assay(sce_mnn, "corrected") <-
    as(assay(sce_mnn, "corrected"), "dgCMatrix")
  
  seu_mnn <- as.Seurat(sce_mnn, counts = NULL, data = "corrected")
  seu_mnn[["mnn_pca"]] <- CreateDimReducObject(
    embeddings = reducedDim(sce_mnn, "mnn_pca"),
    key        = "mnnPC_",
    assay      = "RNA"
  )
  seu_mnn[["umap_mnn"]] <- CreateDimReducObject(
    embeddings = reducedDim(sce_mnn, "UMAP"),
    key        = "umapMNN_",
    assay      = "RNA"
  )
  
  # quick view
  DimPlot(seu_mnn, reduction = "umap_mnn", group.by = batch_col) +
    ggtitle("MNN – batch")
}

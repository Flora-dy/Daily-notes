#!/usr/bin/env Rscript
# ============================================================
# Exos QC v2 — 严格按公司质控报告参数重跑
# 参数来源: X101SC251213857-Z03-J001-B1-1.html (诺禾致源)
#
# 样本特异 high_nGene:
#   Exos1=8000, Exos2=8500, Exos3=8000
#   Exos_A1=7000, Exos_A2=8000, Exos_A3=7500
# 统一: low_nGene=200, pct.mito<5%, pct.HB<5%
# ============================================================

.libPaths(c("/Users/zzp/Code/tianhairui/.rlib", .libPaths()))
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(cowplot)
})

BASE_MATRIX <- paste0(
  "/Volumes/Untitled/$RECYCLE.BIN/",
  "S-1-5-21-1373272656-3585098094-2095149322-500/",
  "$R42GI11.20260204/",
  "Result-X101SC251213857-Z03-J001-B1-1/",
  "2.Summary"
)

EXOS_OUT <- "/Users/zzp/Code/tianhairui/outputs0409/seurat5_exos_v2"
OUT_DIR   <- file.path(EXOS_OUT, "00_qc")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(EXOS_OUT, "rds"), showWarnings=FALSE)

SEED <- 2024
HB_PATTERN <- "^Hba|^Hbb"

# 样本特异 QC 参数（公司报告）
sample_info <- data.frame(
  sample_id  = c("Exos1","Exos2","Exos3","Exos_A1","Exos_A2","Exos_A3"),
  group      = c("Rescue","Rescue","Rescue","Aging","Aging","Aging"),
  low_nGene  = 200,
  high_nGene = c(8000, 8500, 8000, 7000, 8000, 7500),
  max_mito   = 5,
  max_hb     = 5,
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_info$sample_id

N_HVG   <- 3000
N_PCS   <- 30
N_NEIGH <- 20
RES_LIST <- c(0.3, 0.5, 0.8, 1.0)

# ══════════════════════════════════════════════════════════════
# STEP 1: 读取 + QC
# ══════════════════════════════════════════════════════════════
cat("===== STEP 1: Loading & QC (company thresholds) =====\n\n")

qc_list  <- list()
obj_list <- list()

for (sid in sample_info$sample_id) {
  cat(sprintf("--- %s ---\n", sid))
  mat_path <- file.path(BASE_MATRIX, sid, "filter_matrix")
  mat <- Read10X(data.dir=mat_path, gene.column=1)
  if (is.list(mat)) mat <- mat[["Gene Expression"]]

  obj <- CreateSeuratObject(counts=mat, project=sid,
                            min.cells=3, min.features=100)
  obj$sample_id <- sid
  obj$group     <- sample_info[sid, "group"]

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern="^mt-")
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern=HB_PATTERN)

  n_raw <- ncol(obj)
  hi <- sample_info[sid, "high_nGene"]
  lo <- sample_info[sid, "low_nGene"]

  obj_f <- subset(obj,
    subset = nFeature_RNA >= lo &
             nFeature_RNA <= hi &
             percent.mt   <  5  &
             percent.hb   <  5
  )
  n_qc <- ncol(obj_f)

  qc_list[[sid]] <- data.frame(
    sample        = sid,
    group         = sample_info[sid,"group"],
    cells_raw     = n_raw,
    cells_qc      = n_qc,
    pct_removed   = round((1-n_qc/n_raw)*100, 1),
    high_nGene    = hi,
    median_nGene  = round(median(obj_f$nFeature_RNA)),
    median_nCount = round(median(obj_f$nCount_RNA)),
    median_mt     = round(median(obj_f$percent.mt), 3),
    median_hb     = round(median(obj_f$percent.hb), 3)
  )
  cat(sprintf("  Raw: %d → QC: %d (removed %.1f%%)\n",
              n_raw, n_qc, (1-n_qc/n_raw)*100))
  cat(sprintf("  Median nGene=%d, mito=%.3f%%, HB=%.3f%%\n",
              round(median(obj_f$nFeature_RNA)),
              median(obj_f$percent.mt),
              median(obj_f$percent.hb)))

  obj_list[[sid]] <- obj_f
  rm(obj, mat); gc()
}

qc_df <- do.call(rbind, qc_list)
rownames(qc_df) <- NULL
write.csv(qc_df, file.path(OUT_DIR,"qc_summary.csv"), row.names=FALSE)
cat("\n===== QC Summary =====\n")
print(qc_df)

# QC 散点图
scatter_list <- lapply(names(obj_list), function(sid) {
  o <- obj_list[[sid]]
  p1 <- ggplot(o@meta.data, aes(x=nFeature_RNA, y=percent.mt)) +
    geom_point(size=0.3, alpha=0.4, color="#2980B9") +
    geom_hline(yintercept=5, color="red", linetype="dashed") +
    geom_vline(xintercept=c(200, sample_info[sid,"high_nGene"]),
               color="red", linetype="dashed") +
    labs(title=paste0(sid," — mito"), x="nGene", y="% mito") +
    theme_classic(base_size=9) +
    theme(plot.title=element_text(hjust=0.5),
          plot.background=element_rect(fill="white"))
  p2 <- ggplot(o@meta.data, aes(x=nFeature_RNA, y=percent.hb)) +
    geom_point(size=0.3, alpha=0.4, color="#E74C3C") +
    geom_hline(yintercept=5, color="red", linetype="dashed") +
    geom_vline(xintercept=c(200, sample_info[sid,"high_nGene"]),
               color="red", linetype="dashed") +
    labs(title=paste0(sid," — HB"), x="nGene", y="% HB") +
    theme_classic(base_size=9) +
    theme(plot.title=element_text(hjust=0.5),
          plot.background=element_rect(fill="white"))
  p1 | p2
})
p_scatter <- wrap_plots(scatter_list, ncol=2)
ggsave(file.path(OUT_DIR,"qc_scatter.png"), p_scatter,
       width=16, height=3.5*length(obj_list), dpi=180)
cat("Saved: qc_scatter.png\n")

# QC 小提琴图
vln_list <- lapply(names(obj_list), function(sid) {
  VlnPlot(obj_list[[sid]],
          features=c("nFeature_RNA","nCount_RNA","percent.mt","percent.hb"),
          pt.size=0, ncol=4) +
    plot_annotation(title=sid,
                    theme=theme(plot.title=element_text(hjust=0.5,face="bold")))
})
p_vln <- wrap_plots(vln_list, ncol=1)
ggsave(file.path(OUT_DIR,"qc_violin.png"), p_vln,
       width=18, height=4.5*length(obj_list), dpi=150)
cat("Saved: qc_violin.png\n")

# ══════════════════════════════════════════════════════════════
# STEP 2: 合并 + Harmony + 聚类
# ══════════════════════════════════════════════════════════════
cat("\n===== STEP 2: Merge + Harmony + Clustering =====\n")
obj <- merge(obj_list[[1]], y=obj_list[-1],
             add.cell.ids=names(obj_list),
             project="Exos_Testis_v2")
obj <- JoinLayers(obj)
rm(obj_list); gc()

cat(sprintf("Merged: %d cells, %d genes\n", ncol(obj), nrow(obj)))
cat(sprintf("Groups: Aging=%d, Rescue=%d\n",
            sum(obj$group=="Aging"), sum(obj$group=="Rescue")))

obj <- NormalizeData(obj, verbose=FALSE)
obj <- FindVariableFeatures(obj, nfeatures=N_HVG, verbose=FALSE)
obj <- ScaleData(obj, verbose=FALSE)

set.seed(SEED)
obj <- RunPCA(obj, npcs=N_PCS, seed.use=SEED, verbose=FALSE)
obj <- RunHarmony(obj, group.by.vars="sample_id",
                  reduction.use="pca", dims.use=1:N_PCS, verbose=FALSE)
obj <- RunUMAP(obj, reduction="harmony", dims=1:N_PCS, seed.use=SEED, verbose=FALSE)
obj <- FindNeighbors(obj, reduction="harmony", dims=1:N_PCS,
                     k.param=N_NEIGH, verbose=FALSE)
for (res in RES_LIST) {
  obj <- FindClusters(obj, resolution=res, random.seed=SEED, verbose=FALSE)
  cat(sprintf("  res=%.1f → %d clusters\n", res,
              length(unique(obj$seurat_clusters))))
}
obj <- FindClusters(obj, resolution=0.8, random.seed=SEED, verbose=FALSE)
n_cl <- length(unique(obj$seurat_clusters))
cat(sprintf("Using res=0.8: %d clusters\n", n_cl))
print(table(obj$seurat_clusters))

# ══════════════════════════════════════════════════════════════
# STEP 3: 可视化
# ══════════════════════════════════════════════════════════════
cat("\n===== STEP 3: Visualization =====\n")

p_cl <- DimPlot(obj, reduction="umap", label=TRUE, label.size=4,
                pt.size=0.2) +
  labs(title=sprintf("Exos v2 — %d clusters (res=0.8)", n_cl)) +
  theme_cowplot(font_size=12) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white"), legend.position="none")

p_grp <- DimPlot(obj, reduction="umap", group.by="group",
                 cols=c("Rescue"="#4DBBD5","Aging"="#E64B35"), pt.size=0.2) +
  labs(title="Aging vs Rescue") +
  theme_cowplot(font_size=12) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white"))

p_smp <- DimPlot(obj, reduction="umap", group.by="sample_id", pt.size=0.12) +
  labs(title="By Sample") +
  theme_cowplot(font_size=10) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white"))

p_mt <- FeaturePlot(obj, features="percent.mt", reduction="umap", pt.size=0.12) +
  scale_color_viridis_c(option="B") + labs(title="percent.mt") +
  theme_cowplot(font_size=9) +
  theme(plot.background=element_rect(fill="white"), axis.title=element_blank())
p_hb <- FeaturePlot(obj, features="percent.hb", reduction="umap", pt.size=0.12) +
  scale_color_viridis_c(option="A") + labs(title="percent.hb") +
  theme_cowplot(font_size=9) +
  theme(plot.background=element_rect(fill="white"), axis.title=element_blank())
p_nf <- FeaturePlot(obj, features="nFeature_RNA", reduction="umap", pt.size=0.12) +
  scale_color_viridis_c(option="C") + labs(title="nFeature_RNA") +
  theme_cowplot(font_size=9) +
  theme(plot.background=element_rect(fill="white"), axis.title=element_blank())

key_m <- c("Ddx4","Sycp3","Zbtb16","Kit","Acrv1","Tnp1","Sox9","Star","Pecam1","Cd14")
key_m <- key_m[key_m %in% rownames(obj)]
feat_plots <- lapply(key_m, function(g) {
  FeaturePlot(obj, features=g, reduction="umap", pt.size=0.12, order=TRUE) +
    scale_color_viridis_c(option="D") +
    theme_cowplot(font_size=8) +
    theme(plot.background=element_rect(fill="white"),
          plot.title=element_text(hjust=0.5,face="bold.italic"),
          axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())
})
p_feat <- wrap_plots(feat_plots, nrow=2)

p_main <- (p_cl | p_grp | p_smp) / (p_nf | p_mt | p_hb) +
  plot_layout(heights=c(1.2,1)) +
  plot_annotation(
    title="Exos Testis scRNA-seq v2 (company QC: mito<5%, HB<5%)",
    theme=theme(plot.title=element_text(size=13,face="bold",hjust=0.5),
                plot.background=element_rect(fill="white"))
  )

ggsave(file.path(OUT_DIR,"overview.pdf"),  p_main, width=18, height=12)
ggsave(file.path(OUT_DIR,"overview.png"),  p_main, width=18, height=12, dpi=200)
ggsave(file.path(OUT_DIR,"feature_markers.png"), p_feat, width=20, height=8, dpi=200)
cat("  Saved: overview.png, feature_markers.png\n")

# ══════════════════════════════════════════════════════════════
# STEP 4: 保存
# ══════════════════════════════════════════════════════════════
cat("\n===== STEP 4: Saving =====\n")
saveRDS(obj, file.path(EXOS_OUT,"rds","seurat_exos_clustered_v2.rds"))
write.csv(obj@meta.data, file.path(OUT_DIR,"exos_metadata_v2.csv"))
cat("  Saved: seurat_exos_clustered_v2.rds\n")

cat("\n===== DONE =====\n")
cat(sprintf("Total cells after QC: %d  (Aging=%d, Rescue=%d)\n",
            ncol(obj), sum(obj$group=="Aging"), sum(obj$group=="Rescue")))

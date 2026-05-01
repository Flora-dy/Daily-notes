#!/usr/bin/env Rscript
# ============================================================
# EPS QC v2 — 严格按公司质控报告参数重跑
# 参数来源: X101SC251213857-Z03-J001-B1-1.html (诺禾致源)
#
# 样本特异 high_nGene:
#   EPS1=7500, EPS2=7500, EPS3=7500
#   EPS_A1=8000, EPS_A2=6500, EPS_A3=7500
# 统一:
#   low_nGene=200, pct.mito<5%, pct.HB<5%
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

EPS_OUT  <- "/Users/zzp/Code/tianhairui/outputs0409/seurat5_eps"
OUT_DIR  <- file.path(EPS_OUT, "00_qc_v2")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(EPS_OUT, "rds"), showWarnings=FALSE)

SEED <- 2024

# ── 样本特异 QC 参数（来自公司报告）──────────────────────────
sample_info <- data.frame(
  sample_id  = c("EPS1","EPS2","EPS3","EPS_A1","EPS_A2","EPS_A3"),
  group      = c("Rescue","Rescue","Rescue","Aging","Aging","Aging"),
  low_nGene  = 200,          # 统一
  high_nGene = c(7500, 7500, 7500, 8000, 6500, 7500),  # 样本特异
  max_mito   = 5,            # 统一 5%
  max_hb     = 5,            # 统一 5%（血红蛋白）
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_info$sample_id

# 小鼠血红蛋白基因 pattern（用于 PercentageFeatureSet）
HB_PATTERN <- "^Hba|^Hbb"

# 整合 / 聚类参数
N_HVG   <- 3000
N_PCS   <- 30
N_NEIGH <- 20
RES_LIST <- c(0.3, 0.5, 0.8, 1.0)

# ══════════════════════════════════════════════════════════════
# STEP 1: 读取 + QC（样本特异阈值）
# ══════════════════════════════════════════════════════════════
cat("===== STEP 1: Loading & QC (company thresholds) =====\n\n")

qc_list <- list()
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

  # 线粒体基因（小鼠：mt-）
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern="^mt-")
  # 血红蛋白基因
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern=HB_PATTERN)

  n_raw <- ncol(obj)

  # 公司参数过滤
  hi <- sample_info[sid, "high_nGene"]
  lo <- sample_info[sid, "low_nGene"]
  mt <- sample_info[sid, "max_mito"]
  hb <- sample_info[sid, "max_hb"]

  obj_f <- subset(obj,
    subset = nFeature_RNA >= lo &
             nFeature_RNA <= hi &
             percent.mt   <  mt &
             percent.hb   <  hb
  )
  n_qc <- ncol(obj_f)

  qc_list[[sid]] <- data.frame(
    sample         = sid,
    group          = sample_info[sid,"group"],
    cells_raw      = n_raw,
    cells_qc       = n_qc,
    pct_removed    = round((1 - n_qc/n_raw)*100, 1),
    low_nGene      = lo,
    high_nGene     = hi,
    max_mito_pct   = mt,
    max_hb_pct     = hb,
    median_nGene   = round(median(obj_f$nFeature_RNA)),
    median_nCount  = round(median(obj_f$nCount_RNA)),
    median_mt_pct  = round(median(obj_f$percent.mt), 3),
    median_hb_pct  = round(median(obj_f$percent.hb), 3)
  )
  cat(sprintf("  Raw: %d → QC: %d (removed %.1f%%)\n",
              n_raw, n_qc, (1-n_qc/n_raw)*100))
  cat(sprintf("  nGene: [%d, %d], mito<%.0f%%, HB<%.0f%%\n",
              lo, hi, mt, hb))
  cat(sprintf("  Median nGene=%d, median mito=%.3f%%, median HB=%.3f%%\n",
              round(median(obj_f$nFeature_RNA)),
              median(obj_f$percent.mt),
              median(obj_f$percent.hb)))

  obj_list[[sid]] <- obj_f
  rm(obj, mat); gc()
}

qc_df <- do.call(rbind, qc_list)
rownames(qc_df) <- NULL
write.csv(qc_df, file.path(OUT_DIR, "qc_summary_v2.csv"), row.names=FALSE)
cat("\n===== QC Summary =====\n")
print(qc_df[, c("sample","group","cells_raw","cells_qc","pct_removed",
                "high_nGene","median_nGene","median_mt_pct","median_hb_pct")])

# QC 小提琴图
cat("\n=== QC violin plots ===\n")
vln_list <- lapply(names(obj_list), function(sid) {
  o <- obj_list[[sid]]
  VlnPlot(o, features=c("nFeature_RNA","nCount_RNA","percent.mt","percent.hb"),
          pt.size=0, ncol=4) +
    plot_annotation(title=sid,
                    theme=theme(plot.title=element_text(hjust=0.5,face="bold")))
})
p_vln <- wrap_plots(vln_list, ncol=1)
ggsave(file.path(OUT_DIR, "qc_violin_v2.pdf"),  p_vln, width=18, height=4.5*length(obj_list))
ggsave(file.path(OUT_DIR, "qc_violin_v2.png"),  p_vln, width=18, height=4.5*length(obj_list), dpi=150)
cat("  Saved: qc_violin_v2.png\n")

# QC 散点图（nGene vs pct.mt）
scatter_list <- lapply(names(obj_list), function(sid) {
  o <- obj_list[[sid]]
  ggplot(o@meta.data, aes(x=nFeature_RNA, y=percent.mt)) +
    geom_point(size=0.3, alpha=0.4, color="#2980B9") +
    geom_hline(yintercept=5, color="red", linetype="dashed", linewidth=0.7) +
    geom_vline(xintercept=c(200, sample_info[sid,"high_nGene"]),
               color="red", linetype="dashed", linewidth=0.7) +
    labs(title=sid, x="nGene", y="% mito") +
    theme_classic(base_size=10) +
    theme(plot.title=element_text(hjust=0.5,face="bold"),
          plot.background=element_rect(fill="white"))
})
p_scatter <- wrap_plots(scatter_list, ncol=3)
ggsave(file.path(OUT_DIR, "qc_scatter_v2.png"), p_scatter, width=15, height=9, dpi=200)
cat("  Saved: qc_scatter_v2.png\n")

# ══════════════════════════════════════════════════════════════
# STEP 2: 合并 + Harmony + 聚类
# ══════════════════════════════════════════════════════════════
cat("\n===== STEP 2: Merge + Harmony + Clustering =====\n")
obj <- merge(obj_list[[1]], y=obj_list[-1],
             add.cell.ids=names(obj_list),
             project="EPS_Testis_v2")
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
cat("  PCA done\n")

obj <- RunHarmony(obj, group.by.vars="sample_id",
                  reduction.use="pca", dims.use=1:N_PCS, verbose=FALSE)
cat("  Harmony done\n")

obj <- RunUMAP(obj, reduction="harmony", dims=1:N_PCS, seed.use=SEED, verbose=FALSE)
cat("  UMAP done\n")

obj <- FindNeighbors(obj, reduction="harmony", dims=1:N_PCS,
                     k.param=N_NEIGH, verbose=FALSE)
for (res in RES_LIST) {
  obj <- FindClusters(obj, resolution=res, random.seed=SEED, verbose=FALSE)
  cat(sprintf("  res=%.1f → %d clusters\n", res,
              length(unique(obj$seurat_clusters))))
}
obj <- FindClusters(obj, resolution=0.8, random.seed=SEED, verbose=FALSE)
n_cl <- length(unique(obj$seurat_clusters))
cat(sprintf("\nUsing res=0.8: %d clusters\n", n_cl))
print(table(obj$seurat_clusters))

# ══════════════════════════════════════════════════════════════
# STEP 3: 可视化
# ══════════════════════════════════════════════════════════════
cat("\n===== STEP 3: Visualization =====\n")

p_cl <- DimPlot(obj, reduction="umap", label=TRUE, label.size=4,
                pt.size=0.2, repel=FALSE) +
  labs(title=sprintf("EPS v2 — %d clusters (res=0.8)", n_cl)) +
  theme_cowplot(font_size=12) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white"),
        legend.position="none")

p_grp <- DimPlot(obj, reduction="umap", group.by="group",
                 cols=c("Rescue"="#4DBBD5","Aging"="#E64B35"),
                 pt.size=0.2) +
  labs(title="Aging vs Rescue") +
  theme_cowplot(font_size=12) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white"))

p_smp <- DimPlot(obj, reduction="umap", group.by="sample_id",
                 pt.size=0.12) +
  labs(title="By Sample") +
  theme_cowplot(font_size=10) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white"))

# QC 指标在 UMAP 上
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

# Key markers
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
  plot_layout(heights=c(1.2, 1)) +
  plot_annotation(
    title="EPS Testis scRNA-seq v2 (company QC thresholds)",
    theme=theme(plot.title=element_text(size=13,face="bold",hjust=0.5),
                plot.background=element_rect(fill="white"))
  )

ggsave(file.path(OUT_DIR,"overview_v2.pdf"),  p_main, width=18, height=12)
ggsave(file.path(OUT_DIR,"overview_v2.png"),  p_main, width=18, height=12, dpi=200)
cat("  Saved: overview_v2.png\n")

ggsave(file.path(OUT_DIR,"feature_markers_v2.png"), p_feat, width=20, height=8, dpi=200)
cat("  Saved: feature_markers_v2.png\n")

# ══════════════════════════════════════════════════════════════
# STEP 4: 保存
# ══════════════════════════════════════════════════════════════
cat("\n===== STEP 4: Saving =====\n")
saveRDS(obj, file.path(EPS_OUT, "rds/seurat_eps_clustered_v2.rds"))
write.csv(obj@meta.data, file.path(OUT_DIR, "eps_metadata_v2.csv"))
cat("  Saved: seurat_eps_clustered_v2.rds\n")

cat("\n===== DONE =====\n")
cat("Output:", OUT_DIR, "\n")
cat("  qc_summary_v2.csv\n")
cat("  qc_violin_v2.png\n")
cat("  qc_scatter_v2.png\n")
cat("  overview_v2.png\n")
cat("  rds/seurat_eps_clustered_v2.rds\n")

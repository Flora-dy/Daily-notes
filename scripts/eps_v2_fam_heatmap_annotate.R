#!/usr/bin/env Rscript
# ============================================================
# EPS v2 — FindAllMarkers + 聚类热图 + 细胞类型注释
# 输入: seurat_eps_clustered_v2.rds
#       (来自 eps_seurat5_qc_v2.R，公司 QC 参数)
# 输出目录: seurat5_eps/
#   01_heatmap/  → FAM 热图 + DotPlot
#   02_annotation/ → 宽泛注释 UMAP + 生殖细胞亚群分期
# ============================================================

.libPaths(c("/Users/zzp/Code/tianhairui/.rlib", .libPaths()))
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(cowplot)
  library(viridis)
})

EPS_OUT <- "/Users/zzp/Code/tianhairui/outputs0409/seurat5_eps"
H_DIR   <- file.path(EPS_OUT, "01_heatmap")   # 复用已有路径
A_DIR   <- file.path(EPS_OUT, "02_annotation") # 复用已有路径
V2_DIR  <- file.path(EPS_OUT, "v2")            # v2 专属输出
dir.create(H_DIR,   recursive=TRUE, showWarnings=FALSE)
dir.create(A_DIR,   recursive=TRUE, showWarnings=FALSE)
dir.create(V2_DIR,  recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(EPS_OUT,"rds"), showWarnings=FALSE)

# ══════════════════════════════════════════════════════════════
# STEP 1: 加载 v2 聚类结果
# ══════════════════════════════════════════════════════════════
cat("=== Loading EPS v2 clustered object ===\n")
rds_v2 <- file.path(EPS_OUT, "rds/seurat_eps_clustered_v2.rds")
obj <- readRDS(rds_v2)
n_cl <- length(unique(obj$seurat_clusters))
cat(sprintf("Cells: %d, Clusters: %d\n", ncol(obj), n_cl))
Idents(obj) <- "seurat_clusters"

# ══════════════════════════════════════════════════════════════
# STEP 2: FindAllMarkers
# ══════════════════════════════════════════════════════════════
cat("\n=== FindAllMarkers ===\n")

fam_file <- file.path(V2_DIR, "eps_v2_findallmarkers.csv")
if (file.exists(fam_file)) {
  cat("  Loading cached FAM...\n")
  fam <- read.csv(fam_file, stringsAsFactors=FALSE)
  fam$cluster <- factor(fam$cluster, levels=as.character(sort(unique(as.integer(fam$cluster)))))
} else {
  fam <- FindAllMarkers(
    obj,
    only.pos           = TRUE,
    min.pct            = 0.20,
    logfc.threshold    = 0.30,
    max.cells.per.ident = 400,
    test.use           = "wilcox",
    verbose            = FALSE
  )
  fam <- fam[order(fam$cluster, -fam$avg_log2FC), ]
  write.csv(fam, fam_file, row.names=FALSE)
}
cat(sprintf("  FAM: %d markers across %d clusters\n", nrow(fam), n_cl))

# Top-3 打印
cat("\n=== Top-3 markers per cluster ===\n")
top3 <- fam %>%
  group_by(cluster) %>%
  slice_max(order_by=avg_log2FC, n=3, with_ties=FALSE) %>%
  summarise(top3=paste(gene, collapse=", "), .groups="drop")
print(as.data.frame(top3), row.names=FALSE)

# ══════════════════════════════════════════════════════════════
# STEP 3: 聚类热图
# ══════════════════════════════════════════════════════════════
cat("\n=== Cluster heatmap ===\n")

top5 <- fam %>%
  group_by(cluster) %>%
  slice_max(order_by=avg_log2FC, n=5, with_ties=FALSE) %>%
  ungroup()
top_genes <- unique(top5$gene)
cat(sprintf("  Top-5 genes: %d unique\n", length(top_genes)))

set.seed(2024)
cells_sub <- unlist(lapply(levels(obj$seurat_clusters), function(cl) {
  cells <- colnames(obj)[obj$seurat_clusters == cl]
  if (length(cells) == 0) return(character(0))
  sample(cells, min(80, length(cells)))
}))
obj_sub <- obj[, cells_sub]
obj_sub <- ScaleData(obj_sub, features=top_genes, verbose=FALSE)

p_heat <- DoHeatmap(
  obj_sub, features=top_genes,
  group.by="seurat_clusters",
  raster=TRUE, label=TRUE, size=3, angle=0, hjust=0.5
) +
  scale_fill_gradientn(
    colours  = c("#1A237E","#1565C0","#42A5F5","white","#EF5350","#B71C1C"),
    values   = c(0, 0.25, 0.45, 0.55, 0.75, 1),
    na.value = "white", name="Scaled\nExpr"
  ) +
  theme(axis.text.y=element_text(size=5.5),
        plot.title=element_text(hjust=0.5, face="bold", size=12),
        plot.background=element_rect(fill="white", color=NA)) +
  labs(title=sprintf("EPS v2 — Top5 markers per cluster (%d clusters)", n_cl))

h_fig <- max(16, length(top_genes)*0.10 + 5)
ggsave(file.path(V2_DIR,"eps_v2_heatmap.pdf"),  p_heat, width=26, height=h_fig)
ggsave(file.path(V2_DIR,"eps_v2_heatmap.png"),  p_heat, width=26, height=h_fig, dpi=180)
cat("  Saved: eps_v2_heatmap.png\n")
rm(obj_sub); gc()

# ══════════════════════════════════════════════════════════════
# STEP 4: 体细胞 cluster 识别（pct.1 法）
# ══════════════════════════════════════════════════════════════
cat("\n=== Somatic cluster annotation ===\n")

somatic_markers <- list(
  Sertoli     = c("Sox9","Wt1","Amh","Ar","Aldh1a1","Prnd","Cldn11"),
  Leydig      = c("Star","Insl3","Nr2f2","Cyp11a1","Hsd3b1"),
  Myoid       = c("Acta2","Myh11","Tagln","Cnn1","Des"),
  Endothelial = c("Pecam1","Vwf","Cdh5","Tie1","Esam"),
  Macrophage  = c("Cd14","C1qb","Csf1r","Cx3cr1","Adgre1")
)

cluster_ids <- as.character(sort(unique(as.integer(fam$cluster))))
THRESHOLD   <- 0.40

somatic_score <- sapply(names(somatic_markers), function(ct) {
  gs <- somatic_markers[[ct]]
  sapply(cluster_ids, function(cl) {
    sub <- fam[fam$cluster == cl & fam$gene %in% gs, ]
    if (nrow(sub)==0) return(0)
    max(sub$pct.1, na.rm=TRUE)
  })
})
rownames(somatic_score) <- cluster_ids

cat("Somatic pct.1 scores:\n")
print(round(somatic_score, 2))

cluster_ann <- sapply(cluster_ids, function(cl) {
  sc <- somatic_score[cl, ]
  if (max(sc) >= THRESHOLD) return(names(which.max(sc)))
  return("Germ_cells")
})

cat("\nCluster → Type:\n")
print(data.frame(cluster=cluster_ids, type=cluster_ann), row.names=FALSE)

obj$celltype_broad <- unname(cluster_ann[as.character(obj$seurat_clusters)])
cat("\nCell counts by broad type:\n")
print(sort(table(obj$celltype_broad), decreasing=TRUE))

# ── DotPlot（已知 marker）──────────────────────────────────
marker_gene <- c(
  "Ddx4","Ret","Uchl1","Stra8","Sycp3","Tnp1",
  "Zbtb16","Id4","Gfra1","Gfra2",
  "Amh","Wt1","Sox9","Ar","Aldh1a1","Prnd",
  "Acta2","Myh11",
  "Cd14","Cd4",
  "Pecam1","Vwf",
  "Star","Insl3","Nr2f2"
)
marker_found <- marker_gene[marker_gene %in% rownames(obj)]

dp_data <- DotPlot(obj, features=marker_found, group.by="seurat_clusters")$data
dp_data <- dp_data %>%
  group_by(features.plot) %>%
  mutate(exp_norm={
    rng <- range(avg.exp.scaled, na.rm=TRUE)
    if (diff(rng)<1e-6) rep(0,n())
    else (avg.exp.scaled-rng[1])/(rng[2]-rng[1])
  }) %>% ungroup()

dp_data$celltype <- cluster_ann[as.character(dp_data$id)]

p_dot <- ggplot(dp_data, aes(x=features.plot, y=id,
                              size=pct.exp, color=exp_norm)) +
  geom_point() +
  scale_color_gradientn(
    colours=c("grey93","grey80","#AED6F1","#2E86C1","#1A5276"),
    values=c(0,0.2,0.5,0.75,1), limits=c(0,1),
    name="Relative\nExpression"
  ) +
  scale_size_continuous(range=c(0.3,8), name="% Expressing",
                        breaks=c(10,25,50,75)) +
  scale_y_discrete(position="right") +
  labs(x=NULL, y=NULL, title="EPS v2 — Marker DotPlot (all clusters)") +
  theme_classic(base_size=10) +
  theme(axis.text.x=element_text(angle=45,hjust=1,face="italic",size=9),
        axis.text.y=element_text(size=8),
        axis.line=element_blank(),
        panel.border=element_rect(color="grey50",fill=NA,linewidth=0.6),
        plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white",color=NA))

ggsave(file.path(V2_DIR,"eps_v2_dotplot.pdf"), p_dot, width=18, height=14)
ggsave(file.path(V2_DIR,"eps_v2_dotplot.png"), p_dot, width=18, height=14, dpi=200)
cat("  Saved: eps_v2_dotplot.png\n")

write.csv(
  data.frame(cluster=cluster_ids, celltype=cluster_ann),
  file.path(V2_DIR,"broad_cluster_annotation.csv"), row.names=FALSE
)

# ══════════════════════════════════════════════════════════════
# STEP 5: 生殖细胞提取 + 亚群分期
# ══════════════════════════════════════════════════════════════
cat("\n=== Germ cell sub-clustering ===\n")
germ <- subset(obj, celltype_broad=="Germ_cells")
cat(sprintf("Germ cells: %d\n", ncol(germ)))

germ_markers_raw <- c(
  "NANOS3","PAX7","RET","TAF4B","CDH1","ZBTB16","GFRA1","UCHL1","CTCFL",
  "KIT","STRA8",
  "REC8","MEIOC","SYCP3","SYCP2","SPO11","MLH3","SYCP1","PIWIL1",
  "ACRV1","SPACA1","TNP1","PRM2"
)
to_mouse <- function(g) paste0(toupper(substr(g,1,1)), tolower(substr(g,2,nchar(g))))
germ_markers     <- sapply(germ_markers_raw, to_mouse)
germ_markers_found <- germ_markers[germ_markers %in% rownames(germ)]
cat(sprintf("Germ markers: %d / %d found\n",
            length(germ_markers_found), length(germ_markers)))

# 去除 ambient RNA 基因（防止扭曲 PCA）
ambient_remove <- c("Tnp1","Prm2","Prm1","Akap4","Ldhc","Crisp2",
                    "Spata3","Tuba3a","Ace","Rnase10","Spam1","H2al3")

germ <- NormalizeData(germ, verbose=FALSE)
germ <- FindVariableFeatures(germ, nfeatures=3000, verbose=FALSE)
hvg  <- VariableFeatures(germ)
hvg  <- hvg[!hvg %in% ambient_remove]
feat <- unique(c(germ_markers_found, hvg))
feat <- feat[feat %in% rownames(germ)]

germ <- ScaleData(germ, features=feat, verbose=FALSE)
germ <- RunPCA(germ, features=feat, npcs=30, seed.use=2024, verbose=FALSE)
germ <- RunHarmony(germ, group.by.vars="sample_id",
                   reduction.use="pca", dims.use=1:20, verbose=FALSE)
germ <- RunUMAP(germ, reduction="harmony", dims=1:20, seed.use=2024, verbose=FALSE)
germ <- FindNeighbors(germ, reduction="harmony", dims=1:20, verbose=FALSE)

for (res in c(0.3, 0.5, 0.8)) {
  germ <- FindClusters(germ, resolution=res, random.seed=2024, verbose=FALSE)
  cat(sprintf("  res=%.1f → %d germ clusters\n", res,
              length(unique(germ$seurat_clusters))))
}
germ <- FindClusters(germ, resolution=0.8, random.seed=2024, verbose=FALSE)
n_g  <- length(unique(germ$seurat_clusters))
cat(sprintf("Using res=0.8: %d germ clusters\n", n_g))

# ── 发育阶段注释 ──────────────────────────────────────────
counts_g <- GetAssayData(germ, assay="RNA", layer="counts")
clusters  <- as.character(sort(unique(germ$seurat_clusters)))

ann_markers <- list(
  "未分化精原细胞" = c("Nanos3","Ret","Zbtb16","Gfra1","Uchl1"),
  "分化型精原细胞" = c("Kit","Stra8"),
  "精母细胞"       = c("Sycp3","Sycp1","Rec8","Meioc"),
  "精细胞"         = c("Acrv1","Spaca1")
)

stage_pct <- sapply(clusters, function(cl) {
  cells <- colnames(germ)[germ$seurat_clusters == cl]
  sapply(names(ann_markers), function(st) {
    gs <- ann_markers[[st]]; gs <- gs[gs %in% rownames(counts_g)]
    if (length(gs)==0) return(0)
    mean(sapply(gs, function(g) mean(counts_g[g, cells] > 0))) * 100
  })
})
cat("\nStage scores:\n"); print(round(t(stage_pct), 1))

assign_stage <- function(cl) {
  u <- stage_pct["未分化精原细胞", cl]
  d <- stage_pct["分化型精原细胞", cl]
  m <- stage_pct["精母细胞",       cl]
  s <- stage_pct["精细胞",         cl]
  n <- sum(germ$seurat_clusters == cl)
  if (u >= 15)  return("未分化精原细胞")
  if (d >= 12)  return("分化型精原细胞")
  if (s >= 65)  return("精细胞")
  if (n <= 300 && u < 5 && d < 5 && m < 25 && s < 40) return("成熟精子")
  return("精母细胞")
}
germ_stage <- setNames(sapply(clusters, assign_stage), clusters)
cat("\nCluster → Stage:\n"); print(germ_stage)

stage_order <- c("未分化精原细胞","分化型精原细胞","精母细胞","精细胞","成熟精子")
germ$stage_label <- factor(
  unname(germ_stage[as.character(germ$seurat_clusters)]),
  levels = stage_order
)
cat("\nStage counts:\n"); print(table(germ$stage_label))

# ══════════════════════════════════════════════════════════════
# STEP 6: 绘图
# ══════════════════════════════════════════════════════════════
cat("\n=== Plotting ===\n")

stage_colors <- c(
  "未分化精原细胞"="#D62728","分化型精原细胞"="#FF7F0E",
  "精母细胞"="#1F77B4","精细胞"="#2CA02C","成熟精子"="#9467BD"
)

# A: 生殖细胞 UMAP
p_s <- DimPlot(germ, reduction="umap", group.by="stage_label",
               cols=stage_colors, pt.size=0.2) +
  labs(title="EPS v2 — Germ cells by stage") +
  theme_cowplot(font_size=12) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white")) +
  guides(color=guide_legend(override.aes=list(size=4)))

p_c <- DimPlot(germ, reduction="umap", label=TRUE, label.size=3.5,
               pt.size=0.15) +
  labs(title=sprintf("Germ: %d clusters", n_g)) +
  theme_cowplot(font_size=11) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white"),
        legend.position="none")

p_g <- DimPlot(germ, reduction="umap", group.by="group",
               cols=c("Rescue"="#4DBBD5","Aging"="#E64B35"), pt.size=0.2) +
  labs(title="Aging vs Rescue") +
  theme_cowplot(font_size=12) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white"))

key_f <- c("Zbtb16","Kit","Sycp3","Acrv1","Tnp1")
key_f <- key_f[key_f %in% rownames(germ)]
p_feats <- lapply(key_f, function(g) {
  FeaturePlot(germ, features=g, reduction="umap", pt.size=0.12, order=TRUE) +
    scale_color_viridis_c(option="D") +
    theme_cowplot(font_size=9) +
    theme(plot.background=element_rect(fill="white"),
          plot.title=element_text(hjust=0.5,face="bold.italic"),
          axis.title=element_blank(),axis.text=element_blank(),
          axis.ticks=element_blank())
})

p_germ_all <- (p_s | p_c | p_g) / wrap_plots(p_feats, nrow=1) +
  plot_layout(heights=c(1.2, 0.9)) +
  plot_annotation(
    title="EPS v2 — Germ Cell Sub-clustering & Stage Annotation",
    theme=theme(plot.title=element_text(size=13,face="bold",hjust=0.5),
                plot.background=element_rect(fill="white"))
  )

ggsave(file.path(V2_DIR,"eps_v2_germcell_overview.pdf"), p_germ_all, width=18, height=12)
ggsave(file.path(V2_DIR,"eps_v2_germcell_overview.png"), p_germ_all, width=18, height=12, dpi=200)
cat("  Saved: eps_v2_germcell_overview.png\n")

# B: DotPlot（生殖细胞发育阶段）
gene_list_f <- list(
  "未分化精原细胞"=c("Nanos3","Pax7","Ret","Taf4b","Cdh1","Zbtb16","Gfra1","Uchl1","Ctcfl"),
  "分化型精原细胞"=c("Kit","Stra8"),
  "精母细胞"=c("Rec8","Meioc","Sycp3","Sycp2","Spo11","Mlh3","Sycp1","Piwil1"),
  "精细胞"=c("Acrv1","Spaca1"),
  "成熟精子"=c("Tnp1","Prm2")
)
gene_list_f <- lapply(gene_list_f, function(g) g[g %in% rownames(germ)])
gene_order  <- unlist(gene_list_f)

Idents(germ) <- "seurat_clusters"
dp_raw <- DotPlot(germ, features=gene_order, group.by="seurat_clusters")$data
dp_raw$stage <- factor(unname(germ_stage[as.character(dp_raw$id)]), levels=stage_order)
dp_raw$features.plot <- factor(dp_raw$features.plot, levels=gene_order)

dp_st <- dp_raw %>%
  group_by(stage, features.plot) %>%
  summarise(
    avg.exp.scaled = mean(avg.exp.scaled, na.rm=TRUE),
    pct.exp        = mean(pct.exp,        na.rm=TRUE),
    .groups="drop"
  ) %>%
  group_by(features.plot) %>%
  mutate(exp_norm = {
    rng <- range(avg.exp.scaled, na.rm=TRUE)
    if (diff(rng) < 1e-6) rep(0, n())
    else (avg.exp.scaled - rng[1]) / (rng[2] - rng[1])
  }) %>% ungroup()

dp_st$stage <- factor(dp_st$stage, levels=rev(stage_order))

sep_after <- which(diff(as.integer(factor(
  rep(names(gene_list_f), sapply(gene_list_f, length)),
  levels=names(gene_list_f)))) != 0)
sep_x <- sep_after + 0.5

p_gdot <- ggplot(dp_st, aes(x=features.plot, y=stage, size=pct.exp, color=exp_norm)) +
  geom_vline(xintercept=sep_x, color="grey72", linewidth=0.4, linetype="dashed") +
  geom_point() +
  scale_color_gradientn(
    colours=c("grey92","grey80","#AED6F1","#2E86C1","#1A5276"),
    values=c(0,0.2,0.5,0.75,1), limits=c(0,1),
    name="Relative\nExpression\n(per gene)"
  ) +
  scale_size_continuous(range=c(0.5,9), name="% Expressing",
                        breaks=c(10,25,50,75,100)) +
  scale_y_discrete(position="right") +
  labs(x=NULL, y=NULL, title="EPS v2 — Germ cells by stage") +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=45,hjust=1,face="italic",size=10),
        axis.text.y=element_text(size=11,face="bold"),
        axis.line=element_blank(),
        panel.border=element_rect(color="grey40",fill=NA,linewidth=0.7),
        panel.grid.major.y=element_line(color="grey93",linewidth=0.3),
        plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white",color=NA))

# 顶部色条
sec_colors <- c("未分化精原细胞"="#D62728","分化型精原细胞"="#FF7F0E",
                "精母细胞"="#1F77B4","精细胞"="#2CA02C","成熟精子"="#9467BD")
sec_df <- data.frame(
  gene=gene_order,
  grp =rep(names(gene_list_f), sapply(gene_list_f, length)),
  x   =seq_along(gene_order)
)
sec_df$grp <- factor(sec_df$grp, levels=names(gene_list_f))
sec_ctr <- sec_df %>% group_by(grp) %>% summarise(cx=mean(x), .groups="drop")

p_topbar <- ggplot(sec_df, aes(x=x, y=1, fill=grp)) +
  geom_tile(width=0.95, height=0.9) +
  geom_text(data=sec_ctr, aes(x=cx, y=1, label=grp, fill=NULL),
            size=2.8, fontface="bold", color="white") +
  scale_fill_manual(values=sec_colors) +
  scale_x_continuous(expand=c(0,0.5)) +
  theme_void() + theme(legend.position="none",
                       plot.background=element_rect(fill="white",color=NA))

p_gdot_full <- p_topbar / p_gdot +
  plot_layout(heights=c(1,10)) +
  plot_annotation(theme=theme(plot.background=element_rect(fill="white")))

ggsave(file.path(V2_DIR,"eps_v2_germcell_dotplot.pdf"), p_gdot_full, width=15, height=7)
ggsave(file.path(V2_DIR,"eps_v2_germcell_dotplot.png"), p_gdot_full, width=15, height=7, dpi=200)
cat("  Saved: eps_v2_germcell_dotplot.png\n")

# C: 全局 UMAP（含体细胞）
p_umap_broad <- DimPlot(obj, reduction="umap", label=TRUE, label.size=3.5,
                         pt.size=0.15) +
  labs(title=sprintf("EPS v2 — %d clusters", n_cl)) +
  theme_cowplot(font_size=11) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white"),
        legend.position="none")

p_umap_grp <- DimPlot(obj, reduction="umap", group.by="group",
                       cols=c("Rescue"="#4DBBD5","Aging"="#E64B35"),
                       pt.size=0.15) +
  labs(title="Aging vs Rescue") +
  theme_cowplot(font_size=11) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white"))

p_umap_all <- (p_umap_broad | p_umap_grp) +
  plot_annotation(title="EPS v2 — All cells",
                  theme=theme(plot.title=element_text(size=13,face="bold",hjust=0.5),
                              plot.background=element_rect(fill="white")))

ggsave(file.path(V2_DIR,"eps_v2_umap_all.png"), p_umap_all, width=14, height=6, dpi=200)
cat("  Saved: eps_v2_umap_all.png\n")

# ══════════════════════════════════════════════════════════════
# STEP 7: 保存
# ══════════════════════════════════════════════════════════════
cat("\n=== Saving ===\n")
saveRDS(germ, file.path(EPS_OUT, "rds/seurat_eps_germcells_v2.rds"))
write.csv(
  data.frame(cluster=names(germ_stage), stage=germ_stage),
  file.path(V2_DIR,"germ_stage_annotation.csv"), row.names=FALSE
)
cat("  Saved: seurat_eps_germcells_v2.rds\n")

cat("\n===== DONE =====\n")
cat(sprintf("Output dir: %s/v2/\n", EPS_OUT))
cat("  eps_v2_heatmap.png\n")
cat("  eps_v2_dotplot.png\n")
cat("  eps_v2_germcell_overview.png\n")
cat("  eps_v2_germcell_dotplot.png\n")
cat("  rds/seurat_eps_germcells_v2.rds\n")

#!/usr/bin/env Rscript
# ============================================================
# Exos v2 — 8 cell-type DotPlot (全细胞类型)
# 分类: elongating / innate_Lymphoid / Leydig / Macrophage /
#       round_spermatid / Sertoli / spermatocyte / spermatogonia
#
# 输入:
#   seurat_exos_clustered_v2.rds   — 全对象 (95,331 cells, 29 clusters)
#   seurat_exos_germcells_v2.rds   — 生殖细胞子对象 (73,532 cells, 22 sub-clusters)
#   broad_cluster_annotation.csv   — 全对象 cluster → 体细胞/生殖细胞
#   germ_stage_annotation_v2.csv   — 生殖细胞 sub-cluster → 6 阶段
# ============================================================

.libPaths(c("/Users/zzp/Code/tianhairui/.rlib", .libPaths()))
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(cowplot)
})

EXOS_OUT <- "/Users/zzp/Code/tianhairui/outputs0409/seurat5_exos_v2"
A_DIR    <- file.path(EXOS_OUT, "02_annotation")
dir.create(A_DIR, recursive=TRUE, showWarnings=FALSE)

# ══════════════════════════════════════════════════════════════
# 1. 加载对象
# ══════════════════════════════════════════════════════════════
cat("=== Loading objects ===\n")
obj  <- readRDS(file.path(EXOS_OUT, "rds/seurat_exos_clustered_v2.rds"))
germ <- readRDS(file.path(EXOS_OUT, "rds/seurat_exos_germcells_v2.rds"))
cat(sprintf("Full: %d cells, %d clusters\n", ncol(obj), length(unique(obj$seurat_clusters))))
cat(sprintf("Germ: %d cells, %d sub-clusters\n", ncol(germ), length(unique(germ$seurat_clusters))))

# ══════════════════════════════════════════════════════════════
# 2. 构建 8-type 注释
# ══════════════════════════════════════════════════════════════
cat("\n=== Building 8-type annotation ===\n")

# 读取注释表
broad_ann  <- read.csv(file.path(EXOS_OUT, "01_heatmap/broad_cluster_annotation.csv"),
                       stringsAsFactors=FALSE)
germ_stage <- read.csv(file.path(A_DIR, "germ_stage_annotation_v2.csv"),
                       stringsAsFactors=FALSE)

# 生殖细胞 sub-cluster → 6 阶段映射
germ_cl_to_stage <- setNames(germ_stage$stage, as.character(germ_stage$cluster))
germ_cells        <- germ$seurat_clusters
names(germ_cells) <- colnames(germ)

# 6 阶段 → 8 type 映射
stage_to_8 <- c(
  "ProSPG"              = "spermatogonia",
  "Undiff SG"           = "spermatogonia",
  "Diff SG"             = "spermatogonia",
  "Spermatocyte"        = "spermatocyte",
  "Round spermatid"     = "round_spermatid",
  "Elongated spermatid" = "elongating"
)

# 全对象 cluster → 体细胞 type 映射
somatic_to_8 <- c(
  "Sertoli"     = "Sertoli",
  "Leydig"      = "Leydig",
  "Macrophage"  = "Macrophage",
  "Endothelial" = "innate_Lymphoid",   # testis ILC/NK cells
  "Myoid"       = "innate_Lymphoid"
)

# 全对象 cluster → broad type 快速查找
broad_map <- setNames(broad_ann$celltype, as.character(broad_ann$cluster))

# 为每个细胞分配 8-type
cat("Assigning 8-type labels...\n")
cl_full  <- as.character(obj$seurat_clusters)
type_vec <- character(ncol(obj))
names(type_vec) <- colnames(obj)

# 先处理体细胞
for (i in seq_along(type_vec)) {
  cell  <- names(type_vec)[i]
  broad <- broad_map[cl_full[i]]
  if (is.na(broad) || broad == "Germ_cells") next   # 后面处理
  type_vec[i] <- somatic_to_8[broad]
}

# 处理生殖细胞（用 germ sub-cluster 的 stage）
germ_common <- intersect(names(germ_cells), colnames(obj))
cat(sprintf("Germ cells matched: %d / %d\n", length(germ_common), ncol(germ)))

sub_cl   <- as.character(germ_cells[germ_common])
stage_6  <- unname(germ_cl_to_stage[sub_cl])
# stage_6 may be NA for un-annotated sub-clusters (e.g. cluster 23)
stage_6[is.na(stage_6)] <- "Spermatocyte"   # fallback
type_8_g <- unname(stage_to_8[stage_6])
type_8_g[is.na(type_8_g)] <- "spermatocyte" # safety fallback
type_vec[germ_common] <- type_8_g

# 有些 broad=Germ_cells 的细胞不在 germ 子对象里（可能被过滤）→ 默认 spermatocyte
remaining_germ <- names(type_vec)[type_vec == "" & !is.na(broad_map[cl_full]) & broad_map[cl_full] == "Germ_cells"]
if (length(remaining_germ) > 0) {
  type_vec[remaining_germ] <- "spermatocyte"
  cat(sprintf("  %d germ cells not in sub-obj → labelled spermatocyte\n",
              length(remaining_germ)))
}

# 未分配细胞（cluster 不在 broad_map 中）→ 默认 spermatocyte
unassigned <- names(type_vec)[type_vec == "" | is.na(type_vec)]
if (length(unassigned) > 0) {
  type_vec[unassigned] <- "spermatocyte"
  cat(sprintf("  %d unassigned cells → labelled spermatocyte\n", length(unassigned)))
}

obj$celltype_8 <- factor(type_vec, levels=c(
  "elongating","innate_Lymphoid","Leydig","Macrophage",
  "round_spermatid","Sertoli","spermatocyte","spermatogonia"
))

cat("\n8-type counts:\n")
print(sort(table(obj$celltype_8), decreasing=TRUE))

# ══════════════════════════════════════════════════════════════
# 3. Marker 基因
# ══════════════════════════════════════════════════════════════
cat("\n=== Defining marker genes ===\n")

gene_list <- list(
  "spermatogonia"  = c("Id1","Id2","Klf4","Zbtb16","Uchl1","Gfra1",
                       "Nanos3","Sall4","Etv5","Egr4","Kit","Stra8"),
  "spermatocyte"   = c("Sycp3","Hormad1","H2afx","Piwil1","Spo11",
                       "Meioc","Rec8","Mei1","Piwil2","Mns1"),
  "round_spermatid"= c("Acrv1","Spaca1","Tssk1","Pgk2","Tsga8","Hspa1l"),
  "elongating"     = c("Prm1","Prm2","Tnp1","Tnp2","Akap4","Spem1","Odf1"),
  "Sertoli"        = c("Sox9","Wt1","Amh","Ar","Aldh1a1","Prnd","Cldn11"),
  "Leydig"         = c("Star","Insl3","Nr2f2","Cyp11a1","Hsd3b1"),
  "Macrophage"     = c("Cd14","C1qb","Csf1r","Cx3cr1","Adgre1"),
  "innate_Lymphoid"= c("Ncr1","Nkg7","Gata3","Il7r","Klrb1c","Thy1")
)

# 过滤不存在的基因
gene_list <- lapply(gene_list, function(gs) unique(gs[gs %in% rownames(obj)]))
gene_list <- gene_list[sapply(gene_list, length) > 0]

cat("Genes found per type:\n")
print(sapply(gene_list, length))

# 移除在 obj 中完全不存在的空组
gene_order <- unlist(gene_list)
cat(sprintf("Total marker genes: %d\n", length(gene_order)))

# ══════════════════════════════════════════════════════════════
# 4. DotPlot 数据
# ══════════════════════════════════════════════════════════════
cat("\n=== Computing DotPlot ===\n")

# 按 celltype_8 汇总（Seurat DotPlot 需要 Idents 或 group.by）
Idents(obj) <- "celltype_8"
dp_data <- DotPlot(obj, features=gene_order, group.by="celltype_8")$data
dp_data$features.plot <- factor(dp_data$features.plot, levels=gene_order)

# Per-gene rescaling (消除 ambient RNA 尺度问题)
dp_data <- dp_data %>%
  group_by(features.plot) %>%
  mutate(exp_norm = {
    rng <- range(avg.exp.scaled, na.rm=TRUE)
    if (diff(rng) < 1e-6) rep(0, n())
    else (avg.exp.scaled - rng[1]) / (rng[2] - rng[1])
  }) %>% ungroup()

# Y 轴顺序（用户指定，从上到下）
type_order <- c("elongating","innate_Lymphoid","Leydig","Macrophage",
                "round_spermatid","Sertoli","spermatocyte","spermatogonia")
dp_data$id <- factor(dp_data$id, levels=rev(type_order))   # rev → 第一个在最上

# 分隔线位置
seg_grp <- factor(
  rep(names(gene_list), sapply(gene_list, length)),
  levels=names(gene_list)
)
sep_after <- which(diff(as.integer(seg_grp)) != 0)
sep_x     <- sep_after + 0.5

# 配色
type_colors <- c(
  "spermatogonia"   = "#D62728",
  "spermatocyte"    = "#1F77B4",
  "round_spermatid" = "#2CA02C",
  "elongating"      = "#9467BD",
  "Sertoli"         = "#8C564B",
  "Leydig"          = "#E377C2",
  "Macrophage"      = "#7F7F7F",
  "innate_Lymphoid" = "#17BECF"
)

# ══════════════════════════════════════════════════════════════
# 5. 绘图
# ══════════════════════════════════════════════════════════════
cat("\n=== Plotting ===\n")

# 主 DotPlot
p_dot <- ggplot(dp_data,
                aes(x=features.plot, y=id, size=pct.exp, color=exp_norm)) +
  geom_vline(xintercept=sep_x, color="grey72", linewidth=0.4, linetype="dashed") +
  geom_point() +
  scale_color_gradientn(
    colours = c("grey92","grey82","#AED6F1","#2E86C1","#1A5276"),
    values  = c(0, 0.2, 0.5, 0.75, 1),
    limits  = c(0, 1),
    name    = "Relative\nExpression\n(per gene)"
  ) +
  scale_size_continuous(
    range  = c(0.4, 9),
    name   = "% Expressing",
    breaks = c(10, 25, 50, 75, 100)
  ) +
  scale_y_discrete(position="right") +
  labs(x=NULL, y=NULL,
       title="Exos v2 — Cell type marker DotPlot (8 types)") +
  theme_classic(base_size=12) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, face="italic", size=9.5),
    axis.text.y = element_text(size=11, face="bold.italic"),
    axis.line   = element_blank(),
    panel.border= element_rect(color="grey40", fill=NA, linewidth=0.7),
    panel.grid.major.y = element_line(color="grey93", linewidth=0.3),
    plot.title  = element_text(hjust=0.5, face="bold", size=13),
    plot.background = element_rect(fill="white", color=NA),
    legend.title = element_text(size=9),
    legend.text  = element_text(size=8.5)
  )

# 顶部色条
sec_df <- data.frame(
  gene = gene_order,
  grp  = seg_grp,
  x    = seq_along(gene_order)
)
sec_ctr <- sec_df %>%
  group_by(grp) %>%
  summarise(cx=mean(x), .groups="drop")

p_bar <- ggplot(sec_df, aes(x=x, y=1, fill=grp)) +
  geom_tile(width=0.95, height=0.9) +
  geom_text(data=sec_ctr, aes(x=cx, y=1, label=grp, fill=NULL),
            size=2.5, fontface="bold", color="white") +
  scale_fill_manual(values=type_colors[levels(sec_df$grp)]) +
  scale_x_continuous(expand=c(0, 0.5)) +
  theme_void() +
  theme(legend.position="none",
        plot.background=element_rect(fill="white", color=NA))

p_final <- p_bar / p_dot +
  plot_layout(heights=c(1, 10)) +
  plot_annotation(
    theme=theme(plot.background=element_rect(fill="white", color=NA))
  )

# 保存
n_genes <- length(gene_order)
w <- max(16, n_genes * 0.33 + 4)
ggsave(file.path(A_DIR,"exos_v2_8type_dotplot.pdf"), p_final, width=w, height=7)
ggsave(file.path(A_DIR,"exos_v2_8type_dotplot.png"), p_final, width=w, height=7, dpi=220)
cat(sprintf("Saved: exos_v2_8type_dotplot.png  (%.1f x 7 in)\n", w))

# ── UMAP with 8 types ───────────────────────────────────────
p_umap <- DimPlot(obj, reduction="umap", group.by="celltype_8",
                  cols=type_colors, pt.size=0.15, label=TRUE, label.size=3) +
  labs(title="Exos v2 — 8 cell types") +
  theme_cowplot(font_size=11) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white")) +
  guides(color=guide_legend(override.aes=list(size=4)))

p_umap_grp <- DimPlot(obj, reduction="umap", group.by="group",
                       cols=c("Rescue"="#4DBBD5","Aging"="#E64B35"), pt.size=0.15) +
  labs(title="Aging vs Rescue") +
  theme_cowplot(font_size=11) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white"))

p_umap_all <- p_umap | p_umap_grp
ggsave(file.path(A_DIR,"exos_v2_8type_umap.pdf"),  p_umap_all, width=16, height=6.5)
ggsave(file.path(A_DIR,"exos_v2_8type_umap.png"),  p_umap_all, width=16, height=6.5, dpi=200)
cat("Saved: exos_v2_8type_umap.png\n")

cat("\n===== DONE =====\n")
print(table(obj$celltype_8))

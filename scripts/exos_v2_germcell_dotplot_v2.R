#!/usr/bin/env Rscript
# ============================================================
# Exos v2 — Germ cell DotPlot v2
# 6 stages: ProSPG / Undiff SG / Diff SG / Spermatocyte /
#           Round spermatid / Elongated spermatid
# 输入: seurat_exos_germcells_v2.rds
# 输出: exos_v2_germcell_dotplot_v2.{pdf,png}
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
# 1. 加载
# ══════════════════════════════════════════════════════════════
cat("=== Loading germ cell object ===\n")
germ <- readRDS(file.path(EXOS_OUT, "rds/seurat_exos_germcells_v2.rds"))
cat(sprintf("Cells: %d, Clusters: %d\n", ncol(germ), length(unique(germ$seurat_clusters))))

counts_g <- GetAssayData(germ, assay="RNA", layer="counts")
clusters  <- as.character(sort(unique(as.integer(germ$seurat_clusters))))
clusters  <- clusters[clusters != "23"]   # 空 cluster

# ══════════════════════════════════════════════════════════════
# 2. 6 阶段注释
#
# 关键观察（预分析）：
#   - Akap4/Spem1/Tnp2/Odf1/Pgk2/Tssk1 均为 ambient RNA (60-100%)
#   - Kit/Stra8/Dmrtb1 也是 ambient (20-57%，不能用均值)
#   - 必须用 INDIVIDUAL 基因阈值，而不是多基因均值
#
# 可靠标记（各 cluster 真实差异大）：
#   Id2    — ProSPG:  cluster 20 = 55%，其余 < 4%
#   Zbtb16 — Undiff:  cluster 14 = 27.5%，其余 < 5%
#   Kit    — Diff SG: cluster 14 = 38%，cluster 16 = 8.7%，其余 < 6%
#   Sycp3  — Spctye:  cluster 12 = 95.5%，cluster 6/22 = 75-82%，其余 38-70%
#   Tssk1  — Round:   clusters 4/9/11/15 = 96-98%，cluster 8 = 89%
# ══════════════════════════════════════════════════════════════
cat("\n=== 6-stage annotation (individual gene thresholds) ===\n")

# 单基因 pct_detect（用原始 counts，不平均）
pct1 <- function(cl, gene) {
  if (!gene %in% rownames(counts_g)) return(0)
  cells <- colnames(germ)[germ$seurat_clusters == cl]
  mean(counts_g[gene, cells] > 0) * 100
}

# 打印关键个体基因值
cat("\nKey individual gene % per cluster:\n")
cat(sprintf("%-4s %5s %6s %4s %6s %6s %5s %5s\n",
    "Cl","n","Id2","Zbtb16","Kit","Sycp3","Tssk1","Stage_hint"))
for (cl in clusters) {
  n <- sum(germ$seurat_clusters == cl)
  id2    <- pct1(cl,"Id2")
  zbtb16 <- pct1(cl,"Zbtb16")
  kit    <- pct1(cl,"Kit")
  sycp3  <- pct1(cl,"Sycp3")
  tssk1  <- pct1(cl,"Tssk1")
  hint <- if (id2 >= 25) "ProSPG"
          else if (zbtb16 >= 15) "Undiff SG"
          else if (kit >= 7)  "Diff SG"
          else if (sycp3 >= 70) "Spctye"
          else if (tssk1 >= 80) "Round"
          else if (sycp3 >= 40) "Spctye*"
          else "Elongated"
  cat(sprintf("%-4s %5d %6.1f %6.1f %4.1f %6.1f %5.1f  %s\n",
      cl, n, id2, zbtb16, kit, sycp3, tssk1, hint))
}

# 阶段分配（基于单基因阈值，优先级顺序）
assign_stage <- function(cl) {
  id2    <- pct1(cl,"Id2")
  zbtb16 <- pct1(cl,"Zbtb16")
  kit    <- pct1(cl,"Kit")
  sycp3  <- pct1(cl,"Sycp3")
  tssk1  <- pct1(cl,"Tssk1")

  # 严格阈值（基于预分析知道各 marker 的真实分布）
  if (id2    >= 25) return("ProSPG")          # 仅 cluster 20 (Id2=55%)
  if (zbtb16 >= 15) return("Undiff SG")       # 仅 cluster 14 (Zbtb16=27.5%)
  if (kit    >=  7) return("Diff SG")         # cluster 16 (Kit=8.7%)
  if (sycp3  >= 70) return("Spermatocyte")    # clusters 12, 6, 22
  if (tssk1  >= 80) return("Round spermatid") # clusters 4, 8, 9, 11, 15, 21
  # 中等 Sycp3 → 精母细胞（大簇，原始注释的 精母细胞）
  if (sycp3  >= 40) return("Spermatocyte")
  # 其余 spermatid-like clusters → 长形精子
  return("Elongated spermatid")
}

germ_stage <- setNames(sapply(clusters, assign_stage), clusters)
cat("\nCluster → Stage:\n")
print(germ_stage)

stage_order <- c("ProSPG","Undiff SG","Diff SG",
                 "Spermatocyte","Round spermatid","Elongated spermatid")

# 填充空 cluster（如 cluster 23）
all_cl <- as.character(germ$seurat_clusters)
stage_vec <- sapply(all_cl, function(cl) {
  if (cl %in% names(germ_stage)) germ_stage[cl] else "Spermatocyte"
})
germ$stage_label <- factor(unname(stage_vec), levels=stage_order)
cat("\nStage counts:\n")
print(table(germ$stage_label))

# ══════════════════════════════════════════════════════════════
# 3. DotPlot — 更新 marker 基因列表
# ══════════════════════════════════════════════════════════════
cat("\n=== Building DotPlot ===\n")

gene_list_f <- list(
  "ProSPG"              = c("Id1","Id2","Klf4"),
  "Undiff SG"           = c("Zbtb16","Uchl1","Id4","Cdh1","Sall4",
                             "Gfra1","Etv5","Nanos3","Egr4"),
  "Diff SG"             = c("Dmrtb1","Kit","Rhox13","Stra8",
                             "Dnmt3b","Prdm9","Sohlh1"),
  "Spermatocyte"        = c("Tex101","Piwil1","Spo11","Mei1","Piwil2",
                             "Spag6","Sycp3","Tdrd5","Tbpl1","Hormad1",
                             "H2afx","Mns1"),
  "Round spermatid"     = c("Hspa1l","Acrv1","Tssk1","Pgk2",
                             "Spaca1","Tsga8"),
  "Elongated spermatid" = c("Prm1","Prm2","Tnp1","Spem1","Akap4","Odf1")
)

# 过滤不存在的基因
gene_list_f <- lapply(gene_list_f, function(gs) {
  unique(gs[gs %in% rownames(germ)])
})
gene_list_f <- gene_list_f[sapply(gene_list_f, length) > 0]

cat("Genes per stage:\n")
print(sapply(gene_list_f, length))

gene_order <- unlist(gene_list_f)

# 获取 DotPlot 数据
Idents(germ) <- "seurat_clusters"
dp_raw <- DotPlot(germ, features=gene_order, group.by="seurat_clusters")$data
dp_raw$stage <- factor(
  unname(stage_vec)[match(as.character(dp_raw$id),
                          as.character(germ$seurat_clusters[!duplicated(germ$seurat_clusters)]))],
  levels=stage_order
)

# 更可靠的 stage 映射
cl_stage_map <- setNames(
  sapply(as.character(unique(germ$seurat_clusters)), function(cl) {
    if (cl %in% names(germ_stage)) germ_stage[cl] else "Spermatocyte"
  }),
  as.character(unique(germ$seurat_clusters))
)
dp_raw$stage <- factor(
  unname(cl_stage_map[as.character(dp_raw$id)]),
  levels=stage_order
)
dp_raw$features.plot <- factor(dp_raw$features.plot, levels=gene_order)

# 按 stage 聚合
dp_st <- dp_raw %>%
  filter(!is.na(stage)) %>%
  group_by(stage, features.plot) %>%
  summarise(
    avg.exp.scaled = mean(avg.exp.scaled, na.rm=TRUE),
    pct.exp        = mean(pct.exp,        na.rm=TRUE),
    .groups = "drop"
  ) %>%
  # Per-gene rescaling（消除 ambient RNA 尺度问题）
  group_by(features.plot) %>%
  mutate(exp_norm = {
    rng <- range(avg.exp.scaled, na.rm=TRUE)
    if (diff(rng) < 1e-6) rep(0, n())
    else (avg.exp.scaled - rng[1]) / (rng[2] - rng[1])
  }) %>% ungroup()

# y 轴：从上到下 ProSPG → Undiff SG → Diff SG → Spermatocyte → Round → Elongated
dp_st$stage <- factor(dp_st$stage, levels=rev(stage_order))

# 分组分隔线位置
seg_grp <- rep(names(gene_list_f), sapply(gene_list_f, length))
seg_grp <- factor(seg_grp, levels=names(gene_list_f))
sep_after <- which(diff(as.integer(seg_grp)) != 0)
sep_x     <- sep_after + 0.5

# 配色
stage_colors <- c(
  "ProSPG"              = "#7B2D8B",
  "Undiff SG"           = "#D62728",
  "Diff SG"             = "#FF7F0E",
  "Spermatocyte"        = "#1F77B4",
  "Round spermatid"     = "#2CA02C",
  "Elongated spermatid" = "#9467BD"
)

# ── 主 DotPlot ──────────────────────────────────────────────
p_dot <- ggplot(dp_st,
                aes(x=features.plot, y=stage, size=pct.exp, color=exp_norm)) +
  geom_vline(xintercept=sep_x, color="grey72", linewidth=0.4, linetype="dashed") +
  geom_point() +
  scale_color_gradientn(
    colours = c("grey92","grey82","#AED6F1","#2E86C1","#1A5276"),
    values  = c(0, 0.2, 0.5, 0.75, 1),
    limits  = c(0, 1),
    name    = "Relative\nExpression\n(per gene)"
  ) +
  scale_size_continuous(
    range  = c(0.5, 9),
    name   = "% Expressing",
    breaks = c(10, 25, 50, 75, 100)
  ) +
  scale_y_discrete(position="right") +
  labs(x=NULL, y=NULL, title="Exos v2 — Germ cells by developmental stage") +
  theme_classic(base_size=12) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, face="italic", size=9.5),
    axis.text.y = element_text(size=11, face="bold"),
    axis.line   = element_blank(),
    panel.border= element_rect(color="grey40", fill=NA, linewidth=0.7),
    panel.grid.major.y = element_line(color="grey93", linewidth=0.3),
    plot.title  = element_text(hjust=0.5, face="bold", size=13),
    plot.background = element_rect(fill="white", color=NA),
    legend.title = element_text(size=9),
    legend.text  = element_text(size=8.5)
  )

# ── 顶部色条（stage label bar）──────────────────────────────
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
            size=2.6, fontface="bold", color="white") +
  scale_fill_manual(values=stage_colors[levels(sec_df$grp)]) +
  scale_x_continuous(expand=c(0, 0.5)) +
  theme_void() +
  theme(legend.position="none",
        plot.background=element_rect(fill="white", color=NA))

# ── 合并 ────────────────────────────────────────────────────
p_final <- p_bar / p_dot +
  plot_layout(heights=c(1, 10)) +
  plot_annotation(
    theme=theme(plot.background=element_rect(fill="white", color=NA))
  )

# ══════════════════════════════════════════════════════════════
# 4. 保存
# ══════════════════════════════════════════════════════════════
n_genes <- length(gene_order)
w <- max(15, n_genes * 0.32 + 4)

ggsave(file.path(A_DIR,"exos_v2_germcell_dotplot_v2.pdf"),
       p_final, width=w, height=7.5)
ggsave(file.path(A_DIR,"exos_v2_germcell_dotplot_v2.png"),
       p_final, width=w, height=7.5, dpi=220)
cat(sprintf("\nSaved: exos_v2_germcell_dotplot_v2.png  (%.1f × 7.5 in)\n", w))

# ── UMAP with 6 stages ──────────────────────────────────────
p_umap <- DimPlot(germ, reduction="umap", group.by="stage_label",
                  cols=stage_colors, pt.size=0.2) +
  labs(title="Exos v2 — Germ cells (6 stages)") +
  theme_cowplot(font_size=12) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white")) +
  guides(color=guide_legend(override.aes=list(size=4)))

p_umap_grp <- DimPlot(germ, reduction="umap", group.by="group",
                       cols=c("Rescue"="#4DBBD5","Aging"="#E64B35"),
                       pt.size=0.2) +
  labs(title="Aging vs Rescue") +
  theme_cowplot(font_size=12) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white"))

p_combined <- p_umap | p_umap_grp
ggsave(file.path(A_DIR,"exos_v2_germcell_umap_6stage.pdf"),
       p_combined, width=15, height=6.5)
ggsave(file.path(A_DIR,"exos_v2_germcell_umap_6stage.png"),
       p_combined, width=15, height=6.5, dpi=200)
cat("Saved: exos_v2_germcell_umap_6stage.png\n")

# 保存注释表
write.csv(
  data.frame(cluster=names(germ_stage), stage=germ_stage),
  file.path(A_DIR,"germ_stage_annotation_v2.csv"), row.names=FALSE
)

cat("\n===== DONE =====\n")
cat(sprintf("Stage distribution:\n"))
print(table(germ$stage_label))

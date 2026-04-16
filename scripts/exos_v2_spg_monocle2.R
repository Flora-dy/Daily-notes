#!/usr/bin/env Rscript
# ============================================================
# Exos v2 — Spermatogonia Monocle 2 Pseudotime (Fig 5F/G/H style)
# ProSPG → Undiff SPG → Diff SPG
# DDRTree trajectory with pseudotime gradient + subtype + group panels
# ============================================================

.libPaths(c("/Users/zzp/Code/tianhairui/.rlib", .libPaths()))
suppressPackageStartupMessages({
  library(monocle)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(cowplot)
  library(viridis)
})

EXOS_OUT <- "/Users/zzp/Code/tianhairui/outputs0409/seurat5_exos_v2"
OUT_DIR  <- file.path(EXOS_OUT, "05_pseudotime")
RDS_DIR  <- file.path(EXOS_OUT, "rds")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

SEED <- 202409; set.seed(SEED)

GROUP_COLS <- c("Rescue"="#4DBBD5", "Aging"="#E64B35")
spg_colors <- c(
  "Prospermatogonium"               = "#E8735A",   # salmon-red (ProSG)
  "Undifferentiated Spermatogonium" = "#9B59B6",   # purple (Undiff SG)
  "Differentiating Spermatogonium"  = "#3498DB"    # teal-blue (Diff SG)
)

# ══════════════════════════════════════════════════════════════
# 1. 加载 SPG 子对象
# ══════════════════════════════════════════════════════════════
cat("=== Loading SPG sub-object ===\n")
spg <- readRDS(file.path(RDS_DIR, "seurat_exos_spg_subclustered.rds"))
cat(sprintf("SPG cells: %d\n", ncol(spg)))
print(table(spg$spg_subtype))

# ══════════════════════════════════════════════════════════════
# 2. 构建 Monocle 2 CellDataSet
# ══════════════════════════════════════════════════════════════
cat("\n=== Building Monocle2 CellDataSet ===\n")

# 使用 log-normalized 矩阵 + gaussianff() 避免 estimateDispersions 中的
# dplyr::group_by_() 废弃函数问题
expr_mat <- GetAssayData(spg, assay="RNA", layer="data")   # log-normalized
expr_mat <- as(expr_mat, "sparseMatrix")

# Feature data
gene_df <- data.frame(
  gene_short_name = rownames(expr_mat),
  row.names       = rownames(expr_mat)
)

# Cell data (phenotype)
cell_df <- spg@meta.data[, c("orig.ident","group","spg_subtype"), drop=FALSE]
colnames(cell_df)[1] <- "sample_id"

# 创建 CDS (gaussianff + log-normalized data, 跳过 estimateDispersions)
fd  <- new("AnnotatedDataFrame", data=gene_df)
pd  <- new("AnnotatedDataFrame", data=cell_df)
cds <- newCellDataSet(
  expr_mat,
  phenoData            = pd,
  featureData          = fd,
  expressionFamily     = gaussianff(),
  lowerDetectionLimit  = 0.1
)

cds <- estimateSizeFactors(cds)
# estimateDispersions 仅 negbinomial 需要，gaussianff 跳过
cat("CDS built.\n")

# ══════════════════════════════════════════════════════════════
# 3. 选择 ordering genes（差异基因 + 已知 SPG marker）
# ══════════════════════════════════════════════════════════════
cat("\n=== Selecting ordering genes ===\n")

# 方法1: differentialGeneTest（按亚群）
Sys.time() -> t0
diff_test <- differentialGeneTest(
  cds,
  fullModelFormulaStr = "~spg_subtype",
  cores = 1
)
cat(sprintf("differentialGeneTest done in %.1f s\n", as.numeric(Sys.time()-t0)))

# 取 q<0.01 的 top 基因
ordering_genes <- row.names(diff_test)[diff_test$qval < 0.01]
cat(sprintf("Ordering genes (q<0.01): %d\n", length(ordering_genes)))

# 补充已知 SPG marker（确保关键基因在内）
spg_markers <- c("Id2","Id4","Klf4","Gfra1","Ret","Thy1",
                 "Zbtb16","Uchl1","Sall4","Cdh1","Nanos3","Etv5","Egr4",
                 "Kit","Stra8","Dmrtb1","Sohlh1","Dnmt3b","Prdm9")
spg_markers <- spg_markers[spg_markers %in% rownames(cds)]
ordering_genes <- unique(c(ordering_genes, spg_markers))
cat(sprintf("Total ordering genes: %d\n", length(ordering_genes)))

cds <- setOrderingFilter(cds, ordering_genes)

# ══════════════════════════════════════════════════════════════
# 4. 降维 + 排序（DDRTree）
# ══════════════════════════════════════════════════════════════
cat("\n=== DDRTree dimensionality reduction ===\n")
cds <- reduceDimension(
  cds,
  max_components = 2,
  method         = "DDRTree",
  norm_method    = "none"    # required for gaussianff/uninormal family
)

cat("Ordering cells...\n")
cds <- orderCells(cds)

# 根据 ProSPG 细胞的拟时序值决定是否需要翻转方向
# ProSPG 应该是起点（最小拟时序值）
prospe_pt <- pData(cds)$Pseudotime[pData(cds)$spg_subtype == "Prospermatogonium"]
diff_pt   <- pData(cds)$Pseudotime[pData(cds)$spg_subtype == "Differentiating Spermatogonium"]

if (median(prospe_pt, na.rm=TRUE) > median(diff_pt, na.rm=TRUE)) {
  cat("Reversing pseudotime direction (ProSPG should be start)...\n")
  cds <- orderCells(cds, reverse=TRUE)
}

cat("Pseudotime range:", range(pData(cds)$Pseudotime), "\n")
cat("\nMean pseudotime per subtype:\n")
print(tapply(pData(cds)$Pseudotime, pData(cds)$spg_subtype, mean))

# 保存 CDS
saveRDS(cds, file.path(RDS_DIR, "monocle2_spg_cds.rds"))
cat("Saved: monocle2_spg_cds.rds\n")

# ══════════════════════════════════════════════════════════════
# 5. 绘图 (Fig 5F / G / H 风格)
# ══════════════════════════════════════════════════════════════
cat("\n=== Plotting ===\n")

# 公共主题
traj_theme <- theme(
  panel.background  = element_rect(fill="white", color=NA),
  plot.background   = element_rect(fill="white", color=NA),
  panel.grid        = element_blank(),
  axis.line         = element_line(color="grey30", linewidth=0.5),
  axis.text         = element_text(size=10),
  axis.title        = element_text(size=11),
  plot.title        = element_text(hjust=0.5, face="bold", size=12),
  legend.background = element_rect(fill="white", color=NA),
  legend.key        = element_rect(fill="white", color=NA),
  legend.title      = element_text(size=9),
  legend.text       = element_text(size=9)
)

# ── Fig 5F style: pseudotime gradient ────────────────────────
p_F <- plot_cell_trajectory(
  cds,
  color_by   = "Pseudotime",
  cell_size  = 0.8,
  show_branch_points = TRUE
) +
  scale_color_viridis_c(
    option = "D",          # dark purple → yellow (similar to paper)
    name   = "Pseudotime",
    breaks = c(0, 20, 40)
  ) +
  labs(title="Pseudotime", x="Component 1", y="Component 2") +
  traj_theme

# ── Fig 5G style: subtype colors ─────────────────────────────
p_G <- plot_cell_trajectory(
  cds,
  color_by   = "spg_subtype",
  cell_size  = 0.8,
  show_branch_points = TRUE
) +
  scale_color_manual(
    values = spg_colors,
    name   = NULL,
    guide  = guide_legend(override.aes=list(size=3))
  ) +
  labs(title="Subtype", x="Component 1", y="Component 2") +
  traj_theme +
  theme(legend.position = c(0.02, 0.15),
        legend.justification = c(0, 0))

# ── Fig 5H style: Aging vs Rescue ────────────────────────────
p_H <- plot_cell_trajectory(
  cds,
  color_by   = "group",
  cell_size  = 0.8,
  show_branch_points = TRUE
) +
  scale_color_manual(
    values = GROUP_COLS,
    name   = NULL,
    guide  = guide_legend(override.aes=list(size=3))
  ) +
  labs(title="Group", x="Component 1", y="Component 2") +
  traj_theme +
  theme(legend.position = c(0.02, 0.85),
        legend.justification = c(0, 1))

# ── 三图合并 (Fig 5F+G+H) ─────────────────────────────────────
p_FGH <- p_F | p_G | p_H
ggsave(file.path(OUT_DIR,"spg_monocle2_FGH.pdf"), p_FGH, width=18, height=6)
ggsave(file.path(OUT_DIR,"spg_monocle2_FGH.png"), p_FGH, width=18, height=6, dpi=220)
cat("Saved: spg_monocle2_FGH\n")

# 单图单独保存
ggsave(file.path(OUT_DIR,"spg_monocle2_F_pseudotime.pdf"), p_F, width=6, height=6)
ggsave(file.path(OUT_DIR,"spg_monocle2_F_pseudotime.png"), p_F, width=6, height=6, dpi=220)
ggsave(file.path(OUT_DIR,"spg_monocle2_G_subtype.pdf"),    p_G, width=6, height=6)
ggsave(file.path(OUT_DIR,"spg_monocle2_G_subtype.png"),    p_G, width=6, height=6, dpi=220)
ggsave(file.path(OUT_DIR,"spg_monocle2_H_group.pdf"),      p_H, width=6, height=6)
ggsave(file.path(OUT_DIR,"spg_monocle2_H_group.png"),      p_H, width=6, height=6, dpi=220)
cat("Saved individual panels.\n")

# ══════════════════════════════════════════════════════════════
# 6. 拟时序统计（与 Slingshot 版本一致，raw p-value）
# ══════════════════════════════════════════════════════════════
cat("\n=== Pseudotime stats (raw Wilcoxon, no correction) ===\n")
pt_df <- pData(cds) %>%
  as.data.frame() %>%
  mutate(spg_subtype=factor(spg_subtype, levels=names(spg_colors)))

sig_label <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "ns",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ sprintf("p=%.3f", p)
  )
}

stats_df <- pt_df %>%
  group_by(spg_subtype) %>%
  summarise(
    n_Aging    = sum(group=="Aging"),
    n_Rescue   = sum(group=="Rescue"),
    mean_Aging = mean(Pseudotime[group=="Aging"],  na.rm=TRUE),
    mean_Rescue= mean(Pseudotime[group=="Rescue"], na.rm=TRUE),
    delta_mean = mean_Aging - mean_Rescue,
    p_wilcox   = tryCatch(
      wilcox.test(Pseudotime[group=="Aging"],
                  Pseudotime[group=="Rescue"], exact=FALSE)$p.value,
      error=function(e) NA_real_),
    .groups="drop"
  ) %>%
  mutate(sig=sig_label(p_wilcox))

print(as.data.frame(stats_df))
write.csv(stats_df, file.path(OUT_DIR,"spg_monocle2_stats.csv"), row.names=FALSE)

# ── 箱线图 + p 值标注 ─────────────────────────────────────────
y_pos <- pt_df %>%
  group_by(spg_subtype) %>%
  summarise(y_max=quantile(Pseudotime, 0.97, na.rm=TRUE)*1.07, .groups="drop")
annot <- left_join(stats_df, y_pos, by="spg_subtype")

p_box <- ggplot(pt_df, aes(x=group, y=Pseudotime, fill=group)) +
  geom_boxplot(outlier.size=0.4, outlier.alpha=0.4, alpha=0.75, width=0.55) +
  geom_jitter(width=0.18, size=0.3, alpha=0.2) +
  geom_segment(data=annot,
               aes(x=1, xend=2, y=y_max, yend=y_max),
               inherit.aes=FALSE, linewidth=0.5, color="grey30") +
  geom_text(data=annot,
            aes(x=1.5, y=y_max*1.05, label=sig),
            inherit.aes=FALSE, size=4.5, fontface="bold", color="grey20") +
  facet_wrap(~spg_subtype, scales="free_y", ncol=3) +
  scale_fill_manual(values=GROUP_COLS) +
  theme_cowplot(font_size=11) +
  theme(axis.text.x=element_text(angle=30, hjust=1),
        legend.position="none",
        strip.text=element_text(face="bold", size=10),
        plot.background=element_rect(fill="white", color=NA)) +
  labs(title="SPG Pseudotime (Monocle2): Aging vs Rescue  [raw p, no correction]",
       x=NULL, y="Pseudotime")

ggsave(file.path(OUT_DIR,"spg_monocle2_boxplot.pdf"), p_box, width=12, height=5)
ggsave(file.path(OUT_DIR,"spg_monocle2_boxplot.png"), p_box, width=12, height=5, dpi=220)
cat("Saved: spg_monocle2_boxplot\n")

# ══════════════════════════════════════════════════════════════
# 7. 细胞比例柱状图（Fig 5I style）
# ══════════════════════════════════════════════════════════════
prop_df <- pt_df %>%
  group_by(group, spg_subtype) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(group) %>%
  mutate(pct=n/sum(n)*100) %>%
  ungroup()

p_prop <- ggplot(prop_df, aes(x=group, y=pct, fill=spg_subtype)) +
  geom_col(width=0.6) +
  scale_fill_manual(values=spg_colors, name=NULL) +
  scale_y_continuous(expand=c(0,0), limits=c(0,105)) +
  theme_cowplot(font_size=12) +
  theme(plot.background=element_rect(fill="white", color=NA),
        legend.position="right") +
  labs(title="Subtype Proportion", x=NULL, y="Fraction of cells (%)")

ggsave(file.path(OUT_DIR,"spg_monocle2_proportion.pdf"), p_prop, width=5, height=5)
ggsave(file.path(OUT_DIR,"spg_monocle2_proportion.png"), p_prop, width=5, height=5, dpi=220)
cat("Saved: spg_monocle2_proportion\n")

cat("\n===== DONE =====\n")
cat("Output:", OUT_DIR, "\n")
cat("\nFinal stats:\n")
print(stats_df[,c("spg_subtype","n_Aging","n_Rescue","mean_Aging","mean_Rescue","delta_mean","p_wilcox","sig")])

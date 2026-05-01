#!/usr/bin/env Rscript
# ============================================================
# Exos v2 — Spermatogonia Monocle 2 Pseudotime (Fig 5F/G/H style)
# ProSPG → Undiff SPG → Diff SPG
# DDRTree trajectory; ggplot2 manual plotting (避免 orderCells/
# plot_cell_trajectory 的 igraph 2.x 兼容性问题)
# Pseudotime coloring: Slingshot pseudotime (already computed)
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
  library(slingshot)
  library(SingleCellExperiment)
  library(igraph)
})

EXOS_OUT <- "/Users/zzp/Code/tianhairui/outputs0409/seurat5_exos_v2"
OUT_DIR  <- file.path(EXOS_OUT, "05_pseudotime")
RDS_DIR  <- file.path(EXOS_OUT, "rds")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

SEED <- 202409; set.seed(SEED)

GROUP_COLS <- c("Rescue"="#4DBBD5", "Aging"="#E64B35")
spg_colors <- c(
  "Prospermatogonium"               = "#E8735A",
  "Undifferentiated Spermatogonium" = "#9B59B6",
  "Differentiating Spermatogonium"  = "#3498DB"
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
# (gaussianff + log-normalized; 跳过 estimateDispersions)
# ══════════════════════════════════════════════════════════════
cat("\n=== Building Monocle2 CellDataSet ===\n")

expr_mat <- GetAssayData(spg, assay="RNA", layer="data")
expr_mat <- as(expr_mat, "sparseMatrix")

gene_df <- data.frame(gene_short_name=rownames(expr_mat),
                      row.names=rownames(expr_mat))
cell_df <- spg@meta.data[, c("orig.ident","group","spg_subtype"), drop=FALSE]
colnames(cell_df)[1] <- "sample_id"

fd  <- new("AnnotatedDataFrame", data=gene_df)
pd  <- new("AnnotatedDataFrame", data=cell_df)
cds <- newCellDataSet(
  expr_mat,
  phenoData           = pd,
  featureData         = fd,
  expressionFamily    = gaussianff(),
  lowerDetectionLimit = 0.1
)
cds <- estimateSizeFactors(cds)
cat("CDS built.\n")

# ══════════════════════════════════════════════════════════════
# 3. 选择 ordering genes（Seurat FindAllMarkers + known SPG markers）
# ══════════════════════════════════════════════════════════════
cat("\n=== Selecting ordering genes (Seurat FindAllMarkers) ===\n")
t0 <- Sys.time()

Idents(spg) <- "spg_subtype"
markers_all <- FindAllMarkers(
  spg, assay="RNA", only.pos=TRUE,
  min.pct=0.1, logfc.threshold=0.25, test.use="wilcox"
)
cat(sprintf("FindAllMarkers done in %.1f s\n", as.numeric(Sys.time()-t0)))

ordering_genes <- markers_all %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n=400) %>%
  pull(gene) %>%
  unique()
cat(sprintf("Ordering genes from markers: %d\n", length(ordering_genes)))

spg_markers <- c("Id2","Id4","Klf4","Gfra1","Ret","Thy1",
                 "Zbtb16","Uchl1","Sall4","Cdh1","Nanos3","Etv5","Egr4",
                 "Kit","Stra8","Dmrtb1","Sohlh1","Dnmt3b","Prdm9")
spg_markers <- spg_markers[spg_markers %in% featureNames(cds)]
ordering_genes <- unique(c(ordering_genes, spg_markers))
ordering_genes <- ordering_genes[ordering_genes %in% featureNames(cds)]
cat(sprintf("Total ordering genes: %d\n", length(ordering_genes)))

cds <- setOrderingFilter(cds, ordering_genes)

# ══════════════════════════════════════════════════════════════
# 4. DDRTree 降维
# 注：orderCells() 和 plot_cell_trajectory() 依赖已废弃的
# igraph API (nei/neimode)，跳过这两个函数，改用手动坐标提取
# ══════════════════════════════════════════════════════════════
cat("\n=== DDRTree dimensionality reduction ===\n")
cds <- reduceDimension(
  cds,
  max_components = 2,
  method         = "DDRTree",
  norm_method    = "none"
)
cat("DDRTree done.\n")

# 提取细胞在 DDRTree 空间的坐标 (2 x N_cells → N_cells x 2)
cell_pos <- t(monocle::reducedDimS(cds))
colnames(cell_pos) <- c("Component1", "Component2")

# 提取 backbone 节点坐标 (2 x K → K x 2)
tree_nodes <- t(monocle::reducedDimK(cds))
colnames(tree_nodes) <- c("C1", "C2")

# 从 minSpanningTree 提取骨架边 (igraph 现代 API)
mst <- monocle::minSpanningTree(cds)
el  <- igraph::as_edgelist(mst, names=FALSE)
edges_df <- data.frame(
  x    = tree_nodes[el[,1], "C1"],
  y    = tree_nodes[el[,1], "C2"],
  xend = tree_nodes[el[,2], "C1"],
  yend = tree_nodes[el[,2], "C2"]
)

# 保存 CDS（含 DDRTree 降维结果）
saveRDS(cds, file.path(RDS_DIR, "monocle2_spg_cds.rds"))
cat("Saved: monocle2_spg_cds.rds\n")

# ══════════════════════════════════════════════════════════════
# 5. 加载 Slingshot 拟时序（用于 Fig 5F 着色）
# ══════════════════════════════════════════════════════════════
cat("\n=== Loading Slingshot pseudotime for Fig 5F coloring ===\n")
sds    <- readRDS(file.path(RDS_DIR, "slingshot_spg_v2.rds"))
pt_mat <- slingPseudotime(sds)
pt_vec <- rowMeans(pt_mat, na.rm=TRUE)
cat(sprintf("Slingshot pseudotime range: [%.2f, %.2f]\n",
            min(pt_vec, na.rm=TRUE), max(pt_vec, na.rm=TRUE)))

# 构建绘图数据框（按细胞顺序对齐）
cell_names <- rownames(cell_pos)
plot_df <- data.frame(
  x           = cell_pos[, "Component1"],
  y           = cell_pos[, "Component2"],
  spg_subtype = pData(cds)$spg_subtype,
  group       = pData(cds)$group,
  pseudotime  = pt_vec[cell_names],
  row.names   = cell_names
)
plot_df$spg_subtype <- factor(plot_df$spg_subtype, levels=names(spg_colors))

cat("\nMean pseudotime per subtype:\n")
print(tapply(plot_df$pseudotime, plot_df$spg_subtype, mean, na.rm=TRUE))

# ══════════════════════════════════════════════════════════════
# 6. 绘图 (Fig 5F / G / H 风格) — 纯 ggplot2
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

# 骨架层（共用）
backbone_layer <- geom_segment(
  data = edges_df,
  aes(x=x, y=y, xend=xend, yend=yend),
  inherit.aes = FALSE,
  color = "grey35", linewidth = 1.0, alpha = 0.9
)

# ── Fig 5F: pseudotime gradient ──────────────────────────────
p_F <- ggplot(plot_df, aes(x=x, y=y, color=pseudotime)) +
  geom_point(size=0.6, alpha=0.7) +
  backbone_layer +
  scale_color_viridis_c(option="D", name="Pseudotime",
                        na.value="grey80") +
  labs(title="Pseudotime", x="Component 1", y="Component 2") +
  traj_theme

# ── Fig 5G: subtype colors ────────────────────────────────────
p_G <- ggplot(plot_df, aes(x=x, y=y, color=spg_subtype)) +
  geom_point(size=0.6, alpha=0.7) +
  backbone_layer +
  scale_color_manual(values=spg_colors, name=NULL,
                     guide=guide_legend(override.aes=list(size=3))) +
  labs(title="Subtype", x="Component 1", y="Component 2") +
  traj_theme +
  theme(legend.position=c(0.02, 0.15), legend.justification=c(0,0))

# ── Fig 5H: Aging vs Rescue ───────────────────────────────────
p_H <- ggplot(plot_df, aes(x=x, y=y, color=group)) +
  geom_point(size=0.6, alpha=0.7) +
  backbone_layer +
  scale_color_manual(values=GROUP_COLS, name=NULL,
                     guide=guide_legend(override.aes=list(size=3))) +
  labs(title="Group", x="Component 1", y="Component 2") +
  traj_theme +
  theme(legend.position=c(0.02, 0.85), legend.justification=c(0,1))

# ── 三图合并 ──────────────────────────────────────────────────
p_FGH <- p_F | p_G | p_H
ggsave(file.path(OUT_DIR,"spg_monocle2_FGH.pdf"), p_FGH, width=18, height=6)
ggsave(file.path(OUT_DIR,"spg_monocle2_FGH.png"), p_FGH, width=18, height=6, dpi=220)
cat("Saved: spg_monocle2_FGH\n")

ggsave(file.path(OUT_DIR,"spg_monocle2_F_pseudotime.pdf"), p_F, width=6, height=6)
ggsave(file.path(OUT_DIR,"spg_monocle2_F_pseudotime.png"), p_F, width=6, height=6, dpi=220)
ggsave(file.path(OUT_DIR,"spg_monocle2_G_subtype.pdf"),    p_G, width=6, height=6)
ggsave(file.path(OUT_DIR,"spg_monocle2_G_subtype.png"),    p_G, width=6, height=6, dpi=220)
ggsave(file.path(OUT_DIR,"spg_monocle2_H_group.pdf"),      p_H, width=6, height=6)
ggsave(file.path(OUT_DIR,"spg_monocle2_H_group.png"),      p_H, width=6, height=6, dpi=220)
cat("Saved: individual panels F/G/H\n")

# ══════════════════════════════════════════════════════════════
# 7. 拟时序统计（Slingshot pseudotime，raw Wilcoxon，no correction）
# ══════════════════════════════════════════════════════════════
cat("\n=== Pseudotime stats (Slingshot, raw Wilcoxon, no correction) ===\n")

sig_label <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "ns",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ sprintf("p=%.3f", p)
  )
}

pt_valid <- plot_df %>% filter(!is.na(pseudotime))

stats_df <- pt_valid %>%
  group_by(spg_subtype) %>%
  summarise(
    n_Aging     = sum(group=="Aging"),
    n_Rescue    = sum(group=="Rescue"),
    mean_Aging  = mean(pseudotime[group=="Aging"],  na.rm=TRUE),
    mean_Rescue = mean(pseudotime[group=="Rescue"], na.rm=TRUE),
    delta_mean  = mean_Aging - mean_Rescue,
    p_wilcox    = tryCatch(
      wilcox.test(pseudotime[group=="Aging"],
                  pseudotime[group=="Rescue"], exact=FALSE)$p.value,
      error=function(e) NA_real_),
    .groups="drop"
  ) %>%
  mutate(sig=sig_label(p_wilcox))

print(as.data.frame(stats_df))
write.csv(stats_df, file.path(OUT_DIR,"spg_monocle2_stats.csv"), row.names=FALSE)

# ── 箱线图 + p 值标注 ─────────────────────────────────────────
y_pos <- pt_valid %>%
  group_by(spg_subtype) %>%
  summarise(y_max=quantile(pseudotime, 0.97, na.rm=TRUE)*1.07, .groups="drop")
annot <- left_join(stats_df, y_pos, by="spg_subtype")

p_box <- ggplot(pt_valid, aes(x=group, y=pseudotime, fill=group)) +
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
  labs(title="SPG Pseudotime (DDRTree layout, Slingshot PT): Aging vs Rescue  [raw p]",
       x=NULL, y="Pseudotime")

ggsave(file.path(OUT_DIR,"spg_monocle2_boxplot.pdf"), p_box, width=12, height=5)
ggsave(file.path(OUT_DIR,"spg_monocle2_boxplot.png"), p_box, width=12, height=5, dpi=220)
cat("Saved: spg_monocle2_boxplot\n")

# ══════════════════════════════════════════════════════════════
# 8. 细胞比例柱状图（Fig 5I style）
# ══════════════════════════════════════════════════════════════
prop_df <- plot_df %>%
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
print(stats_df[,c("spg_subtype","n_Aging","n_Rescue",
                  "mean_Aging","mean_Rescue","delta_mean","p_wilcox","sig")])

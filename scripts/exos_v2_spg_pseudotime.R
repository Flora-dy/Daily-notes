#!/usr/bin/env Rscript
# ============================================================
# Exos v2 — Spermatogonia Pseudotime (3 subtypes only)
# ProSPG → Undifferentiated SPG → Differentiating SPG
# Stats: raw Wilcoxon p-value, no multiple-testing adjustment
# ============================================================

.libPaths(c("/Users/zzp/Code/tianhairui/.rlib", .libPaths()))
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(cowplot)
  library(viridis)
  library(slingshot)
  library(SingleCellExperiment)
})

EXOS_OUT  <- "/Users/zzp/Code/tianhairui/outputs0409/seurat5_exos_v2"
OUT_DIR   <- file.path(EXOS_OUT, "05_pseudotime")
RDS_DIR   <- file.path(EXOS_OUT, "rds")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

SEED <- 202409; set.seed(SEED)
GROUP_COLS <- c("Rescue"="#4DBBD5", "Aging"="#E64B35")

spg_colors <- c(
  "Prospermatogonium"               = "#7B2D8B",
  "Undifferentiated Spermatogonium" = "#D62728",
  "Differentiating Spermatogonium"  = "#FF7F0E"
)
spg_order <- names(spg_colors)

# ══════════════════════════════════════════════════════════════
# 1. 加载 SPG 子对象
# ══════════════════════════════════════════════════════════════
cat("=== Loading SPG sub-object ===\n")
spg <- readRDS(file.path(RDS_DIR, "seurat_exos_spg_subclustered.rds"))
cat(sprintf("SPG cells: %d\n", ncol(spg)))
cat("Subtype distribution:\n"); print(table(spg$spg_subtype))
cat("Group distribution:\n");   print(table(spg$group))

# ══════════════════════════════════════════════════════════════
# 2. Slingshot pseudotime
# ══════════════════════════════════════════════════════════════
cat("\n=== Running Slingshot on 3 SPG subtypes ===\n")

umap_emb    <- Embeddings(spg, "umap")
stage_labels <- as.character(spg$spg_subtype)

cat(sprintf("Start cluster: %s\n", spg_order[1]))
cat(sprintf("End cluster:   %s\n", spg_order[3]))

sce <- SingleCellExperiment(
  assays  = list(counts = matrix(0L, nrow=1, ncol=ncol(spg))),
  colData = data.frame(
    spg_subtype = stage_labels,
    group       = as.character(spg$group),
    row.names   = colnames(spg)
  )
)
reducedDim(sce, "UMAP") <- umap_emb

sds <- tryCatch(
  slingshot(sce,
            clusterLabels = "spg_subtype",
            reducedDim    = "UMAP",
            start.clus    = "Prospermatogonium",
            end.clus      = "Differentiating Spermatogonium"),
  error = function(e) {
    cat("Slingshot error:", conditionMessage(e), "\n")
    NULL
  }
)

if (is.null(sds)) stop("Slingshot failed.")

pt_mat  <- slingPseudotime(sds)
pt_mean <- rowMeans(pt_mat, na.rm=TRUE)
cat(sprintf("Pseudotime range: [%.2f, %.2f]\n", min(pt_mean, na.rm=TRUE),
            max(pt_mean, na.rm=TRUE)))

pt_df <- data.frame(
  cell        = colnames(spg),
  sample_id   = as.character(spg$orig.ident),
  group       = as.character(spg$group),
  spg_subtype = stage_labels,
  pseudotime  = pt_mean,
  stringsAsFactors = FALSE
)
pt_df$spg_subtype <- factor(pt_df$spg_subtype, levels=spg_order)

write.csv(pt_df, file.path(OUT_DIR, "spg_pseudotime_cells.csv"), row.names=FALSE)
cat("Saved: spg_pseudotime_cells.csv\n")

# ── Stats: raw Wilcoxon, no p-adjustment ─────────────────────
cat("\n=== Aging vs Rescue comparison (raw p-value, no correction) ===\n")
pt_stats <- pt_df %>%
  filter(!is.na(pseudotime)) %>%
  group_by(spg_subtype) %>%
  summarise(
    n_Aging    = sum(group == "Aging"),
    n_Rescue   = sum(group == "Rescue"),
    mean_Aging = mean(pseudotime[group == "Aging"],  na.rm=TRUE),
    mean_Rescue= mean(pseudotime[group == "Rescue"], na.rm=TRUE),
    median_Aging = median(pseudotime[group == "Aging"],  na.rm=TRUE),
    median_Rescue= median(pseudotime[group == "Rescue"], na.rm=TRUE),
    delta_mean = mean_Aging - mean_Rescue,
    p_wilcox   = tryCatch(
      wilcox.test(pseudotime[group == "Aging"],
                  pseudotime[group == "Rescue"],
                  exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

# 不做多重校正，直接展示原始 p 值
cat("\nSPG pseudotime stats (raw p-value, no BH correction):\n")
print(as.data.frame(pt_stats))
write.csv(pt_stats, file.path(OUT_DIR, "spg_pseudotime_stats.csv"), row.names=FALSE)

# 显著性标签（仅用原始 p 值）
sig_label <- function(p) {
  case_when(
    is.na(p)   ~ "ns",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ sprintf("p=%.3f", p)
  )
}
pt_stats$sig <- sig_label(pt_stats$p_wilcox)

# ══════════════════════════════════════════════════════════════
# 3. UMAP 图
# ══════════════════════════════════════════════════════════════
spg$pseudotime <- pt_df$pseudotime[match(colnames(spg), pt_df$cell)]

p_pt_umap <- FeaturePlot(spg, features="pseudotime",
                         min.cutoff="q5", max.cutoff="q95",
                         pt.size=0.6) +
  scale_color_viridis_c(option="plasma", name="Pseudotime") +
  labs(title="Spermatogonia Pseudotime") +
  theme_cowplot(font_size=11) +
  theme(plot.title=element_text(hjust=0.5, face="bold"),
        plot.background=element_rect(fill="white", color=NA)) +
  NoAxes()

p_subtype_umap <- DimPlot(spg, reduction="umap", group.by="spg_subtype",
                           cols=spg_colors, pt.size=0.6, label=TRUE,
                           label.size=4, repel=TRUE) +
  labs(title="SPG Subtypes") +
  theme_cowplot(font_size=11) +
  theme(plot.title=element_text(hjust=0.5, face="bold"),
        plot.background=element_rect(fill="white")) +
  guides(color=guide_legend(override.aes=list(size=4)))

p_umap_row <- p_subtype_umap | p_pt_umap
ggsave(file.path(OUT_DIR,"spg_pseudotime_umap.pdf"), p_umap_row, width=14, height=6)
ggsave(file.path(OUT_DIR,"spg_pseudotime_umap.png"), p_umap_row, width=14, height=6, dpi=200)
cat("Saved: spg_pseudotime_umap\n")

# ══════════════════════════════════════════════════════════════
# 4. 密度图：Aging vs Rescue（整体 + 分亚群）
# ══════════════════════════════════════════════════════════════
pt_valid <- pt_df[!is.na(pt_df$pseudotime), ]

# 整体密度
p_dens_all <- ggplot(pt_valid, aes(x=pseudotime, fill=group, color=group)) +
  geom_density(alpha=0.4, linewidth=0.9) +
  scale_fill_manual(values=GROUP_COLS) +
  scale_color_manual(values=GROUP_COLS) +
  theme_cowplot(font_size=12) +
  theme(plot.background=element_rect(fill="white", color=NA),
        legend.position=c(0.75, 0.85)) +
  labs(title="SPG Pseudotime Distribution (All subtypes)",
       x="Pseudotime", y="Density", fill=NULL, color=NULL)

# 分亚群密度
p_dens_facet <- ggplot(pt_valid, aes(x=pseudotime, fill=group, color=group)) +
  geom_density(alpha=0.4, linewidth=0.8) +
  facet_wrap(~spg_subtype, scales="free_y", ncol=1) +
  scale_fill_manual(values=GROUP_COLS) +
  scale_color_manual(values=GROUP_COLS) +
  theme_cowplot(font_size=11) +
  theme(plot.background=element_rect(fill="white", color=NA),
        strip.text=element_text(face="bold", size=10)) +
  labs(title="SPG Pseudotime by Subtype",
       x="Pseudotime", y="Density", fill=NULL, color=NULL)

p_dens_combined <- p_dens_all | p_dens_facet
ggsave(file.path(OUT_DIR,"spg_pseudotime_density.pdf"), p_dens_combined, width=14, height=7)
ggsave(file.path(OUT_DIR,"spg_pseudotime_density.png"), p_dens_combined, width=14, height=7, dpi=200)
cat("Saved: spg_pseudotime_density\n")

# ══════════════════════════════════════════════════════════════
# 5. 箱线图 + 原始 p 值标注
# ══════════════════════════════════════════════════════════════
# 计算每组每亚群的 y 位置（用于标注）
y_pos <- pt_valid %>%
  group_by(spg_subtype) %>%
  summarise(y_max = quantile(pseudotime, 0.97, na.rm=TRUE) * 1.05, .groups="drop")

annot_df <- left_join(pt_stats, y_pos, by="spg_subtype")

p_box <- ggplot(pt_valid, aes(x=group, y=pseudotime, fill=group)) +
  geom_boxplot(outlier.size=0.4, outlier.alpha=0.5, alpha=0.75, width=0.55) +
  geom_jitter(width=0.18, size=0.3, alpha=0.25) +
  # p 值标注
  geom_segment(data=annot_df,
               aes(x=1, xend=2, y=y_max, yend=y_max),
               inherit.aes=FALSE, linewidth=0.5, color="grey30") +
  geom_text(data=annot_df,
            aes(x=1.5, y=y_max * 1.04, label=sig),
            inherit.aes=FALSE, size=4.5, fontface="bold", color="grey20") +
  facet_wrap(~spg_subtype, scales="free_y", ncol=3) +
  scale_fill_manual(values=GROUP_COLS) +
  theme_cowplot(font_size=11) +
  theme(axis.text.x=element_text(angle=30, hjust=1),
        legend.position="none",
        strip.text=element_text(face="bold", size=10),
        plot.background=element_rect(fill="white", color=NA)) +
  labs(title="SPG Pseudotime: Aging vs Rescue  (raw Wilcoxon p, no correction)",
       x=NULL, y="Pseudotime")

ggsave(file.path(OUT_DIR,"spg_pseudotime_boxplot.pdf"), p_box, width=12, height=5)
ggsave(file.path(OUT_DIR,"spg_pseudotime_boxplot.png"), p_box, width=12, height=5, dpi=200)
cat("Saved: spg_pseudotime_boxplot\n")

# ══════════════════════════════════════════════════════════════
# 6. Violin + p 值（补充图）
# ══════════════════════════════════════════════════════════════
p_vln <- ggplot(pt_valid, aes(x=group, y=pseudotime, fill=group)) +
  geom_violin(trim=FALSE, alpha=0.7, scale="width") +
  geom_boxplot(width=0.12, outlier.shape=NA, fill="white", alpha=0.8) +
  geom_segment(data=annot_df,
               aes(x=1, xend=2, y=y_max, yend=y_max),
               inherit.aes=FALSE, linewidth=0.5, color="grey30") +
  geom_text(data=annot_df,
            aes(x=1.5, y=y_max * 1.04, label=sig),
            inherit.aes=FALSE, size=4.5, fontface="bold", color="grey20") +
  facet_wrap(~spg_subtype, scales="free_y", ncol=3) +
  scale_fill_manual(values=GROUP_COLS) +
  theme_cowplot(font_size=11) +
  theme(axis.text.x=element_text(angle=30, hjust=1),
        legend.position="none",
        strip.text=element_text(face="bold", size=10),
        plot.background=element_rect(fill="white", color=NA)) +
  labs(title="SPG Pseudotime: Aging vs Rescue  (raw Wilcoxon p, no correction)",
       x=NULL, y="Pseudotime")

ggsave(file.path(OUT_DIR,"spg_pseudotime_violin.pdf"), p_vln, width=12, height=5)
ggsave(file.path(OUT_DIR,"spg_pseudotime_violin.png"), p_vln, width=12, height=5, dpi=200)
cat("Saved: spg_pseudotime_violin\n")

# ══════════════════════════════════════════════════════════════
# 7. 保存 Slingshot 对象
# ══════════════════════════════════════════════════════════════
saveRDS(sds, file.path(RDS_DIR,"slingshot_spg_v2.rds"))
cat("Saved: slingshot_spg_v2.rds\n")

cat("\n===== DONE =====\n")
cat("Output:", OUT_DIR, "\n\n")
cat("Final stats summary:\n")
print(pt_stats[, c("spg_subtype","n_Aging","n_Rescue",
                   "mean_Aging","mean_Rescue","delta_mean","p_wilcox","sig")])

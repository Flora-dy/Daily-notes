#!/usr/bin/env Rscript
# ============================================================
# Exos v2 — Cell Composition + DEG + Pseudotime
# 输入: seurat_exos_clustered_v2.rds / seurat_exos_germcells_v2.rds
# 输出: 03_composition / 04_deg / 05_pseudotime
# ============================================================

.libPaths(c("/Users/zzp/Code/tianhairui/.rlib", .libPaths()))
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(cowplot)
  library(scales)
  library(ggrepel)
  library(viridis)
  library(slingshot)
  library(SingleCellExperiment)
})

EXOS_OUT  <- "/Users/zzp/Code/tianhairui/outputs0409/seurat5_exos_v2"
COMP_DIR  <- file.path(EXOS_OUT, "03_composition")
DEG_DIR   <- file.path(EXOS_OUT, "04_deg")
PSEUDO_DIR<- file.path(EXOS_OUT, "05_pseudotime")
RDS_DIR   <- file.path(EXOS_OUT, "rds")
for (d in c(COMP_DIR, DEG_DIR, PSEUDO_DIR)) dir.create(d, recursive=TRUE, showWarnings=FALSE)

SEED <- 202409; set.seed(SEED)
GROUP_COLS <- c("Rescue"="#4DBBD5", "Aging"="#E64B35")

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

stage_colors <- c(
  "ProSPG"              = "#7B2D8B",
  "Undiff SG"           = "#D62728",
  "Diff SG"             = "#FF7F0E",
  "Spermatocyte"        = "#1F77B4",
  "Round spermatid"     = "#2CA02C",
  "Elongated spermatid" = "#9467BD"
)

type_order  <- c("elongating","innate_Lymphoid","Leydig","Macrophage",
                 "round_spermatid","Sertoli","spermatocyte","spermatogonia")
stage_order <- c("ProSPG","Undiff SG","Diff SG",
                 "Spermatocyte","Round spermatid","Elongated spermatid")

# ══════════════════════════════════════════════════════════════
# 1. 加载并重建 celltype_8 注释
# ══════════════════════════════════════════════════════════════
cat("=== Loading v2 objects ===\n")
obj  <- readRDS(file.path(RDS_DIR, "seurat_exos_clustered_v2.rds"))
germ <- readRDS(file.path(RDS_DIR, "seurat_exos_germcells_v2.rds"))
cat(sprintf("Full obj: %d cells | Germ obj: %d cells\n", ncol(obj), ncol(germ)))

broad_ann  <- read.csv(file.path(EXOS_OUT, "01_heatmap/broad_cluster_annotation.csv"),
                       stringsAsFactors=FALSE)
germ_stage <- read.csv(file.path(EXOS_OUT, "02_annotation/germ_stage_annotation_v2.csv"),
                       stringsAsFactors=FALSE)

germ_cl_to_stage <- setNames(germ_stage$stage, as.character(germ_stage$cluster))
germ_cells_vec   <- germ$seurat_clusters
names(germ_cells_vec) <- colnames(germ)

stage_to_8 <- c(
  "ProSPG"="spermatogonia", "Undiff SG"="spermatogonia", "Diff SG"="spermatogonia",
  "Spermatocyte"="spermatocyte", "Round spermatid"="round_spermatid",
  "Elongated spermatid"="elongating"
)
somatic_to_8 <- c(
  "Sertoli"="Sertoli", "Leydig"="Leydig", "Macrophage"="Macrophage",
  "Endothelial"="innate_Lymphoid", "Myoid"="innate_Lymphoid"
)
broad_map <- setNames(broad_ann$celltype, as.character(broad_ann$cluster))

cl_full  <- as.character(obj$seurat_clusters)
type_vec <- character(ncol(obj))
names(type_vec) <- colnames(obj)

for (i in seq_along(type_vec)) {
  broad <- broad_map[cl_full[i]]
  if (is.na(broad) || broad == "Germ_cells") next
  v <- somatic_to_8[broad]
  if (!is.na(v)) type_vec[i] <- v
}

germ_common <- intersect(names(germ_cells_vec), colnames(obj))
sub_cl <- as.character(germ_cells_vec[germ_common])
stage_6 <- unname(germ_cl_to_stage[sub_cl])
stage_6[is.na(stage_6)] <- "Spermatocyte"
type_8_g <- unname(stage_to_8[stage_6])
type_8_g[is.na(type_8_g)] <- "spermatocyte"
type_vec[germ_common] <- type_8_g

unassigned <- names(type_vec)[type_vec == "" | is.na(type_vec)]
if (length(unassigned) > 0) type_vec[unassigned] <- "spermatocyte"

obj$celltype_8 <- factor(type_vec, levels=type_order)

# 始终从 v2 注释 CSV 重建 stage_label（覆盖旧版中文标签）
germ_cl_all <- as.character(germ$seurat_clusters)
sv <- sapply(germ_cl_all, function(cl) {
  if (cl %in% names(germ_cl_to_stage)) germ_cl_to_stage[cl] else "Spermatocyte"
})
germ$stage_label <- factor(unname(sv), levels=stage_order)

cat("8-type counts:\n"); print(sort(table(obj$celltype_8), decreasing=TRUE))
cat("Germ stage counts:\n"); print(table(germ$stage_label))

# 保存含注释的对象（供后续步骤使用）
saveRDS(obj,  file.path(RDS_DIR, "seurat_exos_annotated_v2.rds"))
saveRDS(germ, file.path(RDS_DIR, "seurat_exos_germcells_staged_v2.rds"))
cat("Saved annotated objects.\n")

# ══════════════════════════════════════════════════════════════
# 2. 细胞组成分析 (Cell Composition)
# ══════════════════════════════════════════════════════════════
cat("\n=== Cell Composition ===\n")

meta <- obj@meta.data
if (!"sample_id" %in% colnames(meta)) {
  meta$sample_id <- meta$orig.ident
}
if (!"group" %in% colnames(meta)) {
  meta$group <- ifelse(grepl("A[0-9]", meta$sample_id), "Aging", "Rescue")
}

comp_df <- meta %>%
  group_by(sample_id, group, celltype_8) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(sample_id) %>%
  mutate(total=sum(n), proportion=n/total) %>%
  ungroup() %>%
  filter(!is.na(celltype_8))

write.csv(comp_df, file.path(COMP_DIR, "celltype_proportion_by_sample.csv"), row.names=FALSE)

# 堆叠柱形图
p_bar <- ggplot(comp_df, aes(x=sample_id, y=proportion, fill=celltype_8)) +
  geom_col() +
  scale_y_continuous(labels=percent_format()) +
  scale_fill_manual(values=type_colors) +
  facet_wrap(~group, scales="free_x") +
  theme_cowplot(font_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.background=element_rect(fill="white", color=NA)) +
  labs(x=NULL, y="Proportion", fill="Cell type",
       title="Exos v2 — Cell type composition")
ggsave(file.path(COMP_DIR,"celltype_barplot.pdf"), p_bar, width=10, height=6)
ggsave(file.path(COMP_DIR,"celltype_barplot.png"), p_bar, width=10, height=6, dpi=200)
cat("Saved: celltype_barplot\n")

# 箱线图 + 统计（每种细胞类型）
box_df <- comp_df %>% filter(!is.na(celltype_8))

p_box <- ggplot(box_df, aes(x=group, y=proportion, fill=group)) +
  geom_boxplot(outlier.shape=NA, alpha=0.7) +
  geom_jitter(width=0.2, size=1.5, alpha=0.6) +
  scale_fill_manual(values=GROUP_COLS) +
  facet_wrap(~celltype_8, scales="free_y", ncol=4) +
  theme_cowplot(font_size=10) +
  theme(axis.text.x=element_text(angle=30, hjust=1),
        legend.position="none",
        plot.background=element_rect(fill="white", color=NA)) +
  labs(x=NULL, y="Proportion", title="Cell type proportion: Aging vs Rescue")
ggsave(file.path(COMP_DIR,"celltype_boxplot_by_group.pdf"), p_box, width=14, height=8)
ggsave(file.path(COMP_DIR,"celltype_boxplot_by_group.png"), p_box, width=14, height=8, dpi=200)
cat("Saved: celltype_boxplot_by_group\n")

# Wilcoxon 检验（每种细胞类型）
stat_res <- box_df %>%
  group_by(celltype_8) %>%
  summarise(
    n_Aging   = sum(group=="Aging"),
    n_Rescue  = sum(group=="Rescue"),
    mean_Aging= mean(proportion[group=="Aging"]),
    mean_Rescue=mean(proportion[group=="Rescue"]),
    log2FC    = log2((mean(proportion[group=="Aging"])+1e-6)/
                    (mean(proportion[group=="Rescue"])+1e-6)),
    p_wilcox  = tryCatch(
      wilcox.test(proportion[group=="Aging"],
                  proportion[group=="Rescue"])$p.value,
      error=function(e) NA_real_
    ),
    .groups="drop"
  ) %>%
  mutate(p_adj=p.adjust(p_wilcox, "BH"))

write.csv(stat_res, file.path(COMP_DIR,"celltype_composition_stats.csv"), row.names=FALSE)
cat("\nComposition stats:\n"); print(stat_res)

# 生殖细胞分期组成
germ_meta <- germ@meta.data
if (!"sample_id" %in% colnames(germ_meta)) germ_meta$sample_id <- germ_meta$orig.ident
if (!"group" %in% colnames(germ_meta))     germ_meta$group <- ifelse(grepl("A[0-9]", germ_meta$sample_id), "Aging", "Rescue")

germ_comp <- germ_meta %>%
  group_by(sample_id, group, stage_label) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(sample_id) %>%
  mutate(total=sum(n), proportion=n/total) %>%
  ungroup() %>%
  filter(!is.na(stage_label))

write.csv(germ_comp, file.path(COMP_DIR,"germcell_stage_proportion.csv"), row.names=FALSE)

p_germ_bar <- ggplot(germ_comp, aes(x=sample_id, y=proportion, fill=stage_label)) +
  geom_col() +
  scale_y_continuous(labels=percent_format()) +
  scale_fill_manual(values=stage_colors) +
  facet_wrap(~group, scales="free_x") +
  theme_cowplot(font_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.background=element_rect(fill="white", color=NA)) +
  labs(x=NULL, y="Proportion", fill="Stage", title="Germ cell stage composition")
ggsave(file.path(COMP_DIR,"germcell_stage_barplot.pdf"), p_germ_bar, width=10, height=6)
ggsave(file.path(COMP_DIR,"germcell_stage_barplot.png"), p_germ_bar, width=10, height=6, dpi=200)
cat("Saved: germcell_stage_barplot\n")

# ══════════════════════════════════════════════════════════════
# 3. DEG 分析 (Pseudobulk DESeq2)
# ══════════════════════════════════════════════════════════════
cat("\n=== DEG analysis (pseudobulk DESeq2) ===\n")

suppressPackageStartupMessages(library(DESeq2))

sample_info <- data.frame(
  sample_id = c("Exos1","Exos2","Exos3","Exos_A1","Exos_A2","Exos_A3"),
  group     = c("Rescue","Rescue","Rescue","Aging","Aging","Aging"),
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_info$sample_id

run_pseudobulk_deg <- function(cells_idx, ct_name, meta_all, counts_mat,
                               out_csv, out_pdf_prefix) {
  cells_ct <- colnames(counts_mat)[cells_idx]
  if (length(cells_ct) < 30) {
    cat(sprintf("  Skipping %s: too few cells (%d)\n", ct_name, length(cells_ct)))
    return(NULL)
  }
  cat(sprintf("  %s: %d cells\n", ct_name, length(cells_ct)))

  meta_ct <- meta_all[cells_ct, , drop=FALSE]
  counts_ct <- counts_mat[, cells_ct, drop=FALSE]

  # Pseudobulk aggregation per sample
  samples_present <- intersect(unique(meta_ct$sample_id), rownames(sample_info))
  pb_list <- lapply(samples_present, function(s) {
    sc <- rownames(meta_ct)[meta_ct$sample_id == s]
    if (length(sc) < 3) return(NULL)
    Matrix::rowSums(counts_ct[, sc, drop=FALSE])
  })
  names(pb_list) <- samples_present
  pb_list <- Filter(Negate(is.null), pb_list)

  if (length(pb_list) < 4) {
    cat(sprintf("  %s: <4 valid samples, skipping\n", ct_name))
    return(NULL)
  }

  pb_mat  <- do.call(cbind, pb_list)
  pb_meta <- data.frame(
    sample_id = colnames(pb_mat),
    group     = sample_info[colnames(pb_mat), "group"],
    row.names = colnames(pb_mat)
  )
  pb_meta$group <- factor(pb_meta$group, levels=c("Rescue","Aging"))

  # Filter lowly expressed genes
  keep_genes <- rowSums(pb_mat >= 1) >= 3
  pb_mat <- pb_mat[keep_genes, ]

  res_df <- tryCatch({
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData=pb_mat, colData=pb_meta, design=~group)
    dds <- DESeq2::estimateSizeFactors(dds)
    dds <- DESeq2::DESeq(dds, quiet=TRUE)
    r   <- DESeq2::results(dds, contrast=c("group","Aging","Rescue"),
                           independentFiltering=TRUE, alpha=0.05)
    data.frame(
      gene      = rownames(r),
      baseMean  = as.numeric(r$baseMean),
      log2FC    = as.numeric(r$log2FoldChange),
      lfcSE     = as.numeric(r$lfcSE),
      stat      = as.numeric(r$stat),
      pvalue    = as.numeric(r$pvalue),
      padj      = as.numeric(r$padj),
      celltype  = ct_name,
      direction = ifelse(!is.na(r$log2FoldChange) & r$log2FoldChange > 0,
                         "Aging_up", "Rescue_up"),
      stringsAsFactors=FALSE
    )
  }, error=function(e) {
    cat("  DESeq2 error:", conditionMessage(e), "\n")
    # Fallback: log-CPM t-test
    cpm_mat <- sweep(pb_mat, 2, colSums(pb_mat), "/")*1e6
    lcpm    <- log1p(cpm_mat)
    aging_s  <- rownames(pb_meta)[pb_meta$group=="Aging"]
    rescue_s <- rownames(pb_meta)[pb_meta$group=="Rescue"]
    log2fc   <- rowMeans(lcpm[,aging_s,drop=FALSE]) - rowMeans(lcpm[,rescue_s,drop=FALSE])
    pval     <- apply(lcpm, 1, function(x) {
      tryCatch(t.test(x[aging_s], x[rescue_s])$p.value, error=function(e) NA)
    })
    data.frame(gene=rownames(lcpm), baseMean=rowMeans(pb_mat),
               log2FC=log2fc, lfcSE=NA, stat=NA,
               pvalue=pval, padj=p.adjust(pval,"BH"),
               celltype=ct_name,
               direction=ifelse(log2fc>0,"Aging_up","Rescue_up"),
               stringsAsFactors=FALSE)
  })

  if (!is.null(res_df)) {
    write.csv(res_df, out_csv, row.names=FALSE)
    cat(sprintf("  %s: %d genes tested, %d sig (padj<0.05 & |FC|>1)\n",
                ct_name, nrow(res_df),
                sum(!is.na(res_df$padj) & res_df$padj<0.05 & abs(res_df$log2FC)>1)))

    # Volcano plot
    vdf <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FC), ]
    vdf$sig <- vdf$padj < 0.05 & abs(vdf$log2FC) > 1
    top_genes <- vdf %>% filter(sig) %>%
      arrange(padj) %>% head(20) %>% pull(gene)
    vdf$label <- ifelse(vdf$gene %in% top_genes, vdf$gene, "")
    vdf$color <- case_when(
      vdf$padj < 0.05 & vdf$log2FC > 1  ~ "Aging_up",
      vdf$padj < 0.05 & vdf$log2FC < -1 ~ "Rescue_up",
      TRUE ~ "NS"
    )
    p_vol <- ggplot(vdf, aes(x=log2FC, y=-log10(pmax(padj,1e-300)), color=color)) +
      geom_point(alpha=0.5, size=1) +
      geom_text_repel(aes(label=label), size=2.5, max.overlaps=20,
                      box.padding=0.3, min.segment.length=0) +
      scale_color_manual(values=c("Aging_up"="#E64B35","Rescue_up"="#4DBBD5","NS"="grey70"),
                         guide=guide_legend(override.aes=list(size=3))) +
      geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey50") +
      geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50") +
      theme_cowplot(font_size=11) +
      theme(plot.background=element_rect(fill="white",color=NA)) +
      labs(title=paste0(ct_name, " — Aging vs Rescue"),
           x="log2FC (Aging/Rescue)", y="-log10(adj.p)", color=NULL)
    ggsave(paste0(out_pdf_prefix,"_volcano.pdf"), p_vol, width=7, height=6)
    ggsave(paste0(out_pdf_prefix,"_volcano.png"), p_vol, width=7, height=6, dpi=180)
  }
  return(res_df)
}

# 获取 counts 矩阵和 meta
counts_full <- GetAssayData(obj, assay="RNA", layer="counts")
meta_full   <- obj@meta.data
if (!"sample_id" %in% colnames(meta_full)) meta_full$sample_id <- meta_full$orig.ident
if (!"group"     %in% colnames(meta_full)) {
  meta_full$group <- ifelse(grepl("A[0-9]", meta_full$sample_id), "Aging", "Rescue")
}

all_deg_list <- list()
for (ct in type_order) {
  ct_safe <- gsub("[/ ]","_", ct)
  out_csv <- file.path(DEG_DIR, paste0(ct_safe, "_Aging_vs_Rescue_DEG.csv"))
  out_pfx <- file.path(DEG_DIR, ct_safe)

  if (file.exists(out_csv)) {
    cat(sprintf("  %s: DEG file exists, loading...\n", ct))
    all_deg_list[[ct]] <- read.csv(out_csv, stringsAsFactors=FALSE)
    next
  }

  cells_idx <- which(obj$celltype_8 == ct)
  res <- run_pseudobulk_deg(
    cells_idx=cells_idx, ct_name=ct,
    meta_all=meta_full, counts_mat=counts_full,
    out_csv=out_csv, out_pdf_prefix=out_pfx
  )
  if (!is.null(res)) all_deg_list[[ct]] <- res
  gc()
}

cat("DEG analysis done.\n")
rm(counts_full); gc()

# 合并 top DEGs 热图
all_deg_combined <- do.call(rbind, all_deg_list)
write.csv(all_deg_combined, file.path(DEG_DIR,"all_celltypes_DEG.csv"), row.names=FALSE)

# Top-20 sig DEG summary
top_deg <- all_deg_combined %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FC) > 1) %>%
  group_by(celltype, direction) %>%
  slice_min(padj, n=10, with_ties=FALSE) %>%
  arrange(celltype, direction, padj)
write.csv(top_deg, file.path(DEG_DIR,"top_DEG_summary.csv"), row.names=FALSE)
cat(sprintf("Significant DEGs (padj<0.05, |log2FC|>1): %d total\n",
            sum(!is.na(all_deg_combined$padj) & all_deg_combined$padj<0.05 &
                abs(all_deg_combined$log2FC)>1)))

# ══════════════════════════════════════════════════════════════
# 4. Pseudotime — Slingshot on germ cells
# ══════════════════════════════════════════════════════════════
cat("\n=== Slingshot pseudotime (germ cells) ===\n")

# 最多 25k 细胞
MAX_CELLS <- 25000
germ_sub  <- germ
if (ncol(germ_sub) > MAX_CELLS) {
  set.seed(SEED)
  keep <- sample(colnames(germ_sub), MAX_CELLS)
  germ_sub <- subset(germ_sub, cells=keep)
  cat(sprintf("Downsampled germ cells: %d → %d\n", ncol(germ), ncol(germ_sub)))
}

umap_emb    <- Embeddings(germ_sub, "umap")
stage_labels <- as.character(germ_sub$stage_label)
available_stages <- stage_order[stage_order %in% unique(stage_labels)]
start_stage <- available_stages[1]
cat(sprintf("Start stage: %s\n", start_stage))
cat(sprintf("Stages available: %s\n", paste(available_stages, collapse=", ")))

sce <- SingleCellExperiment(
  assays  = list(counts=matrix(0L, nrow=1, ncol=ncol(germ_sub))),
  colData = data.frame(
    stage_label = stage_labels,
    group       = as.character(germ_sub$group),
    row.names   = colnames(germ_sub)
  )
)
reducedDim(sce, "UMAP") <- umap_emb

sds <- tryCatch(
  slingshot(sce, clusterLabels="stage_label", reducedDim="UMAP",
            start.clus=start_stage),
  error=function(e) { cat("Slingshot error:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(sds)) {
  pt_mat <- slingPseudotime(sds)
  pt_mean <- rowMeans(pt_mat, na.rm=TRUE)

  pt_df <- data.frame(
    cell        = colnames(germ_sub),
    sample_id   = germ_sub$orig.ident,
    group       = as.character(germ_sub$group),
    stage_label = stage_labels,
    pseudotime  = pt_mean,
    stringsAsFactors=FALSE
  )
  write.csv(pt_df, file.path(PSEUDO_DIR,"pseudotime_cells.csv"), row.names=FALSE)

  # UMAP colored by pseudotime
  germ_sub$pseudotime <- pt_df$pseudotime[match(colnames(germ_sub), pt_df$cell)]
  p_pt <- FeaturePlot(germ_sub, features="pseudotime",
                      min.cutoff="q5", max.cutoff="q95", pt.size=0.15) +
    scale_color_viridis_c(option="plasma", name="Pseudotime") +
    labs(title="Exos v2 — Pseudotime (Slingshot)") +
    theme_cowplot(font_size=11) +
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          plot.background=element_rect(fill="white", color=NA)) +
    NoAxes()
  ggsave(file.path(PSEUDO_DIR,"umap_pseudotime.pdf"), p_pt, width=7, height=6)
  ggsave(file.path(PSEUDO_DIR,"umap_pseudotime.png"), p_pt, width=7, height=6, dpi=200)

  # UMAP colored by stage
  p_stage <- DimPlot(germ_sub, group.by="stage_label",
                     cols=stage_colors, pt.size=0.15, label=TRUE, label.size=3) +
    labs(title="Germ cell stages") +
    theme_cowplot(font_size=11) +
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          plot.background=element_rect(fill="white", color=NA)) +
    guides(color=guide_legend(override.aes=list(size=3)))
  ggsave(file.path(PSEUDO_DIR,"umap_germ_stages.pdf"), p_stage, width=8, height=6)
  ggsave(file.path(PSEUDO_DIR,"umap_germ_stages.png"), p_stage, width=8, height=6, dpi=200)

  # 密度图：Aging vs Rescue
  p_dens <- ggplot(pt_df[!is.na(pt_df$pseudotime),],
                   aes(x=pseudotime, fill=group, color=group)) +
    geom_density(alpha=0.4, linewidth=0.8) +
    scale_fill_manual(values=GROUP_COLS) +
    scale_color_manual(values=GROUP_COLS) +
    theme_cowplot(font_size=12) +
    theme(plot.background=element_rect(fill="white", color=NA)) +
    labs(title="Pseudotime distribution: Aging vs Rescue",
         x="Pseudotime", y="Density", fill=NULL, color=NULL)
  ggsave(file.path(PSEUDO_DIR,"pseudotime_density.pdf"), p_dens, width=7, height=5)
  ggsave(file.path(PSEUDO_DIR,"pseudotime_density.png"), p_dens, width=7, height=5, dpi=200)

  # 箱线图 by group + stage
  p_box_pt <- ggplot(pt_df[!is.na(pt_df$pseudotime),],
                     aes(x=group, y=pseudotime, fill=group)) +
    geom_boxplot(alpha=0.7, outlier.size=0.5) +
    geom_jitter(width=0.2, size=0.3, alpha=0.3) +
    scale_fill_manual(values=GROUP_COLS) +
    facet_wrap(~stage_label, ncol=3) +
    theme_cowplot(font_size=10) +
    theme(axis.text.x=element_text(angle=30, hjust=1),
          legend.position="none",
          plot.background=element_rect(fill="white", color=NA)) +
    labs(title="Pseudotime by stage and group", x=NULL, y="Pseudotime")
  ggsave(file.path(PSEUDO_DIR,"pseudotime_boxplot_by_stage.pdf"), p_box_pt, width=10, height=8)
  ggsave(file.path(PSEUDO_DIR,"pseudotime_boxplot_by_stage.png"), p_box_pt, width=10, height=8, dpi=200)

  # 每阶段 Aging vs Rescue Wilcoxon
  pt_stat <- pt_df %>%
    filter(!is.na(pseudotime)) %>%
    group_by(stage_label) %>%
    summarise(
      n_Aging   = sum(group=="Aging"),
      n_Rescue  = sum(group=="Rescue"),
      mean_Aging= mean(pseudotime[group=="Aging"],  na.rm=TRUE),
      mean_Rescue=mean(pseudotime[group=="Rescue"], na.rm=TRUE),
      p_wilcox  = tryCatch(
        wilcox.test(pseudotime[group=="Aging"],
                    pseudotime[group=="Rescue"])$p.value,
        error=function(e) NA_real_
      ),
      .groups="drop"
    ) %>% mutate(p_adj=p.adjust(p_wilcox,"BH"))
  write.csv(pt_stat, file.path(PSEUDO_DIR,"pseudotime_by_stage_stats.csv"), row.names=FALSE)
  cat("\nPseudotime stage stats:\n"); print(pt_stat)

  saveRDS(sds, file.path(RDS_DIR,"slingshot_sds_v2.rds"))
  cat("Saved: slingshot_sds_v2.rds\n")
} else {
  cat("Slingshot failed, skipping pseudotime plots.\n")
}

cat("\n===== DONE =====\n")
cat(sprintf("Outputs:\n  %s\n  %s\n  %s\n", COMP_DIR, DEG_DIR, PSEUDO_DIR))

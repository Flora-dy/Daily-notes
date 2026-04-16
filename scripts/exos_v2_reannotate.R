#!/usr/bin/env Rscript
# ============================================================
# Exos v2 — Full re-annotation + Spermatogonia sub-clustering
#
# Part 1: Full 8-type annotation
#   Spermatogonia / Spermatocyte / Round Spermatid /
#   Elongated Spermatid / Sertoli Cell / Leydig Cell /
#   Macrophage / Endothelial Cell
#   "Cluster sync": UMAP shows both cluster numbers and cell type labels
#
# Part 2: Spermatogonia 3-subtype sub-clustering
#   Extract SPG cells → re-run UMAP → 3 subtypes:
#   Prospermatogonium / Undifferentiated Spermatogonium /
#   Differentiating Spermatogonium
#   "Cluster sync": sub-cluster IDs labeled with subtype names
#
# Part 3: Final 10-type annotation on full object
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
RDS_DIR  <- file.path(EXOS_OUT, "rds")
dir.create(A_DIR,   recursive=TRUE, showWarnings=FALSE)
dir.create(RDS_DIR, recursive=TRUE, showWarnings=FALSE)

SEED <- 202409; set.seed(SEED)
GROUP_COLS <- c("Rescue"="#4DBBD5", "Aging"="#E64B35")

# ─── Color palettes ───────────────────────────────────────────
celltype_colors <- c(
  "Spermatogonia"      = "#D62728",
  "Spermatocyte"       = "#1F77B4",
  "Round Spermatid"    = "#2CA02C",
  "Elongated Spermatid"= "#9467BD",
  "Sertoli Cell"       = "#8C564B",
  "Leydig Cell"        = "#E377C2",
  "Macrophage"         = "#7F7F7F",
  "Endothelial Cell"   = "#17BECF"
)
celltype_order <- names(celltype_colors)

spg_colors <- c(
  "Prospermatogonium"               = "#7B2D8B",
  "Undifferentiated Spermatogonium" = "#D62728",
  "Differentiating Spermatogonium"  = "#FF7F0E"
)
spg_order <- names(spg_colors)

final_colors <- c(
  "Prospermatogonium"               = "#7B2D8B",
  "Undifferentiated Spermatogonium" = "#D62728",
  "Differentiating Spermatogonium"  = "#FF7F0E",
  "Spermatocyte"                    = "#1F77B4",
  "Round Spermatid"                 = "#2CA02C",
  "Elongated Spermatid"             = "#9467BD",
  "Sertoli Cell"                    = "#8C564B",
  "Leydig Cell"                     = "#E377C2",
  "Macrophage"                      = "#7F7F7F",
  "Endothelial Cell"                = "#17BECF"
)
final_order <- names(final_colors)

# ══════════════════════════════════════════════════════════════
# Helper: make clean DotPlot from dp_data (pre-scaled)
# ══════════════════════════════════════════════════════════════
make_dotplot <- function(dp_data, gene_list, color_map, id_order, title="") {
  gene_order <- unname(unlist(gene_list))
  dp_data$features.plot <- factor(as.character(dp_data$features.plot), levels=gene_order)
  dp_data$id <- factor(as.character(dp_data$id),
                       levels=rev(id_order[id_order %in% as.character(dp_data$id)]))

  dp_data <- dp_data %>%
    group_by(features.plot) %>%
    mutate(exp_norm = {
      rng <- range(avg.exp.scaled, na.rm=TRUE)
      if (diff(rng) < 1e-6) rep(0, n())
      else (avg.exp.scaled - rng[1]) / diff(rng)
    }) %>% ungroup()

  seg_grp  <- factor(rep(names(gene_list), sapply(gene_list, length)),
                     levels=names(gene_list))
  sep_x    <- which(diff(as.integer(seg_grp)) != 0) + 0.5
  sec_df   <- data.frame(gene=gene_order, grp=seg_grp, x=seq_along(gene_order))
  sec_ctr  <- sec_df %>% group_by(grp) %>% summarise(cx=mean(x), .groups="drop")

  p_bar <- ggplot(sec_df, aes(x=x, y=1, fill=grp)) +
    geom_tile(width=0.95, height=0.9) +
    geom_text(data=sec_ctr, aes(x=cx, y=1, label=grp, fill=NULL),
              size=2.4, fontface="bold", color="white") +
    scale_fill_manual(values=color_map[levels(sec_df$grp)]) +
    scale_x_continuous(expand=c(0,0.5)) +
    theme_void() + theme(legend.position="none",
                         plot.background=element_rect(fill="white",color=NA))

  p_dot <- ggplot(dp_data,
                  aes(x=features.plot, y=id, size=pct.exp, color=exp_norm)) +
    geom_vline(xintercept=sep_x, color="grey72", linewidth=0.4, linetype="dashed") +
    geom_point() +
    scale_color_gradientn(
      colours=c("grey92","grey82","#AED6F1","#2E86C1","#1A5276"),
      values=c(0,0.2,0.5,0.75,1), limits=c(0,1),
      name="Relative\nExpression") +
    scale_size_continuous(range=c(0.4,9), name="% Expressing",
                          breaks=c(10,25,50,75,100)) +
    scale_y_discrete(position="right") +
    labs(x=NULL, y=NULL, title=title) +
    theme_classic(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,face="italic",size=9.5),
          axis.text.y=element_text(size=11,face="bold"),
          axis.line=element_blank(),
          panel.border=element_rect(color="grey40",fill=NA,linewidth=0.7),
          panel.grid.major.y=element_line(color="grey93",linewidth=0.3),
          plot.title=element_text(hjust=0.5,face="bold",size=13),
          plot.background=element_rect(fill="white",color=NA),
          legend.title=element_text(size=9), legend.text=element_text(size=8.5))

  ratio <- max(1, length(gene_list)) + 1
  p_bar / p_dot + plot_layout(heights=c(1, ratio)) +
    plot_annotation(theme=theme(plot.background=element_rect(fill="white",color=NA)))
}

# ══════════════════════════════════════════════════════════════
# STEP 0. Load data
# ══════════════════════════════════════════════════════════════
cat("=== Loading objects ===\n")
obj  <- readRDS(file.path(RDS_DIR, "seurat_exos_clustered_v2.rds"))
germ <- readRDS(file.path(RDS_DIR, "seurat_exos_germcells_v2.rds"))
cat(sprintf("Full obj: %d cells, %d clusters\n", ncol(obj), length(unique(obj$seurat_clusters))))
cat(sprintf("Germ obj: %d cells, %d sub-clusters\n", ncol(germ), length(unique(germ$seurat_clusters))))

broad_ann  <- read.csv(file.path(EXOS_OUT, "01_heatmap/broad_cluster_annotation.csv"),
                       stringsAsFactors=FALSE)
germ_stage <- read.csv(file.path(A_DIR, "germ_stage_annotation_v2.csv"),
                       stringsAsFactors=FALSE)

broad_map        <- setNames(broad_ann$celltype, as.character(broad_ann$cluster))
germ_cl_to_stage <- setNames(germ_stage$stage,  as.character(germ_stage$cluster))

stage_to_celltype <- c(
  "ProSPG"             = "Spermatogonia",
  "Undiff SG"          = "Spermatogonia",
  "Diff SG"            = "Spermatogonia",
  "Spermatocyte"       = "Spermatocyte",
  "Round spermatid"    = "Round Spermatid",
  "Elongated spermatid"= "Elongated Spermatid"
)
somatic_map <- c(
  "Sertoli"    = "Sertoli Cell",
  "Leydig"     = "Leydig Cell",
  "Macrophage" = "Macrophage",
  "Endothelial"= "Endothelial Cell",
  "Myoid"      = "Endothelial Cell"
)
germ_cl_vec <- setNames(as.character(germ$seurat_clusters), colnames(germ))

# ══════════════════════════════════════════════════════════════
# PART 1 — Full 8-type annotation
# ══════════════════════════════════════════════════════════════
cat("\n=== PART 1: 8-type annotation ===\n")

cl_full <- as.character(obj$seurat_clusters)
ct_vec  <- rep("Spermatocyte", ncol(obj))   # default
names(ct_vec) <- colnames(obj)

# Somatic cells
for (i in seq_along(ct_vec)) {
  broad <- broad_map[cl_full[i]]
  if (is.na(broad) || broad == "Germ_cells") next
  v <- somatic_map[broad]
  if (!is.na(v)) ct_vec[i] <- v
}

# Germ cells via germ sub-cluster stage
germ_in_full <- intersect(names(germ_cl_vec), colnames(obj))
sub_cl  <- germ_cl_vec[germ_in_full]
stage_6 <- germ_cl_to_stage[sub_cl]
stage_6[is.na(stage_6)] <- "Spermatocyte"
ct_germ <- stage_to_celltype[stage_6]
ct_germ[is.na(ct_germ)] <- "Spermatocyte"
ct_vec[germ_in_full] <- ct_germ

obj$celltype <- factor(ct_vec, levels=celltype_order)
cat("Cell type counts:\n")
print(sort(table(obj$celltype), decreasing=TRUE))

# Cluster-to-celltype mapping (cluster sync table)
cl_ct_tbl <- obj@meta.data %>%
  group_by(seurat_clusters, celltype) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(seurat_clusters) %>%
  slice_max(n, n=1, with_ties=FALSE) %>%
  arrange(as.integer(as.character(seurat_clusters)))
write.csv(cl_ct_tbl, file.path(A_DIR,"cluster_celltype_mapping.csv"), row.names=FALSE)
cat("\nCluster -> Cell type (cluster sync):\n")
print(as.data.frame(cl_ct_tbl[,c("seurat_clusters","celltype","n")]))

# Add cluster sync label to metadata (positional assignment — no name mismatch)
cl_ct_vec <- as.character(cl_ct_tbl$celltype)
names(cl_ct_vec) <- as.character(cl_ct_tbl$seurat_clusters)
obj$cluster_sync <- paste0("C", cl_full, ":", cl_ct_vec[cl_full])

Idents(obj) <- "celltype"

# ── Marker genes ─────────────────────────────────────────────
mk8 <- list(
  "Spermatogonia"      = c("Id4","Zbtb16","Uchl1","Gfra1","Nanos3","Kit","Stra8"),
  "Spermatocyte"       = c("Sycp3","Hormad1","H2afx","Piwil1","Spo11","Meioc","Rec8"),
  "Round Spermatid"    = c("Acrv1","Spaca1","Tssk1","Pgk2","Hspa1l"),
  "Elongated Spermatid"= c("Prm1","Prm2","Tnp1","Tnp2","Akap4"),
  "Sertoli Cell"       = c("Sox9","Wt1","Amh","Cldn11","Ar"),
  "Leydig Cell"        = c("Star","Insl3","Cyp11a1","Hsd3b1"),
  "Macrophage"         = c("Cd14","C1qb","Csf1r","Cx3cr1","Adgre1"),
  "Endothelial Cell"   = c("Pecam1","Vwf","Cdh5","Tek","Ncr1")
)
mk8 <- lapply(mk8, function(gs) gs[gs %in% rownames(obj)])
mk8 <- mk8[sapply(mk8, length) > 0]
gene8 <- unname(unlist(mk8))
cat(sprintf("\nMarker genes: %d total\n", length(gene8)))

cat("Computing 8-type DotPlot...\n")
dp8 <- DotPlot(obj, features=gene8, group.by="celltype")$data
p8 <- make_dotplot(dp8, mk8, celltype_colors, celltype_order,
                   title="Exos v2 - Cell Type Marker Expression")
n8 <- length(gene8)
w8 <- max(14, n8 * 0.38 + 4)
ggsave(file.path(A_DIR,"exos_v2_celltype_dotplot.pdf"), p8, width=w8, height=7)
ggsave(file.path(A_DIR,"exos_v2_celltype_dotplot.png"), p8, width=w8, height=7, dpi=220)
cat(sprintf("Saved: exos_v2_celltype_dotplot (%.1f x 7 in)\n", w8))

# ── UMAPs ────────────────────────────────────────────────────
p_ct <- DimPlot(obj, reduction="umap", group.by="celltype",
                cols=celltype_colors, pt.size=0.12, label=TRUE,
                label.size=3.5, repel=TRUE) +
  labs(title="Cell Types") +
  theme_cowplot(font_size=11) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white")) +
  guides(color=guide_legend(override.aes=list(size=4)))

p_cl <- DimPlot(obj, reduction="umap", group.by="seurat_clusters",
                label=TRUE, label.size=3, pt.size=0.12) +
  labs(title="Original Clusters (Cluster Sync)") +
  theme_cowplot(font_size=11) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white")) + NoLegend()

p_grp <- DimPlot(obj, reduction="umap", group.by="group",
                 cols=GROUP_COLS, pt.size=0.12) +
  labs(title="Aging vs Rescue") +
  theme_cowplot(font_size=11) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white"))

ggsave(file.path(A_DIR,"exos_v2_celltype_umap.pdf"),
       p_ct | p_cl | p_grp, width=22, height=6.5)
ggsave(file.path(A_DIR,"exos_v2_celltype_umap.png"),
       p_ct | p_cl | p_grp, width=22, height=6.5, dpi=200)
cat("Saved: exos_v2_celltype_umap\n")

# ══════════════════════════════════════════════════════════════
# PART 2 — Spermatogonia 3-subtype sub-clustering
# ══════════════════════════════════════════════════════════════
cat("\n=== PART 2: Spermatogonia sub-clustering ===\n")

# Extract SPG cells (germ clusters 14=Undiff SG, 16=Diff SG, 20=ProSPG)
spg_stages   <- c("ProSPG","Undiff SG","Diff SG")
spg_clusters <- names(germ_cl_to_stage)[germ_cl_to_stage %in% spg_stages]
cat(sprintf("SPG germ sub-clusters: %s\n", paste(spg_clusters, collapse=", ")))

spg_cells <- colnames(germ)[as.character(germ$seurat_clusters) %in% spg_clusters]
cat(sprintf("Total SPG cells: %d\n", length(spg_cells)))
if (length(spg_cells) == 0) stop("No SPG cells found")

# Subtype label for each SPG cell (from existing validated annotation)
spg_stage_vec   <- germ_cl_to_stage[as.character(germ$seurat_clusters[spg_cells])]
names(spg_stage_vec) <- spg_cells
cat("SPG stage distribution:\n")
print(table(spg_stage_vec))

# Professional subtype names
subtype_prof <- c(
  "ProSPG"   = "Prospermatogonium",
  "Undiff SG"= "Undifferentiated Spermatogonium",
  "Diff SG"  = "Differentiating Spermatogonium"
)

# Create SPG sub-object
spg_obj <- subset(germ, cells=spg_cells)
spg_obj$spg_subtype <- factor(
  unname(subtype_prof[as.character(spg_stage_vec[colnames(spg_obj)])]),
  levels=spg_order
)
cat("\nSPG subtype counts:\n")
print(table(spg_obj$spg_subtype))

# Re-process for a fresh SPG-specific UMAP
cat("Re-processing SPG subset for sub-clustering UMAP...\n")
DefaultAssay(spg_obj) <- "RNA"
spg_obj <- NormalizeData(spg_obj, verbose=FALSE)
spg_obj <- FindVariableFeatures(spg_obj, nfeatures=2000, verbose=FALSE)
spg_obj <- ScaleData(spg_obj, verbose=FALSE)
spg_obj <- RunPCA(spg_obj, npcs=20, verbose=FALSE)
spg_obj <- RunUMAP(spg_obj, dims=1:15, seed.use=SEED, verbose=FALSE)
spg_obj <- FindNeighbors(spg_obj, dims=1:15, verbose=FALSE)

# FindClusters at multiple resolutions — find the one closest to 3 clusters
best_res <- 0.1; best_diff <- Inf; best_ncl <- 0
for (res in c(0.04, 0.06, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25)) {
  tmp <- FindClusters(spg_obj, resolution=res, random.seed=SEED, verbose=FALSE)
  n   <- length(unique(tmp$seurat_clusters))
  cat(sprintf("  res %.2f -> %d clusters\n", res, n))
  if (abs(n - 3) < best_diff ||
      (abs(n - 3) == best_diff && n >= 3)) {
    best_diff <- abs(n - 3); best_res <- res; best_ncl <- n; spg_obj <- tmp
    if (n == 3) break
  }
}
cat(sprintf("Selected resolution %.2f (%d clusters)\n", best_res, best_ncl))

# Add sub-cluster sync: show cluster ID alongside subtype name
# (map each new sub-cluster to the dominant spg_subtype)
sc_to_subtype <- spg_obj@meta.data %>%
  group_by(seurat_clusters, spg_subtype) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(seurat_clusters) %>%
  slice_max(n, n=1, with_ties=FALSE)
sc_sync_map <- setNames(as.character(sc_to_subtype$spg_subtype),
                        as.character(sc_to_subtype$seurat_clusters))
spg_obj$cluster_sync <- paste0(
  "SC", as.character(spg_obj$seurat_clusters), ":",
  sc_sync_map[as.character(spg_obj$seurat_clusters)]
)

cat("\nSub-cluster -> Subtype (cluster sync):\n")
print(as.data.frame(sc_to_subtype))
write.csv(sc_to_subtype, file.path(A_DIR,"spg_subcluster_annotation.csv"),
          row.names=FALSE)

Idents(spg_obj) <- "spg_subtype"

# ── SPG DotPlot ───────────────────────────────────────────────
mk_spg <- list(
  "Prospermatogonium"               = c("Id1","Id2","Id4","Klf4","Gfra1","Ret","Thy1"),
  "Undifferentiated Spermatogonium" = c("Zbtb16","Uchl1","Sall4","Cdh1",
                                         "Nanos3","Etv5","Egr4","Pax7"),
  "Differentiating Spermatogonium"  = c("Kit","Stra8","Dmrtb1","Sohlh1",
                                         "Dnmt3b","Prdm9","Rhox13")
)
mk_spg <- lapply(mk_spg, function(gs) gs[gs %in% rownames(spg_obj)])
mk_spg <- mk_spg[sapply(mk_spg, length) > 0]
gene_spg <- unname(unlist(mk_spg))
cat(sprintf("\nSPG marker genes: %d total\n", length(gene_spg)))

cat("Computing SPG DotPlot...\n")
dp_spg <- DotPlot(spg_obj, features=gene_spg, group.by="spg_subtype")$data
p_spg_dot <- make_dotplot(dp_spg, mk_spg, spg_colors, spg_order,
                          title="Spermatogonia Subtype Markers")
n_sg <- length(gene_spg)
w_sg <- max(12, n_sg * 0.45 + 3)
ggsave(file.path(A_DIR,"exos_v2_spg_dotplot.pdf"), p_spg_dot, width=w_sg, height=6)
ggsave(file.path(A_DIR,"exos_v2_spg_dotplot.png"), p_spg_dot, width=w_sg, height=6, dpi=220)
cat(sprintf("Saved: exos_v2_spg_dotplot (%.1f x 6 in)\n", w_sg))

# ── SPG UMAPs (3-panel) ───────────────────────────────────────
p_spg_type <- DimPlot(spg_obj, reduction="umap", group.by="spg_subtype",
                      cols=spg_colors, pt.size=0.5, label=TRUE,
                      label.size=4, repel=TRUE) +
  labs(title="Spermatogonia Subtypes") +
  theme_cowplot(font_size=11) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white")) +
  guides(color=guide_legend(override.aes=list(size=4)))

p_spg_cl <- DimPlot(spg_obj, reduction="umap", group.by="seurat_clusters",
                    label=TRUE, label.size=4, pt.size=0.5) +
  labs(title="SPG Sub-clusters (Cluster Sync)") +
  theme_cowplot(font_size=11) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white")) + NoLegend()

p_spg_grp <- DimPlot(spg_obj, reduction="umap", group.by="group",
                     cols=GROUP_COLS, pt.size=0.5) +
  labs(title="Aging vs Rescue") +
  theme_cowplot(font_size=11) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white"))

ggsave(file.path(A_DIR,"exos_v2_spg_umap.pdf"),
       p_spg_type | p_spg_cl | p_spg_grp, width=20, height=6)
ggsave(file.path(A_DIR,"exos_v2_spg_umap.png"),
       p_spg_type | p_spg_cl | p_spg_grp, width=20, height=6, dpi=200)
cat("Saved: exos_v2_spg_umap\n")

# ── Violin plots ──────────────────────────────────────────────
vln_g <- c("Id4","Zbtb16","Gfra1","Kit","Stra8")
vln_g <- vln_g[vln_g %in% rownames(spg_obj)]
if (length(vln_g) > 0) {
  p_vln <- VlnPlot(spg_obj, features=vln_g, group.by="spg_subtype",
                   cols=spg_colors, pt.size=0, combine=TRUE) &
    theme(axis.text.x=element_text(angle=30,hjust=1,size=9),
          plot.title=element_text(size=11,face="italic"),
          axis.title.x=element_blank(),
          plot.background=element_rect(fill="white",color=NA))
  ggsave(file.path(A_DIR,"exos_v2_spg_violin.pdf"), p_vln,
         width=length(vln_g)*3.2, height=5)
  ggsave(file.path(A_DIR,"exos_v2_spg_violin.png"), p_vln,
         width=length(vln_g)*3.2, height=5, dpi=200)
  cat("Saved: exos_v2_spg_violin\n")
}

# ══════════════════════════════════════════════════════════════
# PART 3 — Final 10-type annotation on full object
# ══════════════════════════════════════════════════════════════
cat("\n=== PART 3: Final 10-type annotation ===\n")

# Map SPG cells to their specific subtype
spg_in_obj <- intersect(colnames(spg_obj), colnames(obj))
spg_sub_vec <- as.character(spg_obj$spg_subtype)[
  match(spg_in_obj, colnames(spg_obj))]

ct_final <- as.character(obj$celltype)
names(ct_final) <- colnames(obj)
ct_final[spg_in_obj] <- spg_sub_vec

# Any residual Spermatogonia not in SPG sub-obj
still_spg <- names(ct_final)[ct_final == "Spermatogonia"]
if (length(still_spg) > 0) {
  ct_final[still_spg] <- "Undifferentiated Spermatogonium"
  cat(sprintf("  %d SPG cells not sub-typed -> Undifferentiated Spermatogonium\n",
              length(still_spg)))
}

obj$celltype_final <- factor(ct_final, levels=final_order)
Idents(obj) <- "celltype_final"

cat("Final cell type counts:\n")
print(sort(table(obj$celltype_final), decreasing=TRUE))

# Cluster sync for final annotation
cl_final_tbl <- obj@meta.data %>%
  group_by(seurat_clusters, celltype_final) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(seurat_clusters) %>%
  slice_max(n, n=1, with_ties=FALSE) %>%
  arrange(as.integer(as.character(seurat_clusters)))
write.csv(cl_final_tbl, file.path(A_DIR,"cluster_final_celltype_mapping.csv"),
          row.names=FALSE)
cat("\nFinal cluster sync:\n")
print(as.data.frame(cl_final_tbl[,c("seurat_clusters","celltype_final","n")]))

# ── Final UMAPs ──────────────────────────────────────────────
p_fin <- DimPlot(obj, reduction="umap", group.by="celltype_final",
                 cols=final_colors, pt.size=0.12, label=TRUE,
                 label.size=2.8, repel=TRUE) +
  labs(title="Exos v2 - Final Cell Type Annotation") +
  theme_cowplot(font_size=10) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white")) +
  guides(color=guide_legend(override.aes=list(size=4), ncol=1))

p_fin_grp <- DimPlot(obj, reduction="umap", group.by="group",
                     cols=GROUP_COLS, pt.size=0.12) +
  labs(title="Aging vs Rescue") +
  theme_cowplot(font_size=11) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.background=element_rect(fill="white"))

ggsave(file.path(A_DIR,"exos_v2_final_umap.pdf"),
       p_fin | p_fin_grp, width=18, height=6.5)
ggsave(file.path(A_DIR,"exos_v2_final_umap.png"),
       p_fin | p_fin_grp, width=18, height=6.5, dpi=200)
cat("Saved: exos_v2_final_umap\n")

# ── Final DotPlot (10 types) ──────────────────────────────────
mk10 <- list(
  "Prospermatogonium"               = c("Id2","Id4","Klf4","Gfra1","Ret"),
  "Undifferentiated Spermatogonium" = c("Zbtb16","Uchl1","Sall4","Nanos3","Egr4"),
  "Differentiating Spermatogonium"  = c("Kit","Stra8","Dmrtb1","Sohlh1"),
  "Spermatocyte"                    = c("Sycp3","Hormad1","H2afx","Piwil1","Rec8"),
  "Round Spermatid"                 = c("Acrv1","Spaca1","Tssk1","Pgk2"),
  "Elongated Spermatid"             = c("Prm1","Prm2","Tnp1","Akap4"),
  "Sertoli Cell"                    = c("Sox9","Wt1","Amh","Cldn11"),
  "Leydig Cell"                     = c("Star","Insl3","Cyp11a1","Hsd3b1"),
  "Macrophage"                      = c("Cd14","C1qb","Csf1r","Adgre1"),
  "Endothelial Cell"                = c("Pecam1","Cdh5","Tek","Ncr1")
)
mk10 <- lapply(mk10, function(gs) gs[gs %in% rownames(obj)])
mk10 <- mk10[sapply(mk10, length) > 0]
gene10 <- unname(unlist(mk10))
cat(sprintf("\nFinal marker genes: %d total\n", length(gene10)))

cat("Computing final DotPlot...\n")
dp10 <- DotPlot(obj, features=gene10, group.by="celltype_final")$data
p10  <- make_dotplot(dp10, mk10, final_colors, final_order,
                     title="Exos v2 - Final Annotation (10 Cell Types)")
n10 <- length(gene10)
w10 <- max(17, n10 * 0.40 + 4)
ggsave(file.path(A_DIR,"exos_v2_final_dotplot.pdf"), p10, width=w10, height=9)
ggsave(file.path(A_DIR,"exos_v2_final_dotplot.png"), p10, width=w10, height=9, dpi=220)
cat(sprintf("Saved: exos_v2_final_dotplot (%.1f x 9 in)\n", w10))

# ══════════════════════════════════════════════════════════════
# Save objects
# ══════════════════════════════════════════════════════════════
cat("\n=== Saving annotated objects ===\n")
saveRDS(obj,     file.path(RDS_DIR,"seurat_exos_annotated_final.rds"))
saveRDS(spg_obj, file.path(RDS_DIR,"seurat_exos_spg_subclustered.rds"))
cat("Saved: seurat_exos_annotated_final.rds\n")
cat("Saved: seurat_exos_spg_subclustered.rds\n")

# Final summary
ann_sum <- obj@meta.data %>%
  group_by(celltype_final, seurat_clusters) %>%
  summarise(n=n(), .groups="drop") %>%
  arrange(celltype_final, as.integer(as.character(seurat_clusters)))
write.csv(ann_sum, file.path(A_DIR,"final_annotation_summary.csv"), row.names=FALSE)

cat("\n===== DONE =====\n")
cat("Output:", A_DIR, "\n\n")
cat("Final cell type distribution:\n")
print(table(obj$celltype_final))

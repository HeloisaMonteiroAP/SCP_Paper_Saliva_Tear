#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(monocle3)
  library(SingleCellExperiment)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggVennDiagram)
  library(UpSetR)
})

set.seed(123)

# ============================================================
# 0) PATHS:
# ============================================================

expr_file   <- "/home/nilson.coimbra/saliva/Patients_saliva_SCP_baseclean.tsv"
gene_file   <- "/home/nilson.coimbra/saliva/annotation.tsv"
cell_file   <- "/home/nilson.coimbra/saliva/metadata.tsv"
umap_file   <- "/home/nilson.coimbra/saliva/Saliva_patients_UMAP_coordinates_with_leiden.tsv"

choi_file   <- "/home/nilson.coimbra/saliva/choi-CANL_deseq2_DE.tsv"
tissue_file <- "/home/nilson.coimbra/saliva/tissue-oscc-results.ALL.tsv"

outdir <- "/home/nilson.coimbra/saliva/trajectory_analysis_FDR_latest"
prefix <- "SALIVA_integrated"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
for (d in c("umap","traj","genes","overlap")) dir.create(file.path(outdir, d), showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1) PARAMETERS:
# ============================================================

condition_col  <- "Condition"
sample_col     <- "Sample"
root_condition <- "Healthy_donor"

SALIVA_TOP_N    <- 230
SALIVA_TOP_PLOT <- 20
MIN_DETECT_FREQ <- 0.05
SALIVA_FDR <- 0.10  

CHOI_LOG2FC  <- 1
CHOI_FDR     <- 0.10
CHOI_P       <- 0.05

TISSUE_LOG2FC <- 1
TISSUE_FDR    <- 0.10

# ============================================================
# 2) UTILS:
# ============================================================

msg <- function(...) cat(sprintf(...), "\n")

save_plot <- function(p, path, w = 8, h = 6, dpi = 300) {
  ext <- tolower(tools::file_ext(path))
  
  if (ext %in% c("svg")) {
    if (!requireNamespace("svglite", quietly = TRUE)) {
      stop("Falta o pacote 'svglite'. Instale com install.packages('svglite').")
    }
    ggsave(
      filename = path,
      plot = p,
      width = w, height = h,
      device = svglite::svglite,
      bg = "white"
    )
  } else if (ext %in% c("png")) {
    ggsave(
      filename = path,
      plot = p,
      width = w, height = h,
      dpi = dpi,
      bg = "white"
    )
  } else {
    stop("Extensão não suportada: ", ext)
  }
}

std_gene <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- toupper(x)
  x[x == "" | is.na(x)] <- NA_character_
  x
}

stop_if_missing_cols <- function(df, cols, label) {
  miss <- setdiff(cols, colnames(df))
  if (length(miss) > 0) stop(label, " missing columns: ", paste(miss, collapse=", "))
}

# ============================================================
# 3) MONOCLE3 INPUTS:
# ============================================================

prepare_monocle_inputs <- function(expr_file, gene_file, cell_file,
                                   gene_id_col="Genes", gene_name_col="GeneName",
                                   sample_col="Sample") {
  
  expr_df <- read.csv(expr_file, sep="\t", check.names=FALSE)
  stop_if_missing_cols(expr_df, c("Protein.Group","Protein.Names","Genes"), "expr")
  
  rownames(expr_df) <- expr_df[["Genes"]]
  expr_cols <- setdiff(colnames(expr_df), c("Protein.Group","Protein.Names","Genes"))
  expr_mat <- as.matrix(expr_df[, expr_cols, drop=FALSE])
  mode(expr_mat) <- "numeric"
  expr_mat[is.na(expr_mat)] <- 0
  
  gene_df <- read.csv(gene_file, sep="\t", check.names=FALSE)
  if (!gene_id_col %in% colnames(gene_df)) stop("annotation missing: ", gene_id_col)
  rownames(gene_df) <- gene_df[[gene_id_col]]
  
  if (!"gene_short_name" %in% colnames(gene_df)) {
    if (!is.null(gene_name_col) && gene_name_col %in% colnames(gene_df)) {
      gene_df$gene_short_name <- ifelse(is.na(gene_df[[gene_name_col]]) | gene_df[[gene_name_col]]=="",
                                        gene_df[[gene_id_col]], gene_df[[gene_name_col]])
    } else {
      gene_df$gene_short_name <- gene_df[[gene_id_col]]
    }
  }
  
  gene_df <- gene_df[rownames(expr_mat), , drop=FALSE]
  
  cell_df <- read.csv(cell_file, sep="\t", check.names=FALSE)
  if (!sample_col %in% colnames(cell_df)) stop("metadata missing: ", sample_col)
  
  ids_meta <- as.character(cell_df[[sample_col]])
  common_ids <- intersect(colnames(expr_mat), ids_meta)
  if (length(common_ids) == 0) stop("No common cell IDs between expr and metadata.")
  
  expr_mat <- expr_mat[, common_ids, drop=FALSE]
  cell_df  <- cell_df[match(common_ids, ids_meta), , drop=FALSE]
  rownames(cell_df) <- common_ids
  
  stopifnot(identical(colnames(expr_mat), rownames(cell_df)))
  
  list(expr_mat=expr_mat, gene_df=gene_df, cell_df=cell_df)
}

msg("==== [1] Loading expr/gene/meta ====")
parsed <- prepare_monocle_inputs(expr_file, gene_file, cell_file,
                                 gene_id_col="Genes", gene_name_col="GeneName",
                                 sample_col=sample_col)
expr_mat <- parsed$expr_mat
gene_df  <- parsed$gene_df
cell_df  <- parsed$cell_df
msg("Expression: %d features x %d cells", nrow(expr_mat), ncol(expr_mat))

msg("==== [2] Loading external UMAP ====")
umap_ext <- read.delim(umap_file, sep="\t", check.names=FALSE)
stop_if_missing_cols(umap_ext, c("UMAP1","UMAP2"), "umap_file")
stopifnot(nrow(umap_ext) == ncol(expr_mat))
umap_mat <- as.matrix(umap_ext[, c("UMAP1","UMAP2")])

msg("==== [3] Building cds + graph ====")
cds <- new_cell_data_set(expr_mat, cell_metadata=cell_df, gene_metadata=gene_df)

# bring external cols if exist
if ("Leiden" %in% colnames(umap_ext)) colData(cds)$leiden_ext <- umap_ext$Leiden
if ("Condition" %in% colnames(umap_ext) && !(condition_col %in% colnames(colData(cds)))) {
  colData(cds)[[condition_col]] <- umap_ext$Condition
}

cds <- preprocess_cds(cds, num_dim=30)
SingleCellExperiment::reducedDims(cds)$UMAP <- umap_mat
cds <- cluster_cells(cds, reduction_method="UMAP")
cds <- learn_graph(cds, use_partition=TRUE)

if (!(condition_col %in% colnames(colData(cds)))) stop("Condition column not found in cds: ", condition_col)

root_cells <- colnames(cds)[colData(cds)[[condition_col]] == root_condition]
cds <- tryCatch(
  order_cells(cds, reduction_method="UMAP", root_cells=root_cells),
  error=function(e) {
    warning("order_cells(root_cells) failed: ", e$message, " -> running automatic order_cells()")
    order_cells(cds, reduction_method="UMAP")
  }
)
pt <- tryCatch(monocle3::pseudotime(cds), error=function(e) NULL)
if (is.null(pt) || length(pt) != ncol(cds)) stop("Pseudotime failed.")
colData(cds)$pseudotime <- pt

saveRDS(cds, file.path(outdir, paste0(prefix, "_cds.rds")))
msg("Saved cds: %s", file.path(outdir, paste0(prefix, "_cds.rds")))

# ============================================================
# 4) MONOCLE CLASSIC PLOTS:
# ============================================================

msg("==== [4] Monocle plots ====")
colData(cds)[[condition_col]] <- as.factor(colData(cds)[[condition_col]])


p_traj_cond <- plot_cells(
  cds,
  color_cells_by = condition_col,
  group_cells_by = condition_col,
  show_trajectory_graph = TRUE,
  label_branch_points = TRUE,
  label_roots = TRUE,
  label_leaves = TRUE,
  cell_size = 1.6
)+ theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

save_plot(p_traj_cond, file.path(outdir, "umap", paste0(prefix, "_traj_condition.svg")), 8, 6)

p_traj_pt <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_branch_points = TRUE,
  label_roots = TRUE,
  label_leaves = TRUE,
  cell_size = 1.6
) + theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

save_plot(p_traj_pt, file.path(outdir, "umap", paste0(prefix, "_traj_pseudotime.svg")), 8, 6)

p_cells_only_cond <- plot_cells(
  cds,
  color_cells_by = condition_col,
  show_trajectory_graph = FALSE,
  cell_size = 1.6
) + theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
save_plot(p_cells_only_cond, file.path(outdir, "umap", paste0(prefix, "_cells_only_condition.svg")), 8, 6)

# ============================================================
# 5) PSEUDOTIME MAKERS:
# ============================================================

msg("==== [5] Spearman correlation vs pseudotime (+ FDR) ====")

do_corr_vs_pseudotime <- function(cds, min_detect=0.05) {
  pt <- colData(cds)$pseudotime
  
  expr_cells_genes <- Matrix::t(monocle3::normalized_counts(cds)) # cells x genes
  expr_cells_genes <- as.matrix(expr_cells_genes)
  
  det_freq <- colMeans(expr_cells_genes > 0)
  keep <- which(det_freq >= min_detect)
  expr_cells_genes <- expr_cells_genes[, keep, drop=FALSE]
  genes_kept <- colnames(expr_cells_genes)
  
  rho <- numeric(length(genes_kept))
  pval <- numeric(length(genes_kept))
  
  for (i in seq_along(genes_kept)) {
    g <- expr_cells_genes[, i]
    suppressWarnings({
      ct <- suppressWarnings(cor.test(g, pt, method="spearman"))
    })
    rho[i] <- unname(ct$estimate)
    pval[i] <- ct$p.value
  }
  
  res <- tibble(
    feature_id = genes_kept,
    rho = rho,
    rho_abs = abs(rho),
    p_value = pval,
    fdr = p.adjust(pval, method="fdr")
  )
  
  if ("gene_short_name" %in% colnames(rowData(cds))) {
    ginfo <- as.data.frame(rowData(cds)[genes_kept, , drop=FALSE])
    res$gene_short_name <- as.character(ginfo$gene_short_name)
  } else {
    res$gene_short_name <- genes_kept
  }
  
  res %>% arrange(desc(rho_abs))
}

corr_df <- do_corr_vs_pseudotime(cds, min_detect=MIN_DETECT_FREQ)
write_tsv(corr_df, file.path(outdir, "overlap", paste0(prefix, "_genes_vs_pseudotime_spearman.tsv")))

saliva_sig <- corr_df %>%
  filter(!is.na(gene_short_name),
         !is.na(fdr),
         fdr <= SALIVA_FDR) %>%
  arrange(desc(rho_abs))

write_tsv(saliva_sig, file.path(outdir, "overlap", paste0(prefix, "_Saliva_FDR", SALIVA_FDR, "_significant.tsv")))

msg("Saliva corr total: %d genes testados", nrow(corr_df))
msg("Saliva sig (FDR<=%.3f): %d genes", SALIVA_FDR, nrow(saliva_sig))

if (nrow(saliva_sig) >= SALIVA_TOP_N) {
  saliva_topN <- saliva_sig %>% slice_head(n = SALIVA_TOP_N)
  msg("Usando saliva_topN dentro de saliva_sig (FDR<=%.3f): Top%d", SALIVA_FDR, SALIVA_TOP_N)
} else {
  saliva_topN <- corr_df %>% slice_head(n = SALIVA_TOP_N)
  warning(sprintf("Poucos genes em saliva_sig (n=%d) para Top%d. Fallback: usando Top%d por rho_abs sem FDR.",
                  nrow(saliva_sig), SALIVA_TOP_N, SALIVA_TOP_N))
}

write_tsv(saliva_topN, file.path(outdir, "overlap", paste0(prefix, "_Saliva_Top", SALIVA_TOP_N, "_traj_definers.tsv")))

genes_saliva <- std_gene(saliva_topN$gene_short_name) %>% na.omit() %>% unique()
SALIVA_TOP_N_USED <- nrow(saliva_topN)
msg("SalivaTop%d_USED=%d", SALIVA_TOP_N, SALIVA_TOP_N_USED)

# ============================================================
# 6) CHOI (Tissue) + PURAM (Tissue) FILTER:
# ============================================================

msg("==== [6] Load CHOI/Puram and filter ====")

choi <- read_tsv(choi_file, show_col_types = FALSE)
stop_if_missing_cols(choi, c("Gene","log2FoldChange","pvalue"), "choi")

choi2 <- choi %>%
  transmute(
    Gene   = std_gene(Gene),
    log2FC = as.numeric(log2FoldChange),
    p_raw  = as.numeric(pvalue),
    padj   = if ("padj" %in% colnames(choi)) as.numeric(padj) else NA_real_
  )

choi_sig <- if ("padj" %in% colnames(choi) && any(!is.na(choi2$padj))) {
  choi2 %>% filter(!is.na(Gene), !is.na(log2FC), abs(log2FC) >= CHOI_LOG2FC, !is.na(padj), padj <= CHOI_FDR)
} else {
  choi2 %>% filter(!is.na(Gene), !is.na(log2FC), abs(log2FC) >= CHOI_LOG2FC, !is.na(p_raw), p_raw <= CHOI_P)
}
genes_choi <- unique(na.omit(choi_sig$Gene))

tissue <- read_tsv(tissue_file, show_col_types = FALSE)
stop_if_missing_cols(tissue, c("Gene","log2FC_oscc_vs_control","p_control_vs_oscc"), "tissue")

tissue2 <- tissue %>%
  transmute(
    Gene   = std_gene(Gene),
    log2FC = as.numeric(log2FC_oscc_vs_control),
    p_raw  = as.numeric(p_control_vs_oscc)
  ) %>%
  mutate(FDR = p.adjust(p_raw, method="fdr"))

tissue_sig <- tissue2 %>%
  filter(!is.na(Gene), !is.na(log2FC), abs(log2FC) >= TISSUE_LOG2FC, !is.na(FDR), FDR <= TISSUE_FDR)

genes_tissue <- unique(na.omit(tissue_sig$Gene))

write_tsv(choi_sig,   file.path(outdir, "overlap", "CHOI_significant.tsv")) #CHOI
write_tsv(tissue_sig, file.path(outdir, "overlap", "Tissue_significant.tsv")) #Puram

msg("Counts: SalivaTop%d=%d | CHOI_sig=%d | Tissue_sig=%d",
    SALIVA_TOP_N, length(genes_saliva), length(genes_choi), length(genes_tissue))

# ============================================================
# 7) OVERLAP + VENN + TABLES (inclui all_three):
# ============================================================

msg("==== [7] Overlap + Venn ====")

choi_saliva   <- intersect(genes_choi, genes_saliva)
tissue_saliva <- intersect(genes_tissue, genes_saliva)
all_three     <- Reduce(intersect, list(genes_choi, genes_tissue, genes_saliva))

write_tsv(tibble(Gene = choi_saliva),   file.path(outdir, "overlap", "CHOI_SALIVA_overlap.tsv"))
write_tsv(tibble(Gene = tissue_saliva), file.path(outdir, "overlap", "TISSUE_SALIVA_overlap.tsv"))
write_tsv(tibble(Gene = all_three),     file.path(outdir, "overlap", "CHOI_TISSUE_SALIVA_overlap.tsv"))

venn_list <- list(CHOI=genes_choi, Bulk_Saliva=genes_tissue, SCP_Saliva=genes_saliva)


upset_list <- list(
  CHOI   = genes_choi,
  Tissue = genes_tissue,
  Saliva = genes_saliva
)

upset_input <- UpSetR::fromList(upset_list)

svg(file.path(outdir, "overlap", "UpSetR_CHOITissueSaliva.svg"), width = 10, height = 5)
print(UpSetR::upset(upset_input, sets = c("CHOI","Tissue","Saliva"), order.by = "freq"))
dev.off()  

# ============================================================
# 7.1) UPSET PLOT (CHOI / Puram / Saliva):
# ============================================================

all_genes <- sort(unique(c(genes_choi, genes_tissue, genes_saliva)))

upset_df <- tibble::tibble(
  Gene   = all_genes,
  CHOI   = all_genes %in% genes_choi,
  Tissue = all_genes %in% genes_tissue,
  Saliva = all_genes %in% genes_saliva
)

p_upset <- ComplexUpset::upset(
  upset_df,
  intersect = c("CHOI", "Tissue", "Saliva"),
  name = "Genes",
  width_ratio = 0.15
) +
  ggplot2::ggtitle(paste0("UpSet: CHOI vs Puram vs Saliva (Saliva FDR<= ", SALIVA_FDR, ", n=", length(genes_saliva), ")"))

save_plot(p_upset, file.path(outdir, "overlap", "UpSet_CHOITissueSaliva.svg"), 8.5, 5.2)


# ============================================================
# 8) TOP20 PLOTS (MONOCLE3):
# ============================================================

msg("==== [8] Gene panels (Monocle) ====")

gene_to_feature <- function(corr_df, gene_symbols_std) {
  corr_df %>%
    mutate(GENE_STD = std_gene(gene_short_name)) %>%
    filter(GENE_STD %in% gene_symbols_std) %>%
    pull(feature_id) %>% unique()
}

topN_within_set <- function(corr_df, gene_symbols_std, n=20) {
  corr_df %>%
    mutate(GENE_STD = std_gene(gene_short_name)) %>%
    filter(GENE_STD %in% gene_symbols_std) %>%
    arrange(desc(rho_abs)) %>%
    slice_head(n=n)
}

set_saliva_top20 <- saliva_topN %>% slice_head(n = SALIVA_TOP_PLOT)

choi_saliva_top20   <- topN_within_set(corr_df, choi_saliva,   n = SALIVA_TOP_PLOT)
tissue_saliva_top20 <- topN_within_set(corr_df, tissue_saliva, n = SALIVA_TOP_PLOT)
all_three_top20     <- topN_within_set(corr_df, all_three,     n = SALIVA_TOP_PLOT)

write_tsv(set_saliva_top20,       file.path(outdir, "overlap", "Saliva_Top20_traj.tsv"))
write_tsv(choi_saliva_top20,      file.path(outdir, "overlap", "CHOI_Saliva_Top20.tsv"))
write_tsv(tissue_saliva_top20,    file.path(outdir, "overlap", "Puram_Saliva_Top20.tsv"))
write_tsv(all_three_top20,        file.path(outdir, "overlap", "CHOI_Puram_Saliva_Top20.tsv"))

plot_gene_monocle <- function(cds, feature_id, label, out_prefix) {
  if (!feature_id %in% rownames(cds)) return(FALSE)
  
  p1 <- plot_cells(
    cds,
    genes = feature_id,
    show_trajectory_graph = TRUE,
    cell_size = 1.6
  ) + ggtitle(paste0(label, " (UMAP+traj)"))
  
  save_plot(p1, paste0(out_prefix, "_UMAP_traj.svg"), 7.5, 6)
  
  p2 <- plot_genes_in_pseudotime(cds[feature_id, ]) + ggtitle(paste0(label, " (pseudotime)"))
  save_plot(p2, paste0(out_prefix, "_pseudotime.svg"), 7.5, 4.8)
  
  TRUE
}

run_set <- function(set_name, df_top) {
  if (nrow(df_top) == 0) return()
  
  set_dir <- file.path(outdir, "genes", gsub("[^A-Za-z0-9]+","_",set_name))
  dir.create(set_dir, showWarnings=FALSE, recursive=TRUE)
  
  for (i in seq_len(nrow(df_top))) {
    fid <- df_top$feature_id[i]
    gnm <- std_gene(df_top$gene_short_name[i])
    lab <- paste0(set_name, " | ", gnm)
    
    out_pref <- file.path(set_dir, sprintf("%02d_%s", i, gnm))
    plot_gene_monocle(cds, fid, lab, out_pref)
  }
}

run_set(paste0("Saliva_Top", SALIVA_TOP_PLOT), set_saliva_top20)
run_set("CHOI_intersect_Saliva", choi_saliva_top20)
run_set("Tissue_intersect_Saliva", tissue_saliva_top20)
run_set("CHOI_Tissue_Saliva_all3", all_three_top20)
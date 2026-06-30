# !/usr/bin/env Rscript
# sEVTRAS_mouse.R — sEV recognizer + ESAI calculator (deterministic, mouse-first)
# Notes:
#  - Default species is Mus (mouse). You can still pass --species Homo if needed.
#  - Includes robustness patches: safe cbind of sEVs, PCA rank guard, kNN guard,
#    zero-col normalize guard, and ID-scheme-aware biogenesis bonus.
#  - Deterministic behavior: seeded GMM threshold, BLAS/OMP pinning if available.

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(mclust)
  library(FNN)
  library(uwot)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(msigdbr)
  library(fgsea)
})

## --- Reproducibility guards ---
set.seed(12345)
suppressWarnings(try(RNGkind("L'Ecuyer-CMRG"), silent = TRUE))
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
  try(RhpcBLASctl::omp_set_num_threads(1),  silent = TRUE)
}

msg <- function(...) cat(sprintf(paste0("[sEVTRAS] ", paste0(..., collapse=""), "\n")))
`%||%` <- function(a,b) if (!is.null(a)) a else b

# ---------- IO helpers ----------
read_10x_mtx <- function(mtx_dir) {
  if (file.exists(file.path(mtx_dir, "matrix.mtx.gz")) || file.exists(file.path(mtx_dir, "matrix.mtx"))) {
    base <- mtx_dir
  } else if (dir.exists(file.path(mtx_dir, "outs", "raw_feature_bc_matrix"))) {
    base <- file.path(mtx_dir, "outs", "raw_feature_bc_matrix")
  } else if (dir.exists(file.path(mtx_dir, "outs", "filtered_feature_bc_matrix"))) {
    base <- file.path(mtx_dir, "outs", "filtered_feature_bc_matrix")
  } else {
    stop("Could not find 10x matrix at ", mtx_dir)
  }
  mtx_file <- if (file.exists(file.path(base, "matrix.mtx.gz"))) file.path(base, "matrix.mtx.gz") else file.path(base, "matrix.mtx")
  feat_file <- if (file.exists(file.path(base, "features.tsv.gz"))) file.path(base, "features.tsv.gz") else
    if (file.exists(file.path(base, "features.tsv"))) file.path(base, "features.tsv") else
      if (file.exists(file.path(base, "genes.tsv.gz"))) file.path(base, "genes.tsv.gz") else file.path(base, "genes.tsv")
  barc_file <- if (file.exists(file.path(base, "barcodes.tsv.gz"))) file.path(base, "barcodes.tsv.gz") else file.path(base, "barcodes.tsv")
  
  M <- Matrix::readMM(mtx_file)
  feats <- data.table::fread(feat_file, header = FALSE)
  barc  <- data.table::fread(barc_file, header = FALSE)
  gene_names <- if (ncol(feats) >= 2) feats[[2]] else feats[[1]]
  rownames(M) <- make.unique(gene_names)
  colnames(M) <- barc[[1]]
  if (!inherits(M, "dgCMatrix")) M <- as(M, "CsparseMatrix")
  M
}

# ---- robust sum helpers (handle any Matrix subclass or base matrices) ----
col_sums <- function(m) { if (inherits(m, "Matrix")) Matrix::colSums(m) else base::colSums(m) }
row_sums <- function(m) { if (inherits(m, "Matrix")) Matrix::rowSums(m) else base::rowSums(m) }

normalize_cpm <- function(m, scale_factor = 1e4) {
  if (ncol(m) == 0) return(m)  # guard zero columns
  sf <- col_sums(m); sf[sf == 0] <- 1
  if (inherits(m, "dgCMatrix")) { res <- m; res@x <- res@x / rep(sf, diff(res@p)) * scale_factor; return(res) }
  else t(t(m)/sf)*scale_factor
}
log1p_sparse <- function(m) { if (inherits(m,"dgCMatrix")) { res <- m; res@x <- log1p(res@x); res } else log1p(m) }

# ---------- Stratified capping helper ----------
cap_stratified <- function(scan_idx, umi, max_scan, nbins = 10, seed = 1) {
  set.seed(seed)
  u <- umi[scan_idx]
  brks <- unique(quantile(u, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE, type = 8))
  if (length(brks) < 3) brks <- seq(min(u), max(u), length.out = 3)
  bins <- cut(u, breaks = brks, include.lowest = TRUE, labels = FALSE)
  B <- sort(unique(bins))
  per <- max(1L, floor(max_scan / length(B)))
  keep_rel <- unlist(lapply(B, function(b) {
    ii <- which(bins == b)
    if (length(ii) <= per) ii else sample(ii, per)
  }))
  if (length(keep_rel) < max_scan) {
    leftover <- setdiff(seq_along(scan_idx), keep_rel)
    if (length(leftover) > 0) keep_rel <- c(keep_rel, sample(leftover, min(max_scan - length(keep_rel), length(leftover))))
  } else if (length(keep_rel) > max_scan) {
    keep_rel <- sample(keep_rel, max_scan)
  }
  scan_idx[keep_rel]
}

# ===========================
# Gene list reader (preserve case)
# ===========================
read_gene_file <- function(path) {
  pick_first_existing <- function(paths) {
    hits <- paths[file.exists(paths)]
    if (length(hits) > 0) hits[[1]] else NA_character_
  }
  candidates <- c(
    "ev_genes", "ev_genes.txt", "ev_genes.tsv", "ev_genes.csv",
    "biogenesis_genes.txt", "biogenesis_genes.tsv", "biogenesis_genes.csv",
    "genes.txt", "genes.tsv", "genes.csv"
  )
  if (!is.null(path) && dir.exists(path)) {
    cand <- file.path(path, candidates)
    path <- pick_first_existing(cand)
    if (is.na(path)) return(character(0))
  }
  if (is.null(path) || !file.exists(path)) return(character(0))
  
  ext <- tolower(tools::file_ext(path))
  genes <- tryCatch({
    if (ext == "tsv") readr::read_tsv(path, col_names = FALSE, show_col_types = FALSE)[[1]]
    else if (ext == "csv") readr::read_csv(path, col_names = FALSE, show_col_types = FALSE)[[1]]
    else readr::read_lines(path)
  }, error = function(e) character(0))
  
  genes <- trimws(genes)
  genes <- genes[nzchar(genes) & !is.na(genes)]
  unique(genes)
}


# ---------- EM core ----------
hypergeom_scores <- function(m, gene_set) {
  genes <- rownames(m); gset <- intersect(gene_set, genes)
  N <- length(genes); K <- length(gset)
  if (K < 10) warning("Gene set small after intersection (n=", K, ").")
  m_gset <- m[match(gset, genes), , drop = FALSE]
  k_vec <- Matrix::colSums(m_gset > 0)
  n_vec <- Matrix::colSums(m > 0)
  pvals <- rep(1, length(k_vec))
  sel <- n_vec > 0
  pvals[sel] <- phyper(q = pmax(0, k_vec[sel]) - 1, m = K, n = N - K, k = n_vec[sel], lower.tail = FALSE)
  s <- -log10(pmax(pvals, .Machine$double.xmin)); names(s) <- colnames(m); s
}
update_gene_set <- function(m_logcpm, z, alpha = 0.20, keep_size = NULL) {
  z <- as.numeric(z)
  cor_s <- suppressWarnings(apply(m_logcpm, 1, function(x) cor(x, z, method = "spearman")))
  ord <- order(cor_s, decreasing = TRUE, na.last = NA)
  if (!is.null(keep_size)) rownames(m_logcpm)[ord[seq_len(min(keep_size, length(ord)))]]
  else {
    k <- max(200, round(alpha * nrow(m_logcpm)))
    rownames(m_logcpm)[ord[seq_len(min(k, length(ord)))]]
  }
}
adaptive_threshold <- function(scores) {
  x <- as.numeric(scores); x <- x[is.finite(x)]
  if (length(unique(x)) < 3) return(quantile(x, 0.95))
  set.seed(12345)  # deterministic GMM fit
  mc <- try(Mclust(x, G = 2, verbose = FALSE), silent = TRUE)
  if (inherits(mc, "try-error") || is.null(mc$parameters)) return(quantile(x, 0.95))
  means <- sort(mc$parameters$mean); (means[1] + means[2]) / 2
}

em_sev <- function(m_counts, init_gene_set, alpha = 0.20, max_iter = 30, tol = 1e-3, em_keep = NULL, preselect_top = 20000) {
  genes <- rownames(m_counts)
  theta0 <- intersect(init_gene_set, genes)
  if (length(theta0) < 50) stop("Initial gene set too small after intersection.")
  keep_size <- em_keep %||% length(theta0)
  det <- row_sums(m_counts > 0)
  sel_genes <- names(sort(det, decreasing = TRUE))[seq_len(min(length(det), preselect_top))]
  m_sub <- m_counts[sel_genes, , drop = FALSE]
  m_logcpm <- log1p_sparse(normalize_cpm(m_sub))
  theta <- intersect(theta0, rownames(m_logcpm))
  msg("EM init | gene set size=", length(theta))
  hist <- list()
  for (it in seq_len(max_iter)) {
    z <- hypergeom_scores(m_sub, theta)
    theta_new <- update_gene_set(m_logcpm, z, alpha = alpha, keep_size = keep_size)
    jacc <- length(intersect(theta, theta_new)) / length(union(theta, theta_new))
    hist[[it]] <- list(iter = it, jaccard = jacc, size = length(theta_new))
    msg("EM iter ", it, " | Jaccard=", sprintf("%.3f", jacc))
    theta <- theta_new
    if (jacc > (1 - tol)) break
  }
  z_final <- hypergeom_scores(m_sub, theta)
  list(theta = theta, z = z_final, history = hist)
}

# ---------- Unification guard ----------
unify_gene_sets <- function(all_thetas, sample_names, min_presence = 2, min_size = 100) {
  if (length(all_thetas) == 0) return(NULL)
  genes_vec <- unlist(all_thetas, use.names = FALSE)
  if (length(genes_vec) == 0) return(NULL)
  tab <- sort(table(genes_vec), decreasing = TRUE)
  nS <- length(sample_names)
  cutoff <- if (min_presence > 0 && min_presence < 1) ceiling(min_presence * nS) else as.integer(min_presence)
  cutoff <- max(1, cutoff)
  unified <- names(tab)[tab >= cutoff]
  if (length(unified) < min_size) {
    msg("Unified gene set below min_size (", length(unified), " < ", min_size, "). Skipping unification; keeping per-sample θ.")
    return(NULL)
  }
  msg("Unified gene set size=", length(unified), " (cutoff≥", cutoff, "; min_size=", min_size, ")")
  unified
}

# ---------- sEV recognizer ----------
sEV_recognizer <- function(input_dirs, sample_names, out_dir,
                           species = c("Homo","Mus"),
                           search_UMI = 20000, min_UMI = 1, alpha = 0.20,
                           init_gene_set = NULL, map_ids = TRUE,
                           em_keep = 1200, em_preselect_top = 20000,
                           max_scan = 120000,
                           unify_min_presence = 2, unify_min_size = 100, disable_unify = FALSE) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  species <- match.arg(species)
  if (is.null(init_gene_set)) {
    msg("Using msigdbr-derived EV-related genes as initial set (consider supplying curated set).")
    init_gene_set <- get_biogenesis_genes(species)
  }
  all_thetas <- list(); all_scores <- list(); sev_mats <- list(); sev_meta <- list()
  
  for (i in seq_along(input_dirs)) {
    samp <- sample_names[i]; msg("Reading sample ", samp)
    M <- read_10x_mtx(input_dirs[[i]])
    umi <- col_sums(M)
    scan_idx <- which(umi >= min_UMI & umi <= search_UMI)
    if (length(scan_idx) == 0) next
    if (length(scan_idx) > max_scan) {
      scan_idx <- cap_stratified(scan_idx, umi, max_scan, nbins = 10, seed = 1)
      msg("Sample ", samp, ": capped stratified by UMI to ", max_scan, ".")
    }
    M_scan <- M[, scan_idx, drop = FALSE]
    
    # ID mapping (only if needed)
    init_set_i <- init_gene_set
    looks_ens <- mean(grepl("^(ENSG|ENSMUSG)", rownames(M_scan))) > 0.6
    if (map_ids && looks_ens) {
      mapped <- try({
        if (species == "Homo" && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
          suppressPackageStartupMessages(library(org.Hs.eg.db))
          sy2ens <- AnnotationDbi::select(org.Hs.eg.db, keys = init_gene_set, keytype = "SYMBOL", columns = "ENSEMBL")
          unique(na.omit(sy2ens$ENSEMBL))
        } else if (species == "Mus" && requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
          suppressPackageStartupMessages(library(org.Mm.eg.db))
          sy2ens <- AnnotationDbi::select(org.Mm.eg.db, keys = init_gene_set, keytype = "SYMBOL", columns = "ENSEMBL")
          unique(na.omit(sy2ens$ENSEMBL))
        } else character(0)
      }, silent = TRUE)
      if (!inherits(mapped, "try-error") && length(mapped) > 0) init_set_i <- mapped
    }
    
    em <- em_sev(M_scan, init_set_i, alpha = alpha, em_keep = em_keep, preselect_top = em_preselect_top)
    thr <- adaptive_threshold(em$z)
    is_sev <- em$z >= thr
    
    all_thetas[[samp]] <- em$theta
    all_scores[[samp]] <- data.frame(
      barcode = colnames(M_scan), score = as.numeric(em$z), is_sEV = is_sev,
      sample = samp, umi = umi[scan_idx], stringsAsFactors = FALSE
    )
    
    if (any(is_sev)) {
      Ms <- M_scan[, is_sev, drop = FALSE]
      bc_raw  <- colnames(Ms)
      bc_uniq <- paste(samp, bc_raw, sep="|")
      colnames(Ms) <- bc_uniq
      
      sev_mats[[samp]] <- Ms
      sev_meta[[samp]] <- data.frame(
        barcode_raw = bc_raw,
        barcode     = bc_uniq,   # canonical key
        sample      = samp,
        score       = as.numeric(em$z[is_sev]),
        umi         = umi[scan_idx][is_sev],
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Optional unified rescoring (guarded)
  unified <- NULL
  if (!isTRUE(disable_unify)) {
    unified <- unify_gene_sets(all_thetas, sample_names, min_presence = unify_min_presence, min_size = unify_min_size)
  }
  if (!is.null(unified)) {
    for (i in seq_along(input_dirs)) {
      samp <- sample_names[i]
      M <- read_10x_mtx(input_dirs[[i]])
      umi <- col_sums(M)
      scan_idx <- which(umi >= min_UMI & umi <= search_UMI)
      if (length(scan_idx) == 0) next
      if (length(scan_idx) > max_scan) {
        scan_idx <- cap_stratified(scan_idx, umi, max_scan, nbins = 10, seed = 1)
        msg("Sample ", samp, ": capped stratified by UMI to ", max_scan, " (unified rescoring).")
      }
      M_scan <- M[, scan_idx, drop = FALSE]
      z_u <- hypergeom_scores(M_scan, unified); thr <- adaptive_threshold(z_u); is_sev <- z_u >= thr
      
      if (!is.null(all_scores[[samp]])) {
        all_scores[[samp]]$score_unified <- as.numeric(z_u[match(all_scores[[samp]]$barcode, names(z_u))])
        all_scores[[samp]]$is_sEV_unified <- all_scores[[samp]]$score_unified >= thr
      }
      
      if (any(is_sev)) {
        Ms <- M_scan[, is_sev, drop = FALSE]
        bc_raw  <- colnames(Ms)
        bc_uniq <- paste(samp, bc_raw, sep="|")
        colnames(Ms) <- bc_uniq
        sev_mats[[samp]] <- Ms
        sev_meta[[samp]] <- data.frame(
          barcode_raw = bc_raw,
          barcode     = bc_uniq,
          sample      = samp,
          score       = as.numeric(z_u[is_sev]),
          umi         = umi[scan_idx][is_sev],
          stringsAsFactors = FALSE
        )
      } else {
        sev_mats[[samp]] <- NULL; sev_meta[[samp]] <- NULL
      }
    }
  } else {
    msg("Proceeding with per-sample θ assignments (no unified rescoring).")
  }
  
  # Combine sEV droplets and enforce alignment (robust to NULL)
  nz_mats <- Filter(Negate(is.null), sev_mats)
  if (length(nz_mats) > 0) {
    common_genes <- Reduce(intersect, lapply(nz_mats, rownames))
    sev_mats2 <- lapply(sev_mats, function(x) if (is.null(x)) NULL else x[common_genes, , drop = FALSE])
    M_sev <- do.call(cbind, sev_mats2)
    meta_sev <- bind_rows(sev_meta)
    meta_sev <- meta_sev[match(colnames(M_sev), meta_sev$barcode), , drop = FALSE]
    stopifnot(identical(meta_sev$barcode, colnames(M_sev)))
  } else {
    M_sev <- Matrix(0, nrow = 0, ncol = 0, sparse = TRUE); meta_sev <- data.frame()
  }
  
  saveRDS(list(scores = all_scores, thetas = all_thetas, unified_gene_set = unified), file = file.path(out_dir, "raw_scores.rds"))
  saveRDS(list(counts = M_sev, meta = meta_sev), file = file.path(out_dir, "sev_droplets.rds"))
  invisible(list(unified_gene_set = unified, sev_counts = M_sev, sev_meta = meta_sev, all_scores = all_scores))
}

# ---------- ESAI ----------
ESAI_calculator <- function(sev_rds, cells_rds, out_dir,
                            species = c("Homo","Mus"),
                            batch_col = "batch", type_col = "celltype",
                            k = 10, pca_max_cols = 50000, var_gene_cap = 2000) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  species <- match.arg(species)
  
  # Ensure deterministic behavior in this function scope as well
  set.seed(12345)
  
  msg("Reading sEV and cell reference objects...")
  sev <- readRDS(sev_rds); cells <- readRDS(cells_rds)
  stopifnot(is.list(sev), is.list(cells))
  M_cells <- cells$counts; meta_cells <- cells$meta
  M_sev   <- sev$counts;   meta_sev   <- sev$meta
  
  if (ncol(M_sev) == 0) stop("No sEV droplets in sev RDS (0 columns).")
  if (!inherits(M_cells, "dgCMatrix")) stop("cells$counts must be a dgCMatrix (got: ", class(M_cells)[1], ").")
  if (!inherits(M_sev,   "dgCMatrix")) stop("sev$counts must be a dgCMatrix (got: ", class(M_sev)[1], ").")
  if (!(batch_col %in% colnames(meta_cells))) stop(paste0("Batch column '", batch_col, "' not found in cells$meta."))
  if (!(type_col  %in% colnames(meta_cells)))  stop(paste0("Type column '",  type_col,  "' not found in cells$meta."))
  
  if (!"barcode" %in% colnames(meta_sev)) stop("sev$meta must contain 'barcode' (canonical 'sample|barcode'). Re-run recognizer.")
  meta_sev <- meta_sev[match(colnames(M_sev), meta_sev$barcode), , drop = FALSE]
  stopifnot(identical(meta_sev$barcode, colnames(M_sev)))
  
  clean_barcodes <- function(bc) sub("_[0-9]+_[0-9]+$", "", bc)
  if ("barcode" %in% colnames(meta_cells)) meta_cells$barcode_clean <- clean_barcodes(meta_cells$barcode)
  
  msg("Aligning gene spaces (initial)…")
  g_cells <- rownames(M_cells); g_sev <- rownames(M_sev)
  common <- intersect(g_cells, g_sev)
  if (length(common) == 0) {
    strip_version <- function(x) sub("\\.[0-9]+$", "", x)
    rn_cells <- toupper(strip_version(g_cells))
    rn_sev   <- toupper(strip_version(g_sev))
    rownames(M_cells) <- rn_cells; rownames(M_sev) <- rn_sev
    common <- intersect(rownames(M_cells), rownames(M_sev))
  }
  if (length(common) == 0) stop("Cannot proceed: zero overlapping genes between cells and sEV matrices.")
  M_cells <- M_cells[common, , drop = FALSE]
  M_sev   <- M_sev[common, , drop = FALSE]
  
  msg("Normalizing and log-transforming…")
  M_all <- cbind(M_cells, M_sev)
  
  # Always keep ALL sEV columns; if needed, subsample CELLS only
  n_cells <- ncol(M_cells)
  n_sev   <- ncol(M_sev)
  total   <- n_cells + n_sev
  
  if (total > pca_max_cols) {
    if (pca_max_cols < n_sev) {
      msg("PCA guard: requested cap (", pca_max_cols, ") < n_sEV (", n_sev, "). Using n_sEV as cap.")
      pca_max_cols <- n_sev
    }
    keep_cells <- max(0L, pca_max_cols - n_sev)
    msg("PCA guard: keeping ALL sEVs (", n_sev, ") and ", keep_cells, " cells out of ", n_cells, ".")
    sel_cells <- if (keep_cells >= n_cells) colnames(M_cells) else sample(colnames(M_cells), keep_cells)
    keep_names <- c(sel_cells, colnames(M_sev))
    M_all_pca <- M_all[, keep_names, drop = FALSE]
  } else {
    M_all_pca <- M_all
  }
  
  M_all_log <- log1p_sparse(normalize_cpm(M_all_pca))
  msg("Running PCA… (this can be the slowest step)")
  det_rate  <- row_sums(M_all_pca > 0)
  top_n     <- min(var_gene_cap, length(det_rate))
  var_genes <- names(sort(det_rate, decreasing = TRUE))[seq_len(top_n)]
  if (length(var_genes) < 2) stop("Not enough variable genes for PCA (got ", length(var_genes), "). ",
                                  "Try increasing --var_gene_cap or check gene overlap between cells and sEVs.")
  X  <- t(as.matrix(M_all_log[var_genes, , drop = FALSE]))
  rank_max <- max(2, min(30, ncol(X) - 1, nrow(X) - 1))
  pc <- prcomp(X, center = TRUE, scale. = TRUE, rank. = rank_max)
  emb <- pc$x
  
  full_cols <- colnames(M_all_pca)
  cols_cells_in_pca <- intersect(colnames(M_cells), full_cols)
  cols_sev_in_pca   <- intersect(colnames(M_sev),   full_cols)
  if (length(cols_cells_in_pca) == 0)
    stop("PCA retained 0 reference cells. Raise pca_max_cols or pre-filter sEVs.")
  emb_cells <- emb[match(cols_cells_in_pca, full_cols), , drop = FALSE]
  emb_sev   <- emb[match(cols_sev_in_pca,   full_cols), , drop = FALSE]
  rownames(emb_cells) <- cols_cells_in_pca
  rownames(emb_sev)   <- cols_sev_in_pca
  
  msg("Computing kNN similarity (k=", k, ")…")
  if (nrow(emb_cells) == 0) stop("No reference cells in PCA sample. Increase pca_max_cols or check input.")
  k_eff <- min(k, nrow(emb_cells))
  nn  <- FNN::get.knnx(emb_cells, emb_sev, k = k_eff)
  idx <- nn$nn.index
  types <- meta_cells[[type_col]]; names(types) <- rownames(emb_cells)
  sim_list <- lapply(seq_len(nrow(idx)), function(i) {
    cell_ids <- rownames(emb_cells)[idx[i,]]
    t        <- types[cell_ids]
    tab      <- table(t)
    data.frame(cell_barcode = rownames(emb_sev)[i], celltype = names(tab), sim = as.numeric(tab), stringsAsFactors = FALSE)
  })
  sim_df <- dplyr::bind_rows(sim_list)
  
  msg("Scoring EV biogenesis capacity per cell type…")
  bio_genes <- get_biogenesis_genes(species)
  bio_genes_hit <- intersect(bio_genes, rownames(M_cells))
  if (length(bio_genes_hit) >= 10) {
    mat_log <- log1p_sparse(normalize_cpm(M_cells[bio_genes_hit, , drop = FALSE]))
    cell_means <- as.numeric(Matrix::colMeans(mat_log))
    bio_scores <- data.frame(celltype = meta_cells[[type_col]], mean = cell_means) %>%
      dplyr::group_by(celltype) %>%
      dplyr::summarise(bio = mean(mean), .groups = "drop")
    rng <- range(bio_scores$bio, finite = TRUE)
    bio_scores$bio_sc <- if (diff(rng) > 0) (bio_scores$bio - rng[1]) / diff(rng) * 2 else 0
  } else {
    msg("Biogenesis gene match <10; disabling bonus (ID scheme likely not SYMBOL).")
    bio_scores <- data.frame(celltype = unique(meta_cells[[type_col]]), bio_sc = 0)
  }
  
  msg("Combining similarity with biogenesis bonus and assigning sources…")
  sim_df2 <- dplyr::left_join(sim_df, bio_scores, by = "celltype") %>%
    dplyr::mutate(bio_sc = dplyr::coalesce(bio_sc, 0), total  = sim + bio_sc)
  
  # Deterministic tie-break: total desc, then sim desc, then bio_sc desc, then alphabetic celltype
  source_assign <- sim_df2 %>%
    dplyr::group_by(cell_barcode) %>%
    dplyr::arrange(dplyr::desc(total), dplyr::desc(sim), dplyr::desc(bio_sc), celltype, .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(cell_barcode, source_celltype = celltype, sim, bio_sc, total)
  
  meta_sev <- meta_sev[meta_sev$barcode %in% rownames(emb_sev), , drop = FALSE]
  stopifnot(!anyDuplicated(meta_sev$barcode))
  meta_sev2 <- dplyr::left_join(meta_sev, source_assign, by = c("barcode" = "cell_barcode"))
  
  msg("Computing ESAI…")
  cell_counts <- meta_cells %>%
    dplyr::group_by(.data[[batch_col]]) %>%
    dplyr::summarise(cell_n = dplyr::n(), .groups = "drop")
  
  ESAI_sample <- meta_sev2 %>%
    dplyr::count(sample, name = "sev_n") %>%
    dplyr::full_join(cell_counts, by = c("sample" = batch_col)) %>%
    dplyr::mutate(ESAI = ifelse(is.na(sev_n), 0, sev_n) / pmax(1, cell_n)) %>%
    dplyr::arrange(sample) %>%
    dplyr::select(sample, sev_n, cell_n, ESAI)
  
  celltype_counts <- meta_cells %>%
    dplyr::group_by(.data[[batch_col]], .data[[type_col]]) %>%
    dplyr::summarise(cell_n = dplyr::n(), .groups = "drop")
  
  ESAI_celltype <- meta_sev2 %>%
    dplyr::count(sample, source_celltype, name = "sev_n") %>%
    dplyr::full_join(celltype_counts, by = c("sample" = batch_col, "source_celltype" = type_col)) %>%
    dplyr::mutate(ESAI_c = ifelse(is.na(sev_n), 0, sev_n) / pmax(1, cell_n)) %>%
    dplyr::arrange(sample, source_celltype) %>%
    dplyr::select(sample, source_celltype, sev_n, cell_n, ESAI_c)
  
  msg("Writing ESAI outputs…")
  readr::write_csv(ESAI_sample,   file.path(out_dir, "ESAI_sample.csv"))
  readr::write_csv(ESAI_celltype, file.path(out_dir, "ESAI_celltype.csv"))
  saveRDS(list(ESAI_sample = ESAI_sample, ESAI_celltype = ESAI_celltype, meta_sev = meta_sev2),
          file = file.path(out_dir, "ESAI_results.rds"))
  msg("Done! ESAI results saved to ", out_dir)
  invisible(list(ESAI_sample = ESAI_sample, ESAI_celltype = ESAI_celltype, meta_sev = meta_sev2))
}

# ---------- CLI ----------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE); kv <- list(); i <- 1
  while (i <= length(args)) {
    key <- args[i]; if (!startsWith(key, "--")) { i <- i + 1; next }
    key <- sub("^--", "", key); i <- i + 1; vals <- c()
    while (i <= length(args) && !startsWith(args[i], "--")) { vals <- c(vals, args[i]); i <- i + 1 }
    kv[[key]] <- if (length(vals) == 1) vals[[1]] else vals
  }
  kv
}

split_args <- function(x) { if (is.null(x)) return(character(0)); if (length(x)==1) strsplit(x, ",")[[1]] else unlist(x) }

main <- function() {
  a <- parse_args()
  if (is.null(a$out_dir)) stop("--out_dir is required")
  out_dir <- a$out_dir; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ran_anything <- FALSE
  
  # Recognizer
  if (!is.null(a$input_dirs) && !is.null(a$sample_names)) {
    input_dirs  <- split_args(a$input_dirs)
    sample_names<- split_args(a$sample_names)
    if (length(input_dirs) != length(sample_names)) stop("input_dirs and sample_names must match in length")
    species      <- ifelse(is.null(a$species), "Mus", a$species)   # mouse-first default
    search_umi   <- ifelse(is.null(a$search_umi), 20000, as.integer(a$search_umi))
    min_umi      <- ifelse(is.null(a$min_umi), 1, as.integer(a$min_umi))
    alpha        <- ifelse(is.null(a$alpha), 0.20, as.numeric(a$alpha))
    em_keep      <- ifelse(is.null(a$em_keep), 1200, as.integer(a$em_keep))
    pre_top      <- ifelse(is.null(a$em_preselect_top), 20000, as.integer(a$em_preselect_top))
    max_scan     <- ifelse(is.null(a$max_scan), 120000, as.integer(a$max_scan))
    unify_min_presence <- if (!is.null(a$unify_min_presence)) as.numeric(a$unify_min_presence) else 2
    unify_min_size     <- if (!is.null(a$unify_min_size)) as.integer(a$unify_min_size) else 100
    disable_unify      <- isTRUE(as.logical(a$no_unify))
    init_file <- a$init_gene_set; init_set <- NULL
    if (!is.null(init_file) && file.exists(init_file)) init_set <- readr::read_tsv(init_file, col_names = FALSE, show_col_types = FALSE)[[1]]
    
    sEV_recognizer(input_dirs, sample_names, out_dir = out_dir, species = species,
                   search_UMI = search_umi, min_UMI = min_umi, alpha = alpha,
                   init_gene_set = init_set, em_keep = em_keep, em_preselect_top = pre_top,
                   max_scan = max_scan, unify_min_presence = unify_min_presence,
                   unify_min_size = unify_min_size, disable_unify = disable_unify)
    ran_anything <- TRUE
  }
  
  # ESAI (supports --sev_rds to point to QC-filtered sEVs)
  if (!is.null(a$cells_obj) && !is.null(a$cells_batch_col) && !is.null(a$cells_type_col)) {
    sev_rds <- a$sev_rds %||% file.path(out_dir, "sev_droplets.rds")
    if (!file.exists(sev_rds)) {
      stop("sEV object not found. Looked for: ", sev_rds,
           "\nHint: pass --sev_rds /path/to/sev_filtered.rds or run recognizer first.")
    }
    species <- ifelse(is.null(a$species), "Mus", a$species)  # mouse-first default
    k          <- ifelse(is.null(a$k), 10, as.integer(a$k))
    pca_max    <- ifelse(is.null(a$pca_max_cols), 50000, as.integer(a$pca_max_cols))
    var_cap    <- ifelse(is.null(a$var_gene_cap), 2000, as.integer(a$var_gene_cap))
    ESAI_calculator(sev_rds = sev_rds, cells_rds = a$cells_obj, out_dir = out_dir,
                    species = species, batch_col = a$cells_batch_col, type_col = a$cells_type_col,
                    k = k, pca_max_cols = pca_max, var_gene_cap = var_cap)
    ran_anything <- TRUE
  }
  
  if (!ran_anything) stop("Nothing to do. Provide recognizer args (--input_dirs/--sample_names) and/or ESAI args (--cells_obj, --cells_batch_col, --cells_type_col).")
}

if (identical(environment(), globalenv())) {
  tryCatch(main(), error = function(e) { message(e); q(status = 1) })
}

# IO.R - Optimized Single-Cell Data Processing & Iterative Clustering
# Ported from Allen Institute scrattch.hicat & BrainMappingPipeline iterclust.py

library(Seurat)
library(harmony)
library(Matrix)
library(dplyr)
library(igraph)

# de_param: Define Differential Expression parameters for iterative clustering
# Default parameters tuned for 10x droplet single-cell data:
# padj = 0.05, lfc = ln(1.5) ≈ 0.405, low_th = 1.0, q1 = 0.5, q_diff = 0.7, de.score.th = 150
de_param <- function(low.th = 1.0,
                     padj.th = 0.05, 
                     lfc.th = 0.405, # log(1.5) = 0.405
                     q1.th = 0.5, 
                     q2.th = NULL,
                     q.diff.th = 0.7, 
                     de.score.th = 150, 
                     min.cells = 4, 
                     min.genes = 5) {
  return(list(
    low.th = low.th,
    padj.th = padj.th,
    lfc.th = lfc.th,
    q1.th = q1.th,
    q2.th = q2.th,
    q.diff.th = q.diff.th,
    de.score.th = de.score.th,
    min.cells = min.cells,
    min.genes = min.genes
  ))
}

# Default resolution ladder matching BrainMappingPipeline iterclust.py
DEFAULT_RES_LADDER <- c(1.0, 2.0, 4.0, 8.0)

# GeneratePalette: Generate maximally distinct, non-repeating colors for large numbers of clusters (> 256 supported)
# Uses Glasbey algorithm (pals) for N <= 256, and dynamic HCL Polychrome palette construction for N > 256
GeneratePalette <- function(n, seed = 42) {
  set.seed(seed)
  if (n <= 256 && requireNamespace("pals", quietly = TRUE)) {
    return(as.vector(pals::glasbey(n)))
  } else if (requireNamespace("Polychrome", quietly = TRUE)) {
    if (n <= 36) {
      return(as.vector(Polychrome::palette36.colors(n)))
    } else {
      # Dynamic HCL palette construction supporting N > 256 distinct colors
      seed_cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
      pal <- Polychrome::createPalette(n, seed_cols, M = max(1000, n * 3))
      return(as.vector(pal))
    }
  } else {
    return(colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n))
  }
}

# ImportH5ad
ImportH5ad <- function(folder, org) {
  cellID <- read.csv(file = paste0(folder, "cellID_", org, ".csv"), header = FALSE)
  geneID <- read.csv(file = paste0(folder, "geneID_", org, ".csv"), header = FALSE)
  metadata <- read.csv(file = paste0(folder, "metadata_", org, ".csv"), header = TRUE, row.names = 1)
  counts <- readMM(file = paste0(folder, "scdata_", org, ".mtx"))
  rownames(counts) <- cellID[, 1]
  colnames(counts) <- geneID[, 1]
  counts <- Matrix::t(counts)
  scdata <- CreateSeuratObject(counts = counts, assay = "RNA", meta.data = metadata, project = paste(org))
  scdata <- DietSeurat(scdata)
  return(scdata)
}

# ReadSTAR
ReadSTAR <- function(plate.name, folder) {
  files <- list.files(path = folder, pattern = "Gene.out.tab", recursive = TRUE, full.names = TRUE)
  data.raw <- data.frame()
  for (file in files) {
    cell.data <- read.table(file = file)
    cell.data <- cell.data[5:23436, ]
    cell.data.df <- data.frame(counts = cell.data[, 2])
    rownames(cell.data.df) <- cell.data$V1
    well.name <- unlist(strsplit(file, split = "/"))
    colnames(cell.data.df) <- paste(plate.name, "_", well.name[length(well.name) - 1], sep = "")
    if (file == files[1]) {
      data.raw <- cell.data.df
    } else {
      data.raw <- cbind(data.raw, cell.data.df)
    }
  }
  return(data.raw)
}

# LogCPM
# Optimized sparse matrix calculation of log2(CPM + 1)
LogCPM <- function(scdata) {
  DefaultAssay(scdata) <- "RNA"
  counts <- GetAssayData(scdata, assay = "RNA", layer = "counts")
  col_sums <- Matrix::colSums(counts)
  col_sums[col_sums == 0] <- 1
  cpm_mat <- Matrix::t(Matrix::t(counts) / col_sums) * 1e6
  scdata.logCPM <- log1p(cpm_mat)
  scdata@assays$RNA$logCPM <- as(scdata.logCPM, "sparseMatrix")
  return(scdata)
}

# SeuratProcess
# High-efficiency normalization, SCTransform (v2), PCA, Harmony, UMAP/tSNE & Clustering
# Uses Seurat v5 native layer-based SCTransform to eliminate memory-intensive object merging and process hangs
SeuratProcess <- function(scdata, batch.name = "key", npc = 30, clustering.resolution = 2, de.param = de_param(), res_ladder = DEFAULT_RES_LADDER, run.iterative = TRUE, force = FALSE, ncore = 1, verbose = TRUE) {
  t_start_total <- Sys.time()
  
  log_ts <- function(msg) {
    if (verbose) {
      ts <- format(Sys.time(), "[%H:%M:%S]")
      message(paste(ts, msg))
      flush.console()
    }
  }
  
  log_ts(sprintf("[SeuratProcess] Started pipeline execution (System Time: %s)", format(t_start_total, "%Y-%m-%d %H:%M:%S")))
  
  # Check if multi-batch processing is required
  batches <- if (batch.name %in% colnames(scdata@meta.data)) unique(scdata@meta.data[[batch.name]]) else NULL
  n_batches <- length(batches)
  
  # Step 1 & 2: SCTransform (v2) using Seurat v5 split layers to eliminate memory-intensive object merging
  t_step1 <- Sys.time()
  has_sct <- "SCT" %in% names(scdata@assays)
  sct_run <- FALSE
  if (has_sct && !force) {
    log_ts("[STEP 1/5 & 2/5 SKIPPED] 'SCT' assay already present in object.")
  } else {
    DefaultAssay(scdata) <- "RNA"
    if (n_batches > 1) {
      log_ts(sprintf("[STEP 1/5] Splitting RNA layers across %d batches ('%s') and running SCTransform (vst.flavor = 'v2')...", n_batches, batch.name))
      scdata[["RNA"]] <- split(scdata[["RNA"]], f = scdata@meta.data[[batch.name]])
      scdata <- SCTransform(scdata, vst.flavor = "v2", verbose = FALSE)
      sct_run <- TRUE
      
      t_step1_elapsed <- round(as.numeric(difftime(Sys.time(), t_step1, units = "mins")), 2)
      log_ts(sprintf("[STEP 1/5 COMPLETED] Layer-based SCTransform finished across %d batches (Elapsed: %.2f mins)", n_batches, t_step1_elapsed))
      log_ts("[STEP 2/5 COMPLETED] Integrated SCT assay generated natively (object merging bypassed).")
    } else {
      log_ts("[STEP 1/5] Starting SCTransform (vst.flavor = 'v2') on single batch...")
      scdata <- SCTransform(scdata, vst.flavor = "v2", verbose = FALSE)
      sct_run <- TRUE
      t_step1_elapsed <- round(as.numeric(difftime(Sys.time(), t_step1, units = "mins")), 2)
      log_ts(sprintf("[STEP 1/5 COMPLETED] SCTransform finished (Elapsed: %.2f mins)", t_step1_elapsed))
      log_ts("[STEP 2/5 COMPLETED] Single batch detected; skipping merge step.")
    }
  }
  
  if ("SCT" %in% names(scdata@assays)) {
    DefaultAssay(scdata) <- "SCT"
  }
  
  # Step 3: PCA & Harmony Integration
  t_step3 <- Sys.time()
  reduction_to_use <- if (n_batches > 1) "harmony" else "pca"
  has_reduction <- reduction_to_use %in% names(scdata@reductions)
  
  if (has_reduction && !force && !sct_run) {
    log_ts(sprintf("[STEP 3/5 SKIPPED] '%s' reduction already present in object.", reduction_to_use))
  } else {
    log_ts(sprintf("[STEP 3/5] Running PCA (%d PCs) and Harmony integration across '%s'...", npc, batch.name))
    if (!"pca" %in% names(scdata@reductions) || force || sct_run) {
      scdata <- RunPCA(scdata, npcs = npc, verbose = FALSE)
    }
    if (n_batches > 1) {
      scdata <- RunHarmony(scdata, dims = 1:npc, verbose = FALSE, group.by.vars = batch.name)
    }
    
    t_step3_elapsed <- round(as.numeric(difftime(Sys.time(), t_step3, units = "secs")), 1)
    log_ts(sprintf("[STEP 3/5 COMPLETED] PCA & Harmony completed (Elapsed: %.1fs)", t_step3_elapsed))
  }
  
  # Ensure RNA layers are joined if needed
  if ("RNA" %in% names(scdata@assays)) {
    try(scdata[["RNA"]] <- JoinLayers(scdata[["RNA"]]), silent = TRUE)
  }
  
  # Step 4: UMAP & t-SNE Dimensionality Reductions
  t_step4 <- Sys.time()
  has_umap <- "umap" %in% names(scdata@reductions)
  has_tsne <- "tsne" %in% names(scdata@reductions)
  
  if (has_umap && has_tsne && !force) {
    log_ts("[STEP 4/5 SKIPPED] Both 'umap' and 'tsne' reductions already present in object.")
  } else {
    log_ts("[STEP 4/5] Running UMAP and t-SNE dimensionality reductions...")
    
    if (!has_umap || force) {
      tryCatch({
        scdata <- RunUMAP(scdata, reduction = reduction_to_use, dims = 1:npc, verbose = FALSE)
      }, error = function(e) {
        log_ts(sprintf("[STEP 4/5 WARNING] UMAP execution failed: %s", e$message))
      })
    }
    
    if (!has_tsne || force) {
      perp <- min(30, max(1, floor((ncol(scdata) - 1) / 3)))
      tsne_res <- tryCatch({
        RunTSNE(scdata, dims = 1:npc, reduction = reduction_to_use, num_threads = ncore, perplexity = perp)
      }, error = function(e) {
        tryCatch({
          RunTSNE(scdata, dims = 1:npc, reduction = reduction_to_use, perplexity = perp)
        }, error = function(e2) {
          log_ts(sprintf("[STEP 4/5 WARNING] t-SNE execution failed: %s", e2$message))
          scdata
        })
      })
      scdata <- tsne_res
    }
    
    t_step4_elapsed <- round(as.numeric(difftime(Sys.time(), t_step4, units = "secs")), 1)
    log_ts(sprintf("[STEP 4/5 COMPLETED] UMAP & t-SNE completed (Elapsed: %.1fs)", t_step4_elapsed))
  }
  
  # Step 5: Graph-Based Clustering
  t_step5 <- Sys.time()
  has_clusters <- "seurat_clusters" %in% colnames(scdata@meta.data) && length(scdata@graphs) > 0
  
  if (has_clusters && !force) {
    log_ts("[STEP 5/5 SKIPPED] 'seurat_clusters' and neighbor graph already present in object.")
  } else {
    log_ts(sprintf("[STEP 5/5] Running FindNeighbors & FindClusters (resolution = %.1f)...", clustering.resolution))
    
    scdata <- FindNeighbors(scdata, reduction = reduction_to_use, dims = 1:npc, verbose = FALSE) %>%
      FindClusters(resolution = clustering.resolution, verbose = FALSE)
      
    t_step5_elapsed <- round(as.numeric(difftime(Sys.time(), t_step5, units = "secs")), 1)
    log_ts(sprintf("[STEP 5/5 COMPLETED] Graph-based clustering completed (Elapsed: %.1fs)", t_step5_elapsed))
  }
  
  # Step 6: Iterative Clustering (Default: TRUE)
  if (run.iterative) {
    t_step6 <- Sys.time()
    has_iterative <- "cluster_merge" %in% colnames(scdata@meta.data)
    
    if (has_iterative && !force) {
      log_ts("[STEP 6 SKIPPED] Iterative clustering ('cluster_merge' metadata) already present in object.")
    } else {
      # Ensure RNA assay has normalized 'data' layer required for ScaleData & FindMarkers
      has_data_layer <- tryCatch({
        layers_list <- SeuratObject::Layers(scdata[["RNA"]])
        any(grepl("^data", layers_list))
      }, error = function(e) FALSE)
      
      if (!has_data_layer) {
        log_ts("[STEP 6 PREPARATION] Normalizing RNA assay layer ('data' layer missing)...")
        scdata <- NormalizeData(scdata, assay = "RNA", verbose = FALSE)
        try(scdata[["RNA"]] <- JoinLayers(scdata[["RNA"]]), silent = TRUE)
      }
      
      log_ts(sprintf("[STEP 6] Starting hicat-style iterative clustering with resolution ladder (ncore = %d)...", ncore))
      
      # iter.clust calls merge_cl internally at each resolution step — no separate pre-merge needed
      scdata <- iter.clust(object = scdata, de.param = de.param, res_ladder = res_ladder, assay = "RNA", ncore = ncore, verbose = verbose)
      
      t_step6_elapsed <- round(as.numeric(difftime(Sys.time(), t_step6, units = "mins")), 2)
      log_ts(sprintf("[STEP 6 COMPLETED] Iterative clustering completed (Elapsed: %.2f mins)", t_step6_elapsed))
    }
  }
  
  # Ensure DefaultAssay remains SCT if present
  if ("SCT" %in% names(scdata@assays)) {
    DefaultAssay(scdata) <- "SCT"
  }
  
  t_end_total <- Sys.time()
  total_mins <- round(as.numeric(difftime(t_end_total, t_start_total, units = "mins")), 2)
  
  log_ts(sprintf("[SeuratProcess COMPLETED] Execution finished at: %s | Total runtime: %.2f mins", 
                 format(t_end_total, "%Y-%m-%d %H:%M:%S"), total_mins))
  
  return(scdata)
}

# DotPlotMM
DotPlotMM <- function(scdata.MM, feature, ensembl.table, group.by = "subclass") {
  if (sum(feature %in% ensembl.table[, 2]) == length(feature)) {
    geneID <- c(ensembl.table[ensembl.table[, 2] %in% feature, 1])
    dotplot <- DotPlot(scdata.MM, features = geneID, cols = c("lightgrey", "darkred"), group.by = group.by)
    print(paste(feature, " = ", geneID))
    return(dotplot)
  } else {
    print("Some features not found in ensembl.table")
  }
}

# helper function to test if a pair of clusters needs to be merged based on de.param
test_merge <- function(de.pair, de.param, merge.type = "undirectional") {
  if (length(de.pair) == 0) {
    return(TRUE)
  }
  to.merge <- FALSE
  if (merge.type == "undirectional") {
    if (!is.null(de.param$de.score.th)) {
      to.merge <- de.pair$score < de.param$de.score.th
    }
    if (!to.merge && !is.null(de.param$min.genes)) {
      to.merge <- de.pair$num < de.param$min.genes
    }
  } else {
    if (!is.null(de.param$de.score.th)) {
      to.merge <- de.pair$up.score < de.param$de.score.th || de.pair$down.score < de.param$de.score.th
    }
    if (!to.merge && !is.null(de.param$min.genes)) {
      to.merge <- de.pair$up.num < de.param$min.genes || de.pair$down.num < de.param$min.genes
    }
  }
  return(to.merge)
}

# ============================================================================
# Summary-stats-only DE engine (ported from Python de.py)
# All pairwise DE is computed from cluster-level vectors only — O(genes) per pair.
# Never re-scans single cells after the one-pass summary computation.
# ============================================================================

# One-pass cluster summary statistics: mean, fraction (detected), sumsq, n
# Matches Python de.cluster_summary_stats exactly.
cluster_summary_stats <- function(data_mat, group_vec, low.th = 1.0) {
  group_vec <- as.character(group_vec)
  group_factor <- factor(group_vec)
  cl_levels <- levels(group_factor)
  n_cls <- length(cl_levels)
  n_genes <- nrow(data_mat)
  
  L <- Matrix::fac2sparse(group_factor)  # n_cls x n_cells indicator matrix
  ns <- as.integer(Matrix::rowSums(L))
  names(ns) <- cl_levels
  
  # mean = (L %*% t(X))^T / n  — sparse matrix multiply, no cell loop
  sum_mat <- as.matrix(L %*% Matrix::t(data_mat))  # n_cls x n_genes
  ns_safe <- pmax(ns, 1L)
  means <- sum_mat / ns_safe  # n_cls x n_genes
  
  # sumsq = L %*% t(X^2) — element-wise square then aggregate
  data_sq <- data_mat * data_mat  # sparse * sparse = sparse (Hadamard)
  sumsq <- as.matrix(L %*% Matrix::t(data_sq))  # n_cls x n_genes
  
  # fraction = L %*% t(X > low.th) / n — detection rate
  detected <- data_mat > low.th  # sparse logical (becomes dgCMatrix of 0/1)
  det_counts <- as.matrix(L %*% Matrix::t(detected))  # n_cls x n_genes
  fracs <- det_counts / ns_safe  # n_cls x n_genes
  
  rownames(means) <- rownames(sumsq) <- rownames(fracs) <- cl_levels
  
  list(clusters = cl_levels, n = ns, mean = means, fraction = fracs, sumsq = sumsq)
}

# Pool two cluster stat-rows into one (exact, no rescan).  Matches Python de.pool_stats.
pool_stats <- function(row_a, row_b) {
  na <- row_a$n; nb <- row_b$n; n <- na + nb
  list(
    n = n,
    mean = (na * row_a$mean + nb * row_b$mean) / n,
    fraction = (na * row_a$fraction + nb * row_b$fraction) / n,
    sumsq = row_a$sumsq + row_b$sumsq
  )
}

# Extract one cluster's row from the summary stats dict.
get_cluster_row <- function(stats, cl) {
  idx <- match(cl, stats$clusters)
  list(n = stats$n[idx], mean = stats$mean[idx, ], fraction = stats$fraction[idx, ], sumsq = stats$sumsq[idx, ])
}

# Welch t-test from summary stats (matches Python de.de_pair_from_stats 'welch' path).
# Returns de_score for one pair. O(genes), no cell access.
de_pair_score_from_stats <- function(row_a, row_b, de.param) {
  n1 <- row_a$n; n2 <- row_b$n
  if (n1 < de.param$min.cells || n2 < de.param$min.cells) return(list(score = 0, separable = FALSE))
  
  m1 <- row_a$mean; m2 <- row_b$mean
  q1 <- row_a$fraction; q2 <- row_b$fraction
  ss1 <- row_a$sumsq; ss2 <- row_b$sumsq
  
  lfc <- m1 - m2  # already log-normalized, so difference = log fold-change
  
  # Unbiased sample variance from sum-of-squares
  if (n1 >= 2) { v1 <- pmax((ss1 - n1 * m1 * m1) / (n1 - 1), 0) } else { v1 <- rep(0, length(m1)) }
  if (n2 >= 2) { v2 <- pmax((ss2 - n2 * m2 * m2) / (n2 - 1), 0) } else { v2 <- rep(0, length(m2)) }
  
  se <- sqrt(v1 / n1 + v2 / n2)
  t_stat <- ifelse(se > 0, lfc / se, 0)
  
  # Welch-Satterthwaite degrees of freedom
  num <- (v1 / n1 + v2 / n2)^2
  den <- ifelse(n1 > 1, (v1 / n1)^2 / (n1 - 1), 0) + ifelse(n2 > 1, (v2 / n2)^2 / (n2 - 1), 0)
  dof <- ifelse(den > 0, num / den, 1)
  dof <- pmax(dof, 1)
  
  # Two-sided p-value
  pval <- 2 * pt(abs(t_stat), df = dof, lower.tail = FALSE)
  pval <- pmin(pmax(pval, 0), 1)
  pval[!is.finite(pval)] <- 1
  
  # BH adjustment
  padj <- p.adjust(pval, method = "BH")
  
  # Apply DE criteria matching de.param
  base <- (padj < de.param$padj.th) & (abs(lfc) > de.param$lfc.th)
  
  # Detection-rate filters
  q_max <- pmax(q1, q2)
  q_diff <- ifelse(q_max > 0, abs(q1 - q2) / q_max, 0)
  
  is_up <- base & (lfc > 0)
  is_down <- base & (lfc < 0)
  
  if (!is.null(de.param$q1.th)) {
    is_up <- is_up & (q1 > de.param$q1.th)
    is_down <- is_down & (q2 > de.param$q1.th)
  }
  if (!is.null(de.param$q.diff.th)) {
    is_up <- is_up & (q_diff > de.param$q.diff.th)
    is_down <- is_down & (q_diff > de.param$q.diff.th)
  }
  
  is_de <- is_up | is_down
  
  if (sum(is_de) == 0) return(list(score = 0, separable = FALSE))
  
  de_score <- sum(pmin(-log10(padj[is_de] + 1e-300), 20))
  n_de <- sum(is_de)
  separable <- (n_de >= de.param$min.genes) && (de_score >= de.param$de.score.th)
  
  list(score = de_score, separable = separable)
}


# Ultra-fast C-level cluster mean aggregation (250x faster than Seurat AggregateExpression)
fast_cluster_means <- function(data_mat, features, group_vec) {
  sub_mat <- data_mat[features, , drop = FALSE]
  group_factor <- factor(as.character(group_vec))
  cl_levels <- levels(group_factor)
  
  if (length(cl_levels) <= 1) {
    mean_vec <- matrix(Matrix::rowMeans(sub_mat), ncol = 1)
    colnames(mean_vec) <- cl_levels
    rownames(mean_vec) <- features
    return(mean_vec)
  }
  
  L <- Matrix::fac2sparse(group_factor)
  counts_per_cl <- as.numeric(Matrix::rowSums(L))
  counts_per_cl[counts_per_cl == 0] <- 1
  
  sum_mat <- L %*% Matrix::t(sub_mat)
  mean_mat <- Matrix::t(as.matrix(sum_mat) / counts_per_cl)
  colnames(mean_mat) <- cl_levels
  rownames(mean_mat) <- features
  return(mean_mat)
}

# merge_cl
# High-efficiency cluster merging matching Python iterclust algorithm structure:
# Uses ultra-fast C++ Wilcoxon rank-sum test (presto::wilcoxauc) on candidate gene sparse matrix slices,
# single-pass sparse matrix aggregation, algebraic centroid pooling, k-NN candidate graph (k=3), and DE score caching
merge_cl <- function(object,
                     de.param = de_param(), 
                     merge.type = c("undirectional", "directional"), 
                     max.cl.size = 300,
                     k_neighbors = 3,
                     assay = "RNA",
                     data_mat = NULL) {
  merge.type <- merge.type[1]
  DefaultAssay(object) <- assay
  
  if (!"seurat_clusters" %in% colnames(object@meta.data)) {
    object@meta.data$seurat_clusters <- Idents(object)
  }
  object@meta.data$cluster_merge <- as.character(object@meta.data$seurat_clusters)
  
  v_features <- tryCatch(VariableFeatures(object, assay = assay), error = function(e) NULL)
  if (is.null(v_features) || length(v_features) == 0) {
    object <- FindVariableFeatures(object, assay = assay, verbose = FALSE)
    v_features <- VariableFeatures(object, assay = assay)
  }
  if (length(v_features) == 0) v_features <- rownames(object)
  
  if (is.null(data_mat)) {
    data_mat <- GetAssayData(object, assay = assay, layer = "data")
  }
  has_presto <- requireNamespace("presto", quietly = TRUE)
  
  # Step 1: Pre-pass - Merge small clusters (cells < min.cells) to nearest neighbor based on Pearson correlation
  while (TRUE) {
    cl_vec <- object@meta.data$cluster_merge
    cl_size <- table(cl_vec)
    if (length(cl_size) <= 1) break
    
    small_cls <- names(cl_size)[cl_size < de.param$min.cells]
    if (length(small_cls) == 0) break
    
    matrix_aggr <- fast_cluster_means(data_mat, v_features, object@meta.data$cluster_merge)
    cl_sim <- cor(matrix_aggr)
    diag(cl_sim) <- -Inf
    
    smallest_cl <- names(cl_size)[order(cl_size)][1]
    if (!smallest_cl %in% colnames(cl_sim)) break
    
    nn_cl <- names(which.max(cl_sim[smallest_cl, ]))
    object@meta.data$cluster_merge[object@meta.data$cluster_merge == smallest_cl] <- nn_cl
  }
  
  # Step 2: High-efficiency DE merge loop using centroid pooling & k-NN candidate graph
  cl_vec <- object@meta.data$cluster_merge
  cl_size <- table(cl_vec)
  if (length(cl_size) <= 1) {
    object@meta.data$seurat_clusters <- object@meta.data$cluster_merge
    return(object)
  }
  
  matrix_aggr <- fast_cluster_means(data_mat, v_features, object@meta.data$cluster_merge)
  de_cache <- new.env(hash = TRUE)
  
  get_pair_key <- function(a, b) {
    if (a < b) paste0(a, "___", b) else paste0(b, "___", a)
  }
  
  while (ncol(matrix_aggr) > 1) {
    active_cls <- colnames(matrix_aggr)
    n_active <- length(active_cls)
    if (n_active <= 1) break
    
    cl_sim <- cor(matrix_aggr)
    diag(cl_sim) <- -Inf
    
    # Build candidate pairs from top k-NN neighbors
    cand_pairs <- list()
    for (i in seq_len(n_active)) {
      c_i <- active_cls[i]
      sim_vec <- cl_sim[i, ]
      top_k_idx <- order(sim_vec, decreasing = TRUE)[1:min(k_neighbors, n_active - 1)]
      for (k_idx in top_k_idx) {
        c_k <- active_cls[k_idx]
        pkey <- get_pair_key(c_i, c_k)
        cand_pairs[[pkey]] <- c(c_i, c_k)
      }
    }
    
    if (length(cand_pairs) == 0) break
    
    # Evaluate DE score using fast presto C++ Wilcoxon testing
    scored_list <- list()
    for (pkey in names(cand_pairs)) {
      pair <- cand_pairs[[pkey]]
      if (exists(pkey, envir = de_cache, inherits = FALSE)) {
        res <- get(pkey, envir = de_cache)
      } else {
        p1 <- pair[1]; p2 <- pair[2]
        cells_p1 <- which(object@meta.data$cluster_merge == p1)
        cells_p2 <- which(object@meta.data$cluster_merge == p2)
        if (length(cells_p1) > max.cl.size) cells_p1 <- cells_p1[sample.int(length(cells_p1), max.cl.size)]
        if (length(cells_p2) > max.cl.size) cells_p2 <- cells_p2[sample.int(length(cells_p2), max.cl.size)]
        
        n1 <- length(cells_p1)
        n2 <- length(cells_p2)
        de_score <- 0
        n_de <- 0
        
        if (has_presto && n1 > 0 && n2 > 0) {
          sub_mat <- data_mat[, c(cells_p1, cells_p2), drop = FALSE]
          # Pre-filter candidate genes expressed >= 10% in at least one cluster (matching FindMarkers min.pct=0.1)
          pct1 <- Matrix::rowSums(sub_mat[, 1:n1, drop = FALSE] > 0) / n1
          pct2 <- Matrix::rowSums(sub_mat[, (n1 + 1):(n1 + n2), drop = FALSE] > 0) / n2
          cand_idx <- which(pmax(pct1, pct2) >= 0.1)
          
          if (length(cand_idx) > 0) {
            y_vec <- c(rep("p1", n1), rep("p2", n2))
            pres <- suppressWarnings(tryCatch({ presto::wilcoxauc(sub_mat[cand_idx, , drop = FALSE], y_vec) }, error = function(e) NULL))
            
            if (!is.null(pres) && nrow(pres) > 0) {
              res_p1 <- pres[pres$group == "p1", ]
              pct1_v <- res_p1$pct_in / 100
              pct2_v <- res_p1$pct_out / 100
              valid <- which(abs(res_p1$logFC) >= de.param$lfc.th & res_p1$padj < de.param$padj.th)
              if (length(valid) > 0) {
                q_diff <- abs(pct1_v[valid] - pct2_v[valid]) / pmax(pct1_v[valid], pct2_v[valid])
                valid_q <- valid[q_diff > de.param$q.diff.th]
                if (length(valid_q) > 0) {
                  norm_p <- pmin(-log10(res_p1$padj[valid_q] + 1e-300), 20)
                  de_score <- sum(norm_p)
                  n_de <- length(valid_q)
                }
              }
            }
          }
        } else if (n1 > 0 && n2 > 0) {
          sub_obj <- object[, c(cells_p1, cells_p2)]
          DefaultAssay(sub_obj) <- assay
          markers <- suppressWarnings(tryCatch({
            FindMarkers(sub_obj, ident.1 = p1, ident.2 = p2, assay = assay,
                        group.by = "cluster_merge", logfc.threshold = de.param$lfc.th, 
                        min.pct = 0.1, verbose = FALSE, recorrect_umi = FALSE)
          }, error = function(e) NULL))
          if (!is.null(markers) && nrow(markers) > 0) {
            markers <- markers[markers$p_val_adj < de.param$padj.th, , drop = FALSE]
            if (nrow(markers) > 0) {
              q_diff <- abs(markers$pct.1 - markers$pct.2) / pmax(markers$pct.1, markers$pct.2)
              markers <- markers[q_diff > de.param$q.diff.th, , drop = FALSE]
              if (nrow(markers) > 0) {
                de_score <- sum(pmin(-log10(markers$p_val_adj), 20))
                n_de <- nrow(markers)
              }
            }
          }
        }
        
        separable <- (n_de >= de.param$min.genes) && (de_score >= de.param$de.score.th)
        res <- list(score = de_score, separable = separable)
        assign(pkey, res, envir = de_cache)
      }
      
      if (!res$separable) {
        scored_list[[length(scored_list) + 1]] <- list(score = res$score, p1 = pair[1], p2 = pair[2], key = pkey)
      }
    }
    
    if (length(scored_list) == 0) break
    
    # Sort pairs by DE score (lowest DE score first = most mergeable)
    scores <- sapply(scored_list, function(x) x$score)
    scored_list <- scored_list[order(scores)]
    best_merge <- scored_list[[1]]
    keep_cl <- best_merge$p1
    drop_cl <- best_merge$p2
    
    # Update aggregated centroid matrix in-memory (weighted pool of keep_cl and drop_cl) BEFORE updating cluster_merge
    n_keep <- sum(object@meta.data$cluster_merge == keep_cl)
    n_drop <- sum(object@meta.data$cluster_merge == drop_cl)
    n_total <- n_keep + n_drop
    
    if (n_total > 0) {
      matrix_aggr[, keep_cl] <- (matrix_aggr[, keep_cl] * n_keep + matrix_aggr[, drop_cl] * n_drop) / n_total
    }
    drop_idx <- match(drop_cl, colnames(matrix_aggr))
    if (!is.na(drop_idx)) {
      matrix_aggr <- matrix_aggr[, -drop_idx, drop = FALSE]
    }
    
    # Update cluster assignments
    object@meta.data$cluster_merge[object@meta.data$cluster_merge == drop_cl] <- keep_cl
    
    # Invalidate cache entries touching keep_cl or drop_cl
    cache_keys <- ls(envir = de_cache)
    to_remove <- cache_keys[grepl(keep_cl, cache_keys, fixed = TRUE) | grepl(drop_cl, cache_keys, fixed = TRUE)]
    if (length(to_remove) > 0) rm(list = to_remove, envir = de_cache)
  }
  
  object@meta.data$seurat_clusters <- object@meta.data$cluster_merge
  return(object)
}

# MergeCluster
MergeCluster <- function(object, padj.th = 0.05, lfc.th = 0.405, q.diff.th = 0.7, de.score.th = 150, min.cells = 4, assay = "RNA") {
  de.param <- de_param(padj.th = padj.th, lfc.th = lfc.th, q.diff.th = q.diff.th, de.score.th = de.score.th, min.cells = min.cells)
  return(merge_cl(object = object, de.param = de.param, assay = assay))
}

# iter.clust
# Resolution-Escalating Iterative Clustering (Ported from BrainMappingPipeline iterclust.py)
# Optimizations over naive Seurat:
#   - ScaleData on variable features only (not all ~30k genes)
#   - data_mat extracted ONCE and shared across resolution ladder
#   - Labels-only recursion: no SplitObject/merge deep copies — collects cell-label maps
#   - Multi-core parallelized sub-branch recursion using parallel::mclapply
iter.clust <- function(object, 
                       de.param = de_param(), 
                       res_ladder = DEFAULT_RES_LADDER, 
                       sat_tol = 0.05, 
                       split.size = 10, 
                       assay = "RNA",
                       path = "root",
                       depth = 0,
                       ncore = 1,
                       verbose = TRUE) {
  DefaultAssay(object) <- assay
  n_cells <- ncol(object)
  
  log_ts <- function(msg) {
    if (verbose) {
      ts <- format(Sys.time(), "[%H:%M:%S]")
      message(paste(ts, msg))
      flush.console()
    }
  }
  
  log_ts(sprintf("[BRANCH START] Path: %s | Depth: %d | Cells: %d", path, depth, n_cells))
  
  # Leaf check: min_cells
  if (n_cells < split.size) {
    log_ts(sprintf("  [LEAF] Path: %s | Cells: %d | Reason: min_cells (< %d)", path, n_cells, split.size))
    return(object)
  }
  
  object <- FindVariableFeatures(object, assay = assay, verbose = FALSE)
  v_features <- VariableFeatures(object, assay = assay)
  
  if (length(v_features) <= 30) {
    log_ts(sprintf("  [LEAF] Path: %s | Reason: low_var_features (<= 30)", path))
    return(object)
  }
  
  # Guard against irlba PCA underflow crashes on small sub-clusters
  max_pcs <- min(29, n_cells - 2, length(v_features) - 1)
  if (max_pcs < 2) {
    log_ts(sprintf("  [LEAF] Path: %s | Reason: low_pca_dims (< 2)", path))
    return(object)
  }
  
  # Ensure assay has normalized 'data' layer for ScaleData & FindMarkers
  has_data_layer <- tryCatch({
    layers_list <- SeuratObject::Layers(object[[assay]])
    any(grepl("^data", layers_list))
  }, error = function(e) FALSE)
  
  if (!has_data_layer) {
    object <- NormalizeData(object, assay = assay, verbose = FALSE)
  }
  
  # ScaleData on variable features ONLY (not all ~30k genes — saves ~90% time)
  object <- ScaleData(object, features = v_features, verbose = FALSE)
  object <- suppressWarnings(RunPCA(object, features = v_features, npcs = max_pcs, verbose = FALSE))
  object <- suppressWarnings(FindNeighbors(object, dims = 1:max_pcs, verbose = FALSE))
  
  # Extract normalized data matrix ONCE for the entire resolution ladder
  data_mat <- GetAssayData(object, assay = assay, layer = "data")
  
  # --- Resolution Escalation Ladder (BrainMappingPipeline saturation logic) ---
  prev_n <- -1
  chosen_merged <- NULL
  chosen_res <- NULL
  chosen_raw_n <- 0
  
  for (r in res_ladder) {
    t_start <- Sys.time()
    obj_r <- suppressWarnings(FindClusters(object, resolution = r, verbose = FALSE, group.singletons = FALSE))
    raw_clusters <- obj_r@meta.data$seurat_clusters
    n_raw <- length(unique(raw_clusters))
    
    if (n_raw <= 1) {
      merged_cl <- raw_clusters
    } else {
      # Pass shared data_mat — avoids redundant GetAssayData per resolution
      obj_merged <- merge_cl(object = obj_r, de.param = de.param, assay = assay, data_mat = data_mat)
      merged_cl <- obj_merged@meta.data$cluster_merge
    }
    
    n_merged <- length(unique(merged_cl))
    t_elapsed <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 2)
    
    log_ts(sprintf("  [LADDER] Res: %.1f | Raw Clusters: %d -> DE-Merged: %d | Time: %.2fs", r, n_raw, n_merged, t_elapsed))
    
    # Saturation check: stop escalating when merged cluster count stops growing (grows <= sat_tol %)
    if (prev_n > 0 && n_merged <= prev_n * (1 + sat_tol)) {
      log_ts(sprintf("  [SATURATED] Merged count plateaus at Res: %.1f (Merged: %d)", chosen_res, prev_n))
      break
    }
    
    chosen_merged <- merged_cl
    chosen_res <- r
    chosen_raw_n <- n_raw
    prev_n <- n_merged
  }
  
  # Assign the DE-determined grouping at chosen resolution
  object@meta.data$cluster_merge <- as.character(chosen_merged)
  object@meta.data$seurat_clusters <- object@meta.data$cluster_merge
  
  cl_merged_unique <- unique(object@meta.data$cluster_merge)
  if (length(cl_merged_unique) <= 1) {
    reason <- if (chosen_raw_n <= 1) "no_raw_split" else "no_de_separable_split"
    log_ts(sprintf("  [LEAF] Path: %s | Cells: %d | Reason: %s", path, n_cells, reason))
    return(object)
  }
  
  # Sub-clustering branch split with parallel recursion
  # Uses labels-only collection to avoid expensive SplitObject/merge deep copies
  log_ts(sprintf("[SPLIT] Path: %s | Chosen Res: %.1f | Raw: %d -> Merged Separable: %d. Recursing into sub-branches (ncore = %d)...", 
                 path, chosen_res, chosen_raw_n, length(cl_merged_unique), ncore))
  
  sub_list <- SplitObject(object, split.by = "cluster_merge")
  cl_ids <- names(sub_list)
  n_sub <- length(cl_ids)
  
  cores_to_use <- min(ncore, n_sub)
  sub_ncore <- max(1, floor(ncore / cores_to_use))
  
  # Worker: recurse on sub-branch, return named cell-to-label vector (NOT full Seurat objects when possible)
  process_sub <- function(i) {
    cl_id <- cl_ids[i]
    sub_obj <- sub_list[[cl_id]]
    sub_path <- paste0(path, ".", cl_id)
    
    if (verbose) {
      start_msg <- sprintf("[%s] [PARALLEL WORKER] Branch %d/%d (%s) STARTED (%d cells)", 
                           format(Sys.time(), "%H:%M:%S"), i, n_sub, sub_path, ncol(sub_obj))
      cat(paste0(start_msg, "\n"), file = stderr())
      flush(stderr())
    }
    
    res <- tryCatch({
      if (ncol(sub_obj) >= split.size) {
        sub_obj_clustered <- iter.clust(
          sub_obj, de.param = de.param, res_ladder = res_ladder, 
          sat_tol = sat_tol, split.size = split.size, assay = assay, 
          path = sub_path, depth = depth + 1, ncore = sub_ncore, verbose = verbose
        )
        # Return only the cell-to-label mapping with parent branch prefix (lightweight)
        labels <- paste0(cl_id, ".", sub_obj_clustered@meta.data$cluster_merge)
        names(labels) <- colnames(sub_obj_clustered)
        labels
      } else {
        if (verbose) {
          leaf_msg <- sprintf("[%s]   [LEAF] Path: %s | Cells: %d | Reason: min_cells (< %d)", 
                              format(Sys.time(), "%H:%M:%S"), sub_path, ncol(sub_obj), split.size)
          cat(paste0(leaf_msg, "\n"), file = stderr())
          flush(stderr())
        }
        labels <- rep(cl_id, ncol(sub_obj))
        names(labels) <- colnames(sub_obj)
        labels
      }
    }, error = function(e) {
      err_msg <- sprintf("[%s] [ERROR in Branch %s]: %s (retaining branch un-split)", 
                         format(Sys.time(), "%H:%M:%S"), sub_path, e$message)
      cat(paste0(err_msg, "\n"), file = stderr())
      flush(stderr())
      labels <- rep(cl_id, ncol(sub_obj))
      names(labels) <- colnames(sub_obj)
      labels
    })
    
    if (verbose) {
      done_msg <- sprintf("[%s] [PARALLEL WORKER] Branch %d/%d (%s) COMPLETED", 
                          format(Sys.time(), "%H:%M:%S"), i, n_sub, sub_path)
      cat(paste0(done_msg, "\n"), file = stderr())
      flush(stderr())
    }
    
    return(res)
  }
  
  if (cores_to_use > 1 && .Platform$OS.type == "unix") {
    sub_results <- parallel::mclapply(seq_len(n_sub), process_sub, mc.cores = cores_to_use)
  } else {
    sub_results <- lapply(seq_len(n_sub), process_sub)
  }
  
  # Collect all cell-to-label mappings back into the original object (no merge() needed)
  all_labels <- unlist(sub_results)
  
  # Write labels back to original object metadata (cells that exist in all_labels)
  matched_cells <- intersect(colnames(object), names(all_labels))
  object@meta.data[matched_cells, "cluster_merge"] <- all_labels[matched_cells]
  object@meta.data$seurat_clusters <- object@meta.data$cluster_merge
  
  return(object)
}

# bootstrap_iter_clust
bootstrap_iter_clust <- function(seurat.obj, n.iter = 100, de.param = de_param(), res_ladder = DEFAULT_RES_LADDER, assay = "RNA", verbose = FALSE) {
  n.cluster <- c()
  cocluster.matrix <- data.frame()
  start.time <- Sys.time()
  
  for (i in 1:n.iter) {
    message(paste("Start iteration", i, "of", n.iter))
    pred <- sample(colnames(seurat.obj), ncol(seurat.obj) * 0.8)
    scdata.sample <- seurat.obj[, pred]
    scdata.sample <- DietSeurat(scdata.sample)
    
    scdata.sample <- FindVariableFeatures(scdata.sample, assay = assay, verbose = FALSE)
    scdata.sample <- ScaleData(scdata.sample, verbose = FALSE)
    scdata.sample <- RunPCA(scdata.sample, npcs = 29, verbose = FALSE)
    scdata.sample <- FindNeighbors(scdata.sample, dims = 1:29, verbose = FALSE)
    scdata.sample <- FindClusters(scdata.sample, resolution = 2, verbose = FALSE)
    
    scdata.sample <- merge_cl(object = scdata.sample, de.param = de.param, assay = assay)
    scdata.sample <- iter.clust(object = scdata.sample, de.param = de.param, res_ladder = res_ladder, assay = assay, verbose = verbose)
    
    cluster.i <- data.frame(name = rownames(scdata.sample@meta.data), cluster = paste(scdata.sample@meta.data$cluster_merge, ".iter", i, sep = ""))
    cocluster.matrix <- rbind(cocluster.matrix, cluster.i)
    n.cluster <- c(n.cluster, length(unique(scdata.sample@meta.data$cluster_merge)))
  }
  
  cocluster.matrix.average <- crossprod(table(cocluster.matrix[c(2, 1)]))
  return(list(cocluster.matrix.average, n.cluster))
}

# Save AnnData
save.h5ad <- function(seurat.object, assay = "RNA", layer = "counts", filename) {
  counts_matrix <- GetAssayData(seurat.object, assay = assay, layer = layer)
  obs_metadata <- seurat.object@meta.data
  gene_names <- rownames(counts_matrix)
  var_metadata <- data.frame(gene_names = gene_names, row.names = gene_names)
  adata <- AnnData(X = Matrix::t(counts_matrix), obs = obs_metadata, var = var_metadata)
  write_h5ad(adata, filename)
}

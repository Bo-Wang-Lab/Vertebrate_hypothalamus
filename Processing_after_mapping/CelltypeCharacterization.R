# NT.marker.threshold
#
# For a single gene (NT.marker), computes per-cluster expression metrics across
# a Seurat object. To keep computation tractable and balance cluster representation,
# each cluster is downsampled to at most `downsample` cells before any calculation.
#
# Args:
#   NT.marker     - Character. Gene name to evaluate (must exist in the RNA counts layer).
#   scdata        - Seurat object. Must have an "RNA" assay with a "counts" layer and
#                   nCount_RNA in metadata.
#   cluster.label - Character. Name of the metadata column that holds cluster/cell-type
#                   identity (e.g. "seurat_clusters", "celltype").
#   downsample    - Integer (default 100). Maximum number of cells to retain per cluster.
#                   Clusters with fewer cells are kept in full; larger clusters are randomly
#                   sampled without replacement down to this size.
#
# Returns:
#   A data.frame with one row per cluster and three columns:
#     celltype   - Cluster label (character).
#     mean_lnCPM - Mean ln(CPM + 1) across all cells in the (downsampled) cluster,
#                  where CPM = (raw gene counts / total UMIs) * 1e6.
#     pct.expr   - Percentage of cells in the cluster with at least one raw count > 0.

NT.marker.threshold <- function(NT.marker, scdata, cluster.label, downsample = 100) {
  
  # --- 1. Downsample each cluster ---
  # Extracts cluster labels and cell barcodes from metadata, then for each cluster
  # keeps all cells if n <= downsample, or draws a random sample of size `downsample`
  # without replacement if n > downsample. Ensures no cluster dominates the calculation.
  meta_labels <- scdata[[cluster.label]][, 1]
  cell_names <- rownames(scdata@meta.data)
  
  cells_to_keep <- unlist(lapply(unique(meta_labels), function(cluster) {
    cluster_cells <- cell_names[meta_labels == cluster]
    
    if (length(cluster_cells) <= downsample) {
      return(cluster_cells)                              # keep all cells if at or below threshold
    } else {
      return(sample(cluster_cells, downsample, replace = FALSE))  # random draw if above threshold
    }
  }))
  
  scdata.small <- subset(scdata, cells = cells_to_keep)
  
  # --- 2. Extract raw counts (Seurat v5 layer syntax) ---
  # Seurat v5 uses 'layer' instead of the v4 'slot' argument. Raw integer counts
  # are required here — normalised/log-transformed layers would double-transform
  # the data when we compute CPM manually below.
  counts_matrix <- GetAssayData(scdata.small, assay = "RNA", layer = "counts")
  
  if (!NT.marker %in% rownames(counts_matrix)) {
    stop(paste("Error: The gene", NT.marker, "was not found in the RNA assay."))
  }
  
  labels <- as.character(unique(scdata.small[[cluster.label]][, 1]))
  
  # Pre-allocate the output table: one row per cluster
  marker.exp <- data.frame(
    celltype   = labels,
    mean_lnCPM = numeric(length(labels)),
    pct.expr   = numeric(length(labels))
  )
  rownames(marker.exp) <- labels
  
  # --- 3. Per-cluster metric calculation ---
  for (i in seq_along(labels)) {
    current_cluster <- labels[i]
    cell_ids <- rownames(scdata.small@meta.data)[scdata.small[[cluster.label]][, 1] == current_cluster]
    
    gene_counts  <- counts_matrix[NT.marker, cell_ids]   # raw counts for this gene
    total_counts <- scdata.small$nCount_RNA[cell_ids]    # total UMIs per cell (library size)
    
    # CPM normalises each cell by its library size so cells with different sequencing
    # depths are comparable. Adding 1 before log avoids log(0) = -Inf.
    cpm_values    <- (gene_counts / total_counts) * 1e6
    ln_cpm_values <- log(cpm_values + 1)                 # natural log (base e)
    
    marker.exp$mean_lnCPM[i] <- mean(ln_cpm_values)
    marker.exp$pct.expr[i]   <- (sum(gene_counts > 0) / length(gene_counts)) * 100
  }
  
  return(marker.exp)
}




# lookup_gene_markers
#
# Filters a marker data.frame to return all rows matching a single gene.
# Useful for quickly inspecting the cluster-level statistics (avg_log2FC,
# p-value, pct.1, pct.2, etc.) that Seurat's FindAllMarkers() produces for
# one gene of interest.
#
# Args:
#   gene_name  - Character. Gene symbol to look up (must match the 'gene' column exactly).
#   markers_df - Data.frame. Output of FindAllMarkers() or any table with a 'gene' column.
#
# Returns:
#   A subset data.frame containing only rows where gene == gene_name.
#   Returns an empty data.frame (0 rows) if the gene is not present.

lookup_gene_markers <- function(gene_name, markers_df) {
  markers_df[markers_df$gene == gene_name, ]
}

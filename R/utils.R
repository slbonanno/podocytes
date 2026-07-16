# R/utils.R — shared helper functions, sourced by pipeline scripts.

# Canonical human podocyte marker genes, used to identify the podocyte
# cluster after clustering. Standard across the kidney scRNA-seq
# literature (Wilson 2019, Menon 2020, KPMP atlas).
PODOCYTE_MARKERS <- c(
  "NPHS1", "NPHS2", "PODXL", "WT1", "SYNPO",
  "PTPRO", "CLIC5", "THSD7A", "NTNG1"
)

# Standard QC thresholds -- starting point, tune per-dataset after
# looking at real distributions.
qc_defaults <- list(
  min_genes_per_cell = 200,
  max_pct_mito       = 20,
  min_cells_per_gene = 3
)

#' Run standard QC filtering on a Seurat object
qc_filter <- function(obj, qc = qc_defaults) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  subset(
    obj,
    subset = nFeature_RNA > qc$min_genes_per_cell &
      percent.mt < qc$max_pct_mito
  )
}

#' Free memory -- call between processing large datasets
release_memory <- function(...) {
  rm(list = c(...), envir = parent.frame())
  gc(full = TRUE)
}
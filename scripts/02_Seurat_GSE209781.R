library(Seurat)
library(tidyverse)
BiocManager::install("DropletUtils", update = FALSE, ask = FALSE)
library(DropletUtils)

sample_dirs <- c(
  NM01  = "data/raw/GSE209781/NM01",
  NM02  = "data/raw/GSE209781/NM02",
  NM03  = "data/raw/GSE209781/NM03",
  DKD01 = "data/raw/GSE209781/DKD01",
  DKD02 = "data/raw/GSE209781/DKD02",
  DKD03 = "data/raw/GSE209781/DKD03"
)

#############################################################################
# 1. load in raw data and process in Seurat - vet for compatibility issues
#############################################################################

# note that I noticed a 10x file to be 19MB, much bigger than expected.
# this likely means empty droplets are not yet dropped
# pre-filter for this, so Seurat isn't overwhelmed in building the object

library(DropletUtils)
library(Matrix)

process_one_sample <- function(dir, sample_name, condition) {
  counts <- Read10X(data.dir = dir)
  
  # Statistical test against the ambient RNA profile -- proper method,
  # not a hand-picked UMI cutoff. lower= excludes barcodes too sparse to
  # even test reliably (standard default).
  set.seed(100)  # emptyDrops uses Monte Carlo simulation -- reproducible
  ed <- emptyDrops(counts, lower = 100)
  
  # Keep barcodes that are significantly different from ambient RNA
  # (FDR < 0.01 is the field-standard threshold)
  is_cell <- which(ed$FDR < 0.01)
  counts_filtered <- counts[, is_cell]
  
  cat(sprintf("%s: %d raw barcodes -> %d real cells\n",
              sample_name, ncol(counts), length(is_cell)))
  
  obj <- CreateSeuratObject(counts_filtered, project = sample_name, min.cells = 3)
  obj$sample <- sample_name
  obj$condition <- condition
  
  rm(counts, counts_filtered, ed)
  gc(full = TRUE)
  
  obj
}

nm01_obj <- process_one_sample("data/raw/GSE209781/NM01", "NM01", "healthy")
nm01_obj

# this out: NM01: 6794880 raw barcodes -> 26266 real cells
# still much more than the 5k cells the paper reported.
# probably too low threshold on filtering for good cells.
# load into Seurat to look more carefully
counts_nm01 <- Read10X(data.dir = "data/raw/GSE209781/NM01")
ed_nm01 <- emptyDrops(counts_nm01, lower = 100)
is_cell_nm01 <- which(ed_nm01$FDR < 0.01)

subset_mat <- counts_nm01[seq_len(nrow(counts_nm01)), is_cell_nm01]

cell_umi <- colSums(subset_mat)
cell_genes <- colSums(subset_mat > 0)

summary(cell_umi)
summary(cell_genes)
hist(log10(cell_umi + 1), breaks = 50, main = "UMI distribution among emptyDrops-called cells")

# this reveals that there are still a lot of cells retained that aren't to be used in analysis.
# still, Seurat was freezing because we were processing too much data
# the initial drop of empty droplets should be enough data loss to process

# the above was exploratory, to understand this. below is to process all samples

###################################################################################

#############################################################################
# 2. process all the data files, now that we understand pitfalls
# and incorporated dropping empty droplets before loading into Seurat (computational time suck)
#############################################################################

source("R/utils.R")

sample_dirs <- c(
  NM01  = "data/raw/GSE209781/NM01",
  NM02  = "data/raw/GSE209781/NM02",
  NM03  = "data/raw/GSE209781/NM03",
  DKD01 = "data/raw/GSE209781/DKD01",
  DKD02 = "data/raw/GSE209781/DKD02",
  DKD03 = "data/raw/GSE209781/DKD03"
)

load_and_filter_sample <- function(dir, sample_name) {
  cat(sprintf("--- %s ---\n", sample_name))
  
  counts_raw <- Read10X(data.dir = dir)
  
  set.seed(100)
  ed <- emptyDrops(counts_raw, lower = 100)
  is_cell <- which(ed$FDR < 0.01)
  counts_filtered <- counts_raw[seq_len(nrow(counts_raw)), is_cell]
  
  cat(sprintf("  %d raw barcodes -> %d after emptyDrops\n",
              ncol(counts_raw), length(is_cell)))
  
  condition <- ifelse(grepl("^NM", sample_name), "healthy", "DKD")
  
  obj <- CreateSeuratObject(counts_filtered, project = sample_name, min.cells = 3)
  obj$sample <- sample_name
  obj$condition <- condition
  obj <- qc_filter(obj)
  
  cat(sprintf("  %d cells after qc_filter()\n\n", ncol(obj)))
  
  rm(counts_raw, counts_filtered, ed)
  gc(full = TRUE)
  
  obj
}

seurat_list <- imap(sample_dirs, load_and_filter_sample)

# merge into one seurat object, labeled by sample of origin
seurat_209781 <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(sample_dirs))

seurat_209781
table(seurat_209781$condition)
table(seurat_209781$sample)

# save the merged object, this will save ~15m processing to this point
qs2::qs_save(seurat_209781, "data/processed/GSE209781_merged_filtered.qs2")


#############################################################################
# 3. std processing in Seurat - look for clusters
#############################################################################

# finish loading and scaling data
seurat_209781 <- JoinLayers(seurat_209781)

seurat_209781 <- NormalizeData(seurat_209781)
seurat_209781 <- FindVariableFeatures(seurat_209781, nfeatures = 2000)
seurat_209781 <- ScaleData(seurat_209781)

# PCA, decide how many to use
seurat_209781 <- RunPCA(seurat_209781, npcs = 30)
elbow_plot <- ElbowPlot(seurat_209781, ndims = 30)
ggsave("results/GSE209781_PCA_elbow_plot.png", plot = elbow_plot,
       width = 6, height = 4, dpi = 300)

# Clustering
seurat_209781 <- FindNeighbors(seurat_209781, dims = 1:10)
seurat_209781 <- FindClusters(seurat_209781, resolution = 0.5)
seurat_209781 <- RunUMAP(seurat_209781, dims = 1:10)

# first UMAP
umap_plot <- DimPlot(seurat_209781, reduction = "umap", label = TRUE)
ggsave("results/GSE209781_UMAP_clusters.png", plot = umap_plot,
       width = 6, height = 5, dpi = 300)

# overlay podocyte markers on UMAP
# some known podocyte marker genes are in utils.R
library(patchwork)
grid_plot <- FeaturePlot(seurat_209781,
                         features = c("PODXL", "THSD7A", "SYNPO", "ACTB"),
                         ncol = 2)
ggsave("results/GSE209781_podocyte_markers_grid.png", plot = grid_plot,
       width = 10, height = 8, dpi = 300)
# clusters 4, 9, 15 are likely podocytes. these are off to the side of the rest of the UMAP
# if really lucky, they are podocyte cell states

#############################################################################
# 4. map technical params onto clusters (UMAP)
#############################################################################

# pct mt is already in metadata.  add pct ribo
seurat_209781[["percent.ribo"]] <- PercentageFeatureSet(seurat_209781, pattern = "^RP[SL]")

# plot UMAP for tech aspects
p_sample <- DimPlot(seurat_209781, group.by = "sample") + ggtitle("Sample (replicate)")
p_mt     <- FeaturePlot(seurat_209781, features = "percent.mt") + ggtitle("% Mitochondrial")
p_ribo   <- FeaturePlot(seurat_209781, features = "percent.ribo") + ggtitle("% Ribosomal")
p_depth  <- FeaturePlot(seurat_209781, features = "nCount_RNA") + ggtitle("UMI count (depth)")

tech_grid <- wrap_plots(p_sample, p_mt, p_ribo, p_depth, ncol = 2)
ggsave("results/GSE209781_tech_validation_grid.png", plot = tech_grid,
       width = 12, height = 10, dpi = 300)

#############################################################################
# 5. recluster 4, 9, 15 - seem like podocytes by known marker genes
# no obvious technical bias driving these, and good distr. of disease vs healthy samples in them
#############################################################################

podocyte_clusters <- c("4", "9", "15")

seurat_podocytes <- subset(seurat_209781, idents = podocyte_clusters)
seurat_podocytes

# Recluster just this subset, independently
seurat_podocytes <- FindVariableFeatures(seurat_podocytes, nfeatures = 2000)
seurat_podocytes <- ScaleData(seurat_podocytes)
seurat_podocytes <- RunPCA(seurat_podocytes, npcs = 30)
ElbowPlot(seurat_podocytes, ndims = 30)
# seems to flatten ~15 PCs

# redo clustering and UMAP for the original clusters 4, 9, 15 pulled out of initial full dataset UMAP
seurat_podocytes <- FindNeighbors(seurat_podocytes, dims = 1:15)
seurat_podocytes <- FindClusters(seurat_podocytes, resolution = 0.5)
seurat_podocytes <- RunUMAP(seurat_podocytes, dims = 1:15)

podo_umap <- DimPlot(seurat_podocytes, reduction = "umap", label = TRUE)
ggsave("results/GSE209781_podocyte_subcluster_UMAP.png", plot = podo_umap,
       width = 7, height = 6, dpi = 300)

# make technical plots for this reclustering
p_sample <- DimPlot(seurat_podocytes, group.by = "sample") + ggtitle("Sample (replicate)")
p_mt     <- FeaturePlot(seurat_podocytes, features = "percent.mt") + ggtitle("% Mitochondrial")
p_ribo   <- FeaturePlot(seurat_podocytes, features = "percent.ribo") + ggtitle("% Ribosomal")
p_depth  <- FeaturePlot(seurat_podocytes, features = "nCount_RNA") + ggtitle("UMI count (depth)")

tech_grid_podo <- wrap_plots(p_sample, p_mt, p_ribo, p_depth, ncol = 2)
ggsave("results/GSE209781_podocyte_subcluster_tech_validation_grid.png", plot = tech_grid_podo,
       width = 12, height = 10, dpi = 300)


podo_marker_plots <- lapply(PODOCYTE_MARKERS, function(gene) {
  FeaturePlot(seurat_podocytes, features = gene) +
    ggtitle(gene) +
    NoLegend()
})

marker_grid <- wrap_plots(podo_marker_plots, ncol = 3, nrow = 3)
ggsave("results/GSE209781_podocyte_subcluster_all_markers_grid.png", plot = marker_grid,
       width = 15, height = 15, dpi = 300)

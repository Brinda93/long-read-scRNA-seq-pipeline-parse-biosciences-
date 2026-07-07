
# Complete Single-Cell RNA-seq Analysis Pipeline
# Author: Analysis for Parse Biosciences data
# Date: August 2025

# =============================================================================
# 1. LOAD REQUIRED LIBRARIES
# =============================================================================

library(Seurat)
library(dplyr)
library(Matrix)
library(tidyverse)
library(data.table)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# =============================================================================
# 2. SET WORKING DIRECTORY AND LOAD DATA
# =============================================================================

setwd("/gpfs/Labs/Dovat/bxp5423/output_all_combined_resutls/all-sample/DGE_unfiltered")
data_dir <- "/gpfs/Labs/Dovat/bxp5423/output_all_combined_resutls/all-sample/DGE_unfiltered"

# Load and transpose matrix
cat("Loading count matrix...\n")
counts <- t(readMM(file.path(data_dir, "count_matrix.mtx")))

# Load gene metadata
cat("Loading gene metadata...\n")
genes <- fread(file.path(data_dir, "all_genes.csv"), header = TRUE)

# Load cell metadata (adjust filename as needed)
cat("Loading cell metadata...\n")
# Check what files you have first:
cat("Files in directory:\n")
print(list.files(data_dir))

# Adjust this line based on your actual cell metadata file:
cells <- fread(file.path(data_dir, "all_cells.csv"), header = TRUE)  # Adjust filename

# Verify dimensions
cat("Verifying dimensions...\n")
cat("Count matrix dimensions:", dim(counts), "\n")
cat("Number of genes:", nrow(genes), "\n")
cat("Number of cells:", nrow(cells), "\n")

stopifnot(nrow(counts) == nrow(genes))
stopifnot(ncol(counts) == nrow(cells))

# Assign gene names (rows) and cell names (columns)
rownames(counts) <- make.unique(genes$gene_name)
colnames(counts) <- cells$bc_wells  # Adjust column name as needed

# =============================================================================
# 3. CREATE SEURAT OBJECT
# =============================================================================

cat("Creating Seurat object...\n")
seurat_obj <- CreateSeuratObject(
  counts = counts, 
  project = "ParseProject",
  min.cells = 3,      # Filter genes expressed in < 3 cells
  min.features = 200  # Filter cells with < 200 genes
)

# Add metadata
seurat_obj <- AddMetaData(seurat_obj, metadata = cells)

cat("Initial Seurat object created with", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes\n")

# =============================================================================
# 4. QUALITY CONTROL METRICS
# =============================================================================

cat("Calculating QC metrics...\n")

# Calculate mitochondrial gene percentage (mouse: mt-, human: MT-)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# Calculate ribosomal gene percentage
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rpl|^Rps")

# Calculate hemoglobin gene percentage (often high in blood samples)
seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Hb[^(p)]")

# Summary statistics
cat("QC Metrics Summary:\n")
cat("Mitochondrial genes:\n")
print(summary(seurat_obj$percent.mt))
cat("Ribosomal genes:\n")
print(summary(seurat_obj$percent.rb))

# =============================================================================
# 5. VISUALIZE QC METRICS
# =============================================================================

cat("Creating QC plots...\n")

# Violin plots
p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p1)

# Feature scatter plots
p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
p3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(p2 + p3)

# =============================================================================
# 6. QUALITY CONTROL FILTERING
# =============================================================================

cat("Applying quality control filters...\n")

# Show cells before filtering
cat("Cells before filtering:", ncol(seurat_obj), "\n")

# Apply filters - adjust thresholds based on your data
seurat_obj <- subset(seurat_obj, subset = 
  nFeature_RNA > 200 &        # Remove cells with too few genes
  nFeature_RNA < 4000 &       # Remove potential doublets
  nCount_RNA < 25000 &        # Remove cells with excessive UMIs
  percent.mt < 20             # Remove cells with high mitochondrial content
)

cat("Cells after filtering:", ncol(seurat_obj), "\n")

# Re-plot after filtering
p4 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p4)

# =============================================================================
# 7. NORMALIZATION
# =============================================================================

cat("Normalizing data...\n")
seurat_obj <- NormalizeData(
  seurat_obj, 
  normalization.method = "LogNormalize", 
  scale.factor = 10000
)

# =============================================================================
# 8. FIND HIGHLY VARIABLE FEATURES
# =============================================================================

cat("Finding highly variable features...\n")
seurat_obj <- FindVariableFeatures(
  seurat_obj, 
  selection.method = "vst", 
  nfeatures = 2000
)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)
cat("Top 10 most variable features:\n")
print(top10)

# Plot variable features
p5 <- VariableFeaturePlot(seurat_obj)
p6 <- LabelPoints(plot = p5, points = top10, repel = TRUE)
print(p6)

# =============================================================================
# 9. SCALE DATA AND REGRESS OUT UNWANTED VARIATION
# =============================================================================

cat("Scaling data and regressing out unwanted variation...\n")

# Get all genes for scaling (optional: use only variable features)
all.genes <- rownames(seurat_obj)

# Scale data and regress out mitochondrial gene effects and cell cycle
seurat_obj <- ScaleData(
  seurat_obj, 
  features = all.genes,
  vars.to.regress = c("percent.mt", "nCount_RNA")
)

# =============================================================================
# 10. PRINCIPAL COMPONENT ANALYSIS (PCA)
# =============================================================================

cat("Running PCA...\n")
seurat_obj <- RunPCA(
  seurat_obj, 
  features = VariableFeatures(object = seurat_obj)
)

# Visualize PCA results
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)

# PCA plots
p7 <- DimPlot(seurat_obj, reduction = "pca")
p8 <- VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
print(p7 + p8)

# Heatmap of PCs
DimHeatmap(seurat_obj, dims = 1:6, cells = 500, balanced = TRUE)

# =============================================================================
# 11. DETERMINE DIMENSIONALITY
# =============================================================================

cat("Determining dimensionality...\n")

# Elbow plot to determine number of PCs to use
p9 <- ElbowPlot(seurat_obj, ndims = 50)
print(p9)

# JackStraw analysis (more computationally intensive but more precise)
# seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
# seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20)
# JackStrawPlot(seurat_obj, dims = 1:15)

# =============================================================================
# 12. CLUSTERING
# =============================================================================

cat("Performing clustering...\n")

# Find neighbors using first 20 PCs (adjust based on elbow plot)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

# Find clusters at different resolutions
resolutions <- c(0.1, 0.3, 0.5, 0.8, 1.0, 1.2)

for (res in resolutions) {
  seurat_obj <- FindClusters(seurat_obj, resolution = res)
}

# Use resolution 0.5 as default
Idents(seurat_obj) <- "RNA_snn_res.0.5"

cat("Number of clusters at resolution 0.5:", length(levels(seurat_obj)), "\n")

# =============================================================================
# 13. NON-LINEAR DIMENSIONAL REDUCTION (UMAP/tSNE)
# =============================================================================

cat("Running UMAP...\n")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

cat("Running tSNE...\n")
seurat_obj <- RunTSNE(seurat_obj, dims = 1:20)

# =============================================================================
# 14. VISUALIZATION
# =============================================================================

cat("Creating visualization plots...\n")

# UMAP plots
p10 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + 
       NoLegend() + ggtitle("UMAP - Clusters")

p11 <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident", pt.size = 0.5)

print(p10 + p11)

# tSNE plots
p12 <- DimPlot(seurat_obj, reduction = "tsne", label = TRUE, pt.size = 0.5) + 
       NoLegend() + ggtitle("tSNE - Clusters")

print(p12)

# Feature plots for QC metrics
p13 <- FeaturePlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p13)

# =============================================================================
# 15. CLUSTER RESOLUTION COMPARISON
# =============================================================================

cat("Comparing different clustering resolutions...\n")

# Create a plot showing different resolutions
p14 <- DimPlot(seurat_obj, group.by = "RNA_snn_res.0.1", label = TRUE) + ggtitle("Resolution 0.1")
p15 <- DimPlot(seurat_obj, group.by = "RNA_snn_res.0.3", label = TRUE) + ggtitle("Resolution 0.3")
p16 <- DimPlot(seurat_obj, group.by = "RNA_snn_res.0.5", label = TRUE) + ggtitle("Resolution 0.5")
p17 <- DimPlot(seurat_obj, group.by = "RNA_snn_res.0.8", label = TRUE) + ggtitle("Resolution 0.8")

print((p14 + p15) / (p16 + p17))

# =============================================================================
# 16. FIND CLUSTER BIOMARKERS
# =============================================================================

cat("Finding cluster biomarkers...\n")

# Find markers for every cluster compared to all remaining cells
seurat_obj.markers <- FindAllMarkers(
  seurat_obj, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25
)

# Show top markers per cluster
top_markers <- seurat_obj.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

print("Top 2 markers per cluster:")
print(top_markers)

# =============================================================================
# 17. VISUALIZE TOP MARKERS
# =============================================================================

cat("Visualizing top marker genes...\n")

# Get top 5 markers for first few clusters
top5 <- seurat_obj.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Heatmap of top markers
p18 <- DoHeatmap(seurat_obj, features = top5$gene) + NoLegend()
print(p18)

# Feature plots of selected markers
if (nrow(top_markers) >= 4) {
  p19 <- FeaturePlot(seurat_obj, features = top_markers$gene[1:4], ncol = 2)
  print(p19)
}

# Violin plots of selected markers
if (nrow(top_markers) >= 4) {
  p20 <- VlnPlot(seurat_obj, features = top_markers$gene[1:4], ncol = 2)
  print(p20)
}

# =============================================================================
# 18. SAVE RESULTS
# =============================================================================

cat("Saving results...\n")

# Save Seurat object
saveRDS(seurat_obj, file = "seurat_parse_object_complete.rds")

# Save marker genes
write.csv(seurat_obj.markers, file = "cluster_markers.csv", row.names = FALSE)

# Save cluster assignments
cluster_assignments <- data.frame(
  cell_id = colnames(seurat_obj),
  cluster = Idents(seurat_obj),
  UMAP_1 = seurat_obj@reductions$umap@cell.embeddings[,1],
  UMAP_2 = seurat_obj@reductions$umap@cell.embeddings[,2]
)
write.csv(cluster_assignments, file = "cluster_assignments.csv", row.names = FALSE)

# =============================================================================
# 19. SUMMARY STATISTICS
# =============================================================================

cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total cells after QC filtering:", ncol(seurat_obj), "\n")
cat("Total genes:", nrow(seurat_obj), "\n")
cat("Number of clusters (resolution 0.5):", length(levels(seurat_obj)), "\n")
cat("Number of highly variable features:", length(VariableFeatures(seurat_obj)), "\n")

# Cells per cluster
cat("\nCells per cluster:\n")
print(table(Idents(seurat_obj)))

# Top 10 most expressed genes
top_expressed <- sort(rowSums(seurat_obj@assays$RNA@counts), decreasing = TRUE)[1:10]
cat("\nTop 10 most expressed genes:\n")
print(top_expressed)

cat("\nAnalysis complete! Files saved:\n")
cat("- seurat_parse_object_complete.rds\n")
cat("- cluster_markers.csv\n")
cat("- cluster_assignments.csv\n")

# =============================================================================
# 20. OPTIONAL: ADDITIONAL ANALYSES
# =============================================================================

# Uncomment and run these sections for additional analyses:

# # Cell cycle scoring
# cc.genes <- readLines(con = "https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv")
# s.genes <- cc.genes[1:43]
# g2m.genes <- cc.genes[44:97]
# seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# # Doublet detection using DoubletFinder
# # library(DoubletFinder)
# # sweep.res.list <- paramSweep_v3(seurat_obj, PCs = 1:20, sct = FALSE)
# # sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
# # bcmvn <- find.pK(sweep.stats)

# # Trajectory analysis using Monocle3
# # library(monocle3)
# # cds <- as.cell_data_set(seurat_obj)
# # cds <- preprocess_cds(cds, num_dim = 50)
# # cds <- reduce_dimension(cds)

cat("\nPipeline completed successfully!\n")

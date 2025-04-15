# ---
# author: "Laura Jiménez Gracia"
# date: 2025-02-19
# ---
# This R script allows to group cells into different clusters 
# based on the similarity of their gene expression profiles,
# considering different resolution (granularity) levels.


# Pre-processing
## Load packages
library(tidyverse)
library(Seurat)
library(RColorBrewer)

## Paths
path_r_objects_in <- here::here("02_QC/results/R_objects")
path_r_objects_out <- here::here("03_clustering_annotation/results/R_objects")
path_r_figs <- here::here("03_clustering_annotation/results/figs")

# Functions
source(here::here("bin/utils.R"))

# Parameters
## Define resolutions to compute
resolutions_range <- c(0.01, 0.05, 0.1, 0.25, 0.3, 0.4, 0.5, 0.75, 1)
resolution_names <- as.character(purrr::map(resolutions_range, function(res) {
  paste0("RNA_snn_res.", res)}))

color_palette <- Polychrome::createPalette(40, c("#fc6060", "#74f774", "#7c7cfc"))
names(color_palette) <- NULL


# Load Data
## Load Seurat object
seurat_obj <- readRDS(paste0(path_r_objects_in, "/VEIGAEST_harmony_processed.rds"))
DefaultAssay(seurat_obj) <- "RNA"

# Clustering
## Determine the K-nearest neighbor graph
seurat_obj <- FindNeighbors(
  seurat_obj,
  reduction = "harmony",
  dims = 1:30)

## Determine the clusters (Louvain algorithm) for multiple resolutions                                
seurat_obj <- FindClusters(
  seurat_obj,
  resolution = resolutions_range,
  verbose = FALSE)


# Save clustering results
saveRDS(seurat_obj, 
        file = paste0(path_r_objects_out, "/VEIGAEST_clustering_level1_resolutions.rds"))


# Clustering overview
gg_umap_cluster_resolution <- DimPlot(object = seurat_obj,
                                      group.by = resolution_names,
                                      label = TRUE,
                                      label.size = 3,
                                      cols = color_palette,
                                      ncol = 3
                                      ) & NoLegend()

# Save image
ggsave(filename = paste0(path_r_figs, "/VEIGAEST_clustering_level1_resolutions_umap.png"),
       plot = gg_umap_cluster_resolution,
       width = 15,
       height = 15)





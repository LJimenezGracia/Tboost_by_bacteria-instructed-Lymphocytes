# ---
# author: "Laura Jiménez Gracia"
# date: 2025-03-04
# ---
# This R script allows to perform a Differential Expression Analysis (DEA)
# between the different groups associated to a particular condition, and
# to find transcriptional signatures associated to that condition.


# Pre-processing
## Load packages
library(tidyverse)
library(Seurat)
library(ggrepel)

## Parameters
# Paths
path_r_objects_in <- here::here("03_clustering_annotation/results/R_objects")
path_r_results <- here::here("04_GEX_analysis/DEA_and_GSEA/results/DEA_lineages")

# Functions
source(here::here("bin/utils.R"))

## Load data
seurat_obj <- readRDS(paste0(path_r_objects_in, "/VEIGAEST_clustering_annotation_full_cleaned.rds"))
seurat_obj$cell_lineages <- as.character(seurat_obj$cell_lineages)
seurat_obj$cell_lineages[seurat_obj$cell_lineages == "Monocytes/Macrophages"] <- "Mono_Macro"
seurat_obj$cell_lineages <- as.factor(seurat_obj$cell_lineages)

# Define comparison parameters
comparison_name_list <- c("aCD3CD28_BEADS-PBS", "BACB-PBS", "BACT-PBS", "BACB-aCD3CD28_BEADS", "BACT-aCD3CD28_BEADS", "BACB-BACT")
group_1_list <- c("aCD3CD28_BEADS", "BACB", "BACT",  "BACB",  "BACT",  "BACB")
group_2_list <- c("PBS", "PBS", "PBS", "aCD3CD28_BEADS", "aCD3CD28_BEADS",  "BACT")

for (index in 1:length(comparison_name_list)) {
  comparison_name <- comparison_name_list[[index]]
  group_1 <- group_1_list[[index]]
  group_2 <- group_2_list[[index]]
  
  # Define comparison paths
  path_object <- paste0(path_r_results, "/VEIGAEST_DEgenes_lineages_", comparison_name, ".rds")
  path_file <- paste0(path_r_results, "/VEIGAEST_DEgenes_lineages_", comparison_name, ".xlsx")
  path_file_split <- paste0(path_r_results, "/VEIGAEST_DEgenes_lineages_", comparison_name, "_UPDOWN.xlsx")
  
  # DE analysis considering cell-types separately
  ## Subset data
  cells_to_keep <- (seurat_obj$treatment == group_1 | seurat_obj$treatment == group_2)
  seurat_obj_subset <- seurat_obj[, cells_to_keep]
  seurat_obj_subset$cell_lineages <- as.factor(as.character(seurat_obj_subset$cell_lineages))
  
  DE_genes_bycell <- purrr::map(levels(seurat_obj_subset$cell_lineages), function(celltype) {
    seurat_obj_subset <- seurat_obj_subset[, seurat_obj_subset$cell_lineages == celltype]    
    Idents(seurat_obj_subset) <- "treatment"
    
    ## Find DE genes
    DE_genes <- find_DEgenes_celltypes(
      seurat = seurat_obj_subset,
      ident_1 = group_1,
      ident_2 = group_2,
      test_de = "wilcox",
      threshold_pvaladj = 0.05,
      min_cells = 0
    )    
    
    })
  names(DE_genes_bycell) <- levels(seurat_obj_subset$cell_lineages)
  
  # Save results
  ## DE genes as RDS object
  saveRDS(DE_genes_bycell, file = path_object)
  
  ## Export list of DE genes
  openxlsx::write.xlsx(DE_genes_bycell, file = path_file)

  ## SPLIT up- & down-regulated
  DE_genes_split <- purrr::map(names(DE_genes_bycell), function(celltype) {
    # DEA
    DE_genes_up <- DE_genes_bycell[[celltype]] %>% 
      dplyr::filter(is_significant & avg_log2FC > 0.25) %>% 
      dplyr::select(gene, avg_log2FC, p_val_adj)
    
    DE_genes_down <- DE_genes_bycell[[celltype]] %>% 
      dplyr::filter(is_significant & avg_log2FC < -0.25)%>% 
      dplyr::select(gene, avg_log2FC, p_val_adj)
    
    DE_genes_list <- list(DE_genes_up, DE_genes_down)
    names(DE_genes_list) <- c("UP", "DOWN")
    return(DE_genes_list)
    }
    )
  names(DE_genes_split) <- names(DE_genes_bycell)
  DE_genes_split_list <- do.call(c, DE_genes_split)
  DE_genes_split_list <- DE_genes_split_list %>%
    discard(is.null)

  # Save results
  ## Export list of DE genes
  openxlsx::write.xlsx(DE_genes_split_list, file = path_file_split)
}

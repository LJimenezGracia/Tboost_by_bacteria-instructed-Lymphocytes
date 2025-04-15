# ---
# author: "Laura Jiménez Gracia"
# date: 2025-03-03
# ---data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==
# This R script allows to prepare 


# Pre-processing
## Load packages
library(tidyverse)
library(Seurat)

## Paths
path_r_tables_in <- here::here("03_clustering_annotation/results/tables")
path_r_tables_out <- here::here("04_GEX_analysis/compositional_analysis/results")

## Load data
seurat_obj_metadata <- read.csv(paste0(path_r_tables_in, "/VEIGAEST_metadata_full_cleaned.csv"))
seurat_obj_metadata

# Compute number of cells 

## Cell-lineage & sample (treatment)
celllineage_prop_df <- seurat_obj_metadata %>% 
  select(c("treatment", "cell_lineages")) %>%
  dplyr::count(cell_lineages, treatment) %>% # only computing number of cells, not percentages
  reshape2::dcast(treatment~cell_lineages) %>% # change rows to columns
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) # replacing NA by 0
head(celllineage_prop_df)

write_csv(celllineage_prop_df, file = paste0(path_r_tables_out, "/VEIGAEST_celllineages_counts.csv"))

## Cell-types & sample (treatment)
celltypesprop_df <- seurat_obj_metadata %>% 
  select(c("treatment", "cell_types")) %>%
  dplyr::count(cell_types, treatment) %>% # only computing number of cells, not percentages
  reshape2::dcast(treatment~cell_types) %>% # change rows to columns
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) # replacing NA by 0
head(celltypesprop_df)

write_csv(celltypesprop_df, file = paste0(path_r_tables_out, "/VEIGAEST_celltypes_counts.csv"))

# CD8+ T cell activation mediated by bacteria-instructed lymphocytes as antigen-presenting cells capturing engineered bacteria

This is the GitHub repository for the **CD8+ T cell activation mediated by bacteria-instructed lymphocytes as antigen-presenting cells capturing engineered bacteria** project.

The manuscript is currently under a review process.

## Abstract


## Code implementation

The repository is organized into the following folder tree, which contains all the necessary files and scripts to perform the detailed tasks and reproduce all our results.

* **01_cellranger_mapping** --> It includes an overview of the project data information, including which samples and 10X libraries were generated. Also, it contains the scripts needed to create a folder directory to perform the sequencing read mapping to the reference genome. 

* **02_QC** --> R markdown notebooks to perform the quality control and the first pre-processing, including data normalization, scaling, dimensionality reduction and integration.

* **03_clustering_annotation** --> All R markdown notebooks to perform a top-down clustering annotation approach, as well as scripts to find differential expressed markers for each clustering, and to assign a biological-relevant identity to each cluster.

* **04_GEX_analysis** --> All the code used to perform further downstream analysis on the GEX processed data. It includes all the scripts and R markdown notebooks to perform cell composition analysis,  differential expression analysis (DEA) followed by gene set enrichment analysis (GSEA), as well as gene signature evaluation. 

* **05_TCR_analysis** --> All the code used to perform the quality control and the first pre-processing on the TCR data (clonotype calling), as well as the downstream analysis together with the GEX data.


### Package versions

The (most important) packages and versions needed to reproduce the full analysis are listed below:

* CellRanger (v9.0.0) was used to mapped single-cell RNA-seq reads (10X Genomics) to the reference genome.

*--- in R (v4.0.5) ---*
* Seurat (v4.3.0)
* SeuratObject (v4.1.3)
* scRepertoire (v1.8.0)
* Harmony (v1.0)
* fgsea (v1.16)
* Ucell (v2.2.0)

*--- in Python (v3.9.19) ---*
* scanpy (v1.10.3)
* scCODA (v0.1.9)

## Data accessibility

* The complete raw data (FASTQ files) generated in this study, as well as the processed count matrices, have been submitted to the NCBI Gene Expression Omnibus (GEO) under accession number [XXXX](XXXX).
Data will be available upon publication.


## Code accessibility

You can easily download a copy of all the files contained in this repository by:

* Cloning the git repository using the following command in the terminal:

`git clone https://github.com/LJimenezGracia/Tboost_by_bacteria-instructed-Lymphocytes.git`

* Downloading a .ZIP archive [HERE](https://github.com/LJimenezGracia/Tboost_by_bacteria-instructed-Lymphocytes/archive/refs/heads/main.zip).

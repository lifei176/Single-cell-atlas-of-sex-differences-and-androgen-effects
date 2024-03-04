# Sex-differences-orchestrated-by-androgens-at-single-cell-resolution
Codes provided here are to define the DEGs, AASB-DEGs and significantly enriched biological pathways of each cell type across the 17 tissues.
1. Cellranger_pipeline.sh: The code script in Shell for cellranger pipeline which is executed on Linux System.
2. Data_processing.R: The code in R to perform QC and construct seurat object for downstream analysis.
3. DEG_definition.R: The code in R to define the DEGs between two different conditions (p-adjust < 0.05 and |Log2FC| > 0.5).
4. AASB-DEG_definition.R: The code in R to define androgen-associated sex-biased DEGs (AASB-DEGs) for each cell type across the 17 tissues, which comprises positive AASB-DEGs (the expression levels of which are male-biased and positively associated with androgens) and negative AASB-DEGs (the expression levels of which are female-biased and negatively associated with androgens).
5. Pathway_enrichment.R: The code in R to definine the significantly enriched biological pathways based on DEGs and AASB-DEGs respectively (p < 0.01 & q < 0.01).
# Downloading the data
# Requirements

##This project is for Shelby's LKB1 Project
##Project started 16DEC24
##Paper will FASTQ data was mined: DOI: 10.1002/path.6259; Cellular and molecular characteristics of stromal Lkb1 deficiency-induced gastrointestinal polyposis based on single-cell RNA sequencing

# ============================================================
# NOTES TO READER
# ============================================================
# The following code was used by the R. Jones Lab to process and 
# analyze publicly available single-cell RNA-seq data. 
#
# Data Source:
#   FASTQ files were not generated in-house by the R. Jones Lab.
#   They were mined from the following publication:
#   DOI: 10.1002/path.6259
#   "Cellular and molecular characteristics of stromal Lkb1 deficiency-
#   induced gastrointestinal polyposis based on single-cell RNA sequencing"
#
# Purpose:
#   - This code reflects the steps we used to parse and analyze data 
#     for Shelby's LKB1 Project (started 16DEC24).
#   - To the best of our knowledge, analyses were carried out in a manner 
#     consistent with the original authors.
#   - Some analyses shown here were used in the manuscript; other 
#     exploratory code sections were not.
#
# Code Disclaimer:
#   - This is working analysis code, not a polished software package.
#   - Some sections may appear redundant or exploratory; this was part 
#     of the iterative analysis process.
#   - We did not attempt to streamline or extensively clean the script 
#     beyond what was necessary for reproducibility.
#
# Contact:
#   - For questions about this analysis, please contact:
#       Brandon Oswald (brandon.oswald@vai.org)
#       or the senior author Rusty Jones.
#
# Thank you for your understanding.
# ============================================================

##I used the follwoing code to run cellranger and get informtion need to conduct follow on analysis: 
#!/bin/bash
#
#SBATCH --job-name=lkb1_single_cell
#SBATCH --ntasks-per-core=1
#SBATCH --time=10:00:00
#SBATCH --partition=mnp
#SBATCH --nodelist=compute124

cd /varidata/research/projects/rjones/BrandonData/brandon_data_backup/shelbysLKB1singlecell
export PATH=/primary/projects/mnp/brandon/tools/cellranger-7.0.1:$PATH
cellranger count --id=LKB1_WT_vs_HET_GEX_cellranger_output_LKB1_HETS1 --libraries=single_cell_files_to_run_cell_ranger/lkb1_wt_vs_het_libraries.csv --transcriptome=single_cell_files_to_run_cell_ranger/refdata-gex-mm10-2020-A/ --expect-cells=20000 --chemistry=SC3Pv3 --localcores=30 --localmem=256


###In R after running required code above. 
## Load Packages
library(Seurat)
library(SeuratObject)
library(dplyr)
library(patchwork)
library(SeuratData)
library(cowplot)
library(dplyr)
library(glmGamPoi)
library(ggplot2)
library(Nebulosa)
library(devtools)
library(EnhancedVolcano)
library(Pandas)
library(decontX)
library(presto)
library(ComplexHeatmap)
library(org.Mm.eg.db)
library(circlize)
library(ggsignif)
library(ggpubr)
library(clusterProfiler)
library(writexl)
library(scater)
library("scater")
library("scran")
library(ggplot2)

## Start the data processing Single cell Data in terminal for cell ranger
## Cell Ranger output files and took the raw matrix output.
## Ran 3 separate runs for each sample of WT vs LKB1 KO. 
## Will merge post and assign sample name. 


## Set Seed
set.seed(31012023)

# For each sample:
counts_wt <- Read10X(data.dir = "/varidata/research/projects/mnp/brandon/shelby_lkb1_singlecell/LKB1_WT_vs_HET_GEX_cellranger_output_LKB1_WT/outs/raw_feature_bc_matrix")
seurat_wt <- CreateSeuratObject(counts = counts_wt, project = "WT", min.cells = 3, min.features = 200)
seurat_wt$sampleID <- "WT"
seurat_wt

counts_hets1 <- Read10X(data.dir = "/varidata/research/projects/mnp/brandon/shelby_lkb1_singlecell/LKB1_WT_vs_HET_GEX_cellranger_output_LKB1_HETS1/outs/raw_feature_bc_matrix")
seurat_hets1 <- CreateSeuratObject(counts = counts_hets1, project = "HETS1", min.cells = 3, min.features = 200)
seurat_hets1$sampleID <- "HETS1"
seurat_hets1

counts_hets2 <- Read10X(data.dir = "/varidata/research/projects/mnp/brandon/shelby_lkb1_singlecell/LKB1_WT_vs_HET_GEX_cellranger_output_LKB1_HETS2/outs/raw_feature_bc_matrix")
seurat_hets2 <- CreateSeuratObject(counts = counts_hets2, project = "HETS2", min.cells = 3, min.features = 200)
seurat_hets2$sampleID <- "HETS2"
seurat_hets2

# Merge all samples together. 
pbmc <- merge(
  x = seurat_wt, 
  y = list(seurat_hets1, seurat_hets2), 
  add.cell.ids = c("WT", "HETS1", "HETS2")
)

pbmc
## Identify percentage of cells expressing mitochondrial genes, also features (genes), and counts.
## This is used for QC purposes. 

pbmc <- PercentageFeatureSet(pbmc, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, raster = FALSE)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

## Sort out low nCount, low RNA, and keep only singlets with mito percentage less than 20%
pbmc_exp2 <- subset(pbmc, subset = nCount_RNA >=200 & percent.mt < 20)
pbmc_exp2

## SCT Transform RNA data
DefaultAssay(pbmc_exp2) <- 'RNA'
pbmc_exp2 <- SCTransform(pbmc_exp2, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
pbmc_exp2

DefaultAssay(pbmc_exp2) <- 'SCT'
#pbmc_exp2 <- FindVariableFeatures(pbmc_exp2)
pbmc_exp2 <- ScaleData(pbmc_exp2) 
pbmc_exp2 <- RunPCA(pbmc_exp2, dims = 1:40, reduction.name = 'rna_PCA', verbose = FALSE)
pbmc_exp2

#Find Neighbors
my_seurat_object <- FindNeighbors(pbmc_exp2, dims = 1:40, reduction = 'rna_PCA')

#Find clusters
my_seurat_object <- FindClusters(my_seurat_object, resolution = 0.5)

#Run UMAP
my_seurat_object <- RunUMAP(my_seurat_object, dims = 1:40, reduction = 'rna_PCA')
my_seurat_object$orig.ident
DimPlot(my_seurat_object, reduction = "umap", label = TRUE)
FeaturePlot(my_seurat_object, features = "Cd8a", reduction = "umap", split.by = "orig.ident")
FeaturePlot(my_seurat_object, features = "Cd4", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Itk", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Stk11", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Lmo7", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Pdgfra", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Cd79a", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Cd300c2", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Itk", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Tnfrsf17", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Klre1", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Csf2", reduction = "umap")

##Find markers of interest
FindAllMarkers(my_seurat_object, only.pos = TRUE, min.pct = 0.15, logfc.threshold = 0.15)

for (numberset in numbersets){
  Markercluster <- FindMarkers(act_CD8_Tcells, ident.1 = numberset, min.pct = 0.05)
  markers_df <- data.frame(Gene = rownames(Markercluster), Markercluster)
  write_xlsx(markers_df, path=paste0("/data/home/brandon.oswald/act_cd8_naming/ClusterMarkers/",numberset,"clustermarkers.xlsx"))
}
PrepSCTFindMarkers(my_seurat_object, assay = "SCT")
Markercluster <- FindMarkers(my_seurat_object, ident.1 = "10", min.pct = 0.05)
?PrepSCTFindMarkers

DimPlot(
  my_seurat_object, 
  reduction = "umap", 
  group.by = "sampleID",
  cols = c("WT" = "green", "HETS1" = "orange", "HETS2" = "purple")
)

#Subset CD8 and CD4
Tcells <- subset(my_seurat_object, idents = c("10", "3"))
Tcells

## Diet Seurat, get rid of almost everything in DF
DefaultAssay(Tcells) <- 'RNA'
Tcells <-DietSeurat(Tcells, assays = c("RNA"))
Tcells

## SCT Transform RNA data
DefaultAssay(Tcells) <- 'RNA'
Tcells <- SCTransform(Tcells, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
Tcells

DefaultAssay(Tcells) <- 'SCT'
#pbmc_exp2 <- FindVariableFeatures(pbmc_exp2)
Tcells <- ScaleData(Tcells) 
Tcells <- RunPCA(Tcells, dims = 1:40, reduction.name = 'rna_PCA', verbose = FALSE)
Tcells

#Find Neighbors
Tcells <- FindNeighbors(Tcells, dims = 1:40, reduction = 'rna_PCA')

#Find clusters
Tcells <- FindClusters(Tcells, resolution = 0.5)

#Run UMAP
Tcells <- RunUMAP(Tcells, dims = 1:40, reduction = 'rna_PCA')
DimPlot(Tcells, reduction = "umap", label = TRUE)
DimPlot(Tcells, reduction = "umap", label = TRUE, split.by = "orig.ident")
DimPlot(Tcells, reduction = "umap", label = TRUE, cols = "HETS1")

##Feature plot key genes of interest 
FeaturePlot(Tcells, features = "Cd8a", reduction = "umap")
FeaturePlot(Tcells, features = "Cd4", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Itk", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Stk11", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Lmo7", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Pdgfra", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Cd79a", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Cd300c2", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Itk", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Tnfrsf17", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Klre1", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Csf2", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Csf2", reduction = "umap")
FeaturePlot(my_seurat_object, features = "Cd3e", reduction = "umap")
FeaturePlot(Tcells, features = "Il17a", reduction = "umap")
FeaturePlot(Tcells, features = "Il17f", reduction = "umap")
FeaturePlot(Tcells, features = "Rorc", reduction = "umap")
FeaturePlot(Tcells, features = "Il23r", reduction = "umap")
FeaturePlot(Tcells, features = "Il22", reduction = "umap")
FeaturePlot(Tcells, features = "Ifng", reduction = "umap")
FeaturePlot(Tcells, features = "Csf2", reduction = "umap")
FeaturePlot(Tcells, features = "Klre1", reduction = "umap")
FeaturePlot(Tcells, features = "Klrb1c", reduction = "umap")
FeaturePlot(Tcells, features = "Cd3e", reduction = "umap")

#Plot joint density plot
plot_density(Tcells, features = c("Il17a", "Il23r"), reduction = 'umap', joint = TRUE)
summary(Tcells$seurat_clusters)
saveRDS(my_seurat_object, file = "/data/home/brandon.oswald/QC_SC_WT_LKB1_HET.rds")

##Excluding HETS1 from analysis as original analysis (by original authors) did the same. 
my_seurat_object_noHETS1 <- subset(Tcells, subset = sampleID != "HETS1")


## Re-run follwing code with current cells after filtering 
DefaultAssay(my_seurat_object_noHETS1) <- "RNA"
my_seurat_object_noHETS1 <- SCTransform(
  my_seurat_object_noHETS1, 
  method = "glmGamPoi", 
  vars.to.regress = "percent.mt", 
  verbose = FALSE
)

DefaultAssay(my_seurat_object_noHETS1) <- "SCT"

# Scale data (if desired; SCTransform often replaces ScaleData, but it's okay to run if you like)
my_seurat_object_noHETS1 <- ScaleData(my_seurat_object_noHETS1)

# Run PCA
my_seurat_object_noHETS1 <- RunPCA(my_seurat_object_noHETS1, 
                                   dims = 1:40, 
                                   reduction.name = "rna_PCA", 
                                   verbose = FALSE)

# Find neighbors

my_seurat_object_noHETS1 <- FindNeighbors(my_seurat_object_noHETS1, 
                                          dims = 1:40, 
                                          reduction = "rna_PCA")

# Find clusters
my_seurat_object_noHETS1 <- FindClusters(my_seurat_object_noHETS1, resolution = 0.5)

# Run UMAP
my_seurat_object_noHETS1 <- RunUMAP(my_seurat_object_noHETS1, 
                                    dims = 1:40, 
                                    reduction = "rna_PCA")
my_seurat_object_noHETS12 <- my_seurat_object_noHETS1

DimPlot(my_seurat_object_noHETS1, reduction = "umap", label = TRUE, split.by = "orig.ident")
DimPlot(my_seurat_object_noHETS1, 
        reduction = "umap", 
        group.by = "sampleID")

# Example list of genes
# Define your gene list
my_genes <- c("Il17a", "Il17f", "Rorc", "Il22", "Il23r", "Ifng", "Csf2")

# Loop through each gene
for (gene in my_genes) {
  
  # Create the FeaturePlot for the current gene
  p <- FeaturePlot(
    object = my_seurat_object_noHETS1,
    features = gene,
    reduction = "umap"
  ) + 
    ggtitle(gene)  # Optional: Add the gene name as a title
  
  # Save the plot as a PNG file
  ggsave(
    filename = paste0("FeaturePlot_", gene, ".pdf"),
    plot = p,
    device = "pdf",
    width = 5,   # Adjust width as needed
    height = 5,  # Adjust height as needed
    dpi = 300
  )
}


library(Nebulosa)
library(ggplot2)

# Define the gene list
my_genes <- c("Il17a", "Il17f", "Rorc", "Il22", "Il23r", "Ifng", "Csf2""Cd4", "Cd8a", "Foxp3", "Il7r", "Sell", "Tcrg")

# Loop through each gene
for (gene in my_genes) {
  
  # Plot Nebulosa density for the gene on the UMAP
  p <- plot_density(
    my_seurat_object_noHETS1, 
    features = gene, 
    reduction = "umap", 
    joint = FALSE  # TRUE if you want a joint density of multiple features
  ) + 
    ggtitle(gene)
  
  # Save each plot as a PDF
  ggsave(
    filename = paste0("Nebulosa_density_", gene, ".pdf"),
    plot = p,
    device = "pdf",
    width = 5,
    height = 5,
    dpi = 300
  )
}

##Plot more information as needed
plot_density(my_seurat_object_noHETS1, features = "Ifng", reduction = "umap")
plot_density(object = my_seurat_object_noHETS1, features = c("Cd8a", "Cd4", "Il17a"), reduction = "umap")
plot_density(my_seurat_object_noHETS1, features = c("Il17a", "Il22", "Il23r"), joint = TRUE, reduction = "umap")
p <- plot_density(my_seurat_object_noHETS1, c("Cd8a", "Cd4"), joint = TRUE)
p + plot_layout(ncol = 1)


# Plot density for two genes
plot_density(my_seurat_object_noHETS1, features = c("Cd8a", "Cd4", "Cd44"), joint = TRUE)

#Plot umap of each cluster type
DimPlot(my_seurat_object_noHETS1, reduction = "umap", label = TRUE)

# If your clusters are stored as the active identity:
cell_counts <- table(my_seurat_object_noHETS1$sampleID, Idents(my_seurat_object_noHETS1))
print(cell_counts)


my_seurat_object_noHETS1$combined_clusters
library(dplyr)
##For combined clusters
cells_counts_df_combined <- my_seurat_object_noHETS1@meta.data %>%
  group_by(sampleID, combined_clusters) %>%
  summarise(CellCount = n(), .groups = "drop")

##For original clusters
cells_counts_df_orig <- my_seurat_object_noHETS1@meta.data %>%
  group_by(sampleID, seurat_clusters) %>%
  summarise(CellCount = n(), .groups = "drop")
cells_counts_df_orig

# View the data frame
print(cell_counts_df)
cells_counts_df_combined

library(ggplot2)

##Make a bar graph to analyze cell number differences 
ggplot(cell_counts_df, aes(x = seurat_clusters, y = CellCount, fill = sampleID)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Cluster", y = "Number of Cells", fill = "Sample") +
  theme_classic() +
  ggtitle("Cell Counts by Cluster and Sample")

##combined clusters as needed
ggplot(cell_counts_df, aes(x = combined_clusters, y = CellCount, fill = sampleID)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Cluster", y = "Number of Cells", fill = "Sample") +
  theme_classic() +
  ggtitle("Cell Counts by Cluster and Sample")

##plot with cell number on graph (orig)
ggplot(cells_counts_df_orig, aes(x = seurat_clusters, y = CellCount, fill = sampleID)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = CellCount), 
            position = position_dodge(width = 0.9), 
            vjust = -0.25,    # Adjust vertical position of the text
            size = 3) +       # Adjust text size if needed
  labs(x = "Cluster", y = "Number of Cells", fill = "Sample") +
  theme_classic() +
  ggtitle("Cell Counts by Cluster and Sample_Orig")

##plot combined clusters
##plot with cell number on graph (orig)
ggplot(cells_counts_df_combined, aes(x = combined_clusters, y = CellCount, fill = sampleID)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = CellCount), 
            position = position_dodge(width = 0.9), 
            vjust = -0.25,    # Adjust vertical position of the text
            size = 3) +       # Adjust text size if needed
  labs(x = "Cluster", y = "Number of Cells", fill = "Sample") +
  theme_classic() +
  ggtitle("Cell Counts by Cluster and Combined clusters")

# Reorder sampleID factor levels so that WT comes before HETS1
cell_counts_df$sampleID <- factor(cell_counts_df$sampleID, levels = c("WT", "HETS1"))

# Then plot your data
ggplot(cell_counts_df, aes(x = seurat_clusters, y = CellCount, fill = sampleID)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Cluster", y = "Number of Cells", fill = "Sample") +
  theme_classic() +
  ggtitle("Cell Counts by Cluster and Sample")

library(dplyr)
library(ggplot2)

# 1. Compute cell counts by sample and cluster
cell_counts_df <- my_seurat_object_noHETS1@meta.data %>%
  group_by(seurat_clusters, sampleID) %>%  # group by cluster and sample
  summarise(CellCount = n(), .groups = "drop")

# 2. Calculate the percent of cells from each sample within each cluster
percent_df <- cell_counts_df %>%
  group_by(seurat_clusters) %>%
  mutate(TotalCells = sum(CellCount),
         Percent = (CellCount / TotalCells) * 100) %>%
  ungroup()

# Optional: view the resulting data frame
print(percent_df)

# 3. Create a stacked bar graph showing the percent of cells per cluster contributed by each sample
ggplot(percent_df, aes(x = seurat_clusters, y = Percent, fill = sampleID)) +
  geom_bar(stat = "identity") +
  labs(x = "Cluster",
       y = "Percent of Cells",
       fill = "Sample",
       title = "Percent of Cells by Cluster and Sample") +
  theme_classic()

library(dplyr)
library(ggplot2)

# 1. Compute cell counts per sample and cluster
cell_counts_df <- my_seurat_object_noHETS1@meta.data %>%
  group_by(sampleID, seurat_clusters) %>%  # Group by sample and cluster
  summarise(CellCount = n(), .groups = "drop")

# 2. For each sample, compute the percentage contribution of each cluster
cell_counts_df <- cell_counts_df %>%
  group_by(sampleID) %>%
  mutate(TotalCells = sum(CellCount),
         Percent = (CellCount / TotalCells) * 100) %>%
  ungroup()

# (Optional) View the data
print(cell_counts_df)

# 3. Plot a stacked bar graph of cell composition by sample.
#    Each bar (for a given sample) sums to 100%, with segments representing clusters.
ggplot(cell_counts_df, aes(x = sampleID, y = Percent, fill = seurat_clusters)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample",
       y = "Percent of Cells",
       fill = "Cluster",
       title = "Cluster Composition (Percentage) by Sample") +
  theme_minimal()

# Remove cells belonging to cluster "8"
my_seurat_object_noHETS1 <- subset(my_seurat_object_noHETS1, subset = seurat_clusters !=8)
my_seurat_object_noHETS1$seurat_clusters

# Convert current active identity (assumed to be clusters) to character
my_seurat_object_noHETS1$combined_clusters <- as.character(Idents(my_seurat_object_noHETS1))

# For cells in clusters 1 or 7, set the new combined label (here, "1+7")
my_seurat_object_noHETS1$combined_clusters[my_seurat_object_noHETS1$combined_clusters %in% c("1", "7")] <- "1+7"

# Optionally convert to a factor for ordering/control in plotting
my_seurat_object$combined_clusters <- factor(my_seurat_object$combined_clusters)

# Check the new grouping:
table(my_seurat_object$combined_clusters)
# Check that the UMAP remains unchanged for the remaining cells
DimPlot(my_seurat_object_noHETS1, reduction = "umap", label = TRUE)

plot_density(my_seurat_object_noHETS1, features = c("Cd8a", "Cd44"), joint = TRUE)

##plot joint nebulosa plot of the following genes. 
joint_genes <- c("Il17a", "Il17f", "Il22", "Il23r")
library(Nebulosa)
??Nebulosa
library(ggplot2)
p <- plot_density(
  object = my_seurat_object_noHETS1,
  features = joint_genes,
  reduction = "umap",
  joint = TRUE
) +
  ggtitle("Joint Density Plot for Il17a, Il23r, Il22, Il23r")
p[[5]]
print(p)

##Save RDS file for new object
saveRDS(my_seurat_object_noHETS1, file = "/data/home/brandon.oswald/NO_HETS2_QC_SC_WT_LKB1_HET.rds")

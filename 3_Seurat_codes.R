# Loading necessary libraries
library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(ggplot2)

# Merge all 12 doublet-removed sample objects into one Seurat object, adding cell IDs with sample and zygosity information
samples_merged <- merge(x = RR1_SINGLETS, 
                          y = c(RR2_SINGLETS, RR2RD3_SINGLETS, RR3_SINGLETS, RH1_SINGLETS, RH2_SINGLETS, RH2RD3_SINGLETS, RH3_SINGLETS, HH1_SINGLETS, HH2_SINGLETS, HH2RD3_SINGLETS, HH3_SINGLETS),
                          add.cell.id = c("RR1_Normal", "RR2_Normal", "RR2RD3_Normal", "RR3_Normal", "RH1_Heterozygous", "RH2_Heterozygous", "RH2RD3_Heterozygous", "RH3_Heterozygous", "HH1_Homozygous", "HH2_Homozygous", "HH2RD3_Homozygous", "HH3_Homozygous"),  project = 'Day_48')

# Create a new column to store sample information
samples_merged$sample <- rownames(samples_merged@meta.data)

# Split the sample column into separate columns for Sample, Clone_Type, and Barcode
samples_merged@meta.data <- separate(samples_merged@meta.data, col = 'sample', into = c('Sample', 'Clone_Type', 'Barcode'), sep = '_')

# Split the merged object into a list of Seurat objects by sample for subsequent integration
# Normalize each dataset and run PCA for dimensionality reduction
samples_merged.list <- SplitObject(samples_merged, split.by = "Sample")
for (i in 1:length(samples_merged.list)) {
  samples_merged.list[[i]] <- SCTransform(samples_merged.list[[i]], vars.to.regress = c("mitoPercent") , method = "glmGamPoi", vst.flavor = "v2")
}
for (i in 1:length(samples_merged.list)) {
  samples_merged.list[[i]] <- RunPCA(samples_merged.list[[i]], npcs = 50)
}

# Select features for integration and prepare the datasets for integration
features <- SelectIntegrationFeatures(object.list = samples_merged.list, nfeatures = 3000)
samples_merged.list <- PrepSCTIntegration(object.list = samples_merged.list, anchor.features = features)
clones_anchors <- FindIntegrationAnchors(object.list = samples_merged.list, normalization.method = "SCT", anchor.features = features)

# Integrate the datasets into a single Seurat object
samples_merged_integrated <- IntegrateData(anchorset = clones_anchors, normalization.method = "SCT")
# Run PCA and UMAP for visualization, and find neighbors and clusters in the integrated data
samples_merged_integrated <- RunPCA(samples_merged_integrated)
samples_merged_integrated <- RunUMAP(samples_merged_integrated, reduction = "pca", dims = 1:30)
samples_merged_integrated <- FindNeighbors(samples_merged_integrated, reduction = "pca", dims = 1:30)
samples_merged_integrated <- FindClusters(samples_merged_integrated, resolution = c(0.5, 0.2))
samples_merged_integrated <- PrepSCTFindMarkers(samples_merged_integrated, assay = "SCT")

# Visualize clusters using UMAP plots at different resolutionsDimPlot(samples_merged_integrated, reduction = "umap", group.by = "integrated_snn_res.0.5", label = T, pt.size = 0.6)
DimPlot(samples_merged_integrated, reduction = "umap", group.by = "integrated_snn_res.0.2", label = T, pt.size = 0.6)

# Add metadata columns to assign zygosity and sample identity to each cluster

samples_merged_integrated$cell.zygosity <- paste(samples_merged_integrated$Clone_Type, samples_merged_integrated$integrated_snn_res.0.2, 
                                                      sep = "_")

samples_merged_integrated$cell.sample <- paste(samples_merged_integrated$Sample, samples_merged_integrated$integrated_snn_res.0.2, 
                                                   sep = "_")

# Identify cluster-specific markers using the negative binomial test

cluster0_RES0.2 <- FindMarkers(samples_merged_integrated, assay="SCT", slot = "counts", ident.1 = 0, logfc.threshold= 0.25, test.use = "negbinom")
cluster1_RES0.2 <- FindMarkers(samples_merged_integrated, assay="SCT", slot = "counts", ident.1 = 1, logfc.threshold= 0.25, test.use = "negbinom")
cluster2_RES0.2 <- FindMarkers(samples_merged_integrated, assay="SCT", slot = "counts", ident.1 = 2, logfc.threshold= 0.25, test.use = "negbinom")
cluster3_RES0.2 <- FindMarkers(samples_merged_integrated, assay="SCT", slot = "counts", ident.1 = 3, logfc.threshold= 0.25, test.use = "negbinom")
cluster4_RES0.2 <- FindMarkers(samples_merged_integrated, assay="SCT", slot = "counts", ident.1 = 4, logfc.threshold= 0.25, test.use = "negbinom")
cluster5_RES0.2 <- FindMarkers(samples_merged_integrated, assay="SCT", slot = "counts", ident.1 = 5, logfc.threshold= 0.25, test.use = "negbinom")
cluster6_RES0.2 <- FindMarkers(samples_merged_integrated, assay="SCT", slot = "counts", ident.1 = 6, logfc.threshold= 0.25, test.use = "negbinom")

# Annotate clusters with cell type labels based on identified marker genes
samples_merged_integrated$annotation <- NA
samples_merged_integrated$annotation[samples_merged_integrated$seurat_clusters == "0"] <- "SC-alpha"
samples_merged_integrated$annotation[samples_merged_integrated$seurat_clusters == "1"] <- "SC-beta"
samples_merged_integrated$annotation[samples_merged_integrated$seurat_clusters == "2"] <- "SC-EC"
samples_merged_integrated$annotation[samples_merged_integrated$seurat_clusters == "3"] <- "DPP4⁺/ALDH1A1⁺"
samples_merged_integrated$annotation[samples_merged_integrated$seurat_clusters == "4"] <- "SC-ductal"
samples_merged_integrated$annotation[samples_merged_integrated$seurat_clusters == "5"] <- "Proliferating cells"
samples_merged_integrated$annotation[samples_merged_integrated$seurat_clusters == "6"] <- "SC-delta"

# Save the final integrated and annotated object as an RDS file
saveRDS(samples_merged_integrated, file = "samples_merged_integrated.rds")




# Loading necessary libraries
library(SoupX)
library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)
library(DropletUtils)

# SoupX processing for removing ambient RNA from each dataset
# Reading in the matrices
filtered.matrix <- Read10X_h5("filtered_feature_bc_matrix.h5", use.names = T)
raw.matrix  <- Read10X_h5("raw_feature_bc_matrix.h5", use.names = T)

# Creating a “SoupChannel” for SoupX. 
soup.channel  <- SoupChannel(raw.matrix, filtered.matrix)

# Creating a Seurat object from the sparce matrix of filtered cell ranger output:
RR1_DATA  <- CreateSeuratObject(counts =filtered.matrix)

# Generating cluster information for SoupX
RR1_DATA    <- NormalizeData(RR1_DATA) 
RR1_DATA    <- FindVariableFeatures(RR1_DATA, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(RR1_DATA)
RR1_DATA    <- ScaleData(RR1_DATA, features = all.genes)
RR1_DATA    <- RunPCA(RR1_DATA, features = VariableFeatures(object = RR1_DATA))
RR1_DATA    <- FindNeighbors(RR1_DATA, dims = 1:30)
RR1_DATA    <- FindClusters(RR1_DATA)
RR1_DATA    <- RunUMAP(RR1_DATA, dims = 1:30)

# Adding cluster information to the channel using setClusters. setDR is useful for visualizations.
meta    <-RR1_DATA@meta.data
umap    <-RR1_DATA@reductions$umap@cell.embeddings
soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
soup.channel  <- setDR(soup.channel, umap)

# Running SoupX function to calculate ambient RNA profile
soup.channel  <- autoEstCont(soup.channel)

# Adjusting counts
adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)

# Outputting corrected read counts
DropletUtils::write10xCounts("ambientRNA_removed_DATA_RR1", adj.matrix)

# Above steps repeated for the remaining data sets to generate the following counts files
#ambientRNA_removed_DATA_RR2
#ambientRNA_removed_DATA_RR2RD3
#ambientRNA_removed_DATA_RR3
#ambientRNA_removed_DATA_RH1
#ambientRNA_removed_DATA_RH2
#ambientRNA_removed_DATA_RH2RD3
#ambientRNA_removed_DATA_RH3
#ambientRNA_removed_DATA_HH1
#ambientRNA_removed_DATA_HH2
#ambientRNA_removed_DATA_HH2RD3
#ambientRNA_removed_DATA_HH3

# Ambient RNA removed data sets were loaded into R as objects and subjected for doublet removal
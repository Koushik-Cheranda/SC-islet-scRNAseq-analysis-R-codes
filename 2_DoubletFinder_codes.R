# Loading necessary libraries
library(DoubletFinder)
library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)

# Data processed with SoupX is now being used for doublet removal with DoubletFinder

RR1.data <- Read10X(data.dir = "ambientRNA_removed_DATA_RR1")
RR1.data <- CreateSeuratObject(counts = RR1.data, project = "DAY-48", min.cells = 3, min.features = 200)
RR1.data$mitoPercent <- PercentageFeatureSet(object = RR1.data, pattern = "^MT-")

# Generate QC plots before applying any filtering criteria
VlnPlot(RR1.data, features = c("nCount_RNA", "nFeature_RNA", "mitoPercent"), ncol = 3)

# Filter the dataset using defined cutoffs for quality metrics
RR1.data.filt <- subset(
  x = RR1.data,
  subset = nFeature_RNA >= 2500 &
    nFeature_RNA <= 8000 &
    nCount_RNA >= 5000 &
    GAPDH >= 1 &
    ACTB >= 1 &
    mitoPercent <= 25)

# Generate QC plots after filtering for assessment

VlnPlot(RR1.data.filt, features = c("nCount_RNA", "nFeature_RNA", "mitoPercent"), ncol = 3)

# Pre-process the Seurat object using SCTransform
RR1.data.filt <- SCTransform(RR1.data.filt, vars.to.regress = "mitoPercent", method = "glmGamPoi", vst.flavor = "v2")
RR1.data.filt <- RunPCA(RR1.data.filt, npcs = 50)
ElbowPlot(RR1.data.filt)
RR1.data.filt <- RunUMAP(RR1.data.filt, dims = 1:20)
RR1.data.filt <- FindNeighbors(RR1.data.filt, dims = 1:20)
RR1.data.filt <- FindClusters(RR1.data.filt, resolution = 0.2)

# Detect potential doublets for each sample

## Identify optimal pK value (without ground-truth data)
sweep.res.list_RR1 <- paramSweep(RR1.data.filt, PCs = 1:20, sct = TRUE)
sweep.stats_RR1 <- summarizeSweep(sweep.res.list_RR1, GT = FALSE)
bcmvn_RR1 <- find.pK(sweep.stats_RR1)

ggplot(bcmvn_RR1, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

# Select the pK value that maximizes the BCmetric for optimal doublet detection
pK <- bcmvn_RR1 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))


## Estimate the proportion of homotypic doublets 
annotations <- RR1.data.filt@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           

# Calculate the expected number of doublets

nExp_poi <- round(0.075*nrow(RR1.data.filt@meta.data))  ## Assuming 7.5% doublet formation rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder to identify doublets in the dataset
RR1.data.filt <- doubletFinder(RR1.data.filt, 
                                  PCs = 1:20, 
                                  pN = 0.25, 
                                  pK = pK, 
                                  nExp = nExp_poi.adj,
                                  reuse.pANN = FALSE, sct = TRUE)



# Visualize doublets and singlets on UMAP, setting group.by to the appropriate classification
names(RR1.data.filt@meta.data)
table(RR1.data.filt@meta.data$DF.classifications_0.25_0.02_195) 
DimPlot(RR1.data.filt, reduction = 'umap', group.by = "DF.classifications_0.25_0.02_195")

# Subset the Seurat object to separate doublets and singlets for further analysis
RR1_doublets <- subset(x = RR1.data.filt, subset = DF.classifications_0.25_0.02_195 == "Doublet") # Inspect doublets
RR1_singlets <- subset(x = RR1.data.filt, subset = DF.classifications_0.25_0.02_195 == "Singlet") # Use singlets for downstream 
# Apply the same workflow across all datasets to isolate sample-specific singlet populations.

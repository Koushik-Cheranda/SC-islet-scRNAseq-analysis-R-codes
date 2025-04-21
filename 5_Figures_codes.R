# Loading necessary libraries
library(Seurat)
library(ggplot2)
library(gghalves)
library(scCustomize)
library(viridis)

# Load the processed Seurat object and ensure consistent cluster annotation ordering across all figures

my_levels <- c("SC-α", "SC-β", "SC-δ", "SC-EC", "SC-ductal",  "DPP4⁺/ALDH1A1⁺", "Proliferating cells")
samples_merged_integrated@meta.data$annotation <- factor(x = samples_merged_integrated@meta.data$annotation, levels = my_levels)

# Figure 6A: Generate UMAP plot with cells color-coded by cluster to show the seven distinct cell types
Idents(samples_merged_integrated) <- "annotation"
tiff("wholedata_clusters0.2res.tiff", width = 15, height = 10, units = "in", res = 600)
DimPlot(samples_merged_integrated, reduction = "umap", label = F, pt.size = 0.4, cols = c("dodgerblue4", "red1", "blue3","limegreen", "darkorange3", "dodgerblue1", "darkorchid1"))
dev.off()

# Figure 6B: Feature plots displaying the expression levels of key marker genes (log1P SCT normalized counts)

tiff("featureplot_markers.tiff", width = 10, height = 8, units = "in", res = 600)
FeaturePlot(samples_merged_integrated, 
                 reduction = "umap", 
                 features = c("INS", "GCG", "TPH1", "ABCC8", "DPP4", "KRT19", "SST", "ALDH1A1", "MKI67"), 
                 order = TRUE, pt.size = NULL,
                 min.cutoff = 'q20', combine = TRUE, label = FALSE) & theme(
                   plot.title = element_text(size = 10, face = "bold"),
                   axis.text = element_text(size = 10), 
                   axis.title = element_text(size = 10, hjust = 1), 
                   legend.text = element_text(size = 10),
                   axis.line = element_line(size = 0.25)
                 )
dev.off()

# Figure 6C: Heatmap showing scaled expression of top differentially expressed genes for each cluster compared to others

heatmap_genes <- c("INS", "IAPP", "ERO1B", "DLK1", "ACVR1C", "PCSK1", "EEF1A2", "PCDH7", "HADH", "NEFM", "GCG", 
                   "IGFBP2", "TTR", "SERPINA1", "SPINK1", "AGT", "ARX", "SERPINI1", "GLS", "CLU", "DPP4", "GC", 
                   "ALB", "ALDH1A1", "PDK3", "KCTD12", "NR2F1", "SSTR2", "PCDH17", "TPH1", "COL5A2", "DNAJC12", 
                   "FEV", "DDC", "SYT13", "CBLN1", "AC068587.4", "CHGA", "CAMK1D", "SST", "GHRL", "ACSL1", 
                   "HHEX", "SFRP1", "VTN", "CRH", "APOA1", "ST3GAL1", "OLFML3", "COL3A1", "KRT19", "HIST1H4C", 
                   "CENPF", "TOP2A", "TUBA1B", "RRM2", "MKI67", "H2AFZ", "HMGB2", "PTTG1", "TUBB")

tiff("Heatmap_wholedata6*5_2.tiff", height = 6, width = 5, units = "in", res = 600)

DoHeatmap(
  samples_merged_integrated,
  features = heatmap_genes, 
  group.by = "ident",
  group.bar = TRUE, 
  slot = "scale.data",
  assay = "SCT",
  label = FALSE,
  raster = TRUE,
  draw.lines = TRUE,
  group.bar.height = 0.02, 
  angle = 0, 
  group.colors = c("red1", "dodgerblue4", "dodgerblue1", "limegreen", "blue3", "darkorange3",  "darkorchid1")
) + theme(axis.text = element_text(size = 7, face = "bold"))
dev.off()

# Figures 6F-G: Violin plots of core beta and alpha cell markers (log1P SCT normalized counts)

tiff("beta_markers_violinplot.tiff", width = 8, height = 6, units = "in", res = 600)
Stacked_VlnPlot(seurat_object = samples_merged_integrated, features = c("INS", "IAPP","DLK1","HADH","PCSK1","NKX6-1","PDX1"), x_lab_rotate = TRUE, colors_use = c("dodgerblue4", "red1", "blue3","limegreen", "darkorange3", "dodgerblue1", "darkorchid1"),
                plot_spacing = 0.3) + theme(axis.text=element_text(size=7, face = "bold"))
dev.off()

tiff("alpha_markers_violinplot.tiff", width = 8, height = 6, units = "in", res = 600)
Stacked_VlnPlot(seurat_object = samples_merged_integrated,  features = c("GCG", "ARX","ALDH1A1","GC","TTR","DPP4","SERPINE2"), x_lab_rotate = TRUE, colors_use = c("dodgerblue4", "red1", "blue3","limegreen", "darkorange3", "dodgerblue1", "darkorchid1"),
                plot_spacing = 0.3) + theme(axis.text=element_text(size=7, face = "bold"))
dev.off()

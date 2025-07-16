# Loading necessary libraries
library(Seurat)
library(ggplot2)
library(gghalves)
library(scCustomize)


# Load the processed Seurat object and ensure consistent cluster annotation ordering across all figures

my_levels <- c("SC-alpha", "SC-beta",  "SC-delta", "SC-EC", "SC-ductal", "DPP4⁺/ALDH1A1⁺", "Proliferating cells")
samples_merged_integrated@meta.data$annotation <- factor(x = samples_merged_integrated@meta.data$annotation, levels = my_levels)

# Figure 6A: Generate UMAP plot with cells color-coded by cluster to show the seven distinct cell types
Idents(samples_merged_integrated) <- "annotation"
tiff("wholedata_clusters0.2res.tiff", width = 15, height = 10, units = "in", res = 1200)
DimPlot(samples_merged_integrated, reduction = "umap", label = FALSE, pt.size = 0.4, 
        cols = c("dodgerblue4", "red1", "blue3", "limegreen", "darkorange3", "dodgerblue1", "darkorchid1")) & 
        theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 18, hjust = 1),
        axis.line = element_line(size = 0.50))
dev.off()

# Figure 6B: Feature plots displaying the expression levels of key marker genes (log1P SCT normalized counts)
tiff("featureplot_markers.tiff", width = 10, height = 8, units = "in", res = 600)
FeaturePlot(samples_merged_integrated, 
            reduction = "umap", 
            features = c("INS", "GCG", "TPH1", "KRT19", "SST", "MKI67"), 
            order = TRUE, pt.size = NULL,
            min.cutoff = 'q20', combine = TRUE, label = FALSE, ncol = 3) & theme(
            plot.title = element_text(size = 24, face = "plain"),
            axis.text = element_text(size = 0), axis.ticks = element_blank(),
            axis.title = element_text(size = 18, hjust = 1), 
            legend.text = element_text(size = 18),
            axis.line = element_line(size = 0.50))
dev.off()

my_levels <- c("SC-beta", "SC-alpha", "DPP4⁺/ALDH1A1⁺",  "SC-EC", "SC-delta",  "SC-ductal",  "Proliferating cells")
samples_merged_integrated@meta.data$ann <- factor(x = samples_merged_integrated@meta.data$ann, levels = my_levels)
Idents(samples_merged_integrated) <- "annotation"

# Figure 6C: Heatmap showing scaled expression of top differentially expressed genes for each cluster compared to other clusters.
heatmap_genes <- c("INS", "IAPP", "ERO1B", "DLK1", "ACVR1C", "PCSK1", "EEF1A2", "PCDH7", "HADH", "NEFM", "GCG", 
                   "IGFBP2", "TTR", "SERPINA1", "SPINK1", "AGT", "ARX", "SERPINI1", "GLS", "CLU", "DPP4", "GC", 
                   "ALB", "ALDH1A1", "PDK3", "KCTD12", "NR2F1", "SSTR2", "PCDH17", "TPH1", "COL5A2", "DNAJC12", 
                   "FEV", "DDC", "SYT13", "CBLN1", "AC068587.4", "CHGA", "CAMK1D", "SST", "GHRL", "ACSL1", 
                   "HHEX", "SFRP1", "VTN", "CRH", "APOA1", "ST3GAL1", "OLFML3", "COL3A1", "KRT19", "HIST1H4C", 
                   "CENPF", "TOP2A", "TUBA1B", "RRM2", "MKI67", "H2AFZ", "HMGB2", "PTTG1", "TUBB")

tiff("Heatmap_wholedata.tiff", height = 6, width = 5, units = "in", res = 600)

DoHeatmap(samples_merged_integrated,
          features = heatmap_genes, 
          group.by = "ident",
          group.bar = TRUE, cells = 1:30000, 
          slot = "scale.data",
          assay = "SCT",
          label = FALSE,
          raster = TRUE,
          draw.lines = TRUE,
          group.bar.height = 0.02, 
          angle = 0, 
          group.colors = c("red1", "dodgerblue4", "dodgerblue1", "limegreen", "blue3", "darkorange3",  "darkorchid1")) 
+ theme(axis.text = element_text(size = 12, face = "plain"), legend.title = element_text(size = 18), legend.text = element_text(size = 16))
dev.off()

# Figures 6F: Stacked Violin plots of core markers (log1P SCT normalized counts)
P1 <- Stacked_VlnPlot(seurat_object = samples_merged_integrated, idents = c("SC-beta", "SC-alpha"), features = c("NKX6-1","ARX", "PDX1", "GC"), x_lab_rotate = TRUE, colors_use = c("red1", "blue3"),
                      plot_spacing = 0.3) & theme(axis.text = element_text(size = 20, face = "plain"))
P2 <- Stacked_VlnPlot(seurat_object = samples_merged_integrated, idents = c("SC-beta", "SC-EC"), features = c("INS", "TPH1", "LMX1A", "SLC18A1"), x_lab_rotate = TRUE, colors_use = c("red1", "limegreen"),
                      plot_spacing = 0.3) & theme(axis.text = element_text(size = 20, face = "plain"))
P3 <- Stacked_VlnPlot(seurat_object = samples_merged_integrated, idents = c("SC-beta", "SC-delta"), features = c("INS", "SST","GHRL", "HHEX"), x_lab_rotate = TRUE, colors_use = c("red1", "blue3"),
                      plot_spacing = 0.3) & theme(axis.text = element_text(size = 20, face = "plain"))
P4 <- Stacked_VlnPlot(seurat_object = samples_merged_integrated, idents = c("SC-beta", "Proliferating cells"), features = c("INS", "MKI67", "CENPF", "TOP2A"), x_lab_rotate = TRUE, colors_use = c("red1", "darkorchid1"), 
                      plot_spacing = 0.3) & theme(axis.text = element_text(size = 20, face = "plain"))
P5 <- Stacked_VlnPlot(seurat_object = samples_merged_integrated, idents = c("SC-alpha", "SC-EC"), features = c("GCG", "TPH1", "LMX1A", "SLC18A1"), x_lab_rotate = TRUE, colors_use = c("red1", "limegreen"),
                      plot_spacing = 0.3) & theme(axis.text = element_text(size = 20, face = "plain"))
P6 <- Stacked_VlnPlot(seurat_object = samples_merged_integrated, idents = c("SC-alpha", "SC-delta"), features = c("GCG", "SST","GHRL", "HHEX"), x_lab_rotate = TRUE, colors_use = c("red1", "blue3"),
                      plot_spacing = 0.3) & theme(axis.text = element_text(size = 20, face = "plain"))
P7 <- Stacked_VlnPlot(seurat_object = samples_merged_integrated, idents = c("SC-alpha", "Proliferating cells"), features = c("GCG", "MKI67", "CENPF", "TOP2A"), x_lab_rotate = TRUE, colors_use = c("red1", "darkorchid1"),
                      plot_spacing = 0.3) & theme(axis.text = element_text(size = 20, face = "plain"))
P8 <- Stacked_VlnPlot(seurat_object = samples_merged_integrated, idents = c("SC-alpha", "DPP4⁺/ALDH1A1⁺"), features = c("GCG", "GC","DPP4", "ALDH1A1"), x_lab_rotate = TRUE, colors_use = c("red1", "blue3"),
                      plot_spacing = 0.3) & theme(axis.text = element_text(size = 20, face = "plain"), strip.text = element_text(size = 24))

tiff("stacked_violin_plots.tiff", width = 18, height = 8, units = "in", res = 1200)
P1 | P2 | P3 | P4 | P5 | P6 | P7 | P8
dev.off()

# Loading necessary libraries
library(Seurat)
library(ggplot2)
library(gghalves)
library(scCustomize)
library(viridis)

# Load the processed Seurat object and ensure consistent cluster annotation ordering across all figures

my_levels <- c("SC-α", "SC-β", "SC-δ", "SC-EC", "SC-ductal",  "DPP4⁺/ALDH1A1⁺", "Proliferating cells")
samples_merged_integrated@meta.data$annotation <- factor(x = samples_merged_integrated@meta.data$annotation, levels = my_levels)

# Figure 7A: Generate UMAP plot with cells color-coded by cluster to show the seven distinct cell types
Idents(samples_merged_integrated) <- "annotation"
tiff("wholedata_clusters0.2res.tiff", width = 15, height = 10, units = "in", res = 600)
DimPlot(samples_merged_integrated, reduction = "umap", label = F, pt.size = 0.4, cols = c("dodgerblue4", "red1", "blue3","limegreen", "darkorange3", "dodgerblue1", "darkorchid1"))
dev.off()

# Figure 7B: Feature plots displaying the expression levels of key marker genes (log1P SCT normalized counts)

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

# Figure 7C: Heatmap showing scaled expression of top differentially expressed genes for each cluster compared to others

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

# Figures 7F-G: Violin plots of core beta and alpha cell markers (log1P SCT normalized counts)

tiff("beta_markers_violinplot.tiff", width = 8, height = 6, units = "in", res = 600)
Stacked_VlnPlot(seurat_object = samples_merged_integrated, features = c("INS", "IAPP","DLK1","HADH","PCSK1","NKX6-1","PDX1"), x_lab_rotate = TRUE, colors_use = c("dodgerblue4", "red1", "blue3","limegreen", "darkorange3", "dodgerblue1", "darkorchid1"),
                plot_spacing = 0.3) + theme(axis.text=element_text(size=7, face = "bold"))
dev.off()

tiff("alpha_markers_violinplot.tiff", width = 8, height = 6, units = "in", res = 600)
Stacked_VlnPlot(seurat_object = samples_merged_integrated,  features = c("GCG", "ARX","ALDH1A1","GC","TTR","DPP4","SERPINE2"), x_lab_rotate = TRUE, colors_use = c("dodgerblue4", "red1", "blue3","limegreen", "darkorange3", "dodgerblue1", "darkorchid1"),
                plot_spacing = 0.3) + theme(axis.text=element_text(size=7, face = "bold"))
dev.off()

#Supplementary figures
#Figure S12 [A-F], the below codes were used to generate plots for both the unfiltered and filtered data, the seurat objects were switch accordingly.
temp_labels <- samples_merged_integrated@meta.data %>%
  group_by(SAMPLE_ID) %>%
  tally()

#UMI counts
tiff("UMI_count_unfiltered.tiff", width = 12, height = 8, units = "in", res = 600)
ggplot() +
  geom_half_violin(
    data = samples_merged_integrated@meta.data, aes(SAMPLE_ID, nCount_RNA, fill = SAMPLE_ID),
    side = 'l', show.legend = FALSE, trim = FALSE
  ) +
  geom_half_boxplot(
    data = samples_merged_integrated@meta.data, aes(SAMPLE_ID, nCount_RNA, fill = SAMPLE_ID),
    side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
  ) +
  geom_text(
    data = temp_labels,
    aes(x = SAMPLE_ID, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 3.5
  ) +
  scale_color_manual(
    values = c("#D980FA", "#ED4C67", "#FFC312", "#B53471", "#F79F1F", "#C4E538", 
               "#009432", "#1289A7", "#FDA7DF", "#EE5A24", "#A3CB38", "#12CBC4")
  ) +
  scale_fill_manual(
    values = c("#D980FA", "#ED4C67", "#FFC312", "#B53471", "#F79F1F", "#C4E538", 
               "#009432", "#1289A7", "#FDA7DF", "#EE5A24", "#A3CB38", "#12CBC4")
  ) +
  scale_y_continuous(name = 'Number of transcripts', labels = scales::comma, expand = c(0.08, 0)) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold"), 
    axis.text.y = element_text(size = 10, face = "bold"), 
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold")
  )

dev.off()

#Gene counts
tiff("gene_count_samples_RAW.tiff", width = 12, height = 8, units = "in", res = 600)

ggplot() +
  geom_half_violin(
    data = samples_merged_integrated@meta.data, aes(SAMPLE_ID, nFeature_RNA, fill = SAMPLE_ID),
    side = 'l', show.legend = FALSE, trim = FALSE
  ) +
  geom_half_boxplot(
    data = samples_merged_integrated@meta.data, aes(SAMPLE_ID, nFeature_RNA, fill = SAMPLE_ID),
    side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
  ) +
  geom_text(
    data = temp_labels,
    aes(x = SAMPLE_ID, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 3.5
  ) +
  scale_color_manual(
    values = c("#D980FA", "#ED4C67", "#FFC312", "#B53471", "#F79F1F", "#C4E538", 
               "#009432", "#1289A7", "#FDA7DF", "#EE5A24", "#A3CB38", "#12CBC4")
  ) +
  scale_fill_manual(
    values = c("#D980FA", "#ED4C67", "#FFC312", "#B53471", "#F79F1F", "#C4E538", 
               "#009432", "#1289A7", "#FDA7DF", "#EE5A24", "#A3CB38", "#12CBC4")
  ) +
  scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma, expand = c(0.08, 0)) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold"), 
    axis.text.y = element_text(size = 10, face = "bold"), 
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold")
  )

dev.off()

#Mitochondrial percentage
tiff("MT_count_samples_RAW.tiff", width = 12, height = 8, units = "in", res = 600)

ggplot() +
  geom_half_violin(
    data = samples_merged_integrated@meta.data, aes(SAMPLE_ID, mitoPercent, fill = SAMPLE_ID),
    side = 'l', show.legend = FALSE, trim = FALSE
  ) +
  geom_half_boxplot(
    data = samples_merged_integrated@meta.data, aes(SAMPLE_ID, mitoPercent, fill = SAMPLE_ID),
    side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
  ) +
  geom_text(
    data = temp_labels,
    aes(x = SAMPLE_ID, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 3.5
  ) +
  scale_color_manual(
    values = c("#D980FA", "#ED4C67", "#FFC312", "#B53471", "#F79F1F", "#C4E538", 
               "#009432", "#1289A7", "#FDA7DF", "#EE5A24", "#A3CB38", "#12CBC4")
  ) +
  scale_fill_manual(
    values = c("#D980FA", "#ED4C67", "#FFC312", "#B53471", "#F79F1F", "#C4E538", 
               "#009432", "#1289A7", "#FDA7DF", "#EE5A24", "#A3CB38", "#12CBC4")
  ) +
  scale_y_continuous(name = 'Mitochondrial percentage', labels = scales::comma, expand = c(0.08, 0)) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold"), 
    axis.text.y = element_text(size = 10, face = "bold"), 
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold")
  )

dev.off()

#Figure S13B
Idents(samples_merged_integrated) <- "annotation"
alpha_beta <- subset(samples_merged_integrated, ident = c("SC-α", "SC-β"))
tiff("alpha_beta_cluster.tiff", width = 12, height = 8, units = "in", res = 600)
DimPlot(alpha_beta, reduction = "umap", label = F, pt.size = 0.4)
dev.off()

#Figure S13C
#Cell IDs of polyhormonal cells were extracted from the seurta objects and added to the metadata column. 
polyhormonal_cell_IDs <- readLines("polyhormone_cell_info.txt")

if (length(polyhormonal_cell_IDs) == ncol(alpha_beta)) {
  alpha_beta@meta.data$CellType <- polyhormonal_cell_IDs
  print("CellType column added successfully to metadata.")
} else {
  stop("The number of cell types in the file does not match the number of cells in the Seurat object.")
}

Idents(alpha_beta) <- "CellType"
tiff("polyhormonal.tiff", width = 12, height = 8, units = "in", res = 600)
DimPlot(alpha_beta, reduction = "umap", label = F, pt.size = 0.4, cols = c("gray", "red1"))
dev.off()

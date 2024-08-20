# Subsetting SC-α (alpha cells) and SC-β (beta cells) clusters and removing polyhormonal cells

Idents(samples_merged_integrated) <- "integrated_snn_res.0.2"
beta <- subset(samples_merged_integrated, GCG < 4.5, slot = "data", ident = "SC-β")
beta_polyhormone <- subset(samples_merged_integrated, GCG > 4.5, slot = "data", ident = "SC-β")

alpha <- subset(samples_merged_integrated, INS < 5.4, slot = "data", ident = "SC-α")
alpha_polyhormone <- subset(samples_merged_integrated, INS > 5.4, slot = "data", ident = "SC-α")

# Differential expression analysis on the SCT assay after removing polyhormonal cells
# Using the negative binomial test with "Sample" as a latent variable

# Differential expression analysis by zygosity within the beta cluster
Idents(beta) <- "cell.zygosity"

beta_heterozygous_normal <- FindMarkers(beta, assay="SCT", slot = "counts", logfc.threshold= -Inf, min.pct = 0.05, ident.1 = "Heterozygous_1", ident.2 = "Normal_1", test.use = "negbinom", recorrect_umi = FALSE, latent.vars = c("Sample"))
beta_homozygous_normal <- FindMarkers(beta, assay="SCT", slot = "counts", logfc.threshold= -Inf, min.pct = 0.05, ident.1 = "Homozygous_1", ident.2 = "Normal_1", test.use = "negbinom", recorrect_umi = FALSE, latent.vars = c("Sample"))
beta_homozygous_heterozygous <- FindMarkers(beta, assay="SCT", slot = "counts", logfc.threshold= -Inf, min.pct = 0.05, ident.1 = "Homozygous_1", ident.2 = "Heterozygous_1", test.use = "negbinom", recorrect_umi = FALSE, latent.vars = c("Sample"))

beta_heterozygous_normal <- FindMarkers(beta, assay="SCT", slot = "counts", logfc.threshold= -Inf, min.pct = 0.05, ident.1 = "Heterozygous_1", ident.2 = "Normal_1", recorrect_umi = FALSE, test.use = "negbinom")
beta_homozygous_normal <- FindMarkers(beta, assay="SCT", slot = "counts", logfc.threshold= -Inf, min.pct = 0.05, ident.1 = "Homozygous_1", ident.2 = "Normal_1", recorrect_umi = FALSE, test.use = "negbinom")
beta_homozygous_heterozygous <- FindMarkers(beta, assay="SCT", slot = "counts", logfc.threshold= -Inf, min.pct = 0.05, ident.1 = "Homozygous_1", ident.2 = "Heterozygous_1", recorrect_umi = FALSE, test.use = "negbinom")

# Differential expression analysis by zygosity within the alpha cluster
Idents(alpha) <- "cell.zygosity"

alpha_heterozygous_normal <- FindMarkers(alpha, assay="SCT", slot = "counts", logfc.threshold= -Inf, min.pct = 0.05, ident.1 = "Heterozygous_0", ident.2 = "Normal_0", test.use = "negbinom", recorrect_umi = FALSE, latent.vars = c("Sample"))
alpha_homozygous_normal <- FindMarkers(alpha, assay="SCT", slot = "counts", logfc.threshold= -Inf, min.pct = 0.05, ident.1 = "Homozygous_0", ident.2 = "Normal_0", test.use = "negbinom", recorrect_umi = FALSE,latent.vars = c("Sample"))
alpha_homozygous_heterozygous <- FindMarkers(alpha, assay="SCT", slot = "counts", logfc.threshold= -Inf, min.pct = 0.05, ident.1 = "Homozygous_0", ident.2 = "Heterozygous_0", test.use = "negbinom", recorrect_umi = FALSE, latent.vars = c("Sample"))

alpha_heterozygous_normal <- FindMarkers(alpha, assay="SCT", slot = "counts", logfc.threshold= -Inf, min.pct = 0.05, ident.1 = "Heterozygous_0", ident.2 = "Normal_0", recorrect_umi = FALSE, test.use = "negbinom")
alpha_homozygous_normal <- FindMarkers(alpha, assay="SCT", slot = "counts", logfc.threshold= -Inf, min.pct = 0.05, ident.1 = "Homozygous_0", ident.2 = "Normal_0", recorrect_umi = FALSE, test.use = "negbinom")
alpha_homozygous_heterozygous <- FindMarkers(alpha, assay="SCT", slot = "counts", logfc.threshold= -Inf, min.pct = 0.05, ident.1 = "Homozygous_0", ident.2 = "Heterozygous_0", recorrect_umi = FALSE, test.use = "negbinom")


# Differential expression analysis by differentiation state within the beta cluster

Idents(beta) <- "cell.sample"
beta_heterozygous_normal_differentiation_1 <- FindMarkers(beta, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "RH1_1", ident.2 = "RR2_1", recorrect_umi = FALSE, test.use = "negbinom")
beta_homozygous_normal_differentiation_1 <- FindMarkers(beta, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "HH1_1", ident.2 = "RR2_1", recorrect_umi = FALSE, test.use = "negbinom")
beta_homozygous_heterozygous_differentiation_1 <- FindMarkers(beta, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "HH1_1", ident.2 = "RH1_1", recorrect_umi = FALSE, test.use = "negbinom")

beta_heterozygous_normal_differentiation_2 <- FindMarkers(beta, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "RH2_1", ident.2 = "RR1_1", recorrect_umi = FALSE, test.use = "negbinom")
beta_homozygous_normal_differentiation_2 <- FindMarkers(beta, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "HH2_1", ident.2 = "RR1_1", recorrect_umi = FALSE, test.use = "negbinom")
beta_homozygous_heterozygous_differentiation_2 <- FindMarkers(beta, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "HH2_1", ident.2 = "RH2_1", recorrect_umi = FALSE, test.use = "negbinom")

beta_heterozygous_normal_differentiation_3 <- FindMarkers(beta, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "RH3_1", ident.2 = "RR3_1", recorrect_umi = FALSE, test.use = "negbinom")
beta_homozygous_normal_differentiation_3 <- FindMarkers(beta, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "HH3_1", ident.2 = "RR3_1", recorrect_umi = FALSE, test.use = "negbinom")
beta_homozygous_heterozygous_differentiation_3 <- FindMarkers(beta, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "HH3_1", ident.2 = "RH3_1", recorrect_umi = FALSE, test.use = "negbinom")

beta_heterozygous_normal_differentiation_4 <- FindMarkers(beta, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "RH2RD3_1", ident.2 = "RR2RD3_1", recorrect_umi = FALSE, test.use = "negbinom")
beta_homozygous_normal_differentiation_4 <- FindMarkers(beta, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "HH2RD3_1", ident.2 = "RR2RD3_1", recorrect_umi = FALSE, test.use = "negbinom")
beta_homozygous_heterozygous_differentiation_4 <- FindMarkers(beta, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "HH2RD3_1", ident.2 = "RH2RD3_1", recorrect_umi = FALSE, test.use = "negbinom")

# Differential expression analysis by differentiation state within the alpha cluster
Idents(alpha) <- "cell.sample"

alpha_heterozygous_normal_differentiation_1  <- FindMarkers(alpha, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "RH1_0", ident.2 = "RR2_0", recorrect_umi = FALSE, test.use = "negbinom")
alpha_homozygous_normal_differentiation_1 <- FindMarkers(alpha, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "HH1_0", ident.2 = "RR2_0", recorrect_umi = FALSE, test.use = "negbinom")
alpha_homozygous_heterozygous_differentiation_1 <- FindMarkers(alpha, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "HH1_0", ident.2 = "RH1_0", recorrect_umi = FALSE, test.use = "negbinom")

alpha_heterozygous_normal_differentiation_2 <- FindMarkers(alpha, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "RH2_0", ident.2 = "RR1_0", recorrect_umi = FALSE, test.use = "negbinom")
alpha_homozygous_normal_differentiation_2 <- FindMarkers(alpha, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "HH2_0", ident.2 = "RR1_0", recorrect_umi = FALSE, test.use = "negbinom")
alpha_homozygous_heterozygous_differentiation_2 <- FindMarkers(alpha, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "HH2_0", ident.2 = "RH2_0", recorrect_umi = FALSE, test.use = "negbinom")

alpha_heterozygous_normal_differentiation_3 <- FindMarkers(alpha, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "RH3_0", ident.2 = "RR3_0", recorrect_umi = FALSE, test.use = "negbinom")
alpha_homozygous_normal_differentiation_3 <- FindMarkers(alpha, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "HH3_0", ident.2 = "RR3_0", recorrect_umi = FALSE, test.use = "negbinom")
alpha_homozygous_heterozygous_differentiation_3 <- FindMarkers(alpha, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "HH3_0", ident.2 = "RH3_0", recorrect_umi = FALSE, test.use = "negbinom")

alpha_heterozygous_normal_differentiation_4 <- FindMarkers(alpha, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "RH2RD3_0", ident.2 = "RR2RD3_0", recorrect_umi = FALSE, test.use = "negbinom")
alpha_homozygous_normal_differentiation_4 <- FindMarkers(alpha, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "HH2RD3_0", ident.2 = "RR2RD3_0", recorrect_umi = FALSE, test.use = "negbinom")
alpha_homozygous_heterozygous_differentiation_4 <- FindMarkers(alpha, assay="SCT", slot = "counts", min.pct = 0.05, logfc.threshold= -Inf, ident.1 = "HH2RD3_0", ident.2 = "RH2RD3_0", recorrect_umi = FALSE, test.use = "negbinom")

# second vesrion of integration

# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)

library(glmGamPoi)
library(patchwork)
library(clustree)
library(gridExtra)
library(tidyr)
library(dplyr)
library(MAST)

# upload filtered seurat object
Rat1 = readRDS('./Rat1_doubletfinder_singlets.Rds')
DefaultAssay(Rat1) = 'RNA'
Rat1[['SCT']] = NULL
Rat1@meta.data$nCount_SCT = NULL
Rat1@meta.data$nFeature_SCT = NULL
Rat1@meta.data$SCT_snn_res.0.8 = NULL
Rat1@meta.data$seurat_clusters = NULL
Rat1@meta.data$pANN_0.25_0.15_1488 = NULL
Rat1@meta.data$DF.classifications_0.25_0.15_1488 = NULL
head(Rat1@meta.data)

Rat2 = readRDS('./Rat2_doubletfinder_singlets.Rds')
DefaultAssay(Rat2) = 'RNA'
Rat2[['SCT']] = NULL
Rat2@meta.data$nCount_SCT = NULL
Rat2@meta.data$nFeature_SCT = NULL
Rat2@meta.data$SCT_snn_res.0.8 = NULL
Rat2@meta.data$seurat_clusters = NULL
Rat2@meta.data$pANN_0.25_0.01_565 = NULL
Rat2@meta.data$DF.classifications_0.25_0.01_565 = NULL
head(Rat2@meta.data)

Rat3 = readRDS('./Rat3_doubletfinder_singlets.Rds')
DefaultAssay(Rat3) = 'RNA'
Rat3[['SCT']] = NULL
Rat3@meta.data$nCount_SCT = NULL
Rat3@meta.data$nFeature_SCT = NULL
Rat3@meta.data$SCT_snn_res.0.8 = NULL
Rat3@meta.data$seurat_clusters = NULL
Rat3@meta.data$pANN_0.25_0.1_700 = NULL
Rat3@meta.data$DF.classifications_0.25_0.1_700 = NULL
head(Rat3@meta.data)

Rat4 = readRDS('./Rat4_doubletfinder_singlets.Rds')
DefaultAssay(Rat4) = 'RNA'
Rat4[['SCT']] = NULL
Rat4@meta.data$nCount_SCT = NULL
Rat4@meta.data$nFeature_SCT = NULL
Rat4@meta.data$SCT_snn_res.0.8 = NULL
Rat4@meta.data$seurat_clusters = NULL
Rat4@meta.data$pANN_0.25_0.21_439 = NULL
Rat4@meta.data$DF.classifications_0.25_0.21_439 = NULL
head(Rat4@meta.data)

Rat5 = readRDS('./Rat5_doubletfinder_singlets.Rds')
DefaultAssay(Rat5) = 'RNA'
Rat5[['SCT']] = NULL
Rat5@meta.data$nCount_SCT = NULL
Rat5@meta.data$nFeature_SCT = NULL
Rat5@meta.data$SCT_snn_res.0.8 = NULL
Rat5@meta.data$seurat_clusters = NULL
Rat5@meta.data$pANN_0.25_0.03_172 = NULL
Rat5@meta.data$DF.classifications_0.25_0.03_172 = NULL
head(Rat5@meta.data)

Rat8 = readRDS('./Rat8_doubletfinder_singlets.Rds')
DefaultAssay(Rat8) = 'RNA'
Rat8[['SCT']] = NULL
Rat8@meta.data$nCount_SCT = NULL
Rat8@meta.data$nFeature_SCT = NULL
Rat8@meta.data$SCT_snn_res.0.8 = NULL
Rat8@meta.data$seurat_clusters = NULL
Rat8@meta.data$pANN_0.25_0.23_52 = NULL
Rat8@meta.data$DF.classifications_0.25_0.23_52 = NULL
head(Rat8@meta.data)

K2.1 = readRDS('./K2_1_doubletfinder_singlets.Rds')
DefaultAssay(K2.1) = 'RNA'
K2.1[['SCT']] = NULL
K2.1@meta.data$nCount_SCT = NULL
K2.1@meta.data$nFeature_SCT = NULL
K2.1@meta.data$SCT_snn_res.0.8 = NULL
K2.1@meta.data$seurat_clusters = NULL
K2.1@meta.data$pANN_0.25_0.01_166 = NULL
K2.1@meta.data$DF.classifications_0.25_0.01_166 = NULL
head(K2.1@meta.data)

K2.2 = readRDS('./K2_2_doubletfinder_singlets.Rds')
DefaultAssay(K2.2) = 'RNA'
K2.2[['SCT']] = NULL
K2.2@meta.data$nCount_SCT = NULL
K2.2@meta.data$nFeature_SCT = NULL
K2.2@meta.data$SCT_snn_res.0.8 = NULL
K2.2@meta.data$seurat_clusters = NULL
K2.2@meta.data$pANN_0.25_0.19_314 = NULL
K2.2@meta.data$DF.classifications_0.25_0.19_314 = NULL
head(K2.2@meta.data)

K2.3 = readRDS('./K2_3_doubletfinder_singlets.Rds')
DefaultAssay(K2.3) = 'RNA'
K2_3[['SCT']] = NULL
K2.3@meta.data$nCount_SCT = NULL
K2.3@meta.data$nFeature_SCT = NULL
K2.3@meta.data$SCT_snn_res.0.8 = NULL
K2.3@meta.data$seurat_clusters = NULL
K2.3@meta.data$pANN_0.25_0.005_755 = NULL
K2.3@meta.data$DF.classifications_0.25_0.005_755 = NULL
head(K2.3@meta.data)

K2.4 = readRDS('./K2_4_doubletfinder_singlets.Rds')
DefaultAssay(K2.4) = 'RNA'
K2.4[['SCT']] = NULL
K2.4@meta.data$nCount_SCT = NULL
K2.4@meta.data$nFeature_SCT = NULL
K2.4@meta.data$SCT_snn_res.0.8 = NULL
K2.4@meta.data$seurat_clusters = NULL
K2.4@meta.data$pANN_0.25_0.02_1528 = NULL
K2.4@meta.data$DF.classifications_0.25_0.02_1528 = NULL
head(K2.4@meta.data)

K2.5 = readRDS('./K2_5_doubletfinder_singlets.Rds')
DefaultAssay(K2.5) = 'RNA'
K2.5[['SCT']] = NULL
K2.5@meta.data$nCount_SCT = NULL
K2.5@meta.data$nFeature_SCT = NULL
K2.5@meta.data$SCT_snn_res.0.8 = NULL
K2.5@meta.data$seurat_clusters = NULL
K2.5@meta.data$pANN_0.25_0.3_676 = NULL
K2.5@meta.data$DF.classifications_0.25_0.3_676 = NULL
head(K2.5@meta.data)

K2.6 = readRDS('./K2_6_doubletfinder_singlets.Rds')
DefaultAssay(K2.6) = 'RNA'
K2.6[['SCT']] = NULL
K2.6@meta.data$nCount_SCT = NULL
K2.6@meta.data$nFeature_SCT = NULL
K2.6@meta.data$SCT_snn_res.0.8 = NULL
K2.6@meta.data$seurat_clusters = NULL
K2.6@meta.data$pANN_0.25_0.005_413 = NULL
K2.6@meta.data$DF.classifications_0.25_0.005_413 = NULL
head(K2.6@meta.data)

#merge seurat object
ls()
merged_seurat <- merge(K2.1, y = c(K2.2, K2.3, K2.4, K2.5, K2.6, Rat1,
                                   Rat2, Rat3, Rat4, Rat5, Rat8),
                       add.cell.ids = ls(),
                       project = 'agrsv_tamed')

View(merged_seurat@meta.data)
# create column barcode
merged_seurat$sample <- rownames(merged_seurat@meta.data)
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('smpl', 'Barcode'), 
                                    sep = '_')
View(merged_seurat@meta.data)

# QC and filtering --------------
# during detection of doublets I select <20% mt and <1% hb genes (see removing-duoblets-by-doubletfinder.R script)

VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.hb',"percent.rb" ),ncol=5)
FeatureScatter(merged_seurat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'glm')

# standart workflow analysis
merged_seurat <- SCTransform(merged_seurat, method = "glmGamPoi", vars.to.regress = "percent.mt")
merged_seurat <- RunPCA(object = merged_seurat)
ElbowPlot(merged_seurat,ndims = 50)
merged_seurat <- FindNeighbors(object = merged_seurat, dims = 1:40)
merged_seurat <- FindClusters(object = merged_seurat)
merged_seurat <- RunUMAP(object = merged_seurat, dims = 1:40)
# plot
p1 = DimPlot(merged_seurat, reduction = 'umap', group.by = 'smpl')
p2 = DimPlot(merged_seurat, reduction = 'umap', group.by = 'genotype')
grid.arrange(p1, p2, ncol= 2)

# Performing integration on datasets normalized with SCTransform
merged_srat.list <- SplitObject(merged_seurat, split.by = "smpl")
View(merged_srat.list)
merged_srat.list <- lapply(X = merged_srat.list, function(x) SCTransform(x, 
                                                                         vars.to.regress = c("percent.mt"),
                                                                         variable.features.n = 10000,
                                                                         ncells = min(100000, ncol(x)),
                                                                         method = "glmGamPoi",
                                                                         #return.only.var.genes = T, # F in Sankowsky's script
                                                                         verbose = T))
rm(list=ls(pattern="K2"))
rm(list=ls(pattern="Rat"))
gc()
features <- SelectIntegrationFeatures(object.list = merged_srat.list, nfeatures = 3000)
# Perform integration. 
merged_srat.list <- PrepSCTIntegration(object.list = merged_srat.list, anchor.features = features)

merged_srat.anchors <- FindIntegrationAnchors(object.list = merged_srat.list, 
                                              normalization.method = "SCT", 
                                              anchor.features = features, 
                                              dims = 1:45, # influence a lot on layout and cells belonging to exact clusters
                                              verbose = T) # CCA
merged_srat.integrated <- IntegrateData(anchorset = merged_srat.anchors, 
                                      normalization.method = "SCT", 
                                      dims = 1:45, # influence a lot on layout and cells belonging to exact clusters
                                      verbose = T)

# standart workflow analysis
merged_srat.integrated <- RunPCA(object = merged_srat.integrated)
ElbowPlot(merged_srat.integrated,ndims = 50)
DefaultAssay(merged_srat.integrated)

#merged_srat.integrated <- RunUMAP(object = merged_srat.integrated, dims = 1:40)
#merged_srat.integrated <- RunUMAP(object = merged_srat.integrated, dims = 1:40, 
#                                  min.dist = 0.1, n.neighbors = 100)
merged_srat.integrated <- RunUMAP(object = merged_srat.integrated, dims = 1:40, 
                                  min.dist = 0.4, n.neighbors = 100)
p3 = DimPlot(merged_srat.integrated, reduction = 'umap', group.by = 'smpl')
p4 = DimPlot(merged_srat.integrated, reduction = 'umap', group.by = 'genotype')
grid.arrange(p3, p4, ncol = 2)

# quantity batches
bp1 = FeaturePlot(merged_srat.integrated, 'nFeature_RNA', reduction="umap")
bp2 = FeaturePlot(merged_srat.integrated, 'nCount_RNA', reduction="umap")
bp3 = FeaturePlot(merged_srat.integrated, 'percent.mt', reduction="umap")
bp4 = FeaturePlot(merged_srat.integrated, 'percent.hb', reduction="umap")
bp5 = FeaturePlot(merged_srat.integrated, 'percent.rb', reduction="umap")

grid.arrange(p3, bp1,bp2,bp3,bp4, bp5,ncol= 3)


#table(merged_srat.integrated@meta.data$genotype)

#merged_srat.integrated = readRDS('./')
# clustering
merged_srat.integrated <- FindNeighbors(object = merged_srat.integrated, dims = 1:40)
merged_srat.integrated <- FindClusters(object = merged_srat.integrated, resolution = c(0.1, 0.3, 0.5))
View(merged_srat.integrated@meta.data)
p5 = DimPlot(merged_srat.integrated, group.by="integrated_snn_res.0.1", label=T)
p6 = DimPlot(merged_srat.integrated, group.by="integrated_snn_res.0.3", label=T)
p7 = DimPlot(merged_srat.integrated, group.by="integrated_snn_res.0.5", label=T)

grid.arrange(p5,p6,p7, ncol= 3)
clustree(merged_srat.integrated, prefix = "integrated_snn_res.")

saveRDS(object = merged_srat.integrated, file = "cca_integration_v2.Rds")
# MARKERS ANALYSIS
Idents(merged_srat.integrated) = "integrated_snn_res.0.3"
DimPlot(merged_srat.integrated, reduction = 'umap', label=T)
table(merged_srat.integrated@meta.data$integrated_snn_res.0.3)
# find ALL markers
allMarkers <- FindAllMarkers(merged_srat.integrated, max.cells.per.ident = 100, min.pct=0.25,
                             test.use = "MAST", only.pos = T, logfc.threshold = 0.25)
allMarkers$pct.diff <- allMarkers$pct.1 - allMarkers$pct.2
goodMarkers.pct.diff = allMarkers %>% group_by(cluster) %>% top_n(n = 5, wt=pct.diff) %>% pull(gene)
table(goodMarkers.pct.diff$cluster)
goodMarkers.logfc = allMarkers %>% group_by(cluster) %>% top_n(n = 1, wt=avg_log2FC) %>% pull(gene)
table(goodMarkers.logfc$cluster)

DefaultAssay(merged_srat.integrated) = 'RNA'
merged_srat.integrated = NormalizeData(merged_srat.integrated)
# CLUSTERS 0,4,6,15 -- ASTRO
astro_mkrs = c('Agt', 'Slc1a3', 'Gpc5', 'Sparcl1', 'Itih3', 'Gfap')
FeaturePlot(merged_srat.integrated, astro_mkrs, reduction="umap", ncol=3)

# CLUSTERS 1,9, 13,17,24 -- OPC/COP/OG
opc_cop_og_mkrs = c("Pdgfra", "Col5a3", "Tnr", "Gpr17", "Mag", "Gjc2")
FeaturePlot(merged_srat.integrated, opc_cop_og_mkrs, reduction="umap", ncol=3)

# CLUSTERS 2,26,16 -- MG/MG*/macrophages
mg_mgAct_mcrphg = c("Csf1r", "Mpeg1", "Ctss", "Cd83", "Ms4a6bl", "Cd68")
FeaturePlot(merged_srat.integrated, mg_mgAct_mcrphg, reduction="umap", ncol=3)

# CLUSTERS 3,19,21 -- Neurons
neuro_mkrs = c("Tmem130", "Snhg11", 'Kcnip4', "Ppp1r17", 'Col19a1', 'Slc17a6')
FeaturePlot(merged_srat.integrated, neuro_mkrs, reduction="umap", ncol=3)

# CLUSTER 14 -- NEURONS + ?
mkrs14 = FindConservedMarkers(merged_srat.integrated, ident.1=14, grouping.var = 'genotype')
mkrs14$gene = rownames(mkrs14)
mkrs14$agrsv_pct.diff = mkrs14$agrsv_pct.1 - mkrs14$agrsv_pct.2
mkrs14$tamed_pct.diff = mkrs14$tamed_pct.1 - mkrs14$tamed_pct.2
#pct.diff
top14_agrsv_pct.diff = mkrs14 %>% top_n(n = 5, wt=agrsv_pct.diff) %>% pull(gene)
top14_tamed_pct.diff = mkrs14 %>% top_n(n = 5, wt=tamed_pct.diff) %>% pull(gene)
#logfc
top14_agrsv_logfc = mkrs14 %>%top_n(n = 5, wt=agrsv_avg_log2FC) %>% pull(gene)
top14_tamed_logfc = mkrs14 %>% top_n(n = 5, wt=tamed_avg_log2FC) %>% pull(gene)

all14 = c("Syp", "Gaa", "Clstn3" )
FeaturePlot(merged_srat.integrated, all14, reduction="umap", ncol=5)
FeaturePlot(merged_srat.integrated, c(top14_agrsv_pct.diff[1:3], all14) , reduction="umap", ncol=3) # 14 neuro
FeaturePlot(merged_srat.integrated, c('Chgb', 'LOC120102252') , reduction="umap", ncol=2) # ?
FeaturePlot(merged_srat.integrated, c('mt-Rnr1',"mt-Rnr2") , reduction="umap", ncol=2, min.cutoff = 'q20') # ?


# subset neurons
Idents(merged_srat.integrated) = "integrated_snn_res.0.5"
DimPlot(merged_srat.integrated, reduction = 'umap', label=T)
table(merged_srat.integrated@meta.data$integrated_snn_res.0.5)

neurons_cca_integrated_v2 = subset(x = merged_srat.integrated, idents = c('4', '15', '22', '23', '26', '27'))
saveRDS(object = neurons_cca_integrated_v2, file = "neurons_cca_integrated_v2.Rds")  

Idents(merged_srat.integrated) = "integrated_snn_res.0.1"
table(merged_srat.integrated@meta.data$integrated_snn_res.0.1)


###################################################################################################################
########################################## INTEGRATED RES 0.5 #################################################
###################################################################################################################
merged_srat.integrated = readRDS('./cca_integration_v2.Rds')
Idents(merged_srat.integrated) = "integrated_snn_res.0.5"
DimPlot(merged_srat.integrated, reduction = 'umap', label=T)
table(merged_srat.integrated@meta.data$integrated_snn_res.0.5)

DefaultAssay(merged_srat.integrated) = 'integrated'
# find ALL markers
allMarkers <- FindAllMarkers(merged_srat.integrated, max.cells.per.ident = 100, min.pct=0.25,
                             test.use = "MAST", only.pos = T, logfc.threshold = 0.25)
allMarkers$pct.diff <- allMarkers$pct.1 - allMarkers$pct.2
goodMarkers.pct.diff = allMarkers %>% group_by(cluster) %>% top_n(n = 5, wt=pct.diff) %>% pull(gene)
table(goodMarkers.pct.diff$cluster)
goodMarkers.logfc = allMarkers %>% group_by(cluster) %>% top_n(n = 5, wt=avg_log2FC) %>% pull(gene)
table(goodMarkers.logfc$cluster)

DefaultAssay(merged_srat.integrated) = 'RNA'
merged_srat.integrated = NormalizeData(merged_srat.integrated)

# CLUSTERS 0,3,6,16,20,30 -- Astrocutes
astro_mkrs = c("Gpc5", "Slc6a11", "Slc1a3", "Itih3")
FeaturePlot(merged_srat.integrated, astro_mkrs, reduction="umap", ncol=2, min.cutoff = 'q10')
sum(11595,6244,6844,1265,996,225)

# CLUSTERS 1,29 / 21 -- MG/Macrophages
mg_mkrs = c("Selplg", "Ccr5", "Mpeg1", "Cx3cr1")
FeaturePlot(merged_srat.integrated, mg_mkrs, reduction="umap", ncol=2)#, min.cutoff = 'q10')
sum(9093,269)
mphag_mkrs = c("Mrc1", "Clec10a", "Ms4a4a", "Cybb")
FeaturePlot(merged_srat.integrated, mphag_mkrs, reduction="umap", ncol=2)#, min.cutoff = 'q10')

# CLUSTERS 2,9,13,14,19 -- OG/COP/OPC
og_mkrs = c("Pdgfra","Col5a3", "Tnr", "Gpr17", "Mag", "Gjc2")
FeaturePlot(merged_srat.integrated, og_mkrs, reduction="umap", ncol=3)#, min.cutoff = 'q10')
sum(7289,3975,3201,2882,1073)

# CLUSTERS 4,15,22,23,26,27 -- Neurons
neuro_mkrs = c("Snhg11", "Tmem130", "Kcnip4",'Slc17a6')
FeaturePlot(merged_srat.integrated, neuro_mkrs, reduction="umap", ncol=2)#, min.cutoff = 'q10')
sum(6733,1428,665,642,573,431)

# CLUSTERS 5,25  -- Tanicytes
tani_mkrs = c("Ppp1r1b", "Pacrg")
FeaturePlot(merged_srat.integrated, tani_mkrs, reduction="umap", ncol=1)#, min.cutoff = 'q10')
sum(6642, 584)

# CLUSTERS 12  -- Ependymal
epend_mkrs = c("Cdhr4", "Rp1", "Dnah3", "Dnah5")
FeaturePlot(merged_srat.integrated, epend_mkrs, reduction="umap", ncol=2)#, min.cutoff = 'q10')
sum(6642, 584)

# CLUSTERS 8,10,24  -- Pericytes
pery_mkrs = c("Flt1", "Slco1a4" )
FeaturePlot(merged_srat.integrated, pery_mkrs, reduction="umap", ncol=1)#, min.cutoff = 'q10')
sum(4553, 3902,639)

# CLUSTERS 17,31  -- leptomeningeal
lepto_mkrs = c("Col3a1", "Col1a1")
FeaturePlot(merged_srat.integrated, lepto_mkrs, reduction="umap", ncol=1)#, min.cutoff = 'q10')
sum(1149, 75)

# CLUSTERS 7,28  -- ?
mkrs7 = c("Avp", "Nxph4", "Ttn", "Parpbp", "LOC100360573",
               "Lsamp", "S100a1", "Apoc1", "Zeb2", "Ptprd")
FeaturePlot(merged_srat.integrated, mkrs7, reduction="umap", ncol=5)#, min.cutoff = 'q10')

mkrs28 = c("Unc5d", "Calcr", "Snhg11", "Pgr15l", "Arhgap36", 
           "Npy", "Cited1", "Agrp", "Glra2", "Nppa" )
FeaturePlot(merged_srat.integrated, mkrs26, reduction="umap", ncol=5)#, min.cutoff = 'q10')

mkrs7.c = FindConservedMarkers(merged_srat.integrated, ident.1=7, grouping.var = 'genotype')
mkrs7.c$gene = rownames(mkrs7.c)
mkrs7.c$agrsv_pct.diff = mkrs7.c$agrsv_pct.1 - mkrs7.c$agrsv_pct.2
mkrs7.c$tamed_pct.diff = mkrs7.c$tamed_pct.1 - mkrs7.c$tamed_pct.2
#pct.diff
top7_agrsv_pct.diff = mkrs7.c %>% top_n(n = 5, wt=agrsv_pct.diff) %>% pull(gene)
top7_tamed_pct.diff = mkrs7.c %>% top_n(n = 5, wt=tamed_pct.diff) %>% pull(gene)
#logfc
top7_agrsv_logfc = mkrs7.c %>%top_n(n = 5, wt=agrsv_avg_log2FC) %>% pull(gene)
top7_tamed_logfc = mkrs7.c %>% top_n(n = 5, wt=tamed_avg_log2FC) %>% pull(gene)

mkrs7_cnsvrd = c("Apod", "Mbp",  "Cga", "Pmch",
                 "Oxt","Tshb",   "Cartpt", "Mog",
                 "Tmsb4x", "Rps26", "Mt3", "Rpl27a" )
FeaturePlot(merged_srat.integrated, mkrs7_cnsvrd, reduction="umap", ncol=4)#, min.cutoff = 'q10')
FeaturePlot(merged_srat.integrated, c('Apod','Tmsb4x' ), reduction="umap", ncol=2)#, min.cutoff = 'q10')
FeaturePlot(merged_srat.integrated, c('Cga','Tshb' ), reduction="umap", ncol=2)#, min.cutoff = 'q10')

# CLUSTERS 11  -- 
mkrs11 = c("Inava", "Lpar1", "Krt17", "LOC100912227", "LOC120099769", 
          "Ssh2", "Gabra2", "Pcdh9", "Gdpd2", "Coq10b" )
FeaturePlot(merged_srat.integrated, mkrs11, reduction="umap", ncol=5)#, min.cutoff = 'q10')

# CLUSTERS 18  -- 
mkrs18 = c("Etl4", "Trpc5", "LOC102548639", "Gfra1", "Fabp4" ,
           "Peg3", "Mtss2", "LOC120093067", "Net1", "Lgi4" )
FeaturePlot(merged_srat.integrated, mkrs18, reduction="umap", ncol=5)#, min.cutoff = 'q10')

# rename clusters
new.cluster.ids <- c(
  'Astrocytes', 'Microglia','OG/OPC', 'Astrocytes', 'Neurons', 'Tanicytes', 'Astrocytes', 'Unknown_1', 'Pericytes', 'OG/OPC',
  'Pericytes', 'Unknown_2', 'Ependymal', 'OG/OPC', 'OG/OPC', 'Neurons', 'Astrocytes', 'Leptomeningeal', 'Unknown_3', 'OG/OPC',
  'Astrocytes', 'Macrophages', 'Neurons', 'Neurons', 'Pericytes', 'Tanicytes', 'Neurons', 'Neurons', 'Unknown_1', 'Microglia',
  'Astrocytes', 'Leptomeningeal')
names(new.cluster.ids) <- levels(merged_srat.integrated)
merged_srat.integrated <- RenameIdents(merged_srat.integrated, new.cluster.ids)  
DimPlot(merged_srat.integrated, reduction = "umap", label = TRUE) + NoLegend()

# Dot plot

general.mkrs = c("Slc6a11", "Itih3", "Selplg", "Cx3cr1", "Mrc1", "Ms4a4a","Pdgfra","Col5a3", "Tnr", "Gpr17", "Mag", "Gjc2",
                 "Snhg11", "Tmem130","Ppp1r1b", "Pacrg", "Dnah3", "Dnah5", "Flt1", "Slco1a4", "Col3a1", "Col1a1")
DotPlot(merged_srat.integrated, features = general.mkrs) + RotatedAxis()

# save object with umap parameters 0.4/100, resolution 0.5
saveRDS(merged_srat.integrated, 'cca_integrated_v2_final.Rds')

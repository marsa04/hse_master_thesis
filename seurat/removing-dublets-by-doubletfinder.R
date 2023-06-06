library(remotes)
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
#install.packages('Seurat')
#remotes::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(Seurat)
library(ggplot2)
#install.packages('tidyverse')
library(tidyverse)
library(data.table)
library(glmGamPoi)
library(clustree)
library(gridExtra)
################################################################################
setwd('~/novaseq-rats-analysis/singlets-data-after-doubletfinder/')
#DoubletFinder is more sensetive to heterotypic doublets and less sensetive to homotypic doublets




                ########################################################################################
                ################################# Rat1 #################################################
                ########################################################################################


# uplaod raw data
sample = 'Rat1'
counts = Read10X_h5('~/mstool/cellranger-7.0.0/Rat1_cellout/outs/filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)

# create seurat object
srat <- CreateSeuratObject(counts=counts, min.cells = 2, min.features = 200, project = sample)
head(srat@meta.data)
dim(srat)


#upload Rat's mito genes
mitogenes = fread(file = "~/novaseq-rats-analysis/mRatBN7.2_mitoGenes.txt", sep='\t', header = F)
mitogenes = mitogenes$V1[mitogenes$V1 %in% rownames(srat)]
srat[['percent.mt']] = PercentageFeatureSet(srat, features = mitogenes)
# add blood markers
bloodgenes = c("Hba-a1", "Hba-a2", "Hba-a3","Hbb", "Hbb-b1", "Hbb-bs","Hbe1", "Hbq1", "Hbq1b")
bloodgenes = bloodgenes[bloodgenes %in% rownames(srat)]
srat[['percent.hb']] <- PercentageFeatureSet(srat, features = bloodgenes)
sum(srat$percent.hb > 2)
hist(srat$percent.hb, breaks = 100000, xlim = range(1,10))
# add ribogenes
ribogenes = grep(pattern = "^Rp[sl]", x=rownames(x=srat), value = TRUE, ignore.case = T)
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^Rp[sl]")
# plots
VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.hb',"percent.rb" ), ncol=5)
FeatureScatter(srat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'glm')

# QC: filtering
srat <- subset(srat, subset = percent.mt < 20 & percent.hb < 1)
dim(srat)
#removing Hb- and ribogens
#counts <- counts[-(which(rownames(counts) %in% c(bloodgenes, ribogenes))),]
#srat <- subset(srat, features = rownames(counts))
#dim(srat)

counts=0
#SCT
srat <- SCTransform(srat, method = "glmGamPoi", vars.to.regress = "percent.mt")
VariableFeaturePlot(srat) + scale_y_log10()
#PCA
srat <- RunPCA(object = srat)
ElbowPlot(srat,ndims = 50)

srat <- FindNeighbors(object = srat, dims = 1:30)
srat <- FindClusters(object = srat)
#clustree(srat, prefix = "SCT_snn_res.")
#Idents(srat) = "SCT_snn_res.0.05"
srat <- RunUMAP(object = srat, dims = 1:30)

## pK Identification (no ground-truth) 
sweep.res.list_srat <- paramSweep_v3(srat, PCs = 1:30, sct = TRUE)
sweep.stats_srat <- summarizeSweep(sweep.res.list_srat, GT = FALSE)
bcmvn_srat <- find.pK(sweep.stats_srat)

ggplot(bcmvn_srat, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_srat %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- srat@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.1*nrow(srat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
srat <- doubletFinder_v3(srat, 
                          PCs = 1:30, 
                          pN = 0.25, 
                          pK = pK, 
                          nExp = nExp_poi.adj,
                          reuse.pANN = FALSE, sct = TRUE)
colnames(srat@meta.data)
#srat <- doubletFinder_v3(srat, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = "DF.classifications_0.25_0.27_1607", sct = TRUE)

DimPlot(srat, reduction = 'umap', group.by = "DF.classifications_0.25_0.15_1488")
DimPlot(srat, reduction='umap', label = T)
# number of singlets and doublets
table(srat@meta.data$DF.classifications_0.25_0.15_1488)
dim(srat)
# save BC cells-doublets
write.table(rownames(srat@meta.data[srat$DF.classifications_0.25_0.15_1488 == 'Doublet',]), 
            file = str_c(sample, '_doublets.txt'), sep = "\t",
            row.names = TRUE, col.names = NA)

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)
table(srat$SCT_snn_res.0.8)
rp1 = DimPlot(srat, reduction = 'umap', label=T)
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
grid.arrange(rp1,rp2,rp3,rp4,rp5,rp6, ncol= 3)

srat = subset(srat, subset = DF.classifications_0.25_0.15_1488 == 'Singlet')
dim(srat)
srat$genotype = 'tamed'
srat$sex = 'mixed'
srat@meta.data = srat@meta.data[,c(1,2,3,13,14,4,5,6,7,8,9,10,11,12)]
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))
#size = file.size("Rat1_doubletfinder_singlets.Rds")
#print(structure(size, class = "object_size"), units = "Mb")
# cells = names(Idents(srat)[Idents(srat) !=7]) -- убрать кластер 7
# srat = subset(srat, cells = cells)

# find ALL markers
Idents(srat)
DefaultAssay(srat)
allMarkers <- FindAllMarkers(srat, max.cells.per.ident = 100, min.pct=0.25,
                                test.use = "MAST", only.pos = T, logfc.threshold = 0.25)
sum(allMarkers$p_val_adj < 0.05)
allMarkers$pct.diff <- allMarkers$pct.1 - allMarkers$pct.2
goodMarkers.pct.diff = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=pct.diff) %>% pull(gene)
goodMarkers.logfc = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=avg_log2FC) %>% pull(gene)

DefaultAssay(ast.combined) = 'RNA'
ast.combined = NormalizeData(ast.combined)
FeaturePlot(srat, goodMarkers.pct.diff, reduction="umap", ncol=6)
FeaturePlot(srat, goodMarkers.logfc, reduction="umap", ncol=6)

table(srat@meta.data$SCT_snn_res.0.05)
cell_markers <- c(
  # neurons:
  "Tmem130", "Snhg11", "Kcnip4", "Slc17a6", "RGD1566401",
  # oligodendrocytes:
  "Mag", "Vcan",  "Tnr", "Bmp4", 
  # astrocytes:
  "Gpc5", "Agt",
  # pericytes:
  "Flt1", "Rgs5", "Slco1c1", "Ebf1",
  # ependymocytes:
  "Ccdc153",
  # tanycytes:
  "Crym",
  # microglia:
  "Cd83", "Ctss", "Lyn",
  # activated MG or MFs:
  "Mrc1", "Tpt1",
  # leptomeningeal?:
  "Slc47a1", "Igf2"
)


FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)

rp1 = DimPlot(srat, reduction = 'umap', label=T, group.by = 'SCT_snn_res.0.05')
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
grid.arrange(rp1,rp2,rp3,rp4,rp5,rp6, ncol= 3)






########################################################################################
################################# Rat2 #################################################
########################################################################################

# uplaod raw data
sample = 'Rat2'
counts = Read10X_h5('~/mstool/cellranger-7.0.0/Rat2_cellout/outs/filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)

# create seurat object
srat <- CreateSeuratObject(counts=counts, min.cells = 2, min.features = 200, project = sample)
head(srat@meta.data)
dim(srat)
srat$genotype = 'tamed'
srat$sex = 'mixed'

#upload Rat's mito genes
mitogenes = fread(file = "~/novaseq-rats-analysis/mRatBN7.2_mitoGenes.txt", sep='\t', header = F)
mitogenes = mitogenes$V1[mitogenes$V1 %in% rownames(srat)]
srat[['percent.mt']] = PercentageFeatureSet(srat, features = mitogenes)
# add blood markers
bloodgenes = c("Hba-a1", "Hba-a2", "Hba-a3","Hbb", "Hbb-b1", "Hbb-bs","Hbe1", "Hbq1", "Hbq1b")
bloodgenes = bloodgenes[bloodgenes %in% rownames(srat)]
srat[['percent.hb']] <- PercentageFeatureSet(srat, features = bloodgenes)
sum(srat$percent.hb > 1)
hist(srat$percent.hb, breaks = 100000, xlim = range(0,10))
# add ribogenes
ribogenes = grep(pattern = "^Rp[sl]", x=rownames(x=srat), value = TRUE, ignore.case = T)
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^Rp[sl]")

#
VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.hb',"percent.rb" ), ncol=5)
FeatureScatter(srat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'glm')

# QC: filtering
srat <- subset(srat, subset = percent.mt < 20 & percent.hb < 1)
dim(srat)
#removing Hb- and ribogens
counts <- counts[-(which(rownames(counts) %in% c(bloodgenes, ribogenes))),]
srat <- subset(srat, features = rownames(counts))
dim(srat)
counts=0
#SCT
srat <- SCTransform(srat, method = "glmGamPoi", vars.to.regress = "percent.mt")
VariableFeaturePlot(srat) + scale_y_log10()
#PCA
srat <- RunPCA(object = srat)
ElbowPlot(srat,ndims = 50)
srat <- FindNeighbors(object = srat, dims = 1:30)
srat <- FindClusters(object = srat)
#clustree(srat, prefix = "SCT_snn_res.")
#Idents(srat) = "SCT_snn_res.0.05"
srat <- RunUMAP(object = srat, dims = 1:30)
## pK Identification (no ground-truth) 
sweep.res.list_srat <- paramSweep_v3(srat, PCs = 1:30, sct = TRUE)
sweep.stats_srat <- summarizeSweep(sweep.res.list_srat, GT = FALSE)
bcmvn_srat <- find.pK(sweep.stats_srat)

ggplot(bcmvn_srat, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_srat %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- srat@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.069*nrow(srat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
srat <- doubletFinder_v3(srat, 
                         PCs = 1:30, 
                         pN = 0.25, 
                         pK = pK, 
                         nExp = nExp_poi.adj,
                         reuse.pANN = FALSE, sct = TRUE)
colnames(srat@meta.data)

DimPlot(srat, reduction = 'umap', group.by = "DF.classifications_0.25_0.01_565")
DimPlot(srat, reduction = 'umap', label=T)

# save BC cells-doublets
write.table(rownames(srat@meta.data[srat$DF.classifications_0.25_0.01_565 == 'Doublet',]), 
            file = str_c(sample, '_doublets.txt'), sep = "\t",
            row.names = TRUE, col.names = NA)

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)
table(srat$SCT_snn_res.0.8)
rp1 = DimPlot(srat, reduction = 'umap', label=T)
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
grid.arrange(rp1,rp2,rp3,rp4,rp5,rp6, ncol= 3)

srat = subset(srat, subset = DF.classifications_0.25_0.01_565 == 'Singlet')
dim(srat)
head(srat@meta.data)
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))
#  -------------------------------- find ALL markers --------------------------------
Idents(srat)
DefaultAssay(srat)
allMarkers <- FindAllMarkers(srat, max.cells.per.ident = 100, min.pct=0.25,
                             test.use = "MAST", only.pos = T, logfc.threshold = 0.25)
sum(allMarkers$p_val_adj < 0.05)
allMarkers$pct.diff <- allMarkers$pct.1 - allMarkers$pct.2
goodMarkers.pct.diff = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=pct.diff) %>% pull(gene)
table(goodMarkers.pct.diff$cluster)
goodMarkers.logfc = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=avg_log2FC) %>% pull(gene)


#DefaultAssay(ast.combined) = 'RNA'
#ast.combined = NormalizeData(ast.combined)
FeaturePlot(srat, goodMarkers.pct.diff[c(1:4,6:17,19:26)], reduction="umap", ncol=6)
FeaturePlot(srat, goodMarkers.logfc, reduction="umap", ncol=6)

table(srat@meta.data$SCT_snn_res.0.1)
cell_markers <- c(
  # neurons:
  "Tmem130", "Snhg11", "Kcnip4", "Slc17a6",
  # OPC: 
  "Pdgfra", "Fyn",  "Gpr17", "Sox10",
  # OG
  'Mog', 'Sirt2',
  # astrocytes: 
  "Gpc5", "Agt", 'Sox9',
  # pericytes: 
  "Flt1", "Rgs5", "Slco1c1", "Ebf1",
  # ependymocytes: 
  "Ccdc153", 
  # tanycytes:
  "Crym", 
  # microglia: 
  "Cd83", "Ctss", "Lyn",
  # activated MG or MFs: 
  "Mrc1", "Tpt1", 
  # leptomeningeal?: 
  "Slc47a1", "Igf2" 
  )

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)

rp1 = DimPlot(srat, reduction = 'umap', label=T, group.by = 'SCT_snn_res.0.05')
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
grid.arrange(rp1,rp2,rp3,rp4,rp5,rp6, ncol= 3)
# save singlets
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))
dim(srat)






########################################################################################
################################# Rat3 #################################################
########################################################################################

# uplaod raw data
sample = 'Rat3'
counts = Read10X_h5('~/mstool/cellranger-7.0.0/Rat3_cellout/outs/filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)

# create seurat object
srat <- CreateSeuratObject(counts=counts, min.cells = 3, min.features = 200, project = sample)
srat$genotype = 'agrsv'
srat$sex = 'female'
head(srat@meta.data)
dim(srat)
#upload Rat's mito genes
mitogenes = fread(file = "~/novaseq-rats-analysis/mRatBN7.2_mitoGenes.txt", sep='\t', header = F)
mitogenes = mitogenes$V1[mitogenes$V1 %in% rownames(srat)]
srat[['percent.mt']] = PercentageFeatureSet(srat, features = mitogenes)
sum(srat$percent.mt > 10)
# add blood markers
bloodgenes = c("Hba-a1", "Hba-a2", "Hba-a3","Hbb", "Hbb-b1", "Hbb-bs","Hbe1", "Hbq1", "Hbq1b")
bloodgenes = bloodgenes[bloodgenes %in% rownames(srat)]
srat[['percent.hb']] <- PercentageFeatureSet(srat, features = bloodgenes)
sum(srat$percent.hb > 1)
hist(srat$percent.hb, breaks = 100000, xlim = range(0,10))
# add ribogenes
ribogenes = grep(pattern = "^Rp[sl]", x=rownames(x=srat), value = TRUE, ignore.case = T)
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^Rp[sl]")

#
VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.hb',"percent.rb" ), ncol=5)
FeatureScatter(srat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'glm')

# QC: filtering
srat <- subset(srat, subset = percent.mt < 20 & percent.hb < 1)
dim(srat)
#removing Hb- and ribogens
counts <- counts[-(which(rownames(counts) %in% c(bloodgenes, ribogenes))),]
srat <- subset(srat, features = rownames(counts))
dim(srat)
counts=0
#SCT
srat <- SCTransform(srat, method = "glmGamPoi", vars.to.regress = "percent.mt")
VariableFeaturePlot(srat) + scale_y_log10()
#PCA
srat <- RunPCA(object = srat)
ElbowPlot(srat,ndims = 50)
srat <- FindNeighbors(object = srat, dims = 1:30)
srat <- FindClusters(object = srat)
clustree(srat, prefix = "SCT_snn_res.")
Idents(srat) = "SCT_snn_res.0.1"
srat <- RunUMAP(object = srat, dims = 1:30)
## pK Identification (no ground-truth) 
sweep.res.list_srat <- paramSweep_v3(srat, PCs = 1:30, sct = TRUE)
sweep.stats_srat <- summarizeSweep(sweep.res.list_srat, GT = FALSE)
bcmvn_srat <- find.pK(sweep.stats_srat)

ggplot(bcmvn_srat, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_srat %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- srat@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.076*nrow(srat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
srat <- doubletFinder_v3(srat, 
                         PCs = 1:30, 
                         pN = 0.25, 
                         pK = pK, 
                         nExp = nExp_poi.adj,
                         reuse.pANN = FALSE, sct = TRUE)
colnames(srat@meta.data)

DimPlot(srat, reduction = 'umap', group.by = "DF.classifications_0.25_0.1_700")
DimPlot(srat, reduction = 'umap', label=T)

# save BC cells-doublets
write.table(rownames(srat@meta.data[srat$DF.classifications_0.25_0.1_700 == 'Doublet',]), 
            file = str_c(sample, '_doublets.txt'), sep = "\t",
            row.names = TRUE, col.names = NA)

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)
table(srat$SCT_snn_res.0.8)
rp1 = DimPlot(srat, reduction = 'umap', label=T)
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
grid.arrange(rp1,rp2,rp3,rp4,rp5,rp6, ncol= 3)

srat = subset(srat, subset = DF.classifications_0.25_0.1_700 == 'Singlet')
dim(srat)
head(srat@meta.data)
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))
#  -------------------------------- find ALL markers --------------------------------
table(Idents(srat))
DefaultAssay(srat)
allMarkers <- FindAllMarkers(srat, max.cells.per.ident = 100, min.pct=0.25,
                             test.use = "MAST", only.pos = T, logfc.threshold = 0.25)

allMarkers$pct.diff <- allMarkers$pct.1 - allMarkers$pct.2
goodMarkers.pct.diff = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=pct.diff) %>% pull(gene)
goodMarkers.logfc = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=avg_log2FC) %>% pull(gene)

#DefaultAssay(ast.combined) = 'RNA'
#ast.combined = NormalizeData(ast.combined)
FeaturePlot(srat, goodMarkers.pct.diff, reduction="umap", ncol=6)
FeaturePlot(srat, goodMarkers.logfc, reduction="umap", ncol=6)

table(srat@meta.data$SCT_snn_res.0.05)
cell_markers <- c(
  # neurons:
  "Tmem130", "Snhg11", "Kcnip4", "Slc17a6", "RGD1566401",
  # oligodendrocytes: 
  "Mag", "Vcan",  "Tnr", "Bmp4", 
  # astrocytes: 
  "Gpc5", "Agt", 
  # pericytes: 
  "Flt1", "Rgs5", "Slco1c1", "Ebf1",
  # ependymocytes: 
  "Ccdc153", 
  # tanycytes:
  "Crym", 
  # microglia: 
  "Cd83", "Ctss", "Lyn",
  # activated MG or MFs: 
  "Mrc1", "Tpt1", 
  # leptomeningeal?: 
  "Slc47a1", "Igf2" 
)

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)

rp1 = DimPlot(srat, reduction = 'umap', label=T, group.by = 'SCT_snn_res.0.1')
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
grid.arrange(rp1,rp2,rp3,rp4,rp5,rp6, ncol= 3)
# save singlets
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))
dim(srat)








########################################################################################
################################# Rat4 #################################################
########################################################################################

# uplaod raw data
sample = 'Rat4'
counts = Read10X_h5('~/mstool/cellranger-7.0.0/Rat4_cellout/outs/filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)

# create seurat object
srat <- CreateSeuratObject(counts=counts, min.cells = 3, min.features = 200, project = sample)
srat$genotype = 'agrsv'
srat$sex = 'male'
head(srat@meta.data)

dim(srat)
#upload Rat's mito genes
mitogenes = fread(file = "~/novaseq-rats-analysis/mRatBN7.2_mitoGenes.txt", sep='\t', header = F)
mitogenes = mitogenes$V1[mitogenes$V1 %in% rownames(srat)]
srat[['percent.mt']] = PercentageFeatureSet(srat, features = mitogenes)
sum(srat$percent.mt > 10)
# add blood markers
bloodgenes = c("Hba-a1", "Hba-a2", "Hba-a3","Hbb", "Hbb-b1", "Hbb-bs","Hbe1", "Hbq1", "Hbq1b")
bloodgenes = bloodgenes[bloodgenes %in% rownames(srat)]
srat[['percent.hb']] <- PercentageFeatureSet(srat, features = bloodgenes)
sum(srat$percent.hb > 1)
hist(srat$percent.hb, breaks = 100000, xlim = range(0,10))
# add ribogenes
ribogenes = grep(pattern = "^Rp[sl]", x=rownames(x=srat), value = TRUE, ignore.case = T)
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^Rp[sl]")

#
VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.hb',"percent.rb" ), ncol=5)
FeatureScatter(srat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'glm')

# QC: filtering
srat <- subset(srat, subset = percent.mt < 20 & percent.hb < 1)
dim(srat)
#removing Hb- and ribogens
counts <- counts[-(which(rownames(counts) %in% c(bloodgenes, ribogenes))),]
srat <- subset(srat, features = rownames(counts))
dim(srat)
counts=0
#SCT
srat <- SCTransform(srat, method = "glmGamPoi", vars.to.regress = "percent.mt")
VariableFeaturePlot(srat) + scale_y_log10()
#PCA
srat <- RunPCA(object = srat)
ElbowPlot(srat,ndims = 50)
srat <- FindNeighbors(object = srat, dims = 1:30)
srat <- FindClusters(object = srat)
clustree(srat, prefix = "SCT_snn_res.")
Idents(srat) = "SCT_snn_res.0.05"
srat <- RunUMAP(object = srat, dims = 1:30)
## pK Identification (no ground-truth) 
sweep.res.list_srat <- paramSweep_v3(srat, PCs = 1:30, sct = TRUE)
sweep.stats_srat <- summarizeSweep(sweep.res.list_srat, GT = FALSE)
bcmvn_srat <- find.pK(sweep.stats_srat)

ggplot(bcmvn_srat, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_srat %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- srat@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.061*nrow(srat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
srat <- doubletFinder_v3(srat, 
                         PCs = 1:30, 
                         pN = 0.25, 
                         pK = pK, 
                         nExp = nExp_poi.adj,
                         reuse.pANN = FALSE, sct = TRUE)
colnames(srat@meta.data)
#srat <- doubletFinder_v3(srat, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.23_730", sct = TRUE)
# number of singlets and doublets
table(srat@meta.data$DF.classifications_0.25_0.21_439)


DimPlot(srat, reduction = 'umap', group.by = "DF.classifications_0.25_0.21_439")
DimPlot(srat, reduction = 'umap', label=T)

# save BC cells-doublets
write.table(rownames(srat@meta.data[srat$DF.classifications_0.25_0.21_439 == 'Doublet',]), 
            file = str_c(sample, '_doublets.txt'), sep = "\t",
            row.names = TRUE, col.names = NA)

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)
table(srat$SCT_snn_res.0.8)
rp1 = DimPlot(srat, reduction = 'umap', label=T)
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
grid.arrange(rp1,rp2,rp3,rp4,rp5,rp6, ncol= 3)

srat = subset(srat, subset = DF.classifications_0.25_0.21_439 == 'Singlet')
dim(srat)
head(srat@meta.data)
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))
#  -------------------------------- find ALL markers --------------------------------
table(Idents(srat))
DefaultAssay(srat)
allMarkers <- FindAllMarkers(srat, max.cells.per.ident = 100, min.pct=0.25,
                             test.use = "MAST", only.pos = T, logfc.threshold = 0.25)

allMarkers$pct.diff <- allMarkers$pct.1 - allMarkers$pct.2
goodMarkers.pct.diff = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=pct.diff) %>% pull(gene)
goodMarkers.logfc = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=avg_log2FC) %>% pull(gene)

#DefaultAssay(ast.combined) = 'RNA'
#ast.combined = NormalizeData(ast.combined)
FeaturePlot(srat, goodMarkers.pct.diff, reduction="umap", ncol=6)
FeaturePlot(srat, goodMarkers.logfc, reduction="umap", ncol=6)

table(srat@meta.data$SCT_snn_res.0.05)
cell_markers <- c(
  # neurons:
  "Tmem130", "Snhg11", "Kcnip4", "Slc17a6", "RGD1566401",
  # oligodendrocytes: 
  "Mag", "Vcan",  "Tnr", "Bmp4", 
  # astrocytes: 
  "Gpc5", "Agt", 
  # pericytes: 
  "Flt1", "Rgs5", "Slco1c1", "Ebf1",
  # ependymocytes: 
  "Ccdc153", 
  # tanycytes:
  "Crym", 
  # microglia: 
  "Cd83", "Ctss", "Lyn",
  # activated MG or MFs: 
  "Mrc1", "Tpt1", 
  # leptomeningeal?: 
  "Slc47a1", "Igf2" 
)

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)

rp1 = DimPlot(srat, reduction = 'umap', label=T, group.by = 'SCT_snn_res.0.05')
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
grid.arrange(rp1,rp2,rp3,rp4,rp5,rp6, ncol= 3)
# save singlets
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))
dim(srat)










########################################################################################
################################# Rat5 #################################################
########################################################################################

# uplaod raw data
sample = 'Rat5'
counts = Read10X_h5('~/mstool/cellranger-7.0.0/Rat5_cellout/outs/filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)

# create seurat object
srat <- CreateSeuratObject(counts=counts, min.cells = 3, min.features = 200, project = sample)
srat$genotype = 'agrsv'
srat$sex = 'male'
head(srat@meta.data)
dim(srat)
#upload Rat's mito genes
mitogenes = fread(file = "~/novaseq-rats-analysis/mRatBN7.2_mitoGenes.txt", sep='\t', header = F)
mitogenes = mitogenes$V1[mitogenes$V1 %in% rownames(srat)]
srat[['percent.mt']] = PercentageFeatureSet(srat, features = mitogenes)
sum(srat$percent.mt > 10)
# add blood markers
bloodgenes = c("Hba-a1", "Hba-a2", "Hba-a3","Hbb", "Hbb-b1", "Hbb-bs","Hbe1", "Hbq1", "Hbq1b")
bloodgenes = bloodgenes[bloodgenes %in% rownames(srat)]
srat[['percent.hb']] <- PercentageFeatureSet(srat, features = bloodgenes)
sum(srat$percent.hb > 1)
hist(srat$percent.hb, breaks = 100000, xlim = range(0,10))
# add ribogenes
ribogenes = grep(pattern = "^Rp[sl]", x=rownames(x=srat), value = TRUE, ignore.case = T)
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^Rp[sl]")

#
VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.hb',"percent.rb" ), ncol=5)
FeatureScatter(srat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'glm')

# QC: filtering
srat <- subset(srat, subset = percent.mt < 20 & percent.hb < 1)
dim(srat)
#removing Hb- and ribogens
counts <- counts[-(which(rownames(counts) %in% c(bloodgenes, ribogenes))),]
srat <- subset(srat, features = rownames(counts))
dim(srat)
counts=0
#SCT
srat <- SCTransform(srat, method = "glmGamPoi", vars.to.regress = "percent.mt")
VariableFeaturePlot(srat) + scale_y_log10()
#PCA
srat <- RunPCA(object = srat)
ElbowPlot(srat,ndims = 50)
srat <- FindNeighbors(object = srat, dims = 1:30)
srat <- FindClusters(object = srat)
clustree(srat, prefix = "SCT_snn_res.")
Idents(srat) = "SCT_snn_res.0.1"
srat <- RunUMAP(object = srat, dims = 1:30)
## pK Identification (no ground-truth) 
sweep.res.list_srat <- paramSweep_v3(srat, PCs = 1:30, sct = TRUE)
sweep.stats_srat <- summarizeSweep(sweep.res.list_srat, GT = FALSE)
bcmvn_srat <- find.pK(sweep.stats_srat)

ggplot(bcmvn_srat, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_srat %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- srat@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.039*nrow(srat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
srat <- doubletFinder_v3(srat, 
                         PCs = 1:30, 
                         pN = 0.25, 
                         pK = pK, 
                         nExp = nExp_poi.adj,
                         reuse.pANN = FALSE, sct = TRUE)
colnames(srat@meta.data)
#srat <- doubletFinder_v3(srat, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.23_730", sct = TRUE)
# number of singlets and doublets
table(srat@meta.data$DF.classifications_0.25_0.03_172)
#table(srat@meta.data$DF.classifications_0.25_0.09_599)

DimPlot(srat, reduction = 'umap', group.by = "DF.classifications_0.25_0.03_172")
DimPlot(srat, reduction = 'umap', label=T)
# save BC cells-doublets
write.table(rownames(srat@meta.data[srat$DF.classifications_0.25_0.03_172 == 'Doublet',]), 
            file = str_c(sample, '_doublets.txt'), sep = "\t",
            row.names = TRUE, col.names = NA)

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)
table(srat$SCT_snn_res.0.8)
rp1 = DimPlot(srat, reduction = 'umap', label=T)
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
grid.arrange(rp1,rp2,rp3,rp4,rp5,rp6, ncol= 3)

srat = subset(srat, subset = DF.classifications_0.25_0.03_172 == 'Singlet')
dim(srat)
head(srat@meta.data)
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))

#  -------------------------------- find ALL markers --------------------------------
table(Idents(srat))
DefaultAssay(srat)
allMarkers <- FindAllMarkers(srat, max.cells.per.ident = 100, min.pct=0.25,
                             test.use = "MAST", only.pos = T, logfc.threshold = 0.25)

allMarkers$pct.diff <- allMarkers$pct.1 - allMarkers$pct.2
goodMarkers.pct.diff = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=pct.diff) %>% pull(gene)
goodMarkers.logfc = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=avg_log2FC) %>% pull(gene)

#DefaultAssay(ast.combined) = 'RNA'
#ast.combined = NormalizeData(ast.combined)
FeaturePlot(srat, goodMarkers.pct.diff, reduction="umap", ncol=4)
FeaturePlot(srat, goodMarkers.logfc, reduction="umap", ncol=4)

table(srat@meta.data$SCT_snn_res.0.05)
cell_markers <- c(
  # neurons:
  "Tmem130", "Snhg11", "Kcnip4", "Slc17a6", "RGD1566401",
  # oligodendrocytes: 
  "Mag", "Vcan",  "Tnr", "Bmp4", 
  # astrocytes: 
  "Gpc5", "Agt", 
  # pericytes: 
  "Flt1", "Rgs5", "Slco1c1", "Ebf1",
  # ependymocytes: 
  "Ccdc153", 
  # tanycytes:
  "Crym", 
  # microglia: 
  "Cd83", "Ctss", "Lyn",
  # activated MG or MFs: 
  "Mrc1", "Tpt1", 
  # leptomeningeal?: 
  "Slc47a1", "Igf2" 
)

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)

rp1 = DimPlot(srat, reduction = 'umap', label=T, group.by = 'SCT_snn_res.0.1')
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
grid.arrange(rp1,rp2,rp3,rp4,rp5,rp6, ncol= 3)
# save singlets
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))
dim(srat)






########################################################################################
################################# Rat8 #################################################
########################################################################################

# uplaod raw data
sample = 'Rat8'
counts = Read10X_h5('~/mstool/cellranger-7.0.0/Rat8_cellout/outs/filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)

# create seurat object
srat <- CreateSeuratObject(counts=counts, min.cells = 3, min.features = 200, project = sample)
srat$genotype = 'tamed'
srat$sex = 'male'
head(srat@meta.data)
dim(srat)
#upload Rat's mito genes
mitogenes = fread(file = "~/novaseq-rats-analysis/mRatBN7.2_mitoGenes.txt", sep='\t', header = F)
mitogenes = mitogenes$V1[mitogenes$V1 %in% rownames(srat)]
srat[['percent.mt']] = PercentageFeatureSet(srat, features = mitogenes)
sum(srat$percent.mt > 10)
# add blood markers
bloodgenes = c("Hba-a1", "Hba-a2", "Hba-a3","Hbb", "Hbb-b1", "Hbb-bs","Hbe1", "Hbq1", "Hbq1b")
bloodgenes = bloodgenes[bloodgenes %in% rownames(srat)]
srat[['percent.hb']] <- PercentageFeatureSet(srat, features = bloodgenes)
sum(srat$percent.hb > 1)
hist(srat$percent.hb, breaks = 100000, xlim = range(0,10))
# add ribogenes
ribogenes = grep(pattern = "^Rp[sl]", x=rownames(x=srat), value = TRUE, ignore.case = T)
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^Rp[sl]")

#
VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.hb',"percent.rb" ), ncol=5)
FeatureScatter(srat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'glm')

# QC: filtering
srat <- subset(srat, subset = percent.mt < 20 & percent.hb < 1)
dim(srat)
#removing Hb- and ribogens
counts <- counts[-(which(rownames(counts) %in% c(bloodgenes, ribogenes))),]
srat <- subset(srat, features = rownames(counts))
dim(srat)
counts=0
#SCT
srat <- SCTransform(srat, method = "glmGamPoi", vars.to.regress = "percent.mt")
VariableFeaturePlot(srat) + scale_y_log10()
#PCA
srat <- RunPCA(object = srat)
ElbowPlot(srat,ndims = 50)
srat <- FindNeighbors(object = srat, dims = 1:30)
srat <- FindClusters(object = srat)
clustree(srat, prefix = "SCT_snn_res.")
Idents(srat) = "SCT_snn_res.0.3"
srat <- RunUMAP(object = srat, dims = 1:30)
## pK Identification (no ground-truth) 
sweep.res.list_srat <- paramSweep_v3(srat, PCs = 1:30, sct = TRUE)
sweep.stats_srat <- summarizeSweep(sweep.res.list_srat, GT = FALSE)
bcmvn_srat <- find.pK(sweep.stats_srat)

ggplot(bcmvn_srat, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_srat %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- srat@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.025*nrow(srat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
srat <- doubletFinder_v3(srat, 
                         PCs = 1:30, 
                         pN = 0.25, 
                         pK = pK, 
                         nExp = nExp_poi.adj,
                         reuse.pANN = FALSE, sct = TRUE)
colnames(srat@meta.data)
#srat <- doubletFinder_v3(srat, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.23_730", sct = TRUE)
# number of singlets and doublets
table(srat@meta.data$DF.classifications_0.25_0.23_52)
#table(srat@meta.data$DF.classifications_0.25_0.09_599)

DimPlot(srat, reduction = 'umap', group.by = "DF.classifications_0.25_0.23_52")
DimPlot(srat, reduction = 'umap', label=T)

# save BC cells-doublets
write.table(rownames(srat@meta.data[srat$DF.classifications_0.25_0.23_52 == 'Doublet',]), 
            file = str_c(sample, '_doublets.txt'), sep = "\t",
            row.names = TRUE, col.names = NA)

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)
table(srat$SCT_snn_res.0.8)
rp1 = DimPlot(srat, reduction = 'umap', label=T)
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
grid.arrange(rp1,rp2,rp3,rp4,rp5,rp6, ncol= 3)

srat = subset(srat, subset = DF.classifications_0.25_0.23_52 == 'Singlet')
dim(srat)
head(srat@meta.data)
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))

#  -------------------------------- find ALL markers --------------------------------
table(Idents(srat))
DefaultAssay(srat)
allMarkers <- FindAllMarkers(srat, max.cells.per.ident = 100, min.pct=0.25,
                             test.use = "MAST", only.pos = T, logfc.threshold = 0.25)

allMarkers$pct.diff <- allMarkers$pct.1 - allMarkers$pct.2
goodMarkers.pct.diff = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=pct.diff) %>% pull(gene)
goodMarkers.logfc = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=avg_log2FC) %>% pull(gene)

#DefaultAssay(ast.combined) = 'RNA'
#ast.combined = NormalizeData(ast.combined)
FeaturePlot(srat, goodMarkers.pct.diff, reduction="umap", ncol=4)
FeaturePlot(srat, goodMarkers.logfc, reduction="umap", ncol=4)

#table(srat@meta.data$SCT_snn_res.0.05)
cell_markers <- c(
  # neurons:
  "Tmem130", "Snhg11", "Kcnip4", "Slc17a6", "RGD1566401",
  # oligodendrocytes: 
  "Mag", "Vcan",  "Tnr", "Bmp4", 
  # astrocytes: 
  "Gpc5", "Agt", 
  # pericytes: 
  "Flt1", "Rgs5", "Slco1c1", "Ebf1",
  # ependymocytes: 
  "Ccdc153", 
  # tanycytes:
  "Crym", 
  # microglia: 
  "Cd83", "Ctss", "Lyn",
  # activated MG or MFs: 
  "Mrc1", "Tpt1", 
  # leptomeningeal?: 
  "Slc47a1", "Igf2" 
)

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)

rp1 = DimPlot(srat, reduction = 'umap', label=T, group.by = 'SCT_snn_res.0.3')
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
grid.arrange(rp1,rp2,rp3,rp4,rp5,rp6, ncol= 3)
# save singlets
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))
dim(srat)









########################################################################################
################################# K2_1 #################################################
########################################################################################

# uplaod raw data
sample = 'K2_1'
counts = Read10X_h5('~/mstool/cellranger-7.0.0/K2_1_cellout_7000/outs/filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)

# create seurat object
srat <- CreateSeuratObject(counts=counts, min.cells = 2, min.features = 200, project = sample)
srat$genotype = 'tamed'
srat$sex = 'mixed'
head(srat@meta.data)
dim(srat)
#upload Rat's mito genes
mitogenes = fread(file = "~/novaseq-rats-analysis/mRatBN7.2_mitoGenes.txt", sep='\t', header = F)
mitogenes = mitogenes$V1[mitogenes$V1 %in% rownames(srat)]
srat[['percent.mt']] = PercentageFeatureSet(srat, features = mitogenes)
sum(srat$percent.mt > 10)
# add blood markers
bloodgenes = c("Hba-a1", "Hba-a2", "Hba-a3","Hbb", "Hbb-b1", "Hbb-bs","Hbe1", "Hbq1", "Hbq1b")
bloodgenes = bloodgenes[bloodgenes %in% rownames(srat)]
srat[['percent.hb']] <- PercentageFeatureSet(srat, features = bloodgenes)
sum(srat$percent.hb > 1)
hist(srat$percent.hb, breaks = 100000, xlim = range(0,10))
# add ribogenes
ribogenes = grep(pattern = "^Rp[sl]", x=rownames(x=srat), value = TRUE, ignore.case = T)
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^Rp[sl]")

#
VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.hb',"percent.rb" ), ncol=5)
FeatureScatter(srat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'glm')

# QC: filtering
srat <- subset(srat, subset = percent.mt < 20 & percent.hb < 1)
dim(srat)
#removing Hb- and ribogens
counts <- counts[-(which(rownames(counts) %in% c(bloodgenes, ribogenes))),]
srat <- subset(srat, features = rownames(counts))
dim(srat)
counts=0
#SCT
srat <- SCTransform(srat, method = "glmGamPoi", vars.to.regress = "percent.mt")
VariableFeaturePlot(srat) + scale_y_log10()
#PCA
srat <- RunPCA(object = srat)
ElbowPlot(srat,ndims = 50)
srat <- FindNeighbors(object = srat, dims = 1:30)
srat <- FindClusters(object = srat)
clustree(srat, prefix = "SCT_snn_res.")
Idents(srat) = "SCT_snn_res.0.1"
srat <- RunUMAP(object = srat, dims = 1:30)
## pK Identification (no ground-truth) 
sweep.res.list_srat <- paramSweep_v3(srat, PCs = 1:30, sct = TRUE)
sweep.stats_srat <- summarizeSweep(sweep.res.list_srat, GT = FALSE)
bcmvn_srat <- find.pK(sweep.stats_srat)

ggplot(bcmvn_srat, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_srat %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- srat@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.039*nrow(srat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
srat <- doubletFinder_v3(srat, 
                         PCs = 1:30, 
                         pN = 0.25, 
                         pK = pK, 
                         nExp = nExp_poi.adj,
                         reuse.pANN = FALSE, sct = TRUE)
colnames(srat@meta.data)

# number of singlets and doublets
table(srat@meta.data$DF.classifications_0.25_0.01_166)


DimPlot(srat, reduction = 'umap', group.by = "DF.classifications_0.25_0.01_166")
DimPlot(srat, reduction = 'umap', label=T)
# save BC cells-doublets
write.table(rownames(srat@meta.data[srat$DF.classifications_0.25_0.01_166 == 'Doublet',]), 
            file = str_c(sample, '_doublets.txt'), sep = "\t",
            row.names = TRUE, col.names = NA)

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)
table(srat$SCT_snn_res.0.8)
rp1 = DimPlot(srat, reduction = 'umap', label=T)
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
grid.arrange(rp1,rp2,rp3,rp4,rp5,rp6, ncol= 3)

FeatureScatter(srat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'glm')

srat = subset(srat, subset = DF.classifications_0.25_0.01_166 == 'Singlet')
dim(srat)
head(srat@meta.data)
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))
#  -------------------------------- find ALL markers --------------------------------
DefaultAssay(srat) = 'RNA'
#SCT
srat <- SCTransform(srat, method = "glmGamPoi", vars.to.regress = "percent.mt")
VariableFeaturePlot(srat) + scale_y_log10()
#PCA
srat <- RunPCA(object = srat)
ElbowPlot(srat,ndims = 50)
srat <- FindNeighbors(object = srat, dims = 1:30)
srat <- FindClusters(object = srat, resolution = 0.1)
clustree(srat, prefix = "SCT_snn_res.")
Idents(srat) = "SCT_snn_res.0.1"
srat <- RunUMAP(object = srat, dims = 1:30)
table(Idents(srat))
DefaultAssay(srat)

allMarkers <- FindAllMarkers(srat, max.cells.per.ident = 100, min.pct=0.25,
                             test.use = "MAST", only.pos = T, logfc.threshold = 0.25)

allMarkers$pct.diff <- allMarkers$pct.1 - allMarkers$pct.2
allMarkers[allMarkers$p_val_adj<0.05,]
goodMarkers.pct.diff = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=pct.diff) %>% pull(gene)
table(goodMarkers.pct.diff$cluster)
goodMarkers.logfc = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=avg_log2FC) %>% pull(gene)

DimPlot(srat, reduction = 'umap', label=T)
FeaturePlot(srat, goodMarkers.pct.diff, reduction="umap", ncol=4)
FeaturePlot(srat, goodMarkers.logfc, reduction="umap", ncol=4)

#table(srat@meta.data$SCT_snn_res.0.05)
cell_markers <- c(
  # neurons:
  "Tmem130", "Snhg11", "Kcnip4", "Slc17a6", "RGD1566401",
  # oligodendrocytes: 
  "Mag", "Vcan",  "Tnr", "Bmp4", 
  # astrocytes: 
  "Gpc5", "Agt", 
  # pericytes: 
  "Flt1", "Rgs5", "Slco1c1", "Ebf1",
  # ependymocytes: 
  "Ccdc153", 
  # tanycytes:
  "Crym", 
  # microglia: 
  "Cd83", "Ctss", "Lyn",
  # activated MG or MFs: 
  "Mrc1", "Tpt1", 
  # leptomeningeal?: 
  "Slc47a1", "Igf2" 
)

########################################################################################
################################# K2_2 #################################################
########################################################################################
gc()
# uplaod raw data
sample = 'K2_2'
counts = Read10X_h5('~/mstool/cellranger-7.0.0/K2_2_cellout_7000/outs/filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)

# create seurat object
srat <- CreateSeuratObject(counts=counts, min.cells = 2, min.features = 200, project = sample)
srat$genotype = 'agrsv'
srat$sex = 'mixed'
head(srat@meta.data)
dim(srat)
#upload Rat's mito genes
mitogenes = fread(file = "~/novaseq-rats-analysis/mRatBN7.2_mitoGenes.txt", sep='\t', header = F)
mitogenes = mitogenes$V1[mitogenes$V1 %in% rownames(srat)]
srat[['percent.mt']] = PercentageFeatureSet(srat, features = mitogenes)
sum(srat$percent.mt > 10)
# add blood markers
bloodgenes = c("Hba-a1", "Hba-a2", "Hba-a3","Hbb", "Hbb-b1", "Hbb-bs","Hbe1", "Hbq1", "Hbq1b")
bloodgenes = bloodgenes[bloodgenes %in% rownames(srat)]
srat[['percent.hb']] <- PercentageFeatureSet(srat, features = bloodgenes)
sum(srat$percent.hb > 1)
hist(srat$percent.hb, breaks = 100000, xlim = range(0,10))
# add ribogenes
ribogenes = grep(pattern = "^Rp[sl]", x=rownames(x=srat), value = TRUE, ignore.case = T)
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^Rp[sl]")

#
VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.hb',"percent.rb" ), ncol=5)
FeatureScatter(srat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'glm')

# QC: filtering
srat <- subset(srat, subset = percent.mt < 20 & percent.hb < 1)
dim(srat)
#removing Hb- and ribogens
counts <- counts[-(which(rownames(counts) %in% c(bloodgenes, ribogenes))),]
srat <- subset(srat, features = rownames(counts))
dim(srat)
counts=0
#SCT
srat <- SCTransform(srat, method = "glmGamPoi", vars.to.regress = "percent.mt")
VariableFeaturePlot(srat) + scale_y_log10()
#PCA
srat <- RunPCA(object = srat)
ElbowPlot(srat,ndims = 50)
srat <- FindNeighbors(object = srat, dims = 1:30)
srat <- FindClusters(object = srat)
clustree(srat, prefix = "SCT_snn_res.")
Idents(srat) = "SCT_snn_res.0.05"
srat <- RunUMAP(object = srat, dims = 1:30)
## pK Identification (no ground-truth) 
sweep.res.list_srat <- paramSweep_v3(srat, PCs = 1:30, sct = TRUE)
sweep.stats_srat <- summarizeSweep(sweep.res.list_srat, GT = FALSE)
bcmvn_srat <- find.pK(sweep.stats_srat)

ggplot(bcmvn_srat, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_srat %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- srat@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.054*nrow(srat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
srat <- doubletFinder_v3(srat, 
                         PCs = 1:30, 
                         pN = 0.25, 
                         pK = pK, 
                         nExp = nExp_poi.adj,
                         reuse.pANN = FALSE, sct = TRUE)
colnames(srat@meta.data)
#srat <- doubletFinder_v3(srat, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.23_730", sct = TRUE)
# number of singlets and doublets
table(srat@meta.data$DF.classifications_0.25_0.19_314)

DimPlot(srat, reduction = 'umap', group.by = "DF.classifications_0.25_0.19_314")
DimPlot(srat, reduction = 'umap', label=T)

# save BC cells-doublets
write.table(rownames(srat@meta.data[srat$DF.classifications_0.25_0.19_314 == 'Doublet',]), 
            file = str_c(sample, '_doublets.txt'), sep = "\t",
            row.names = TRUE, col.names = NA)

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)
table(srat$SCT_snn_res.0.8)
rp1 = DimPlot(srat, reduction = 'umap', label=T)
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
grid.arrange(rp1,rp2,rp3,rp4,rp5,rp6, ncol= 3)

FeatureScatter(srat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'glm')

srat = subset(srat, subset = DF.classifications_0.25_0.19_314 == 'Singlet')
dim(srat)
head(srat@meta.data)
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))

#  -------------------------------- find ALL markers --------------------------------
table(Idents(srat))
DefaultAssay(srat)
allMarkers <- FindAllMarkers(srat, max.cells.per.ident = 100, min.pct=0.25,
                             test.use = "MAST", only.pos = T, logfc.threshold = 0.25)

allMarkers$pct.diff <- allMarkers$pct.1 - allMarkers$pct.2
goodMarkers.pct.diff = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=pct.diff) %>% pull(gene)
table(goodMarkers.pct.diff$cluster)
goodMarkers.logfc = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=avg_log2FC) %>% pull(gene)

#DefaultAssay(ast.combined) = 'RNA'
#ast.combined = NormalizeData(ast.combined)
FeaturePlot(srat, goodMarkers.pct.diff[-17], reduction="umap", ncol=4)
FeaturePlot(srat, goodMarkers.logfc, reduction="umap", ncol=4)

#table(srat@meta.data$SCT_snn_res.0.05)
cell_markers <- c(
  # neurons:
  "Tmem130", "Snhg11", "Kcnip4", "Slc17a6", "RGD1566401",
  # oligodendrocytes: 
  "Mag", "Vcan",  "Tnr", "Bmp4", 
  # astrocytes: 
  "Gpc5", "Agt", 
  # pericytes: 
  "Flt1", "Rgs5", "Slco1c1", "Ebf1",
  # ependymocytes: 
  "Ccdc153", 
  # tanycytes:
  "Crym", 
  # microglia: 
  "Cd83", "Ctss", "Lyn",
  # activated MG or MFs: 
  "Mrc1", "Tpt1", 
  # leptomeningeal?: 
  "Slc47a1", "Igf2" 
)

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)

rp1 = DimPlot(srat, reduction = 'umap', label=T)
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
grid.arrange(rp1,rp2,rp3,rp4,rp5,rp6, ncol= 3)
# save singlets
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))
dim(srat)






########################################################################################
################################# K2_3 #################################################
########################################################################################

# uplaod raw data
sample = 'K2_3'
counts = Read10X_h5('~/mstool/cellranger-7.0.0/K2_3_cellout/outs/filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)

# create seurat object
srat <- CreateSeuratObject(counts=counts, min.cells = 3, min.features = 200, project = sample)
srat$genotype = 'tamed'
srat$sex = 'mixed'
head(srat@meta.data)
dim(srat)
#upload Rat's mito genes
mitogenes = fread(file = "~/novaseq-rats-analysis/mRatBN7.2_mitoGenes.txt", sep='\t', header = F)
mitogenes = mitogenes$V1[mitogenes$V1 %in% rownames(srat)]
srat[['percent.mt']] = PercentageFeatureSet(srat, features = mitogenes)
sum(srat$percent.mt > 10)
# add blood markers
bloodgenes = c("Hba-a1", "Hba-a2", "Hba-a3","Hbb", "Hbb-b1", "Hbb-bs","Hbe1", "Hbq1", "Hbq1b")
bloodgenes = bloodgenes[bloodgenes %in% rownames(srat)]
srat[['percent.hb']] <- PercentageFeatureSet(srat, features = bloodgenes)
sum(srat$percent.hb > 1)
hist(srat$percent.hb, breaks = 100000, xlim = range(0,10))
# add ribogenes
ribogenes = grep(pattern = "^Rp[sl]", x=rownames(x=srat), value = TRUE, ignore.case = T)
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^Rp[sl]")

#
VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.hb',"percent.rb" ), ncol=5)
FeatureScatter(srat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'glm')

# QC: filtering
srat <- subset(srat, subset = percent.mt < 20 & percent.hb < 1 )
dim(srat)
#removing Hb- and ribogens
counts <- counts[-(which(rownames(counts) %in% c(bloodgenes, ribogenes))),]
srat <- subset(srat, features = rownames(counts))
dim(srat)
counts=0
#SCT
srat <- SCTransform(srat, method = "glmGamPoi", vars.to.regress = "percent.mt")
VariableFeaturePlot(srat) + scale_y_log10()
#PCA
srat <- RunPCA(object = srat)
ElbowPlot(srat,ndims = 50)
srat <- FindNeighbors(object = srat, dims = 1:30)
srat <- FindClusters(object = srat)
clustree(srat, prefix = "SCT_snn_res.")
Idents(srat) = "SCT_snn_res.0.05"
srat <- RunUMAP(object = srat, dims = 1:30)
## pK Identification (no ground-truth) 
sweep.res.list_srat <- paramSweep_v3(srat, PCs = 1:30, sct = TRUE)
sweep.stats_srat <- summarizeSweep(sweep.res.list_srat, GT = FALSE)
bcmvn_srat <- find.pK(sweep.stats_srat)

ggplot(bcmvn_srat, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_srat %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- srat@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.076*nrow(srat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
srat <- doubletFinder_v3(srat, 
                         PCs = 1:30, 
                         pN = 0.25, 
                         pK = pK, 
                         nExp = nExp_poi.adj,
                         reuse.pANN = FALSE, sct = TRUE)
colnames(srat@meta.data)
#srat <- doubletFinder_v3(srat, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.23_730", sct = TRUE)
# number of singlets and doublets
table(srat@meta.data$DF.classifications_0.25_0.005_755)


DimPlot(srat, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_755")
DimPlot(srat, reduction = 'umap', label=T)

# save BC cells-doublets
write.table(rownames(srat@meta.data[srat$DF.classifications_0.25_0.005_755 == 'Doublet',]), 
            file = str_c(sample, '_doublets.txt'), sep = "\t",
            row.names = TRUE, col.names = NA)

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)
table(srat$SCT_snn_res.0.8)
rp1 = DimPlot(srat, reduction = 'umap', label=T)
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
grid.arrange(rp1,rp2,rp3,rp4,rp5,rp6, ncol= 3)

FeatureScatter(srat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'glm')

srat = subset(srat, subset = DF.classifications_0.25_0.005_755 == 'Singlet')
dim(srat)
head(srat@meta.data)
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))
#  -------------------------------- find ALL markers --------------------------------
table(Idents(srat))
DefaultAssay(srat)
allMarkers <- FindAllMarkers(srat, max.cells.per.ident = 100, min.pct=0.25,
                             test.use = "MAST", only.pos = T, logfc.threshold = 0.25)

allMarkers$pct.diff <- allMarkers$pct.1 - allMarkers$pct.2
goodMarkers.pct.diff = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=pct.diff) %>% pull(gene)

goodMarkers.logfc = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=avg_log2FC) %>% pull(gene)

#DefaultAssay(ast.combined) = 'RNA'
#ast.combined = NormalizeData(ast.combined)
FeaturePlot(srat, goodMarkers.pct.diff, reduction="umap", ncol=6)
FeaturePlot(srat, goodMarkers.logfc, reduction="umap", ncol=6)

#table(srat@meta.data$SCT_snn_res.0.05)
cell_markers <- c(
  # neurons:
  "Tmem130", "Snhg11", "Kcnip4", "Slc17a6", "RGD1566401",
  # oligodendrocytes: 
  "Mag", "Vcan",  "Tnr", "Bmp4", 
  # astrocytes: 
  "Gpc5", "Agt", 
  # pericytes: 
  "Flt1", "Rgs5", "Slco1c1", "Ebf1",
  # ependymocytes: 
  "Ccdc153", 
  # tanycytes:
  "Crym", 
  # microglia: 
  "Cd83", "Ctss", "Lyn",
  # activated MG or MFs: 
  "Mrc1", "Tpt1", 
  # leptomeningeal?: 
  "Slc47a1", "Igf2" 
)

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)

rp1 = DimPlot(srat, reduction = 'umap', label=T, group.by = 'SCT_snn_res.0.05')
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
grid.arrange(rp1,rp2,rp3,rp4,rp5,rp6, ncol= 3)
# save singlets
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))
dim(srat)








########################################################################################
################################# K2_4 #################################################
########################################################################################

# uplaod raw data
sample = 'K2_4'
counts = Read10X_h5('~/mstool/cellranger-7.0.0/K2_4_cellout/outs/filtered_feature_bc_matrix.h5',
                    use.names = TRUE, unique.features = TRUE)

# create seurat object
srat <- CreateSeuratObject(counts=counts, min.cells = 2, min.features = 200, project = sample)
srat$genotype = 'agrsv'
srat$sex = 'mixed'
head(srat@meta.data)
dim(srat)
#upload Rat's mito genes
mitogenes = fread(file = "~/novaseq-rats-analysis/mRatBN7.2_mitoGenes.txt", sep='\t', header = F)
mitogenes = mitogenes$V1[mitogenes$V1 %in% rownames(srat)]
srat[['percent.mt']] = PercentageFeatureSet(srat, features = mitogenes)
sum(srat$percent.mt > 10)
# add blood markers
bloodgenes = c("Hba-a1", "Hba-a2", "Hba-a3","Hbb", "Hbb-b1", "Hbb-bs","Hbe1", "Hbq1", "Hbq1b")
bloodgenes = bloodgenes[bloodgenes %in% rownames(srat)]
srat[['percent.hb']] <- PercentageFeatureSet(srat, features = bloodgenes)
sum(srat$percent.hb > 1)
hist(srat$percent.hb, breaks = 100000, xlim = range(0,10))
# add ribogenes
ribogenes = grep(pattern = "^Rp[sl]", x=rownames(x=srat), value = TRUE, ignore.case = T)
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^Rp[sl]")

#
VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.hb',"percent.rb" ), ncol=5)
FeatureScatter(srat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'glm')

# QC: filtering
srat <- subset(srat, subset = percent.mt < 20 & percent.hb < 1 )
dim(srat)
#removing Hb- and ribogens
counts <- counts[-(which(rownames(counts) %in% c(bloodgenes, ribogenes))),]
srat <- subset(srat, features = rownames(counts))
dim(srat)
counts=0
#SCT
srat <- SCTransform(srat, method = "glmGamPoi", vars.to.regress = "percent.mt")
VariableFeaturePlot(srat) + scale_y_log10()
#PCA
srat <- RunPCA(object = srat)
ElbowPlot(srat,ndims = 50)
srat <- FindNeighbors(object = srat, dims = 1:30)
srat <- FindClusters(object = srat)
clustree(srat, prefix = "SCT_snn_res.")
Idents(srat) = "SCT_snn_res.0.1"
srat <- RunUMAP(object = srat, dims = 1:30)
## pK Identification (no ground-truth) 
sweep.res.list_srat <- paramSweep_v3(srat, PCs = 1:30, sct = TRUE)
sweep.stats_srat <- summarizeSweep(sweep.res.list_srat, GT = FALSE)
bcmvn_srat <- find.pK(sweep.stats_srat)

ggplot(bcmvn_srat, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_srat %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- srat@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.1*nrow(srat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
srat <- doubletFinder_v3(srat, 
                         PCs = 1:30, 
                         pN = 0.25, 
                         pK = pK, 
                         nExp = nExp_poi.adj,
                         reuse.pANN = FALSE, sct = TRUE)
colnames(srat@meta.data)
#srat <- doubletFinder_v3(srat, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.23_730", sct = TRUE)
# number of singlets and doublets
table(srat@meta.data$DF.classifications_0.25_0.02_1528)
#table(srat@meta.data$DF.classifications_0.25_0.09_599)

DimPlot(srat, reduction = 'umap', group.by = "DF.classifications_0.25_0.02_1528")
DimPlot(srat, reduction = 'umap', label=T)

# save BC cells-doublets
write.table(rownames(srat@meta.data[srat$DF.classifications_0.25_0.02_1528 == 'Doublet',]), 
            file = str_c(sample, '_doublets.txt'), sep = "\t",
            row.names = TRUE, col.names = NA)

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)
table(srat$SCT_snn_res.0.8)
rp1 = DimPlot(srat, reduction = 'umap', label=T)
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
grid.arrange(rp1,rp2,rp3,rp4,rp5,rp6, ncol= 3)

FeatureScatter(srat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'glm')

srat = subset(srat, subset = DF.classifications_0.25_0.02_1528 == 'Singlet')
dim(srat)
head(srat@meta.data)
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))
#  -------------------------------- find ALL markers --------------------------------
table(Idents(srat))
DefaultAssay(srat)
allMarkers <- FindAllMarkers(srat, max.cells.per.ident = 100, min.pct=0.25,
                             test.use = "MAST", only.pos = T, logfc.threshold = 0.25)

allMarkers$pct.diff <- allMarkers$pct.1 - allMarkers$pct.2
goodMarkers.pct.diff = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=pct.diff) %>% pull(gene)

goodMarkers.logfc = allMarkers %>% group_by(cluster) %>% top_n(n = 2, wt=avg_log2FC) %>% pull(gene)

#DefaultAssay(ast.combined) = 'RNA'
#ast.combined = NormalizeData(ast.combined)
FeaturePlot(srat, goodMarkers.pct.diff, reduction="umap", ncol=6)
FeaturePlot(srat, goodMarkers.logfc, reduction="umap", ncol=6)

#table(srat@meta.data$SCT_snn_res.0.05)
cell_markers <- c(
  # neurons:
  "Tmem130", "Snhg11", "Kcnip4", "Slc17a6", "RGD1566401",
  # oligodendrocytes: 
  "Mag", "Vcan",  "Tnr", "Bmp4", 
  # astrocytes: 
  "Gpc5", "Agt", 
  # pericytes: 
  "Flt1", "Rgs5", "Slco1c1", "Ebf1",
  # ependymocytes: 
  "Ccdc153", 
  # tanycytes:
  "Crym", 
  # microglia: 
  "Cd83", "Ctss", "Lyn",
  # activated MG or MFs: 
  "Mrc1", "Tpt1", 
  # leptomeningeal?: 
  "Slc47a1", "Igf2" 
)

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)

rp1 = DimPlot(srat, reduction = 'umap', label=T, group.by = 'SCT_snn_res.0.05')
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
grid.arrange(rp1,rp2,rp3,rp4,rp5,rp6, ncol= 3)
# save singlets
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))
dim(srat)

















########################################################################################
################################# K2_5 #################################################
########################################################################################

# uplaod raw data
sample = 'K2_5'
counts = Read10X_h5('~/mstool/cellranger-7.0.0/K2_5_cellout/outs/filtered_feature_bc_matrix.h5',
                    use.names = TRUE, unique.features = TRUE)

# create seurat object
srat <- CreateSeuratObject(counts=counts, min.cells = 2, min.features = 200, project = sample)
srat$genotype = 'tamed'
srat$sex = 'mixed'
head(srat@meta.data)
dim(srat)
#upload Rat's mito genes
mitogenes = fread(file = "~/novaseq-rats-analysis/mRatBN7.2_mitoGenes.txt", sep='\t', header = F)
mitogenes = mitogenes$V1[mitogenes$V1 %in% rownames(srat)]
srat[['percent.mt']] = PercentageFeatureSet(srat, features = mitogenes)
sum(srat$percent.mt > 10)
# add blood markers
bloodgenes = c("Hba-a1", "Hba-a2", "Hba-a3","Hbb", "Hbb-b1", "Hbb-bs","Hbe1", "Hbq1", "Hbq1b")
bloodgenes = bloodgenes[bloodgenes %in% rownames(srat)]
srat[['percent.hb']] <- PercentageFeatureSet(srat, features = bloodgenes)
sum(srat$percent.hb > 1)
hist(srat$percent.hb, breaks = 100000, xlim = range(0,10))
# add ribogenes
ribogenes = grep(pattern = "^Rp[sl]", x=rownames(x=srat), value = TRUE, ignore.case = T)
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^Rp[sl]")

#
VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.hb',"percent.rb" ), ncol=5)
FeatureScatter(srat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'glm')

# QC: filtering
srat <- subset(srat, subset = percent.mt < 20 & percent.hb < 1 )
dim(srat)
#removing Hb- and ribogens
counts <- counts[-(which(rownames(counts) %in% c(bloodgenes, ribogenes))),]
srat <- subset(srat, features = rownames(counts))
dim(srat)
counts=0
#SCT
srat <- SCTransform(srat, method = "glmGamPoi", vars.to.regress = "percent.mt")
VariableFeaturePlot(srat) + scale_y_log10()
#PCA
srat <- RunPCA(object = srat)
ElbowPlot(srat,ndims = 50)
srat <- FindNeighbors(object = srat, dims = 1:30)
srat <- FindClusters(object = srat)
clustree(srat, prefix = "SCT_snn_res.")
Idents(srat) = "SCT_snn_res.0.1"
srat <- RunUMAP(object = srat, dims = 1:30)
## pK Identification (no ground-truth) 
sweep.res.list_srat <- paramSweep_v3(srat, PCs = 1:30, sct = TRUE)
sweep.stats_srat <- summarizeSweep(sweep.res.list_srat, GT = FALSE)
bcmvn_srat <- find.pK(sweep.stats_srat)

ggplot(bcmvn_srat, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_srat %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- srat@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.076*nrow(srat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
srat <- doubletFinder_v3(srat, 
                         PCs = 1:30, 
                         pN = 0.25, 
                         pK = pK, 
                         nExp = nExp_poi.adj,
                         reuse.pANN = FALSE, sct = TRUE)
colnames(srat@meta.data)
#srat <- doubletFinder_v3(srat, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.23_730", sct = TRUE)
# number of singlets and doublets
table(srat@meta.data$DF.classifications_0.25_0.3_676)
#table(srat@meta.data$DF.classifications_0.25_0.09_599)

# save BC cells-doublets
write.table(rownames(srat@meta.data[srat$DF.classifications_0.25_0.3_676 == 'Doublet',]), 
            file = str_c(sample, '_doublets.txt'), sep = "\t",
            row.names = TRUE, col.names = NA)

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)
table(srat$SCT_snn_res.0.8)
rp1 = DimPlot(srat, reduction = 'umap', label=T)
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
rp7 = DimPlot(srat, reduction = 'umap', group.by = "DF.classifications_0.25_0.3_676")
grid.arrange(rp1,rp7,rp2,rp3,rp4,rp5,rp6, ncol= 4)

FeatureScatter(srat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'glm')

srat = subset(srat, subset = DF.classifications_0.25_0.3_676 == 'Singlet')
dim(srat)
head(srat@meta.data)
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))

















########################################################################################
################################# K2_6 #################################################
########################################################################################

# uplaod raw data
sample = 'K2_6'
counts = Read10X_h5('~/mstool/cellranger-7.0.0/K2_6_cellout/outs/filtered_feature_bc_matrix.h5',
                    use.names = TRUE, unique.features = TRUE)

# create seurat object
srat <- CreateSeuratObject(counts=counts, min.cells = 2, min.features = 200, project = sample)
srat$genotype = 'agrsv'
srat$sex = 'mixed'
head(srat@meta.data)
dim(srat)
#upload Rat's mito genes
mitogenes = fread(file = "~/novaseq-rats-analysis/mRatBN7.2_mitoGenes.txt", sep='\t', header = F)
mitogenes = mitogenes$V1[mitogenes$V1 %in% rownames(srat)]
srat[['percent.mt']] = PercentageFeatureSet(srat, features = mitogenes)
sum(srat$percent.mt > 10)
# add blood markers
bloodgenes = c("Hba-a1", "Hba-a2", "Hba-a3","Hbb", "Hbb-b1", "Hbb-bs","Hbe1", "Hbq1", "Hbq1b")
bloodgenes = bloodgenes[bloodgenes %in% rownames(srat)]
srat[['percent.hb']] <- PercentageFeatureSet(srat, features = bloodgenes)
sum(srat$percent.hb > 1)
hist(srat$percent.hb, breaks = 100000, xlim = range(0,10))
# add ribogenes
ribogenes = grep(pattern = "^Rp[sl]", x=rownames(x=srat), value = TRUE, ignore.case = T)
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^Rp[sl]")

#
VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.hb',"percent.rb" ), ncol=5)
FeatureScatter(srat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'glm')

# QC: filtering
srat <- subset(srat, subset = percent.mt < 20 & percent.hb < 1 )
dim(srat)
#removing Hb- and ribogens
counts <- counts[-(which(rownames(counts) %in% c(bloodgenes, ribogenes))),]
srat <- subset(srat, features = rownames(counts))
dim(srat)
counts=0
#SCT
srat <- SCTransform(srat, method = "glmGamPoi", vars.to.regress = "percent.mt")
VariableFeaturePlot(srat) + scale_y_log10()
#PCA
srat <- RunPCA(object = srat)
ElbowPlot(srat,ndims = 50)
srat <- FindNeighbors(object = srat, dims = 1:30)
srat <- FindClusters(object = srat)
clustree(srat, prefix = "SCT_snn_res.")
Idents(srat) = "SCT_snn_res.0.1"
srat <- RunUMAP(object = srat, dims = 1:30)
## pK Identification (no ground-truth) 
sweep.res.list_srat <- paramSweep_v3(srat, PCs = 1:30, sct = TRUE)
sweep.stats_srat <- summarizeSweep(sweep.res.list_srat, GT = FALSE)
bcmvn_srat <- find.pK(sweep.stats_srat)

ggplot(bcmvn_srat, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_srat %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- srat@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.061*nrow(srat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
srat <- doubletFinder_v3(srat, 
                         PCs = 1:30, 
                         pN = 0.25, 
                         pK = pK, 
                         nExp = nExp_poi.adj,
                         reuse.pANN = FALSE, sct = TRUE)
colnames(srat@meta.data)
#srat <- doubletFinder_v3(srat, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.23_730", sct = TRUE)
# number of singlets and doublets
table(srat@meta.data$DF.classifications_0.25_0.005_413)
#table(srat@meta.data$DF.classifications_0.25_0.09_599)

# save BC cells-doublets
write.table(rownames(srat@meta.data[srat$DF.classifications_0.25_0.005_413 == 'Doublet',]), 
            file = str_c(sample, '_doublets.txt'), sep = "\t",
            row.names = TRUE, col.names = NA)

FeaturePlot(srat, cell_markers, reduction="umap", ncol=6)
table(srat$SCT_snn_res.0.8)
rp1 = DimPlot(srat, reduction = 'umap', label=T)
rp2 = FeaturePlot(srat, 'nFeature_RNA', reduction="umap")
rp3 = FeaturePlot(srat, 'nCount_RNA', reduction="umap")
rp4 = FeaturePlot(srat, 'percent.mt', reduction="umap")
rp5 = FeaturePlot(srat, 'percent.hb', reduction="umap")
rp6 = FeaturePlot(srat, 'percent.rb', reduction="umap")
rp7 = DimPlot(srat, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_413")
grid.arrange(rp1,rp7,rp2,rp3,rp4,rp5,rp6, ncol= 4)

FeatureScatter(srat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'glm')

srat = subset(srat, subset = DF.classifications_0.25_0.005_413 == 'Singlet')
dim(srat)
head(srat@meta.data)
saveRDS(object = srat, file = str_c(sample, "_doubletfinder_singlets.Rds"))

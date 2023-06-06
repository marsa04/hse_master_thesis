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

# subset scafe.astro
scafe.obj = readRDS("scafe.count.mtx.integrated.celltype.Rds")
DimPlot(scafe.obj, reduction = 'umap', label = T)
dim(scafe.obj)
DefaultAssay(scafe.obj) = 'RNA'
scafe.astro = subset(x = scafe.obj, idents = "Astrocytes")
dim(scafe.astro)
scafe.obj = NULL
gc()
# filtering scafe.astro by barcodes from astrocytes.lastvar.Rds
cr.astro = readRDS('../../../../singlets-data-after-doubletfinder/astrocytes/astrocytes.lastvar.Rds')
dim(cr.astro)
bc = colnames(cr.astro)
sum(colnames(scafe.astro) %in% bc)
sum(bc %in% colnames(scafe.astro))
scafe.astro <- subset(scafe.astro, cells = bc)
dim(scafe.astro)
sum(names(Idents(cr.astro)) == rownames(cr.astro[[]]))
cr.astro[['celltype']] = Idents(cr.astro)

colnames(scafe.astro[[]])
to.remove <- c("nCount_SCT","nFeature_SCT","SCT_snn_res.0.8","seurat_clusters", 
               "integrated_snn_res.0.1", "integrated_snn_res.0.3", "integrated_snn_res.0.5", "celltype")
scafe.astro$nCount_SCT = NULL
scafe.astro$nFeature_SCT = NULL
scafe.astro$SCT_snn_res.0.8 = NULL
scafe.astro$seurat_clusters = NULL
scafe.astro$seurat_clusters = NULL
scafe.astro$integrated_snn_res.0.1 = NULL
scafe.astro$integrated_snn_res.0.3 = NULL
scafe.astro$integrated_snn_res.0.5 = NULL
scafe.astro$celltype = NULL
scafe.astro@assays$SCT = NULL
scafe.astro@assays$integrated = NULL
scafe.astro@reductions$pca = NULL
scafe.astro@reductions$umap = NULL
DimPlot(scafe.astro)
#########################
#   split & integrate   #
#########################
ast = scafe.astro
scafe.astro = NULL
rownames(ast)[1]
DefaultAssay(ast)
# split the dataset into a list of samples
ast.list = SplitObject(ast, split.by = "orig.ident")
rm(ast)
gc()
# SCTransfoem for all samples
ast.list <- lapply(X = ast.list, function(x) SCTransform(x, 
                                                         vars.to.regress = c("nCount_RNA","nFeature_RNA"),
                                                         variable.features.n = 10000,
                                                         ncells = min(100000, ncol(x)),
                                                         method = "glmGamPoi",
                                                         #return.only.var.genes = T, # F in Sankowsky's script
                                                         verbose = T))
gc()
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ast.list, nfeatures = 30000)
# Perform integration. 
ast.list <- PrepSCTIntegration(object.list = ast.list, 
                               anchor.features = features, 
                               verbose = T)
ast.anchors <- FindIntegrationAnchors(object.list = ast.list, 
                                      normalization.method = "SCT", 
                                      anchor.features = features, 
                                      dims = 1:45, # influence a lot on layout and cells belonging to exact clusters
                                      verbose = T) # CCA
ast.integrated <- IntegrateData(anchorset = ast.anchors, 
                                normalization.method = "SCT", 
                                dims = 1:45, # influence a lot on layout and cells belonging to exact clusters
                                verbose = T)
gc()

# standart workflow analysis
ast.integrated <- RunPCA(object = ast.integrated)
ElbowPlot(ast.integrated,ndims = 50)
DefaultAssay(ast.integrated)
ast.integrated <- RunUMAP(object = ast.integrated, dims = 1:20)
DimPlot(ast.integrated)
# batch checking
p1 = DimPlot(ast.integrated, reduction = 'umap', group.by = 'smpl')
p2 = DimPlot(ast.integrated, reduction = 'umap', group.by = 'genotype')
bp0 = DimPlot(ast.integrated, reduction = 'umap', group.by = 'sex')
bp1 = FeaturePlot(ast.integrated, 'nFeature_RNA', reduction="umap")
bp2 = FeaturePlot(ast.integrated, 'nCount_RNA', reduction="umap")
#bp3 = FeaturePlot(ast.integrated, 'percent.mt', reduction="umap")
#bp4 = FeaturePlot(ast.integrated, 'percent.hb', reduction="umap")
#bp5 = FeaturePlot(ast.integrated, 'percent.rb', reduction="umap")
grid.arrange(p1,p2,bp1,bp2, ncol= 2)

ast.scafe = ast.integrated

# UMAP
ast.scafe = RunUMAP(object = ast.scafe, dims = 1:20, min.dist = 0.3, n.neighbors = 30)
DimPlot(ast.scafe)
# clustering
ast.scafe <- FindNeighbors(object = ast.scafe, dims = 1:20)
ast.scafe <- FindClusters(object = ast.scafe, resolution = 0.5)
DimPlot(ast.scafe, reduction = 'umap', label = T)
dim(ast.scafe)
#Transfer cell types
dim(cr.astro)
cr.astro.subset = subset(cr.astro, cells = colnames(ast.scafe))
dim(cr.astro.subset)
dim(ast.scafe)
sum(colnames(ast.scafe) == colnames(cr.astro.subset))
sum(rownames(ast.scafe[[]]) == rownames(cr.astro.subset[[]]))

ast.scafe[['celltype']] = cr.astro.subset$celltype
head(ast.scafe)
ast.scafe@active.ident = ast.scafe$celltype
unique(Idents(ast.scafe))
DimPlot(ast.scafe, reduction = 'umap', label = T)

saveRDS(ast.scafe, file = './astro.scafe/ast.scafe.nfeatures30000.Rds')


# Find dif tCRE between condition ----------------------------------------------
ast.scafe = readRDS('~/novaseq-rats-analysis/scafe/scafe.tool.sc.solo.aggregate.outs/count/integrated.objects/astro.scafe/ast.scafe.nfeatures30000.Rds')

DefaultAssay(ast.scafe) = 'RNA'
ast.scafe$celltype.genotype <- paste0(ast.scafe$celltype,'_', ast.scafe$genotype)
Idents(ast.scafe) <- ast.scafe$celltype.genotype

diftCRE.btwn.genotype = data.frame()
for(i in unique(ast.scafe$celltype)){
  print(i)
  degs = FindMarkers(ast.scafe, 
                     ident.1 = paste0(i,'_','agrsv'), 
                     ident.2 = paste0(i,'_','tamed'),
                     test.use = 'MAST')
  degs$celltype = i
  degs$gene = rownames(degs)
  diftCRE.btwn.genotype = rbind(diftCRE.btwn.genotype, degs)}

diftCRE.btwn.genotype = diftCRE.btwn.genotype %>% dplyr::filter(p_val_adj < 0.05)
for(i in unique(diftCRE.btwn.genotype$celltype)){
  print(i) 
  print(nrow(diftCRE.btwn.genotype[diftCRE.btwn.genotype$celltype == i,]))}

diftCRE.btwn.genotype.top100 = degs.btwn.genotype.rna %>% group_by(celltype) %>% top_n(n = 100, wt=avg_log2FC)


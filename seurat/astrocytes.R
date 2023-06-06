library(Seurat)
library(SeuratData)
library(glmGamPoi)
library(patchwork)
library(clustree)
library(gridExtra)
library(tidyr)
library(dplyr)

ast.raw = readRDS('ast.raw.integration.ccav2.subset.Rds')
DefaultAssay(ast.raw) = 'RNA'
dim(ast.raw) # 27106 27169

integrated.full.data = NULL
gc()
#########################
#   split & integrate   #
#########################
# split the dataset into a list of samples
ast.list = SplitObject(ast, split.by = "orig.ident")
rm(ast)
# SCTransfoem for all samples
ast.list <- lapply(X = ast.list, function(x) SCTransform(x, 
                                                         vars.to.regress = c("percent.mt"),
                                                         variable.features.n = 10000,
                                                         ncells = min(100000, ncol(x)),
                                                         method = "glmGamPoi",
                                                         #return.only.var.genes = T, # F in Sankowsky's script
                                                         verbose = T))
gc()
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ast.list, nfeatures = 3000)
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
ast.integrated <- RunUMAP(object = ast.integrated, dims = 1:30)
# batch checking
p1 = DimPlot(ast.integrated, reduction = 'umap', group.by = 'smpl')
p2 = DimPlot(ast.integrated, reduction = 'umap', group.by = 'genotype')
bp0 = DimPlot(ast.integrated, reduction = 'umap', group.by = 'sex')
bp1 = FeaturePlot(ast.integrated, 'nFeature_RNA', reduction="umap")
bp2 = FeaturePlot(ast.integrated, 'nCount_RNA', reduction="umap")
bp3 = FeaturePlot(ast.integrated, 'percent.mt', reduction="umap")
bp4 = FeaturePlot(ast.integrated, 'percent.hb', reduction="umap")
bp5 = FeaturePlot(ast.integrated, 'percent.rb', reduction="umap")
grid.arrange(p1,p2,bp0, bp1,bp2,bp3,bp4,bp5,ncol= 4)

table(ast.integrated@meta.data$genotype)

# clustering
ast.integrated <- FindNeighbors(object = ast.integrated, dims = 1:30)
ast.integrated <- FindClusters(object = ast.integrated, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
View(ast.integrated@meta.data)
p5 = DimPlot(ast.integrated, group.by="integrated_snn_res.0.1", label=T)
p6 = DimPlot(ast.integrated, group.by="integrated_snn_res.0.3", label=T)
p7 = DimPlot(ast.integrated, group.by="integrated_snn_res.0.5", label=T)
p8 = DimPlot(ast.integrated, group.by="integrated_snn_res.0.7", label=T)
p9 = DimPlot(ast.integrated, group.by="integrated_snn_res.1", label=T)

grid.arrange(p5,p6,p7,p8,p9, ncol= 3)
clustree(ast.integrated, prefix = "integrated_snn_res.")

saveRDS(object = ast.integrated, file = "astrocytes.reintegrated.v1.Rds")

# FIND MARKERS
## res 0.1
Idents(ast.integrated) = "integrated_snn_res.0.1"
str(ast.integrated)
DimPlot(ast.integrated, reduction = 'umap', label=T)
table(ast.integrated@meta.data$integrated_snn_res.0.1)
# find ALL markers
allMarkers <- FindAllMarkers(ast.integrated, max.cells.per.ident = 100, min.pct=0.25,
                                test.use = "MAST", only.pos = T, logfc.threshold = 0.25)
allMarkers$pct.diff <- allMarkers$pct.1 - allMarkers$pct.2
goodMarkers.pct.diff = allMarkers %>% group_by(cluster) %>% top_n(n = 5, wt=pct.diff) %>% pull(gene)

goodMarkers.logfc = allMarkers %>% group_by(cluster) %>% top_n(n = 5, wt=avg_log2FC) %>% pull(gene)


DefaultAssay(ast.integrated) = 'RNA'
ast.integrated = NormalizeData(ast.integrated)


astro_mkrs = c("Gpc5", "Slc6a11", "Slc1a3", "Itih3")
FeaturePlot(ast.integrated, astro_mkrs, reduction="umap", ncol=2)#, min.cutoff = 'q10')

#cluster 6 - PEricytes

cl6 = c("Ldlrad4",   "Sdk1",      "Sgms1",     "Zfp36",     "Mef2c",
        "Ldlrad4", "Tfrc", "Abcb1a", "Kdr", "Slco1a4" )
FeaturePlot(ast.integrated, cl6, reduction="umap", ncol=4)#, min.cutoff = 'q10')

# cluster 7 -- blood??

cl7 = c("Mx2","Rnf213","Mx1","Oasl2", "MGC108823", "Spats2l", "Stat1", "Ifi44", "LOC103690031")
FeaturePlot(ast.integrated, cl7, reduction="umap", ncol=3)#, min.cutoff = 'q10')


# cluster 5 - OG/OPC/Neu
cl5 = c("Nxph1", "Tmeff2", "Nkain2", "St18")
FeaturePlot(ast.integrated, cl5, reduction="umap", ncol=2)#, min.cutoff = 'q10')

general.mkrs = c("Selplg", "Cx3cr1", "Mrc1", "Ms4a4a","Pdgfra","Col5a3", "Tnr", "Gpr17", "Mag", "Gjc2",
                 "Snhg11", "Tmem130","Ppp1r1b", "Pacrg", "Dnah3", "Dnah5", "Flt1", "Slco1a4", "Col3a1", "Col1a1")

FeaturePlot(ast.integrated, general.mkrs, reduction="umap", ncol=5)#, min.cutoff = 'q10')


##########################################################################################
############################## SELECT TRUE ASTRO CLUSTERS ################################
##########################################################################################
ast.raw = readRDS('astrocytes.reintegrated.v1.Rds')
Idents(ast.raw) = "integrated_snn_res.0.1"
ast.raw = subset(x = ast.raw, idents = as.character(0:4))
DefaultAssay(ast.raw) = 'RNA'
dim(ast.raw) # 27106 26399


# subset re-integration 
# split the dataset into a list of samples
ast.list = SplitObject(ast.raw, split.by = "orig.ident")
rm(ast.raw)
# SCTransfoem for all samples
ast.list <- lapply(X = ast.list, function(x) SCTransform(x, 
                                                         vars.to.regress = c("percent.mt"),
                                                         variable.features.n = 10000,
                                                         ncells = min(100000, ncol(x)),
                                                         method = "glmGamPoi",
                                                         #return.only.var.genes = T, # F in Sankowsky's script
                                                         verbose = T))
gc()
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ast.list, nfeatures = 3000)
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
ast.list = NULL
ast.anchors = NULL
gc()

saveRDS(ast.integrated, 'ast.filtered_integrated.Rds')
ast = readRDS('ast.filtered_integrated.Rds')
#ast = ast.integrated

#ast@meta.data$integrated_snn_res.0.5 = NULL

DefaultAssay(ast) = 'integrated'
ast <- RunPCA(object = ast)
ElbowPlot(ast,ndims = 50)
# clustering
ast <- FindNeighbors(object = ast, dims = 1:30)
#ast <- FindClusters(object = ast, resolution = c(0.1, 0.3, 0.5, 0.7,0.85, 1))
ast <- FindClusters(object = ast, resolution = c(0.7, 0.85, 1, 1.2, 1.4, 2, 2.5, 3, 3.5, 4, 4.5,5))
clustree(ast, prefix = "integrated_snn_res.")
# UMAP
Idents(ast) = "integrated_snn_res.0.7"
ast <- RunUMAP(object = ast, dims = 1:30, min.dist = 0.3, n.neighbors = 30)
DimPlot(ast, reduction = 'umap',label=T)
#ast <- RunUMAP(object = ast, dims = 1:30, min.dist = 0.2, n.neighbors = 30)
#DimPlot(ast, reduction = 'umap',label=T)
#ast <- RunUMAP(object = ast, dims = 1:30, min.dist = 0.1, n.neighbors = 30)
#DimPlot(ast, reduction = 'umap',label=T)
###
#ast <- RunUMAP(object = ast, dims = 1:30, min.dist = 0.3, n.neighbors = 60)
#DimPlot(ast, reduction = 'umap',label=T)
#ast <- RunUMAP(object = ast, dims = 1:30, min.dist = 0.2, n.neighbors = 60)
#DimPlot(ast, reduction = 'umap',label=T)
#ast <- RunUMAP(object = ast, dims = 1:30, min.dist = 0.1, n.neighbors = 60)
#DimPlot(ast, reduction = 'umap',label=T)
###
#ast <- RunUMAP(object = ast, dims = 1:30, min.dist = 0.3, n.neighbors = 100)
#DimPlot(ast, reduction = 'umap',label=T)
#ast <- RunUMAP(object = ast, dims = 1:30, min.dist = 0.1, n.neighbors = 100)
#DimPlot(ast, reduction = 'umap',label=T)
#ast <- RunUMAP(object = ast, dims = 1:30, min.dist = 0.2, n.neighbors = 100)
#DimPlot(ast, reduction = 'umap',label=T)

# find markers
table(ast@meta.data$integrated_snn_res.0.7)

allMarkers07 <- FindAllMarkers(ast, max.cells.per.ident = 100, min.pct=0.25,
                             test.use = "MAST", only.pos = T, logfc.threshold = 0.25)
allMarkers07$pct.diff <- allMarkers07$pct.1 - allMarkers07$pct.2

goodMarkers.pct.diff = allMarkers07 %>% group_by(cluster) %>% top_n(n = 5, wt=pct.diff) %>% pull(gene)
#table(goodMarkers.pct.diff$cluster)
goodMarkers.logfc = allMarkers07 %>% group_by(cluster) %>% top_n(n = 5, wt=avg_log2FC) %>% pull(gene)

#find CONSERVED markers
res0.7_mkrs = data.frame()
for(i in 0:13){
  mkrs = FindConservedMarkers(ast, ident.1=i, grouping.var = 'genotype')
  mkrs$cluster = i
  mkrs$gene = rownames(mkrs)
  res0.7_mkrs = rbind(res0.7_mkrs, mkrs)}
unique(res0.7_mkrs$cluster)
View(res0.7_mkrs)

res0.7_mkrs$agrsv_pct.diff = res0.7_mkrs$agrsv_pct.1 - res0.7_mkrs$agrsv_pct.2
res0.7_mkrs$tamed_pct.diff = res0.7_mkrs$tamed_pct.1 - res0.7_mkrs$tamed_pct.2
#pct.diff
top5_res0.7_agrsv_pct.diff = res0.7_mkrs %>% group_by(cluster) %>% top_n(n = 5, wt=agrsv_pct.diff) %>% pull(gene)
table(top5_res0.7_agrsv_pct.diff$cluster)
top5_res0.7_tamed_pct.diff = res0.7_mkrs %>% group_by(cluster) %>% top_n(n = 5, wt=tamed_pct.diff) %>% pull(gene)
#table(top5_res0.1_tamed_pct.diff$cluster)
#logfc
top5_res0.7_agrsv_logfc = res0.7_mkrs %>% group_by(cluster) %>% top_n(n = 5, wt=agrsv_avg_log2FC) %>% pull(gene)
#table(top5_res0.1_agrsv_logfc$cluster)
top5_res0.7_tamed_logfc = res0.7_mkrs %>% group_by(cluster) %>% top_n(n = 5, wt=tamed_avg_log2FC) %>% pull(gene)
#table(top5_res0.1_tamed_logfc$cluster)

DefaultAssay(ast) = 'RNA'
ast = NormalizeData(ast)


#CLUSTER 0 -- 
cells0 = 4748
barplot(table(ast@meta.data[ast@meta.data$integrated_snn_res.0.7 == 0,]$genotype) / cells0)
barplot(table(ast@meta.data[ast@meta.data$integrated_snn_res.0.7 == 0,]$smpl) / cells0)
cl0 = c("Sgcd", "Grid2") 
FeaturePlot(ast, cl0, reduction="umap", ncol=4)#, min.cutoff = 'q10')

#CLUSTER1
cl1 = c("Kcnn2") 
FeaturePlot(ast, 'Nrxn3', reduction="umap", ncol=4, min.cutoff = 'q10')
#CLUSTER2
cl2 = c("Sparc" )
FeaturePlot(ast, cl2, reduction="umap", ncol=5, min.cutoff = 'q10')

#CLUSTER3
cl3 = c("Apod","Cldn11")
FeaturePlot(ast, cl3, reduction="umap", ncol=3)#, min.cutoff = 'q10')

#CLUSTER4
cl4 = c("LOC310926")
FeaturePlot(ast, cl4, reduction="umap", ncol=5)#, min.cutoff = 'q5')

#CLUSTER 5
cl5 = c('Efna5')
FeaturePlot(ast, cl5, reduction="umap", ncol=1)#, min.cutoff = 'q5')
#custer6
cl6 = c("Gpc6", "Lama2","Sgip1") # 
FeaturePlot(ast, cl6, reduction="umap", ncol=3)#, min.cutoff = 'q5')
#cluster7
cl7 = c("Slc38a1", "Aebp1","Gfap","A2m","Lcat") 
FeaturePlot(ast, cl7, reduction="umap", ncol=4)#, min.cutoff = 'q5')
# CLUSTER 8
cl8 = c("Cntnap2")
FeaturePlot(ast, cl8, reduction="umap", ncol=1)#, min.cutoff = 'q5')
#cluster 9
cl9 = c("Mbp", "Mobp")
#cluster10
cl10 = c("Lrp1b")

#cluster11
cl11 = c("Sgcz")
FeaturePlot(ast, cl11, reduction="umap", ncol=4)#, min.cutoff = 'q5')
#CLUSTER12
cl12 = c("'Fos', 'Jun'")
#CLUSTER13
cl13 = c("Kcnb2")
FeaturePlot(ast, cl13, reduction="umap", ncol=4)#, min.cutoff = 'q5')

# final mkrs
finalmkrs0.7 = c("Sgcd", "Grid2",
                 "Kcnn2",
                 "Sparc", 
                 "Apod","Cldn11",
                 "LOC310926", 
                 'Efna5',
                 "Gpc6",
                 "Slc38a1","Gfap",
                 "Cntnap2",
                 "Mbp", "Mobp",
                 "Lrp1b",
                 "Sgcz",
                 "Fos", "Jun",
                 "Kcnb2")

new.cluster.ids <- c("Sgcd+/Grid2+", 
                     "Kcnn2+",
                     "Lrp1b+",
                     "Sgcz+",
                     "Fos+/Jun+",
                     "Kcnb2+",
                     "Sparc+",       
                     "Apod+/Cldn11+", 
                     "LOC310926+",    
                     'Efna5+',        
                     "Gpc6+",         
                     "Slc38a1+/Gfap+",  
                     "Cntnap2+",
                     "Mbp+/Mobp+")
new.cluster.ids = c('0','1','10','11','12','13', '2', '3', '4', '5', '6','7','8','9')
names(new.cluster.ids) <- levels(ast)
ast <- RenameIdents(ast, new.cluster.ids)  

DimPlot(ast, reduction = "umap", label = T)
DotPlot(ast, features = finalmkrs0.7) + RotatedAxis()

# перегруппировать по нормальному кластеры в объекте!
dot_plot = DotPlot(ast, features = x) + RotatedAxis()
x = c('Sgcd', 'Grid2', 'Sgcz', 'Kcnn2', 'Gpc6', 'Fos', 'Jun', 'Sparc', 'Apod', 'Cldn11', 'LOC310926','Lrp1b',
      'Kcnb2', 'Efna5', 'Slc38a1', 'Gfap', 'Cntnap2', 'Mbp', 'Mobp')
dot_plot$data$id <- factor(x = dot_plot$data$id, levels = c("Sgcd+/Grid2+", "Sgcz+", "Kcnn2+", "Gpc6+","Fos+/Jun+",
                                                            "Sparc+", "Apod+/Cldn11+", "LOC310926+", "Lrp1b+",
                                                            "Kcnb2+", 'Efna5+', "Slc38a1+/Gfap+",
                                                            "Cntnap2+", "Mbp+/Mobp+"))
dot_plot

DimPlot(ast, label = T, label.size = 4) + NoLegend()


ast = readRDS('./astrocytes.lastvar.Rds')
SaveH5Seurat(ast, filename = "./ast.lastvar.h5Seurat")

FeaturePlot(ast, cl12, reduction="umap", ncol=2)#, min.cutoff = 'q5')





#### Find degs between genotype -------------------------------------------------
ast = readRDS('~/novaseq-rats-analysis/singlets-data-after-doubletfinder/astrocytes/astrocytes.lastvar.Rds')

DefaultAssay(ast)

# findMarkers between conditions
# add new column to separate cells by markers and genotype
ast$celltype = Idents(ast)
ast$celltype.genotype <- paste0(ast$celltype,'_', ast$genotype)
Idents(ast) <- ast$celltype.genotype
DimPlot(ast, label=T)

# find markers beetwen genotype by SCT
DefaultAssay(ast) = 'SCT'
ast = PrepSCTFindMarkers(ast)
unique(Idents(ast))
degs.btwn.genotype = data.frame()
for(i in unique(ast$celltype)){
  print(i)
  degs = FindMarkers(ast, ident.1 = paste0(i,'_','agrsv'), 
                          ident.2 = paste0(i,'_','tamed'), 
                          assay = "SCT", test.use = 'MAST')
  degs$celltype = i
  degs$gene = rownames(degs)
  degs.btwn.genotype = rbind(degs.btwn.genotype, degs)}

degs.btwn.genotype = degs.btwn.genotype %>% dplyr::filter(p_val_adj < 0.05)
head(sort(degs.btwn.genotype[degs.btwn.genotype$celltype == 'Efna5+',]$avg_log2FC, decreasing = T))

degs.btwn.genotype$pct.diff = degs.btwn.genotype$pct.1 - degs.btwn.genotype$pct.2

##pct.diff
#top5_res0.7_agrsv_pct.diff = res0.7_mkrs %>% group_by(cluster) %>% top_n(n = 5, wt=agrsv_pct.diff) %>% pull(gene)
#table(top5_res0.7_agrsv_pct.diff$cluster)
#top5_res0.7_tamed_pct.diff = res0.7_mkrs %>% group_by(cluster) %>% top_n(n = 5, wt=tamed_pct.diff) %>% pull(gene)
#table(top5_res0.1_tamed_pct.diff$cluster)

##logfc
top.degs.btwn.genotype.logfc = degs.btwn.genotype %>% group_by(celltype) %>% top_n(n = 5, wt=avg_log2FC)# %>% pull(gene)
table(top.degs.btwn.genotype.logfc$celltype)

# find markers between genotype by RNA 
DefaultAssay(ast) = "RNA"
unique(Idents(ast))
degs.btwn.genotype.rna = data.frame()
for(i in unique(ast$celltype)){
  print(i)
  degs = FindMarkers(ast, ident.1 = paste0(i,'_','agrsv'), 
                     ident.2 = paste0(i,'_','tamed'), 
                     test.use = 'MAST')
  degs$celltype = i
  degs$gene = rownames(degs)
  degs.btwn.genotype.rna = rbind(degs.btwn.genotype.rna, degs)}

degs.btwn.genotype.rna = degs.btwn.genotype.rna %>% dplyr::filter(p_val_adj < 0.05)

##logfc
degs.btwn.genotype.rna.logfc = degs.btwn.genotype.rna %>% group_by(celltype) %>% top_n(n = 5, wt=avg_log2FC)# %>% pull(gene)
table(degs.btwn.genotype.rna.logfc$celltype)


# find intersection (Upset plot) -----------------------------------------------
#SCT
# create list of DEGs
degs.list = list()
for(i in unique(degs.btwn.genotype$celltype)){
  v = degs.btwn.genotype[degs.btwn.genotype$celltype == i,]$gene
  degs.list[[i]] = v
}

library(UpSetR)
upset1 =  upset(fromList(degs.list),nsets = length(degs.list), order.by = 'freq')

upset1$New_data
x1 <- unlist(degs.list, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
x1
degs.shared = x1[ rowSums(upset1$New_data) ==  length(degs.list)]

#RNA
degs.list.rna = list()
for(i in unique(degs.btwn.genotype.rna$celltype)){
  v = degs.btwn.genotype.rna[degs.btwn.genotype.rna$celltype == i,]$gene
  degs.list.rna[[i]] = v
}

library(UpSetR)
upset2 =  upset(fromList(degs.list.rna),nsets = length(degs.list.rna), order.by = 'freq')

upset2$New_data
x2 <- unlist(degs.list.rna, use.names = FALSE)
x2 <- x2[ !duplicated(x2) ]
x2
degs.shared = x2[ rowSums(upset2$New_data) ==  length(degs.list.rna)]


grid.arrange(upset1, upset2,ncol=2)
# subset of shared DEGs between all cell types
degs.btwn.genotype.shared = degs.btwn.genotype[degs.btwn.genotype$gene %in% degs.shared,]
top.degs.btwn.genotype.shared = degs.btwn.genotype.shared %>% group_by(celltype) %>% top_n(n = 5, wt=p_val_adj)

sort(table(top.degs.btwn.genotype.shared$gene))

# Vizualization on violin plot
DefaultAssay(ast) = 'RNA'
ast = NormalizeData(ast)

FeaturePlot(ast, features = 'Fos', split.by = 'genotype', cols = c("grey", "red"),min.cutoff = 'q20')

plots = VlnPlot(ast, features = c('Dnah11', 'Apoe', 'Oxt'), split.by = "genotype", group.by = "celltype", 
        pt.size = 0.3, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

degs.btwn.genotype[degs.btwn.genotype$gene == "Gfap",]



plots[[1]]


# Wilcoxon runk sum test between share of cell agrsv vs tamed -----------------------------------------
Idents(ast) = ast$celltype
# two ways how to count share of cells (near and by dplyr)
table(ast@meta.data[ast@meta.data$celltype == "Fos+/Jun+",]$smpl) / table(ast@meta.data$smpl)
# create list of dataframes with proportion of cells for each celltype
props.per.celltype = list()
for(i in unique(ast$celltype)){
  df <- ast@meta.data %>%
    dplyr::select(orig.ident, celltype, genotype) %>%
    group_by(orig.ident) %>%
    dplyr::mutate(total_cells = n()) %>%
    group_by(orig.ident, celltype, genotype) %>%
    dplyr::mutate(total_celltype = n()) %>%
    dplyr::mutate(prop_celltype = total_celltype / total_cells) %>%
    distinct(orig.ident, celltype, genotype, prop_celltype) %>%
    dplyr::filter(celltype == i) %>%
    arrange(genotype)
  props.per.celltype[[i]] = df
  }

props.plots = list()
for(celltype in props.per.celltype){
  p = celltype %>%
    ggplot(aes(x=factor(orig.ident, levels = df$orig.ident),fill=genotype, y=prop_celltype)) +
    geom_bar(stat = "identity",) +
    scale_fill_manual(values=c('red','blue'))+
    xlab("aggresive / tamed") +
    ylab(paste("Proportion", unique(celltype$celltype))) +
    ggtitle(unique(celltype$celltype)) + 
    guides(fill="none") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") + 
    theme_bw()
  props.plots[[unique(celltype$celltype)]] = p
}

grid.arrange(props.plots[[2]], props.plots[[9]], props.plots[[10]], props.plots[[7]], props.plots[[4]],
             props.plots[[5]],props.plots[[13]],props.plots[[8]],props.plots[[1]],props.plots[[6]],
             props.plots[[11]],props.plots[[12]],props.plots[[3]],props.plots[[14]],ncol=4)

sort(sapply(props.per.celltype, function(x) 
  wilcox.test(x$prop_celltype ~ x$genotype,
              alternative="two.sided")$p.value))




head(ast@meta.data)



devtools::install_github("kharchenkolab/sccore", ref="dev")
devtools::install_github('kharchenkolab/cacoa')

#sample.groups: vector with condition labels per sample named with sample ids
#cell.groups: cell type annotation vector named by cell ids
#sample.per.cell: vector with sample labels per cell named with cell ids
#ref.level: id of the condition, corresponding to the reference (i.e. control)
#target.level: id of the condition, corresponding to the target (i.e. case)
library(clusterProfiler)
library(DESeq2)
library(DOSE)
library(EnhancedVolcano) 
library(enrichplot) 
library(fabia) 
library(GOfuncR) 
library(Rgraphviz)
library(sccore)
library(dplyr)
library(cacoa)
setwd("~/novaseq-rats-analysis/singlets-data-after-doubletfinder/astrocytes/")
#
so = readRDS('./astrocytes.lastvar.Rds')

#
sample.groups.df = so@meta.data %>% group_by(orig.ident) %>% filter(row_number()==1)# %>% select(genotype)
sample.groups =  factor(sample.groups.df$genotype)
names(sample.groups) = sample.groups.df$orig.ident
sample.groups
#
so@meta.data$celltype = Idents(so)
cell.groups <- as.character(so@meta.data$celltype)
names(cell.groups) <- rownames(so@meta.data)
cell.groups
#
sample.per.cell <- so@meta.data$orig.ident
names(sample.per.cell) <- rownames(so@meta.data)
sample.per.cell
#
ref.level <- "tamed"
target.level <- "agrsv"

#
cao <- Cacoa$new(so, sample.groups=sample.groups, cell.groups=cell.groups, 
                 sample.per.cell=sample.per.cell, ref.level=ref.level, target.level=target.level, 
                 graph.name='integrated_nn', embedding = so@reductions$umap@cell.embeddings, n.cores = 15)
# set default plot parameters
cao$plot.params <- list(size=0.1, alpha=0.1, font.size=c(2, 3))
cao$plot.theme <- cao$plot.theme + theme(legend.background=element_blank())

# Estimate cluster-based changes
cao$estimateCellLoadings()
cao$estimateExpressionShiftMagnitudes()
#Plot compositional changes
cao$plotCellLoadings(show.pvals=FALSE)
#Plot expression changes
cao$plotExpressionShiftMagnitudes()

# sample_meta is a list or data.frame with metadata per sample
sample_meta = so@meta.data %>% group_by(orig.ident) %>% filter(row_number()==1) %>% dplyr::select(genotype, sex, smpl)
sample_meta = sample_meta[,2:4]
rownames(sample_meta) = sample_meta$smpl
sample_meta
#This function shows warning if metadata significantly affects results:
pvals <- cao$estimateMetadataSeparation(sample.meta=sample_meta)$padjust
## Warning in cao$estimateMetadataSeparation(sample.meta = sample_meta):
## Significant separation by: Sample_Source, Diagnosis, Status
sort(pvals)
# Now we can show sample expression structure colored by one of the significant factors:
cao$plotSampleDistances(
  space="expression.shifts",
  sample.colors=sample_meta$sex, legend.position=c(0, 1), font.size=2
)
#The same can be done for compostional structure of samples:
cao$plotSampleDistances(
  space="coda",
  sample.colors=sample_meta$Sample_Source, legend.position=c(0, 1), font.size=2
)

# Functional interpretation
#Estimate DE and GSEA:
BiocManager::install("org.Rn.eg.db")
library(org.Rn.eg.db)

cao$estimateDEPerCellType(independent.filtering=TRUE, test='DESeq2.Wald', verbose=FALSE)
cao$estimateOntology(type="GSEA", org.db=org.Rn.eg.db::org.Rn.eg.db, verbose=FALSE, n.cores=1)

#A quick way to look how many significant genes there are per cell type is to show a panel of volcano plots:
cao$plotVolcano(xlim=c(-3, 3), ylim=c(0, 3.5), lf.cutoff=1)
# extract DEGs from 
degs.cocoa = data.frame()
for(i in 1:14){
  x = as.data.frame(cao$test.results$de[[i]]$res[cao$test.results$de[[i]]$res$padj <0.05,])
  x$cluster = names(cao$test.results$de[i])
  degs.cocoa = rbind(degs.cocoa, x)
}
rownames(degs.cocoa) = 1:nrow(degs.cocoa)
write.csv(degs.cocoa, './degs.cocoa.csv')

degs.cocoa = data.frame()
for(i in 1:14){
  x = as.data.frame(cao$test.results$de[[i]]$res)
  x$cluster = names(cao$test.results$de[i])
  degs.cocoa = rbind(degs.cocoa, x)
}
rownames(degs.cocoa) = 1:nrow(degs.cocoa)
write.csv(degs.cocoa, './degs.raw.cocoa.csv')

# To better visualize GOs, we have a function that collapses 
# ontologies with similar enriched genes or similar enrichment pattern:
cao$plotOntologyHeatmapCollapsed(
  name="GSEA", genes="up", n=40, clust.method="ward.D", size.range=c(1, 4)
)
#You can also plot GO terms without collapsing them by enrichment patterns. 
#Usually, it takes hundreds of rows, so to make it readable you can focus on the process of interest 
# by filtering them by description:
cao$plotOntologyHeatmap(
  name="GSEA", genes="up", description.regex="extracellular|matrix"
)




##### Cluster-free

# Compositional changes
# Estimate:
cao$estimateCellDensity(method='graph')
cao$estimateDiffCellDensity(type='wilcox')
# PLot
plot_grid(
  cao$plotEmbedding(color.by='cell.groups'), 
  cao$plotDiffCellDensity(legend.position=c(0, 1)), 
  ncol=2
)

# Expression changes
# Estimate:
cao$estimateClusterFreeExpressionShifts(
  gene.selection="expression", min.n.between=3, min.n.within=3
)
# Plot
cao$plotClusterFreeExpressionShifts(legend.position=c(0, 1), font.size=2)

# Cluster-free DE
cao$estimateClusterFreeDE(
  n.top.genes=1000, min.expr.frac=0.01, adjust.pvalues=TRUE, smooth=TRUE, 
  verbose=TRUE
)
#Estimate gene programs:
cao$estimateGenePrograms(method="pam", z.adj=TRUE, smooth=FALSE)
# Plot gene programs:
cao$plotGeneProgramScores(
  legend.position=c(0, 1), plot.na=FALSE, 
  adj.list=theme(legend.key.width=unit(8, "pt"), legend.key.height=unit(12, "pt"))
)

#Plot genes from one program:
plot_grid(plotlist=cao$plotGeneProgramGenes(
  program.id=12, max.genes=9, plot.na=FALSE, legend.position=c(0, 1)
), ncol=3)






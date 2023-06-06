set.seed(1234)

library(monocle3)
#library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(cicero)

# upload genome size
mratbn7.2 = read.table('~/novaseq-rats-analysis/scafe/rat.bn7.2.genome.txt')

# download and unzip annotation
temp <- tempfile()
download.file("http://ftp.ensembl.org/pub/release-109/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.109.gtf.gz", temp)
gene_anno <- rtracklayer::readGFF(temp)
unlink(temp)
# rename some columns to match requirements
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name


# upload cicero astro full -----------------------------------------------------
astro.full.cicero_cds = readRDS('~/novaseq-rats-analysis/scafe/scafe.tool.sc.solo.aggregate.outs/count/integrated.objects/astro.scafe/cicero/scafe.astro.cicero_cds.full.Rds')
astro.full.conns = readRDS("~/novaseq-rats-analysis/scafe/scafe.tool.sc.solo.aggregate.outs/count/integrated.objects/astro.scafe/cicero/astro.full.cons.Rds")
astro.full.ccans = readRDS("~/novaseq-rats-analysis/scafe/scafe.tool.sc.solo.aggregate.outs/count/integrated.objects/astro.scafe/cicero/astro.full.ccans.Rds")

plot_connections(astro.full.conns, "chr6", 105121170-50000, 105124036+50000,
                 gene_model = gene_anno, 
                 coaccess_cutoff = .25, 
                 connection_width = .5, 
                 collapseTranscripts = "longest" )


# cicero_cds AGRSV ------------------------------------------------------------------------
astro.agrsv.cds = readRDS("~/novaseq-rats-analysis/scafe/scafe.tool.sc.solo.aggregate.outs/count/integrated.objects/astro.scafe/cicero/scafe.astro.cds.agrsv.Rds")
# make cicero cds
umap_coords <- reducedDims(astro.agrsv.cds)$UMAP
astro.agrsv.cicero_cds <- make_cicero_cds(astro.agrsv.cds, reduced_coordinates = umap_coords)
# run cicero
astro.agrsv.conns <- run_cicero(astro.agrsv.cicero_cds, mratbn7.2, sample_num = 100)
astro.agrsv.conns$in_tamed = compare_connections(astro.agrsv.conns, astro.tamed.conns)
saveRDS(astro.agrsv.conns, './astro.scafe/cicero/astro.agrsv.conns.Rds')
astro.agrsv.conns = readRDS("~/novaseq-rats-analysis/scafe/scafe.tool.sc.solo.aggregate.outs/count/integrated.objects/astro.scafe/cicero/astro.agrsv.conns.Rds")
# CCAN
astro.agrsv.ccan <- generate_ccans(astro.agrsv.conns)
saveRDS(astro.agrsv.ccan, './astro.scafe/cicero/astro.agrsv.ccan.Rds')

# cicero_cds TAMED ------------------------------------------------------------------------
astro.tamed.cds = readRDS("~/novaseq-rats-analysis/scafe/scafe.tool.sc.solo.aggregate.outs/count/integrated.objects/astro.scafe/cicero/scafe.astro.cds.tamed.Rds")
# make cicero cds
umap_coords <- reducedDims(astro.tamed.cds)$UMAP
astro.tamed.cicero_cds <- make_cicero_cds(astro.tamed.cds, reduced_coordinates = umap_coords)
# run cicero
astro.tamed.conns <- run_cicero(astro.tamed.cicero_cds, mratbn7.2, sample_num = 100)
astro.tamed.conns$in_agrsv <- compare_connections(astro.tamed.conns, astro.agrsv.conns)
saveRDS(astro.tamed.conns, './astro.scafe/cicero/astro.tamed.conns.Rds')
astro.tamed.conns = readRDS("~/novaseq-rats-analysis/scafe/scafe.tool.sc.solo.aggregate.outs/count/integrated.objects/astro.scafe/cicero/astro.tamed.conns.Rds")

# CCAN
astro.tamed.ccan <- generate_ccans(astro.tamed.conns)
saveRDS(astro.tamed.ccan, './astro.scafe/cicero/astro.tamed.ccan.Rds')

# plotting ---------------------------------------------------------------------
gene_anno[gene_anno$symbol %in% 'Dnah11',]
plot_connections(astro.tamed.conns, "chr6", 138839177-1000000, 139155536+1000000,
                 viewpoint = "chr6_138839177_138840177",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0, 
                 connection_width = .5, 
                 collapseTranscripts = "longest")

plot_connections(astro.agrsv.conns, "chr6", 138839177-1000000, 139155536+1000000,
                 viewpoint = "chr6_138839177_138840177",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0, 
                 connection_width = .5, 
                 collapseTranscripts = "longest")

#
plot_connections(astro.agrsv.conns, "chr6", 138839177-50000, 139155536+50000, 
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0,
                 connection_width = .5,
                 comparison_track = astro.tamed.conns,
                 comparison_connection_width = .5,
                 nclude_axis_track = FALSE,
                 nclude_axis_track = F,
                 collapseTranscripts = "longest") 






gene_anno[gene_anno$gene %in% 'ENSRNOG00000065988',]



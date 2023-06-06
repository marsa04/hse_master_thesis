setwd('~/novaseq-rats-analysis/scafe/scafe.tool.sc.solo.aggregate.outs/count/')
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
library(stringr)

# Rat1
sample = 'Rat1'
counts = Read10X('./Rat1/matrix/', unique.features = TRUE)
srat <- CreateSeuratObject(counts=counts, min.cells = 2, min.features = 200, project = sample)
head(srat@assays)
srat[['genotype']] = 'tamed'
srat[['sex']] = 'mixed'
dim(srat)
bc = read.table('~/novaseq-rats-analysis/singlets-data-after-doubletfinder/filtered.barkodes.for.each.sample/Rat1.filtered_barcodes.txt')
sum(colnames(srat) %in% bc$V1)
sum(bc$V1 %in% colnames(srat))
srat <- subset(srat, cells = bc$V1)
dim(srat)
saveRDS(object = srat, file = str_c('./filtered.objects/', sample, ".scafe.count.mtx.filtered.Rds"))

# Rat2
sample = 'Rat2'
counts = Read10X('./Rat2/matrix/', unique.features = TRUE)
srat <- CreateSeuratObject(counts=counts, min.cells = 2, min.features = 200, project = sample)
head(srat@assays)
srat[['genotype']] = 'tamed'
srat[['sex']] = 'mixed'
head(srat)
dim(srat)
bc = read.table('~/novaseq-rats-analysis/singlets-data-after-doubletfinder/filtered.barkodes.for.each.sample/Rat2.filtered_barcodes.txt')
sum(colnames(srat) %in% bc$V1)
sum(bc$V1 %in% colnames(srat))
srat <- subset(srat, cells = bc$V1)
dim(srat)
saveRDS(object = srat, file = str_c('./filtered.objects/', sample, ".scafe.count.mtx.filtered.Rds"))

# Rat3
sample = 'Rat3'
counts = Read10X('./Rat3/matrix/', unique.features = TRUE)
srat <- CreateSeuratObject(counts=counts, min.cells = 2, min.features = 200, project = sample)
head(srat@assays)
srat[['genotype']] = 'agrsv'
srat[['sex']] = 'female'
head(srat)
dim(srat)
bc = read.table('~/novaseq-rats-analysis/singlets-data-after-doubletfinder/filtered.barkodes.for.each.sample/Rat3.filtered_barcodes.txt')
sum(colnames(srat) %in% bc$V1)
sum(bc$V1 %in% colnames(srat))
srat <- subset(srat, cells = bc$V1)
dim(srat)
saveRDS(object = srat, file = str_c('./filtered.objects/', sample, ".scafe.count.mtx.filtered.Rds"))

# Rat4
sample = 'Rat4'
counts = Read10X('./Rat4/matrix/', unique.features = TRUE)
srat <- CreateSeuratObject(counts=counts, min.cells = 2, min.features = 200, project = sample)
head(srat@assays)
srat[['genotype']] = 'agrsv'
srat[['sex']] = 'male'
head(srat)
dim(srat)
bc = read.table('~/novaseq-rats-analysis/singlets-data-after-doubletfinder/filtered.barkodes.for.each.sample/Rat4.filtered_barcodes.txt')
sum(colnames(srat) %in% bc$V1)
sum(bc$V1 %in% colnames(srat))
srat <- subset(srat, cells = bc$V1)
dim(srat)
saveRDS(object = srat, file = str_c('./filtered.objects/', sample, ".scafe.count.mtx.filtered.Rds"))

# Rat5
sample = 'Rat5'
counts = Read10X('./Rat5/matrix/', unique.features = TRUE)
srat <- CreateSeuratObject(counts=counts, min.cells = 2, min.features = 200, project = sample)
head(srat@assays)
srat[['genotype']] = 'agrsv'
srat[['sex']] = 'male'
head(srat)
dim(srat)
bc = read.table('~/novaseq-rats-analysis/singlets-data-after-doubletfinder/filtered.barkodes.for.each.sample/Rat5.filtered_barcodes.txt')
sum(colnames(srat) %in% bc$V1)
sum(bc$V1 %in% colnames(srat))
srat <- subset(srat, cells = bc$V1)
dim(srat)
saveRDS(object = srat, file = str_c('./filtered.objects/', sample, ".scafe.count.mtx.filtered.Rds"))

# Rat8
sample = 'Rat8'
counts = Read10X('./Rat8/matrix/', unique.features = TRUE)
srat <- CreateSeuratObject(counts=counts, min.cells = 2, min.features = 200, project = sample)
head(srat@assays)
srat[['genotype']] = 'tamed'
srat[['sex']] = 'male'
head(srat)
dim(srat)
bc = read.table('~/novaseq-rats-analysis/singlets-data-after-doubletfinder/filtered.barkodes.for.each.sample/Rat8.filtered_barcodes.txt')
sum(colnames(srat) %in% bc$V1)
sum(bc$V1 %in% colnames(srat))
srat <- subset(srat, cells = bc$V1)
dim(srat)
saveRDS(object = srat, file = str_c('./filtered.objects/', sample, ".scafe.count.mtx.filtered.Rds"))

# K2.1
sample = 'K2.1'
counts = Read10X('./K2_1/matrix/', unique.features = TRUE)
srat <- CreateSeuratObject(counts=counts, min.cells = 2, min.features = 200, project = sample)
head(srat@assays)
srat[['genotype']] = 'tamed'
srat[['sex']] = 'mixed'
head(srat)
dim(srat)
bc = read.table('~/novaseq-rats-analysis/singlets-data-after-doubletfinder/filtered.barkodes.for.each.sample/K2_1.filtered_barcodes.txt')
sum(colnames(srat) %in% bc$V1)
sum(bc$V1 %in% colnames(srat))
srat <- subset(srat, cells = bc$V1)
dim(srat)
saveRDS(object = srat, file = str_c('./filtered.objects/', sample, ".scafe.count.mtx.filtered.Rds"))

# K2.2
sample = 'K2.2'
counts = Read10X('./K2_2/matrix/', unique.features = TRUE)
srat <- CreateSeuratObject(counts=counts, min.cells = 2, min.features = 200, project = sample)
head(srat@assays)
srat[['genotype']] = 'agrsv'
srat[['sex']] = 'mixed'
head(srat)
dim(srat)
bc = read.table('~/novaseq-rats-analysis/singlets-data-after-doubletfinder/filtered.barkodes.for.each.sample/K2_2.filtered_barcodes.txt')
sum(colnames(srat) %in% bc$V1)
sum(bc$V1 %in% colnames(srat))
srat <- subset(srat, cells = bc$V1)
dim(srat)
saveRDS(object = srat, file = str_c('./filtered.objects/', sample, ".scafe.count.mtx.filtered.Rds"))

# K2.3
sample = 'K2.3'
counts = Read10X('./K2_3/matrix/', unique.features = TRUE)
srat <- CreateSeuratObject(counts=counts, min.cells = 2, min.features = 200, project = sample)
head(srat@assays)
srat[['genotype']] = 'tamed'
srat[['sex']] = 'mixed'
head(srat)
dim(srat)
bc = read.table('~/novaseq-rats-analysis/singlets-data-after-doubletfinder/filtered.barkodes.for.each.sample/K2_3.filtered_barcodes.txt')
sum(colnames(srat) %in% bc$V1)
sum(bc$V1 %in% colnames(srat))
srat <- subset(srat, cells = bc$V1)
dim(srat)
saveRDS(object = srat, file = str_c('./filtered.objects/', sample, ".scafe.count.mtx.filtered.Rds"))

# K2.4
sample = 'K2.4'
counts = Read10X('./K2_4/matrix/', unique.features = TRUE)
srat <- CreateSeuratObject(counts=counts, min.cells = 2, min.features = 200, project = sample)
head(srat@assays)
srat[['genotype']] = 'agrsv'
srat[['sex']] = 'mixed'
head(srat)
dim(srat)
bc = read.table('~/novaseq-rats-analysis/singlets-data-after-doubletfinder/filtered.barkodes.for.each.sample/K2_4.filtered_barcodes.txt')
sum(colnames(srat) %in% bc$V1)
sum(bc$V1 %in% colnames(srat))
srat <- subset(srat, cells = bc$V1)
dim(srat)
saveRDS(object = srat, file = str_c('./filtered.objects/', sample, ".scafe.count.mtx.filtered.Rds"))

# K2.5
sample = 'K2.5'
counts = Read10X('./K2_5/matrix/', unique.features = TRUE)
srat <- CreateSeuratObject(counts=counts, min.cells = 2, min.features = 200, project = sample)
head(srat@assays)
srat[['genotype']] = 'tamed'
srat[['sex']] = 'mixed'
head(srat)
dim(srat)
bc = read.table('~/novaseq-rats-analysis/singlets-data-after-doubletfinder/filtered.barkodes.for.each.sample/K2_5.filtered_barcodes.txt')
sum(colnames(srat) %in% bc$V1)
sum(bc$V1 %in% colnames(srat))
srat <- subset(srat, cells = bc$V1)
dim(srat)
saveRDS(object = srat, file = str_c('./filtered.objects/', sample, ".scafe.count.mtx.filtered.Rds"))

# K2.5
sample = 'K2.6'
counts = Read10X('./K2_6/matrix/', unique.features = TRUE)
srat <- CreateSeuratObject(counts=counts, min.cells = 2, min.features = 200, project = sample)
head(srat@assays)
srat[['genotype']] = 'agrsv'
srat[['sex']] = 'mixed'
head(srat)
dim(srat)
bc = read.table('~/novaseq-rats-analysis/singlets-data-after-doubletfinder/filtered.barkodes.for.each.sample/K2_6.filtered_barcodes.txt')
sum(colnames(srat) %in% bc$V1)
sum(bc$V1 %in% colnames(srat))
srat <- subset(srat, cells = bc$V1)
dim(srat)
saveRDS(object = srat, file = str_c('./filtered.objects/', sample, ".scafe.count.mtx.filtered.Rds"))

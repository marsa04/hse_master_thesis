library(ChIPseeker)
library(TxDb.Rnorvegicus.UCSC.rn7.refGene)
txdb <- TxDb.Rnorvegicus.UCSC.rn7.refGene
library(clusterProfiler)
library(AnnotationDbi)
library(org.Rn.eg.db)

path = "~/novaseq-rats-analysis/scafe/scafe.tool.sc.solo.aggregate.outs/annotate/agrsv.tamed/bed/agrsv.tamed.CRE.annot.bed.gz"
#peaks = readPeakFile(path)
#peakAnno = annotatePeak(peak=peaks, tssRegion=c(-1000, 1000), TxDb = txdb,  annoDb = "org.Rn.eg.db")

tcre = read.table(path)
library(stringr)

tcre$V1 = str_to_lower(tcre$V1)
mycap <- function(mystr = "") {
  mystr <- tolower(mystr)
  n <- nchar(mystr)
  substr(mystr, n, n) <- toupper(substr(mystr, n, n))
  return(mystr)
}
tcre$V1 = sapply(tcre$V1, mycap)
colnames(tcre) = c('chr', 'start', 'end', 'id', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10', 'V11', 'V12')

tcre.gr = makeGRangesFromDataFrame(tcre,keep.extra.columns = T)
peakHeatmap(tcre.gr, TxDb = txdb)

peakAnno = annotatePeak(peak=tcre.gr, tssRegion=c(-1000, 1000), TxDb = txdb,  annoDb = "org.Rn.eg.db")
saveRDS(peakAnno, '~/novaseq-rats-analysis/scafe/scafe.tool.sc.solo.aggregate.outs/annotate/agrsv.tamed/tCREpeakAnno.Rds')
plotAnnoBar(peakAnno)

# get peaks from proms and distal ----------------------------------------------
# proms
proms.gr = peakAnno@anno[peakAnno@detailGenomicAnnotation$Promoter,]
proms.df = as.data.frame(proms.gr)
proms.df$peak = paste(proms.df$seqnames, proms.df$start, proms.df$end, sep = '-')

#dist
dist.gr = peakAnno@anno[peakAnno@detailGenomicAnnotation$distal_intergenic,]
dist.df = as.data.frame(dist.gr)
dist.df$peak = paste(dist.df$seqnames, dist.df$start, dist.df$end, sep='-')

nrow(dist.df)
sum(startsWith(dist.df$annotation, "Distal"))

#int
int.gr = peakAnno@anno[peakAnno@detailGenomicAnnotation$Intron,]
int.df = as.data.frame(int.gr)
int.df = int.df[startsWith(int.df$annotation, "Intron"),]
int.df$peak = paste(int.df$seqnames, int.df$start, int.df$end, sep='-')


source('func.R')
#рандомно разбить геном на участки 
mock.peaks.all.genome <- GetMockPeaksWholeGenome(txdb = txdb,
                                                 peak.length = median(tcre.gr@ranges@width))
#самплирование последовательностей
mock.peaks.all.genome <- mock.peaks.all.genome[sample(1:nrow(mock.peaks.all.genome), 1000, replace = F)]
#аннтирование последовательностей
mock.peaks.all.genome.anno <- GetPeakAnno(makeGRangesFromDataFrame(mock.peaks.all.genome), 
                                          tssRegion = c(-1000, 1000),
                                          txdb = txdb,
                                          annoDb = "org.Rn.eg.db")

#BarPlotForArticle(anno.list = list(tCRE = peakAnno,
#                                   Genome = mock.peaks.all.genome.anno))

peakAnnoList = list(peakAnno, mock.peaks.all.genome.anno)
names(peakAnnoList) = c('tCRE', 'Genome')
p2 = plotAnnoBar(peakAnnoList)




############### dist to TSS plot ###############
dist_tcre_tss <- peakAnno@anno@elementMetadata@listData$distanceToTSS

# extract peak locations and take first word in string 
peakloc <- peakAnno@anno@elementMetadata@listData$annotation
peakloc <- word(peakloc, 1)

df <- as.data.frame(peakAnno@anno@ranges)
df$seqnames <- as.data.frame(peakAnno@anno@seqnames)
df$dist <- dist_tcre_tss
df$peakloc <- peakloc
df$peakloc <- recode_factor(df$peakloc, 
                            Downstream = "Promoter")
# draw a histogram for dist from TSS
limit <- 2 * sd(df$dist)
p1 <- ggplot(df, aes(x=dist)) + 
  geom_histogram(bins = 50) + 
  scale_x_continuous(limits = c(-limit, limit), breaks = c(-300000,0,300000), labels = c('-3e+05', '0', '3e+05')) +
#  xlim(-limit, limit) +
  theme_bw() +
  xlab("Distance to Nearest Protein-Coding TSS") +
  ylab("Number of tCRE") +
  theme(text = element_text(size = 14))
p1

round(limit)

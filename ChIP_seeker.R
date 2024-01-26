rm(list = ls())
setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing")

# Load Data ---------------------------------------------------------------
library(ChIPseeker)

loadFile_peakCall <- function(){
 file <- choose.files()
 file <- readPeakFile(file, header = FALSE)
 return(file)
}

# iPSC_H3K27ac <- loadFile_peakCall()
# iPSC_H3K9me3 <- loadFile_peakCall()

iPSC20b_H3K27ac <- loadFile_peakCall()
iPSC20b_H3K9me3 <- loadFile_peakCall()

# iPSC_CM_H3K27ac <- loadFile_peakCall()
# iPSC_CM_H3K9me3 <- loadFile_peakCall()

CM_H3K27ac <- loadFile_peakCall()

LV_H3K27ac <- loadFile_peakCall()
LV_H3K9me3 <- loadFile_peakCall()

# Protocol 1: ChIPseeker and epigenomic dataset prep ------------------------------------------

prepGRangeObj <- function(seek_object){
 seek_object$Peaks = seek_object$V4
 seek_object$level = seek_object$V5
 seek_object$V4 = seek_object$V5 = NULL
 return(seek_object)
}

# iPSC_H3K27ac <- prepGRangeObj(iPSC_H3K27ac)
# iPSC_H3K9me3 <- prepGRangeObj(iPSC_H3K9me3)

iPSC20b_H3K27ac <- prepGRangeObj(iPSC20b_H3K27ac)
iPSC20b_H3K9me3 <- prepGRangeObj(iPSC20b_H3K9me3)

# iPSC_CM_H3K27ac <- prepGRangeObj(iPSC_CM_H3K27ac)
# iPSC_CM_H3K9me3 <- prepGRangeObj(iPSC_CM_H3K9me3)

CM_H3K27ac <- prepGRangeObj(CM_H3K27ac)

LV_H3K27ac <- prepGRangeObj(LV_H3K27ac)
LV_H3K9me3 <- prepGRangeObj(LV_H3K9me3)

# Protocol 2: Annotation of Epigenomic datasets----------------------------------------------------------
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")
TxDb_hg38 = TxDb.Hsapiens.UCSC.hg38.knownGene

Epi_list <- GRangesList(iPSC_H3K27ac = iPSC_H3K27ac, iPSC_H3K9me3 = iPSC_H3K9me3,
                        iPSC20b_H3K27ac = iPSC20b_H3K27ac, iPSC20b_H3K9me3 = iPSC20b_H3K9me3,
                        iPSC_CM_H3K27ac = iPSC_CM_H3K27ac, iPSC_CM_H3K9me3 = iPSC_CM_H3K9me3,
                        CM_H3K27ac = CM_H3K27ac,
                        LV_H3K27ac = LV_H3K27ac, LV_H3K9me3 = LV_H3K9me3)

peakAnnoList <- lapply(Epi_list, annotatePeak, tssRegion = c(-2000,2000), 
                       TxDb = TxDb_hg38)
peakAnnoList_genes = lapply(peakAnnoList, function(i)as.data.frame(i)$geneId)

#can directly extract the GeneID column
# iPSC_H3K27ac_anno_genes <- as.data.frame(peakAnnoList$iPSC_H3K27ac)$geneId
# iPSC_H3K9me3_anno_genes <- as.data.frame(peakAnnoList$iPSC_H3K9me3)$geneId

iPSC20b_H3K27ac_anno_genes <- as.data.frame(peakAnnoList$iPSC20b_H3K27ac)$geneId
iPSC20b_H3K9me3_anno_genes <- as.data.frame(peakAnnoList$iPSC20b_H3K9me3)$geneId
# 
# iPSC_CM_H3K27ac_anno_genes <- as.data.frame(peakAnnoList$iPSC_CM_H3K27ac)$geneId
# iPSC_CM_H3K9me3_anno_genes <- as.data.frame(peakAnnoList$iPSC_CM_H3K9me3)$geneId

CM_H3K27ac_anno_genes <- as.data.frame(peakAnnoList$CM_H3K27ac)$geneId

LV_H3K9me3_anno_genes <- as.data.frame(peakAnnoList$LV_H3K9me3)$geneId
LV_H3K27ac_anno_genes <- as.data.frame(peakAnnoList$LV_H3K27ac)$geneId

# Protocol 4: Visualiztion of annotated results ----------------------------------------------
library(ggVennDiagram)
library(ggupset)
library(clusterProfiler)
library(DOSE)

plotAnnoBar(peakAnnoList, main = "Genomic Feature Distribution")

# iPSC_H3K27ac_anno_enrichGO <- enrichGO(gene = iPSC_H3K27ac_anno_genes, OrgDb = "org.Hs.eg.db",
#                                   ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10,
#                                   maxGSSize = 200)
# iPSC_H3K9me3_anno_enrichGO <- enrichGO(gene = iPSC_H3K9me3_anno_genes, OrgDb = "org.Hs.eg.db",
#                                   ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10,
#                                   maxGSSize = 200)

iPSC20b_H3K27ac_anno_enrichGO <- enrichGO(gene = iPSC20b_H3K27ac_anno_genes, 
                                          OrgDb = "org.Hs.eg.db",
                                          ont = "ALL", pvalueCutoff = 0.05, 
                                          minGSSize = 10, maxGSSize = 200)
iPSC20b_H3K9me3_anno_enrichGO <- enrichGO(gene = iPSC20b_H3K9me3_anno_genes, 
                                          OrgDb = "org.Hs.eg.db",
                                          ont = "ALL", pvalueCutoff = 0.05, 
                                          minGSSize = 10, maxGSSize = 200)

# iPSC_CM_H3K27ac_anno_enrichGO <- enrichGO(gene = iPSC_CM_H3K27ac_anno_genes, OrgDb = "org.Hs.eg.db",
#                                    ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10,
#                                    maxGSSize = 200)
# iPSC_CM_H3K9me3_anno_enrichGO <- enrichGO(gene = iPSC_CM_H3K9me3_anno_genes, OrgDb = "org.Hs.eg.db",
#                                    ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10,
#                                    maxGSSize = 200)

CM_H3K27ac_anno_enrichGO <- enrichGO(gene = CM_H3K27ac_anno_genes, OrgDb = "org.Hs.eg.db",
                                     ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10,
                                     maxGSSize = 200)

LV_H3K27ac_anno_enrichGO <- enrichGO(gene = LV_H3K27ac_anno_genes, OrgDb = "org.Hs.eg.db",
                                     ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10,
                                     maxGSSize = 200)
LV_H3K9me3_anno_enrichGO <- enrichGO(gene = LV_H3K9me3_anno_genes, OrgDb = "org.Hs.eg.db",
                                     ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10,
                                     maxGSSize = 200)

#dotplot(Epi_list_enrichGO, x = "p.adjust", title = "GO Enrichment Analysis Epigenetic List",
#        showCategory = 14)

# dotplot(iPSC_H3K27ac_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for iPSC H3K27ac", showCategory = 18)
# dotplot(iPSC_H3K9me3_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for iPSC H3K9me3", showCategory = 18)

dotplot(iPSC20b_H3K27ac_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for iPSC 20b H3K27ac", showCategory = 18)
dotplot(iPSC20b_H3K9me3_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for iPSC 20b H3K9me3", showCategory = 18)
# 
# dotplot(iPSC_CM_H3K27ac_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for iPSC CM H3K27ac", showCategory = 18)
# dotplot(iPSC_CM_H3K9me3_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for iPSC CM H3K9me3", showCategory = 18)

dotplot(CM_H3K27ac_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for CM H3K27ac", showCategory = 18)

dotplot(LV_H3K27ac_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for LV H3K27ac", showCategory = 18)
dotplot(LV_H3K9me3_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for LV H3K9me3", showCategory = 18)

#####Metaplots#####

TSS = getBioRegion(TxDb=TxDb_hg38, upstream=3000, downstream=3000, by = "gene", 
                   type = "start_site")
Epi_list_tagMatrix = lapply(Epi_list, getTagMatrix, windows = TSS)

#might have to plot these separately. in addition to together
plotAvgProf(Epi_list_tagMatrix, xlim=c(-3000, 3000), ylab = "Count Frequency")
plotPeakProf(Epi_list_tagMatrix, facet = "none", conf = 0.95)



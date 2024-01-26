rm(list = ls())
setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing/Practice with Data")

# Load Data ---------------------------------------------------------------
library(ChIPseeker)

loadFile_peakCall <- function(){
 file <- file.choose()
 file <- readPeakFile(file, header = FALSE)
 return(file)
}

iPSC20b_H3K27ac <- loadFile_peakCall()
iPSC20b_H3K9me3 <- loadFile_peakCall()

JohnD0_H3K27ac <- loadFile_peakCall()
JohnD15_H3K27ac <- loadFile_peakCall()
JohnD30_H3K27ac <- loadFile_peakCall()

LV_H3K27ac <- loadFile_peakCall()
LV_H3K9me3 <- loadFile_peakCall()

# Protocol 1: ChIPseeker and epigenomic dataset prep ------------------------------------------

prepGRangeObj <- function(seek_object){
 seek_object$Peaks = seek_object$V4
 seek_object$level = seek_object$V5
 seek_object$V4 = seek_object$V5 = NULL
 return(seek_object)
}

iPSC20b_H3K27ac <- prepGRangeObj(iPSC20b_H3K27ac)
iPSC20b_H3K9me3 <- prepGRangeObj(iPSC20b_H3K9me3)

JohnD0_H3K27ac <- prepGRangeObj(JohnD0_H3K27ac)
JohnD15_H3K27ac <- prepGRangeObj(JohnD15_H3K27ac)
JohnD30_H3K27ac <- prepGRangeObj(JohnD30_H3K27ac)

LV_H3K27ac <- prepGRangeObj(LV_H3K27ac)
LV_H3K9me3 <- prepGRangeObj(LV_H3K9me3)

# Protocol 2: Annotation of Epigenomic datasets----------------------------------------------------------
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")
TxDb_hg38 = TxDb.Hsapiens.UCSC.hg38.knownGene

Epi_list <- GRangesList(D0K27ac = JohnD0_H3K27ac, D15K27ac = JohnD15_H3K27ac, 
                        D30K27ac = JohnD30_H3K27ac, iPSC_K27ac = iPSC20b_H3K27ac, 
                        iPSC_K9me3 = iPSC20b_H3K9me3, LV_K27ac = LV_H3K27ac,
                        LV_K9me3 = LV_H3K9me3)


peakAnnoList <- lapply(Epi_list, annotatePeak, tssRegion = c(-2000,0), TxDb = TxDb_hg38)
peakAnnoList_genes = lapply(peakAnnoList, function(i)as.data.frame(i)$geneId)

#can directly extract the GeneID column
D0K27ac_anno_genes <- as.data.frame(peakAnnoList$D0K27ac)$geneId
D15K27ac_anno_genes <- as.data.frame(peakAnnoList$D15K27ac)$geneId
D30K27ac_anno_genes <- as.data.frame(peakAnnoList$D30K27ac)$geneId

iPSC_K9me3_anno_genes <- as.data.frame(peakAnnoList$iPSC_K9me3)$geneId
iPSC_K27ac_anno_genes <- as.data.frame(peakAnnoList$iPSC_K27ac)$geneId

LV_K9me3_anno_genes <- as.data.frame(peakAnnoList$LV_K9me3)$geneId
LV_K27ac_anno_genes <- as.data.frame(peakAnnoList$LV_K27ac)$geneId

# Protocol 4: Visualiztion of annotated results ----------------------------------------------
library(ggVennDiagram)
library(ggupset)
library(clusterProfiler)
library(DOSE)

plotAnnoBar(peakAnnoList)

iPSC_K27ac_anno_enrichGO <- enrichGO(gene = iPSC_K27ac_anno_genes, OrgDb = "org.Hs.eg.db",
                                     ont = "ALL", pvalueCutoff = 0.01, minGSSize = 10,
                                     maxGSSize = 200)

iPSC_K9me3_anno_enrichGO <- enrichGO(gene = iPSC_K9me3_anno_genes, OrgDb = "org.Hs.eg.db",
                                     ont = "ALL", pvalueCutoff = 0.01, minGSSize = 10,
                                     maxGSSize = 200)

D0K27ac_anno_enrichGO <- enrichGO(gene = D0K27ac_anno_genes, OrgDb = "org.Hs.eg.db",
                                  ont = "ALL", pvalueCutoff = 0.01, minGSSize = 10,
                                  maxGSSize = 200)

D15K27ac_anno_enrichGO <- enrichGO(gene = D15K27ac_anno_genes, OrgDb = "org.Hs.eg.db",
                                   ont = "ALL", pvalueCutoff = 0.01, minGSSize = 10,
                                   maxGSSize = 200)

D30K27ac_anno_enrichGO <- enrichGO(gene = D30K27ac_anno_genes, OrgDb = "org.Hs.eg.db",
                                   ont = "ALL", pvalueCutoff = 0.01, minGSSize = 10,
                                   maxGSSize = 200)

LV_K27ac_anno_enrichGO <- enrichGO(gene = LV_K27ac_anno_genes, OrgDb = "org.Hs.eg.db",
                                   ont = "ALL", pvalueCutoff = 0.01, minGSSize = 10,
                                   maxGSSize = 200)

LV_K9me3_anno_enrichGO <- enrichGO(gene = LV_K9me3_anno_genes, OrgDb = "org.Hs.eg.db",
                                   ont = "ALL", pvalueCutoff = 0.01, minGSSize = 10,
                                   maxGSSize = 200)

#Epi_list_enrichGO = compareCluster(geneCluster = c(D0K27ac_anno_genes,D15K27ac_anno_genes,
#                                  D30K27ac_anno_genes), fun = "enrichGO",
#                                  pvalueCutoff = 0.01, minGSSize = 10, maxGSSize = 200,
#                                  OrgDb = "org.Hs.eg.db")

#Epi_list_enrichGO = compareCluster(geneCluster = peakAnnoList_genes, fun = "enrichGO",
#                                   pvalueCutoff = 0.01, minGSSize = 10, maxGSSize = 200,
#                                   OrgDb = "org.Hs.eg.db")


#dotplot(Epi_list_enrichGO, x = "p.adjust", title = "GO Enrichment Analysis Epigenetic List",
#        showCategory = 14)

dotplot(iPSC_K27ac_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for iPSC H3K27ac", showCategory = 18)
dotplot(iPSC_K9me3_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for iPSC H3K9me3", showCategory = 18)

dotplot(D0K27ac_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for Day 0 H3K27ac", showCategory = 18)
dotplot(D15K27ac_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis Day 15  H3K27ac", showCategory = 18)
dotplot(D30K27ac_anno_enrichGO,  x = 'p.adjust', title = "GO Enrichment Analysis for Day 30 H3K27ac", showCategory = 18)

dotplot(LV_K27ac_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for LV H3K27ac", showCategory = 18)
dotplot(LV_K9me3_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for LV H3K9me3", showCategory = 18)

barplot(D0K27ac_anno_enrichGO, showCategory = 20) 

#####Metaplots#####

TSS = getBioRegion(TxDb=TxDb_hg38, upstream=3000, downstream=3000, by = "gene", 
                   type = "start_site")
Epi_list_tagMatrix = lapply(Epi_list, getTagMatrix, windows = TSS)

#might have to plot these separately. in addition to together
plotAvgProf(Epi_list_tagMatrix, xlim=c(-3000, 3000), ylab = "Count Frequency")
plotPeakProf(Epi_list_tagMatrix, facet = "none", conf = 0.95)



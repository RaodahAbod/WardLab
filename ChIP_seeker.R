rm(list = ls())
setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing")

# Load Data ---------------------------------------------------------------
library(ChIPseeker)

loadFile_peakCall <- function(){
 file <- choose.files()
 file <- readPeakFile(file, header = FALSE)
 return(file)
}

iPSC_H3K27ac <- loadFile_peakCall()

iPSC_CM_H3K27ac <- loadFile_peakCall()

CM_H3K27ac <- loadFile_peakCall()

LVQ3_H3K27ac <- loadFile_peakCall()

LVQ4_H3K27ac <- loadFile_peakCall()

# Protocol 1: ChIPseeker and epigenomic dataset prep ------------------------------------------

prepGRangeObj <- function(seek_object){
 seek_object$Peaks = seek_object$V4
 seek_object$level = seek_object$V5
 seek_object$V4 = seek_object$V5 = NULL
 return(seek_object)
}

iPSC_H3K27ac <- prepGRangeObj(iPSC_H3K27ac)

iPSC_CM_H3K27ac <- prepGRangeObj(iPSC_CM_H3K27ac)

CM_H3K27ac <- prepGRangeObj(CM_H3K27ac)

LVQ3_H3K27ac <- prepGRangeObj(LVQ3_H3K27ac)

LVQ4_H3K27ac <- prepGRangeObj(LVQ4_H3K27ac)

# Protocol 2: Annotation of Epigenomic datasets----------------------------------------------------------
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")
TxDb_hg38 = TxDb.Hsapiens.UCSC.hg38.knownGene

Epi_list <- GRangesList(iPSC_H3K27ac = iPSC_H3K27ac, iPSC_CM_H3K27ac = iPSC_CM_H3K27ac,
                        CM_H3K27ac = CM_H3K27ac, LVQ4_H3K27ac = LVQ4_H3K27ac)

peakAnnoList <- lapply(Epi_list, annotatePeak, tssRegion = c(-2000,2000), 
                       TxDb = TxDb_hg38)
peakAnnoList_genes = lapply(peakAnnoList, function(i)as.data.frame(i)$geneId)
peakAnnoList_genes2 <- peakAnnoList_genes[2:4]

#can directly extract the GeneID column
iPSC_H3K27ac_anno_genes <- as.data.frame(peakAnnoList$iPSC_H3K27ac)$geneId

iPSC_CM_H3K27ac_anno_genes <- as.data.frame(peakAnnoList$iPSC_CM_H3K27ac)$geneId

CM_H3K27ac_anno_genes <- as.data.frame(peakAnnoList$CM_H3K27ac)$geneId

LVQ3_H3K27ac_anno_genes <- as.data.frame(peakAnnoList$LVQ3_H3K27ac)$geneId

LVQ4_H3K27ac_anno_genes <- as.data.frame(peakAnnoList$LVQ4_H3K27ac)$geneId

# Protocol 4: Visualiztion of annotated results ----------------------------------------------
library(ggVennDiagram)
library(ggupset)
library(clusterProfiler)
library(DOSE)

plotAnnoBar(peakAnnoList, main = "Genomic Feature Distribution")

iPSC_H3K27ac_anno_enrichGO <- enrichGO(gene = iPSC_H3K27ac_anno_genes, OrgDb = "org.Hs.eg.db",
                                       ont = "ALL", pvalueCutoff = 0.01, minGSSize = 5,
                                       maxGSSize = 250)

iPSC_CM_H3K27ac_anno_enrichGO <- enrichGO(gene = iPSC_CM_H3K27ac_anno_genes, OrgDb = "org.Hs.eg.db",
                                          ont = "ALL", pvalueCutoff = 0.01, minGSSize = 5,
                                          maxGSSize = 250)

CM_H3K27ac_anno_enrichGO <- enrichGO(gene = CM_H3K27ac_anno_genes, OrgDb = "org.Hs.eg.db",
                                     ont = "ALL", pvalueCutoff = 0.01, minGSSize = 5,
                                     maxGSSize = 250)

LVQ3_H3K27ac_anno_enrichGO <- enrichGO(gene = LVQ3_H3K27ac_anno_genes, OrgDb = "org.Hs.eg.db",
                                       ont = "ALL", pvalueCutoff = 0.01, minGSSize = 5,
                                       maxGSSize = 250)

LVQ4_H3K27ac_anno_enrichGO <- enrichGO(gene = LVQ4_H3K27ac_anno_genes, OrgDb = "org.Hs.eg.db",
                                       ont = "ALL", pvalueCutoff = 0.01, minGSSize = 5,
                                       maxGSSize = 250)

Epi_list_enrichGO <- compareCluster(geneCluster = peakAnnoList_genes, fun = "enrichGO", 
                                    pvalueCutoff = 0.01, OrgDb = "org.Hs.eg.db", minGSSize = 5,
                                    maxGSSize = 250, ont = "ALL")

Epi_list_enrichGO2 <- compareCluster(geneCluster = peakAnnoList_genes2, fun = "enrichGO", 
                                     pvalueCutoff = 0.01, OrgDb = "org.Hs.eg.db", minGSSize = 5,
                                     maxGSSize = 250, ont = "ALL")


dotplot(Epi_list_enrichGO, size = "count", title = "GO Enrichment Analysis Epigenetic List")
dotplot(Epi_list_enrichGO2, size = "count", title = "GO Enrichment Analysis Epigenetic List no iPSC")


dotplot(iPSC_H3K27ac_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for iPSC H3K27ac", showCategory = 18)

dotplot(iPSC_CM_H3K27ac_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for iPSC CM H3K27ac", showCategory = 18)

dotplot(CM_H3K27ac_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for CM H3K27ac", showCategory = 18)

dotplot(LVQ3_H3K27ac_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for LV Q3 H3K27ac", showCategory = 18)

dotplot(LVQ4_H3K27ac_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for LV Q4 H3K27ac", showCategory = 18)

#####Metaplots#####

TSS = getBioRegion(TxDb=TxDb_hg38, upstream=3000, downstream=3000, by = "gene", 
                   type = "start_site")
Epi_list_tagMatrix = lapply(Epi_list, getTagMatrix, windows = TSS)

#might have to plot these separately. in addition to together
plotAvgProf(Epi_list_tagMatrix, xlim=c(-3000, 3000), ylab = "Count Frequency")
plotPeakProf(Epi_list_tagMatrix, facet = "none", conf = 0.95)



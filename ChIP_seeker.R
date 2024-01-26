rm(list = ls())
setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing")

# Load Data ---------------------------------------------------------------
library(ChIPseeker)

loadFile_peakCall <- function(){
 file <- choose.files()
 file <- readPeakFile(file, header = FALSE)
 return(file)
}

DataSetName <- readline(prompt = "Enter the name of working data set: ")

Q1 <- loadFile_peakCall()
Q2 <- loadFile_peakCall()
Q3 <- loadFile_peakCall()
Q4 <- loadFile_peakCall()

# Protocol 1: ChIPseeker and epigenomic dataset prep ------------------------------------------

prepGRangeObj <- function(seek_object){
 seek_object$Peaks = seek_object$V4
 seek_object$level = seek_object$V5
 seek_object$V4 = seek_object$V5 = NULL
 return(seek_object)
}

Q1 <- prepGRangeObj(Q1)
Q2 <- prepGRangeObj(Q2)
Q3 <- prepGRangeObj(Q3)
Q4 <- prepGRangeObj(Q4)

# Protocol 2: Annotation of Epigenomic datasets----------------------------------------------------------
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")
TxDb_hg38 = TxDb.Hsapiens.UCSC.hg38.knownGene

Epi_list <- GRangesList(Q1 = Q1, Q2 = Q2, Q3 = Q3, Q4 = Q4) 

peakAnnoList <- lapply(Epi_list, annotatePeak, tssRegion = c(-2000,2000), 
                       TxDb = TxDb_hg38)
peakAnnoList_genes = lapply(peakAnnoList, function(i)as.data.frame(i)$geneId)

Q1_anno_genes <- as.data.frame(peakAnnoList$Q1)$geneId
Q2_anno_genes <- as.data.frame(peakAnnoList$Q2)$geneId
Q3_anno_genes <- as.data.frame(peakAnnoList$Q3)$geneId
Q4_anno_genes <- as.data.frame(peakAnnoList$Q4)$geneId

# Protocol 4: Visualiztion of annotated results ----------------------------------------------
library(ggVennDiagram)
library(ggupset)
library(clusterProfiler)
library(DOSE)

plotAnnoBar(peakAnnoList, main = "Genomic Feature Distribution")

Q1_anno_enrichGO <- enrichGO(gene = Q1_anno_genes, OrgDb = "org.Hs.eg.db", ont = "ALL", 
                             pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 200)
Q2_anno_enrichGO <- enrichGO(gene = Q2_anno_genes, OrgDb = "org.Hs.eg.db", ont = "ALL", 
                             pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 200)
Q3_anno_enrichGO <- enrichGO(gene = Q3_anno_genes, OrgDb = "org.Hs.eg.db", ont = "ALL", 
                             pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 200)
Q4_anno_enrichGO <- enrichGO(gene = Q4_anno_genes, OrgDb = "org.Hs.eg.db", ont = "ALL", 
                             pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 200)


dotplot(Q1_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for Quartile 1", 
        showCategory = 18)
dotplot(Q2_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for Quartile 2", 
        showCategory = 18)
dotplot(Q3_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for Quartile 3", 
        showCategory = 18)
dotplot(Q4_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for Quartile 4", 
        showCategory = 18)


#####Metaplots#####

TSS = getBioRegion(TxDb=TxDb_hg38, upstream=3000, downstream=3000, by = "gene", 
                   type = "start_site")
Epi_list_tagMatrix = lapply(Epi_list, getTagMatrix, windows = TSS)

#might have to plot these separately. in addition to together
plotAvgProf(Epi_list_tagMatrix, xlim=c(-3000, 3000), ylab = "Count Frequency")


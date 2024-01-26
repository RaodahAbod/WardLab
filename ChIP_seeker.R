rm(list = ls())
setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing")

# Load Data ---------------------------------------------------------------
library(ChIPseeker)

loadFile_peakCall <- function(){
 file <- choose.files()
 file <- readPeakFile(file, header = FALSE)
 return(file)
}

road_or_Encode <- readline(prompt = "Is this Roadmap or ENCODE?: ")

topFourth <- loadFile_peakCall()
topThird <- loadFile_peakCall()
topHalf <- loadFile_peakCall()
topTwoThirds <- loadFile_peakCall()
original <- loadFile_peakCall()

# Protocol 1: ChIPseeker and epigenomic dataset prep ------------------------------------------

prepGRangeObj <- function(seek_object){
 seek_object$Peaks = seek_object$V4
 seek_object$level = seek_object$V5
 seek_object$V4 = seek_object$V5 = NULL
 return(seek_object)
}

topFourth <- prepGRangeObj(topFourth)
topThird <- prepGRangeObj(topThird)
topHalf <- prepGRangeObj(topHalf)
topTwoThirds <- prepGRangeObj(topTwoThirds)
original <- prepGRangeObj(original)

# Protocol 2: Annotation of Epigenomic datasets----------------------------------------------------------
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")
TxDb_hg38 = TxDb.Hsapiens.UCSC.hg38.knownGene

Epi_list <- GRangesList(topFourth = topFourth, topThird = topThird, topHalf = topHalf, 
                        topTwoThirds = topTwoThirds, original = original)

peakAnnoList <- lapply(Epi_list, annotatePeak, tssRegion = c(-2000,2000), 
                       TxDb = TxDb_hg38)

#can directly extract the GeneID column
topFourth_anno_genes <- as.data.frame(peakAnnoList$topFourth)$geneId

topThird_anno_genes <- as.data.frame(peakAnnoList$topThird)$geneId

topHalf_anno_genes <- as.data.frame(peakAnnoList$topHalf)$geneId

topTwoThirds_anno_genes <- as.data.frame(peakAnnoList$topTwoThirds)$geneId

orig_anno_genes <- as.data.frame(peakAnnoList$original)$geneId

list2env(peakAnnoList,envir = .GlobalEnv) #renee!

# Protocol 4: Visualiztion of annotated results ----------------------------------------------
library(ggVennDiagram)
library(ggupset)
library(clusterProfiler)
library(DOSE)

plotAnnoBar(peakAnnoList, main = "Genomic Feature Distribution")

topFourth_anno_enrichGO <- enrichGO(gene = topFourth_anno_genes, OrgDb = "org.Hs.eg.db",
                                    ont = "ALL", pvalueCutoff = 0.01, minGSSize = 5,
                                    maxGSSize = 500)

topThird_anno_enrichGO <- enrichGO(gene = topThird_anno_genes, OrgDb = "org.Hs.eg.db",
                                   ont = "ALL", pvalueCutoff = 0.01, minGSSize = 5,
                                   maxGSSize = 500)

topHalf_anno_enrichGO <- enrichGO(gene = topHalf_anno_genes, OrgDb = "org.Hs.eg.db",
                                  ont = "ALL", pvalueCutoff = 0.01, minGSSize = 5,
                                  maxGSSize = 500)

topTwoThirds_anno_enrichGO <- enrichGO(gene = topTwoThirds_anno_genes, OrgDb = "org.Hs.eg.db",
                                       ont = "ALL", pvalueCutoff = 0.01, minGSSize = 5,
                                       maxGSSize = 500)

# Epi_list_enrichGO <- compareCluster(geneCluster = peakAnnoList_genes, fun = "enrichGO", 
#                                     pvalueCutoff = 0.01, OrgDb = "org.Hs.eg.db", minGSSize = 5,
#                                     maxGSSize = 500, ont = "ALL")

dotplot(Epi_list_enrichGO, size = "count", title = "GO Enrichment Analysis Epigenetic List")

dotplot(topFourth_anno_enrichGO,split = "ONTOLOGY", x = 'p.adjust', title = "GO Enrichment Analysis for LV-E top 1/4th", showCategory = 7) + facet_grid(ONTOLOGY~., scale="free")

dotplot(topThird_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for LV-E top 1/3rd", showCategory = 18)

dotplot(topHalf_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for LV-E top 1/2", showCategory = 18)

dotplot(topTwoThirds_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis LV-E for top 2/3rd", showCategory = 18)

#####Metaplots#####

TSS = getBioRegion(TxDb=TxDb_hg38, upstream=3000, downstream=3000, by = "gene", 
                   type = "start_site")
Epi_list_tagMatrix = lapply(Epi_list, getTagMatrix, windows = TSS)

#might have to plot these separately. in addition to together
plotAvgProf(Epi_list_tagMatrix, xlim=c(-3000, 3000), ylab = "Count Frequency")
plotPeakProf(Epi_list_tagMatrix, facet = "none", conf = 0.95)




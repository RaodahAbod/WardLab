rm(list = ls())

setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing/Practice with Data")

# Load Data ---------------------------------------------------------------
library(ChIPseeker)
filename = readline()
CTCF_seeker <- readPeakFile("CTCF_H1.test.bed", header = FALSE)
head(CTCF_seeker)

seek_D0_K27ac <- readPeakFile("D0_H3K27ac.bed", header = FALSE)
seek_D15_K27ac <- readPeakFile("D15_H3K27ac.bed", header = FALSE)
seek_D0_K9me3 <- readPeakFile("D0_H3K9me3.bed", header = FALSE)
seek_D15_K9me3 <- readPeakFile("D15_H3K9me3.bed", header = FALSE)


# Protocol 1: ChIPseeker and epigenomic dataset prep ------------------------------------------

CTCF_seeker$CTCF_peaks = CTCF_seeker$V4
CTCF_seeker$level = CTCF_seeker$V5
CTCF_seeker$V4 = CTCF_seeker$V5 = NULL
CTCF_seeker

prepGRangeObj <- function(seek_object){
 seek_object$Peaks = seek_object$V4
 seek_object$level = seek_object$V5
 seek_object$V4 = seek_object$V5 = NULL
 return(seek_object)
}

seek_D0_K27ac <- prepGRangeObj(seek_D0_K27ac)
seek_D15_K27ac <- prepGRangeObj(seek_D15_K27ac)
seek_D0_K9me3 <- prepGRangeObj(seek_D0_K9me3)
seek_D15_K9me3 <- prepGRangeObj(seek_D15_K9me3)

# Protocol 2: Annotation of Epigenomic datasets----------------------------------------------------------
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")

TxDb_hg38 = TxDb.Hsapiens.UCSC.hg38.knownGene

CTCF_annotate <- annotatePeak(CTCF_seeker, TxDb = TxDb_hg38)

D0K27ac_anno <- annotatePeak(seek_D0_K27ac, TxDb = TxDb_hg38, annoDb="org.Hs.eg.db")
D15K27ac_anno <- annotatePeak(seek_D15_K27ac, TxDb = TxDb_hg38,annoDb="org.Hs.eg.db")
D0K9me3_anno <- annotatePeak(seek_D0_K9me3, TxDb = TxDb_hg38, annoDb="org.Hs.eg.db")
D15K9me3_anno <- annotatePeak(seek_D15_K9me3, TxDb = TxDb_hg38, annoDb="org.Hs.eg.db")

# Protocol 4: Visualiztion of annotated results ----------------------------------------------
library(ggVennDiagram)
library(ggupset)

plotAnnoPie(CTCF_annotate)
plotAnnoBar(CTCF_annotate)
upsetplot(CTCF_annotate, show.numbers = "yes")
upsetplot(CTCF_annotate, vennpie=TRUE)

visualizeAnno <- function(anno_object) {
 plotAnnoPie(anno_object)
 plotAnnoBar(anno_object)
 upsetplot(anno_object)
 
}

visualizeAnno(D0K27ac_anno)
visualizeAnno(D15K27ac_anno)
visualizeAnno(D0K9me3_anno)
visualizeAnno(D15K9me3_anno)

#supports multiple entries if using the peakAnnotatedList object

#  Protocol 5: functional analysis of epigenomic datasets ------------------------------------------
library(clusterProfiler)

Epi_list <- GRangesList(D0K27ac = seek_D0_K27ac, D15K27ac = seek_D15_K27ac, 
                        D0K9me3 = seek_D0_K9me3, D15K9me3 = seek_D15_K9me3)


peakAnnoList <- lapply(Epi_list, annotatePeak, tssRegion = c(-2000,0), TxDb = TxDb_hg38)
plotAnnoBar(peakAnnoList)
peakAnnoList_genes = lapply(peakAnnoList, function(i)as.data.frame(i)$geneId)
#can directly extract the GeneID column

D0K27ac_anno_genes <- as.data.frame(peakAnnoList$D0K27ac)$geneId
D15K27ac_anno_genes <- as.data.frame(peakAnnoList$D15K27ac)$geneId

D0K9me3_anno_genes <- as.data.frame(peakAnnoList$D0K9me3)$geneId
#some reason my D0K9me3 is not outputting a gene enrichment and IDK WHYYYYY

D15K9me3_anno_genes <- as.data.frame(peakAnnoList$D15K9me3)$geneId

#####Enrichment Analysis#####

D0K27ac_anno_enrichGO <- enrichGO(gene = D0K27ac_anno_genes, OrgDb = "org.Hs.eg.db",
                                  ont = "ALL", pvalueCutoff = 0.01, minGSSize = 10,
                                  maxGSSize = 200)

# D15K27ac_anno_enrichGO <- enrichGO(gene = D15K27ac_anno_genes, OrgDb = "org.Hs.eg.db",
#                                    ont = "MF", pvalueCutoff = 0.05, minGSSize = 10,
#                                    maxGSSize = 200)

D15K27ac_anno_enrichGO <- enrichGO(gene = D15K27ac_anno_genes, OrgDb = "org.Hs.eg.db",
                                   ont = "ALL", pvalueCutoff = 0.01, minGSSize = 10,
                                   maxGSSize = 200)

# D0K9me3_anno_enrichGO <- enrichGO(gene = D0K9me3_anno_genes, OrgDb = "org.Hs.eg.db",
#                                   ont = "ALL", pvalueCutoff = 0.01, minGSSize = 2,
#                                   maxGSSize = 200)

D15K9me3_anno_enrichGO <- enrichGO(gene = D15K9me3_anno_genes, OrgDb = "org.Hs.eg.db",
                                   ont = "ALL", pvalueCutoff = 0.01, minGSSize = 10,
                                   maxGSSize = 200)


Epi_list_enrichGO = compareCluster(geneCluster = peakAnnoList_genes, fun = "enrichGO",
                                   pvalueCutoff = 0.01, minGSSize = 10, maxGSSize = 200,
                                   OrgDb = "org.Hs.eg.db")

#####Visualization#####

dotplot(D0K27ac_anno_enrichGO, x= 'p.adjust', title = "GO Enrichment Analysis for Day 0 H3K27ac",showCategory = 18)
dotplot(D15K27ac_anno_enrichGO, x= 'p.adjust', title = "GO Enrichment Analysis for Day 15 H3K27ac", showCategory = 18)
#dotplot(D0K9me3_anno_enrichGO, title = "GO Enrichment Analysis for Day 0 H3K9me3", showCategory = 20)
dotplot(D15K9me3_anno_enrichGO, x = 'p.adjust', title = "GO Enrichment Analysis for Day 15 H3K9me3", showCategory = 18)

dotplot(Epi_list_enrichGO, title = "GO Enrichment Analysis Epigenetic List ")


#####Extract Gene Enrichment List#####
library(tidyverse)

enrichedTibble <- function(anno_enrich){
 enrichment_table <- anno_enrich %>% as_tibble()
 # Arrange by adjusted p value
 enrichment_table <- enrichment_table %>% arrange(desc(p.adjust))
 # Factor by P.adjusted value
 enrichment_table$Description <- factor(enrichment_table$Description, 
                                        levels = enrichment_table$Description[order(enrichment_table$p.adjust)])
 return(enrichment_table)
}

D0K27ac_enriched_Tibble <- enrichedTibble(D0K27ac_anno_enrichGO)
D15K27ac_enriched_Tibble <- enrichedTibble(D15K27ac_anno_enrichGO)
D0K9me3_enriched_Tibble <- enrichedTibble(D0K9me3_anno_enrichGO)
D15K9me3_enriched_Tibble <- enrichedTibble(D15K9me3_anno_enrichGO)


# Protocol 6: genome wide and locus specific distribution -----------------
library(ggplot2)
library(RColorBrewer)

#covplot(CTCF_seeker, weightCol = "level" , title="Genome-wide distribution of CTCF peaks")

covplot(seek_D0_K27ac, weightCol= NULL , title="Genome-wide distribution of D0K27ac peaks")
covplot(seek_D15_K27ac, weightCol= NULL , title="Genome-wide distribution of D15K27ac peaks")
covplot(seek_D0_K9me3, weightCol= NULL , title="Genome-wide distribution of D0K9me3 peaks")
covplot(seek_D15_K9me3, weightCol= NULL , title="Genome-wide distribution of D15K9me3 peaks")

# Protocol 7: Heatmaps and Metaplots --------------------------------------

#tagHeatmap(D0K27ac_tagMatrix, title = "H3K4me3 peaks around TES")
#idk what this is or what info it provides...

#tagHeatmap(Epi_list_tagMatrix, xlim=c(-3000, 3000), color = brewer.pal(length(Epi_list_tagMatrix), "Dark2"))
#TSS is transcription start sites. TES is transcription end sites

peakHeatmap(seek_D0_K27ac, TxDb = TxDb_hg38, upstream = 2000, downstream = 2000, 
            title = "Day 0 H3K27ac peaks around TSS")
peakHeatmap(seek_D15_K27ac, TxDb = TxDb_hg38, upstream = 2000, downstream = 2000, 
            title = "Day 15 H3K27ac peaks around TSS")
peakHeatmap(seek_D0_K9me3, TxDb = TxDb_hg38, upstream = 2000, downstream = 2000, 
            title = "Day 0 H3K9me3 peaks around TSS")
peakHeatmap(seek_D15_K9me3, TxDb = TxDb_hg38, upstream = 2000, downstream = 2000, 
            title = "Day 15 H3K9me3 peaks around TSS")

#Meta plots - can change the type to start_site or end_site
D0K27ac_tagMatrix = getTagMatrix(seek_D0_K27ac, TxDb=TxDb_hg38, type = "start_site", 
                                 upstream = 3000, downstream = 3000, by = "gene")
D15K27ac_tagMatrix = getTagMatrix(seek_D15_K27ac, TxDb=TxDb_hg38, type = "start_site", 
                                  upstream = 3000, downstream = 3000, by = "gene")
D0K9me3_tagMatrix = getTagMatrix(seek_D0_K9me3, TxDb=TxDb_hg38, type = "start_site", 
                                 upstream = 3000, downstream = 3000, by = "gene")
D15K9me3_tagMatrix = getTagMatrix(seek_D15_K9me3, TxDb=TxDb_hg38, type = "start_site", 
                                  upstream = 3000, downstream = 3000, by = "gene")

plotAvgProf(D0K27ac_tagMatrix, xlim=c(-3000, 3000), ylab = "Count Frequency")
plotPeakProf(D0K27ac_tagMatrix, facet = "none", conf = 0.95)

plotAvgProf(D15K27ac_tagMatrix, xlim=c(-3000, 3000), ylab = "Count Frequency")
plotPeakProf(D15K27ac_tagMatrix, facet = "none", conf = 0.95)

plotAvgProf(D0K9me3_tagMatrix, xlim=c(-3000, 3000), ylab = "Count Frequency")
plotPeakProf(D0K9me3_tagMatrix, facet = "none", conf = 0.95)

plotAvgProf(D15K9me3_tagMatrix, xlim=c(-3000, 3000), ylab = "Count Frequency")
plotPeakProf(D15K9me3_tagMatrix, facet = "none", conf = 0.95)


TSS = getBioRegion(TxDb=TxDb_hg38, upstream=3000, downstream=3000, by = "gene", 
                   type = "start_site")
Epi_list_tagMatrix = lapply(Epi_list, getTagMatrix, windows = TSS)

#might have to plot these separately. in addition to together
plotAvgProf(Epi_list_tagMatrix, xlim=c(-3000, 3000), ylab = "Count Frequency")
plotPeakProf(Epi_list_tagMatrix, facet = "none", conf = 0.95)




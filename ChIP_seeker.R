rm(list = ls())

setwd("C:/Users/raoda/Desktop/R Stuff/CUT&Tag Processing/Practice with Data")

# Load Data ---------------------------------------------------------------
library(ChIPseeker)
#CTCF_demo <- read.table("CTCF_H1.test.bed")

CTCF_seeker <- readPeakFile("CTCF_H1.test.bed", header = TRUE)
CTCF_seeker
test <- readPeakFile("23_h21792_D0_H3K9me3.tabular", header = FALSE)


# Protocol 1: ChIPseeker and epigenomic dataset prep ------------------------------------------
CTCF_seeker$CTCF_peaks = CTCF_seeker$V4
CTCF_seeker$level = CTCF_seeker$V5
CTCF_seeker$V4 = CTCF_seeker$V5 = NULL
CTCF_seeker
head(getGEOgenomeVersion())

# Protocol 2: Annotation of Epigenomic datasets----------------------------------------------------------
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
TxDb_hg38 = TxDb.Hsapiens.UCSC.hg38.knownGene
ChIPseq_CTCF_demo_anno_default <- annotatePeak(CTCF_seeker, TxDb = TxDb_hg38)

test_annotate <- annotatePeak(test, TxDb = TxDb_hg38)

head(as.GRanges(ChIPseq_CTCF_demo_anno_default), 5)

#write.table(as.data.frame(ChIPseq_CTCF_demo_anno_default), file = "/path/to/download/location")

#if you want to have an annotation priority
#use parameter genomicAnnotationPriority -> for example

#ChIPseq_CTCF_demo_anno_change_priority = annotatePeak(ChIPseq_CTCF_demo, TxDb=TxDb_hg19,
#genomicAnnotationPriority = c("Exon", "Intron", "5UTR", "3UTR", "Promoter", 
#"Downstream", "Intergenic"))


#for more personalization of data sets post annotation, you can use options() and
#reference the protocol

#can also identify gene names using the package "org.Hs.eg.db" and parameter annoDb
ChIPseq_CTCF_demo_anno_gene_name <- annotatePeak(ChIPseq_CTCF_demo, tssRegion
                        = c(-2000,0), TxDb = TxDb_hg19, annoDb = "org.Hs.eg.db")
#gives you new columns of ENSEMBL/ENTREZID, SYMBOL, and GENENAME. use head(as.GRanges) to overview


#this would give you a more condensed version of this sections output of annotated
#peaks. outputs a list. you can also have multiple inputs for multiple epigenomic
#data sets. 
Epi_data_list <- GRangesList(CTCF = CTCF_seeker, h3k9me3 = test)
peakAnnoList_user_defined <- lapply(Epi_data_list, annotatePeak, tssRegion 
                                    = c(-2000,0), TxDb=TxDb_hg38)

#can also specify ranges of chromosomes, start stop sites to your liking. 
# Protocol 4: Visualiztion of annotated results ----------------------------------------------
library(ggVennDiagram)
library(ggupset)

plotAnnoPie(ChIPseq_CTCF_demo_anno_default)
plotAnnoBar(ChIPseq_CTCF_demo_anno_default)
vennpie(ChIPseq_CTCF_demo_anno_default)
upsetplot(ChIPseq_CTCF_demo_anno_default)
upsetplot(ChIPseq_CTCF_demo_anno_default, vennpie=TRUE)

#supports multiple entries if using the peakAnnotatedList object


#  Protocol 5: functional analysis of epigenomic datasets ------------------------------------------
library(clusterProfiler)

#can directly extract the GeneID column
ChIPseq_CTCF_demo_anno_default_genes <- as.data.frame(peakAnnoList_user_defined$CTCF)$geneId

#####Enrichment Analysis#####

#pass genes for functional enrichment analysis
ChIPseq_CTCF_demo_anno_default_enrichGO <- enrichGO(gene = 
                    ChIPseq_CTCF_demo_anno_default_genes, OrgDb = "org.Hs.eg.db")

#can also use compareCluster() to perdorm enrichmen analysis on multiple sets of epi regions
Epi_data_list_genes_enrichGO = compareCluster(geneCluster = peakAnnoList_user_defined_genes,
              fun = "enrichGO", pvalueCutoff = 0.05, OrgDb = "org.Hs.eg.db")

#####Visualization#####

dotplot(Epi_data_list_genes_enrichGO, title = "GO Enrichment Analysis")


# Protocol 6: genome wide and locus specific distribution -----------------


# Protocol 7: Heatmaps and Metaplots --------------------------------------



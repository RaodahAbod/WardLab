rm(list = ls())
setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing")

# Load Data ---------------------------------------------------------------
library(ChIPseeker)

loadFile_peakCall <- function(){
 file <- choose.files()
 file <- readPeakFile(file, header = FALSE)
 return(file)
}
day0_Ch <- loadFile_peakCall()
day1_Ch <- loadFile_peakCall()
day2_Ch <- loadFile_peakCall()
day5_Ch <- loadFile_peakCall()
day15_Ch <- loadFile_peakCall()

# Protocol 1: ChIPseeker and epigenomic dataset prep ------------------------------------------

prepGRangeObj <- function(seek_object){
 seek_object$Peaks = seek_object$V4
 seek_object$level = seek_object$V5
 seek_object$V4 = seek_object$V5 = NULL
 return(seek_object)
}

day0_Ch <- prepGRangeObj(day0_Ch)
day1_Ch <- prepGRangeObj(day1_Ch)
day2_Ch <- prepGRangeObj(day2_Ch)
day5_Ch <- prepGRangeObj(day5_Ch)
day15_Ch <- prepGRangeObj(day15_Ch)


# Protocol 2: Annotation of Epigenomic datasets----------------------------------------------------------
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("TxDb.Ptroglodytes.UCSC.panTro4.refGene")
library("org.Hs.eg.db")
TxDb_hg38 = TxDb.Hsapiens.UCSC.hg38.knownGene
Tx_Pantro4 = TxDb.Ptroglodytes.UCSC.panTro4.refGene

Epi_list <- GRangesList(day0_Ch = day0_Ch, day1_Ch = day1_Ch, day2_Ch = day2_Ch,
                        day5_Ch = day5_Ch, day15_Ch = day15_Ch)

peakAnnoList <- lapply(Epi_list, annotatePeak, tssRegion = c(-2000,2000), 
                       TxDb = TxDb_hg38)
peakAnnoList <- lapply(Epi_list, annotatePeak, tssRegion = c(-2000,2000), 
                       TxDb = Tx_Pantro4)

#list2env(peakAnnoList,envir = .GlobalEnv) #renee!

# Protocol 4: Visualiztion of annotated results ----------------------------------------------
library(ggVennDiagram)
library(ggupset)
library(clusterProfiler)
library(DOSE)

plotAnnoBar(peakAnnoList, main = "Genomic Feature Distribution")

#####Metaplots#####

TSS = getBioRegion(TxDb=TxDb_hg38, upstream=3000, downstream=3000, by = "gene", 
                   type = "start_site")
TSS = getBioRegion(TxDb=Tx_Pantro4, upstream=3000, downstream=3000, by = "gene", 
                   type = "start_site")
Epi_list_tagMatrix = lapply(Epi_list, getTagMatrix, windows = TSS)

#might have to plot these separately. in addition to together
plotAvgProf(Epi_list_tagMatrix, xlim=c(-3000, 3000), ylab = "Count Frequency")
#plotPeakProf(Epi_list_tagMatrix, facet = "none", conf = 0.95)




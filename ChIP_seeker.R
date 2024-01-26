rm(list = ls())
setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing")

# Load Data ---------------------------------------------------------------
library(ChIPseeker)

loadFile_peakCall <- function(){
 file <- choose.files()
 file <- readPeakFile(file, header = FALSE)
 return(file)
}
day0_S <- loadFile_peakCall()
day0h_S <- loadFile_peakCall()
day1_S <- loadFile_peakCall()
day2_S <- loadFile_peakCall()
day3_S <- loadFile_peakCall()
day4_S <- loadFile_peakCall()
day5_S <- loadFile_peakCall()
day7_S <- loadFile_peakCall()
day10_S <- loadFile_peakCall()
day15_S <- loadFile_peakCall()
day30_S <- loadFile_peakCall()

# Protocol 1: ChIPseeker and epigenomic dataset prep ------------------------------------------

prepGRangeObj <- function(seek_object){
 seek_object$Peaks = seek_object$V4
 seek_object$level = seek_object$V5
 seek_object$V4 = seek_object$V5 = NULL
 return(seek_object)
}

day0_S <- prepGRangeObj(day0_S)
day0h_S <- prepGRangeObj(day0h_S)
day1_S <- prepGRangeObj(day1_S)
day2_S <- prepGRangeObj(day2_S)
day3_S <- prepGRangeObj(day3_S)
day4_S <- prepGRangeObj(day4_S)
day5_S <- prepGRangeObj(day5_S)
day7_S <- prepGRangeObj(day7_S)
day10_S <- prepGRangeObj(day10_S)
day15_S <- prepGRangeObj(day15_S)
day30_S <- prepGRangeObj(day30_S)

# Protocol 2: Annotation of Epigenomic datasets----------------------------------------------------------
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")
TxDb_hg38 = TxDb.Hsapiens.UCSC.hg38.knownGene

Epi_list <- GRangesList(day0_S = day0_S, day0h_S = day0h_S, day1_S = day1_S, day2_S = day2_S,
                        day3_S = day3_S, day4_S = day4_S, day5_S = day5_S, day7_S = day7_S,
                        day10_S = day10_S, day15_S = day15_S, day30_S = day30_S)

peakAnnoList <- lapply(Epi_list, annotatePeak, tssRegion = c(-2000,2000), 
                       TxDb = TxDb_hg38)

list2env(peakAnnoList,envir = .GlobalEnv) #renee!

# Protocol 4: Visualiztion of annotated results ----------------------------------------------
library(ggVennDiagram)
library(ggupset)
library(clusterProfiler)
library(DOSE)

plotAnnoBar(peakAnnoList, main = "Genomic Feature Distribution")

#####Metaplots#####

TSS = getBioRegion(TxDb=TxDb_hg38, upstream=3000, downstream=3000, by = "gene", 
                   type = "start_site")
Epi_list_tagMatrix = lapply(Epi_list, getTagMatrix, windows = TSS)

#might have to plot these separately. in addition to together
plotAvgProf(Epi_list_tagMatrix, xlim=c(-3000, 3000), ylab = "Count Frequency")
#plotPeakProf(Epi_list_tagMatrix, facet = "none", conf = 0.95)




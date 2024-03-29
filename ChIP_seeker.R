#Comment Annotations - 02/27/24 - RA

#Removes all variables from environment and workspace.
rm(list = ls())

setwd("C:/Users/raoda/Desktop/R Stuff/CUT&Tag Processing")

# Load Data ---------------------------------------------------------------
library(ChIPseeker)

#Function opens file directory. Input is read and returns as a peak files
loadFile_peakCall <- function(){
 file <- choose.files()
 file <- readPeakFile(file, header = FALSE)
 return(file)
}

#calls function loadFile - stores as peak file
sample1 <- loadFile_peakCall()
sample2 <- loadFile_peakCall()
sample3 <- loadFile_peakCall()
sample4 <- loadFile_peakCall()
sample5 <- loadFile_peakCall()
sample6 <- loadFile_peakCall()
sample10 <- loadFile_peakCall()
sample11 <- loadFile_peakCall()

sample5_old <- loadFile_peakCall()
sample6_old <- loadFile_peakCall()
sample10_old <- loadFile_peakCall()
sample11_old <- loadFile_peakCall()

# Protocol 1: ChIPseeker and epigenomic dataset prep ------------------------------------------

#function organizes peak file and returns as a GRanges object compatible 
#with ChIPSeeker functions
prepGRangeObj <- function(seek_object){
 seek_object$Peaks = seek_object$V4
 seek_object$level = seek_object$V5
 seek_object$V4 = seek_object$V5 = NULL
 return(seek_object)
}

#calls the prepGrangeObj function and overwrites original peakfile variable
sample1 <- prepGRangeObj(sample1)
sample2 <- prepGRangeObj(sample2)
sample3 <- prepGRangeObj(sample3)
sample4 <-prepGRangeObj(sample4)
sample5 <- prepGRangeObj(sample5)
sample6 <- prepGRangeObj(sample6)
sample10 <- prepGRangeObj(sample10)
sample11 <- prepGRangeObj(sample11)

sample5_old <- prepGRangeObj(sample5_old)
sample6_old <- prepGRangeObj(sample6_old)
sample10_old <- prepGRangeObj(sample10_old)
sample11_old <- prepGRangeObj(sample11_old)

# Protocol 2: Annotation of Epigenomic datasets----------------------------------------------------------

#loads in reference genomes used to align peak files
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#library("TxDb.Ptroglodytes.UCSC.panTro4.refGene")
library("org.Hs.eg.db")
#library("org.Pt.eg.db")

#pantro <- getGEOInfo(genome="pantro 4", simplify=TRUE)

#creates a variable of the loaded reference genome for use in ChIPseeker functions
TxDb_hg38 = TxDb.Hsapiens.UCSC.hg38.knownGene
#Tx_Pantro4 = TxDb.Ptroglodytes.UCSC.panTro4.refGene

#Stores a concatenated list of peak files. user specified, order and name of peaks will be 
#the new name of peaks in the list and will be used in figure outputs. 
Epi_list <- GRangesList(`78-1 H3K27me3 Accutase` = sample1, 
                        `78-1 H3K27me3 Trypsin` = sample2,
                        `78-1 H3K27me3 Flash Froz` = sample3,
                        `78-1 H3K27me3 Scraped` = sample4,
                        `H21792 H3K27ac Accutase` = sample5,
                        `H21792 H3K27ac Trypsin` = sample6,
                        `78-1 H3K27ac Accutase` = sample10,
                        `78-1 H3K27ac Scraped` = sample11) 

old_Epi_list <- GRangesList(`H21792 H3K27ac Accutase (orig Broad M2 peaks)` = sample5_old,
                            `H21792 H3K27ac Trypsin (orig Broad M2 peaks)` = sample6_old,
                            `78-1 H3K27ac Accutase (orig Broad M2 peaks)` = sample10_old,
                            `78-1 H3K27ac Scraped (orig Broad M2 peaks)` = sample11_old)

peakAnnoList <- lapply(Epi_list, annotatePeak, tssRegion = c(-2000,2000),TxDb = TxDb_hg38)

peakAnnoList2 <- lapply(old_Epi_list, annotatePeak, tssRegion = c(-2000,2000),TxDb = TxDb_hg38)

#list2env(peakAnnoList,envir = .GlobalEnv) #renee!

# Protocol 4: Visualiztion of annotated results ----------------------------------------------
library(ggVennDiagram)
library(ggupset)
library(clusterProfiler)
library(DOSE)

#plotting the annotation bar plot. Input is peakAnnoList
plotAnnoBar(peakAnnoList)
plotAnnoBar(peakAnnoList2)
# Metaplots -------------------------------------------------------------------------

#Function to generate the locations of transcription start sites in your reference genome.
TSS = getBioRegion(TxDb=TxDb_hg38, upstream=3000, downstream=3000, by = "gene", 
                   type = "start_site")

Epi_list_tagMatrix = lapply(Epi_list, getTagMatrix, windows = TSS)

#plotting the transcription Start Site density meta plot.
plotAvgProf(Epi_list_tagMatrix, xlim=c(-3000, 3000), ylab = "Count Frequency")

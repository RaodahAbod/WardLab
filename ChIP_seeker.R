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
chimpFile <- loadFile_peakCall()

Day0_IgG_H21792 <- loadFile_peakCall()

Day0_iPSC_H3K27ac <- loadFile_peakCall()
Day0_iPSC_H3K9me3 <- loadFile_peakCall()
Day15_CM_H3K27ac <- loadFile_peakCall()
Day15_CM_H3K9me3 <- loadFile_peakCall()
Day30_CM_H3K27ac <- loadFile_peakCall()

ENCODE_iPSC20b_H3K27ac <- loadFile_peakCall()
ENCODE_LV_H3K27ac <- loadFile_peakCall()

Roadmap_iPSC20b_H3K27ac <- loadFile_peakCall()
Roadmap_iPSC20b_H3K9me3 <- loadFile_peakCall()
Roadmap_LV_H3K27ac <- loadFile_peakCall()
Roadmap_LV_H3K9me3 <- loadFile_peakCall()


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
chimpFile <- prepGRangeObj(chimpFile)

Day0_IgG_H21792 <- prepGRangeObj(Day0_IgG_H21792)

Day0_iPSC_H3K27ac <- Day0_iPSC_H3K27ac %>% prepGRangeObj()
Day0_iPSC_H3K9me3 <- Day0_iPSC_H3K9me3 %>% prepGRangeObj()
Day15_CM_H3K27ac <- Day15_CM_H3K27ac %>% prepGRangeObj()
Day15_CM_H3K9me3 <- Day15_CM_H3K9me3 %>% prepGRangeObj()
Day30_CM_H3K27ac <- prepGRangeObj(Day30_CM_H3K27ac)

ENCODE_iPSC20b_H3K27ac <- prepGRangeObj(ENCODE_iPSC20b_H3K27ac)
ENCODE_LV_H3K27ac <- prepGRangeObj(ENCODE_LV_H3K27ac)

Roadmap_iPSC20b_H3K27ac <- prepGRangeObj(Roadmap_iPSC20b_H3K27ac)
Roadmap_iPSC20b_H3K9me3 <- prepGRangeObj(Roadmap_iPSC20b_H3K9me3)
Roadmap_LV_H3K27ac <- prepGRangeObj(Roadmap_LV_H3K27ac)
Roadmap_LV_H3K9me3 <- prepGRangeObj(Roadmap_LV_H3K9me3)

# Protocol 2: Annotation of Epigenomic datasets----------------------------------------------------------

#loads in reference genomes used to align peak files
#library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("TxDb.Ptroglodytes.UCSC.panTro4.refGene")
#library("org.Hs.eg.db")
library("org.Pt.eg.db")
pantro <- getGEOInfo(genome="pantro 4", simplify=TRUE)

#creates a variable of the loaded reference genome for use in ChIPseeker functions
#TxDb_hg38 = TxDb.Hsapiens.UCSC.hg38.knownGene
Tx_Pantro4 = TxDb.Ptroglodytes.UCSC.panTro4.refGene

#Stores a concatenated list of peak files. user specified, order and name of peaks will be 
#the new name of peaks in the list and will be used in figure outputs. 
Epi_list <- GRangesList(chimpD15H3K27ac = chimpFile)

#GRanges(originalGRangeObj = newDataName) <-- setup of GRangesList input
Epi_list_k27 <- GRangesList(Day0_iPSC_H3K27ac = Day0_iPSC_H3K27ac, 
                            Day15_CM_H3K27ac = Day15_CM_H3K27ac, 
                            Day30_CM_H3K27ac = Day30_CM_H3K27ac, 
                            ENCODE_iPSC20b_H3K27ac = ENCODE_iPSC20b_H3K27ac,
                            ENCODE_LV_H3K27ac = ENCODE_LV_H3K27ac,
                            Roadmap_iPSC20b_H3K27ac = Roadmap_iPSC20b_H3K27ac,
                            Roadmap_LV_H3K27ac = Roadmap_LV_H3K27ac)

Epi_list_k9 <- GRangesList(Day0_iPSC_H3K9me3 = Day0_iPSC_H3K9me3,
                           Day15_CM_H3K9me3 = Day15_CM_H3K9me3,
                           Roadmap_iPSC20b_H3K9me3 = Roadmap_iPSC20b_H3K9me3,
                           Roadmap_LV_H3K9me3 = Roadmap_LV_H3K9me3)

#peakAnnoList_k27 <- lapply(Epi_list_k27, annotatePeak, tssRegion = c(-2000,2000), 
                       #TxDb = TxDb_hg38)
#peakAnnoList_k9 <- lapply(Epi_list_k9, annotatePeak, tssRegion = c(-2000,2000), 
                      # TxDb = TxDb_hg38)

#function of ChIPseeker to annotate peak files within the list. 
peakAnnoList <- lapply(Epi_list, annotatePeak, tssRegion = c(-2000,2000), 
                       TxDb = pantro)

peakAnnoList <- lapply(Epi_list, annotatePeak, tssRegion = c(-2000,2000), 
                       TxDb = Tx_Pantro4)

#list2env(peakAnnoList,envir = .GlobalEnv) #renee!

# Protocol 4: Visualiztion of annotated results ----------------------------------------------
library(ggVennDiagram)
library(ggupset)
library(clusterProfiler)
library(DOSE)

#plotting the annotation bar plot. Input is peakAnnoList
plotAnnoBar(peakAnnoList_k27)
plotAnnoBar(peakAnnoList_k9)

plotAnnoBar(peakAnnoList)

# Metaplots -------------------------------------------------------------------------


#Function to generate the locations of transcription start sites in your reference genome.
TSS = getBioRegion(TxDb=TxDb_hg38, upstream=3000, downstream=3000, by = "gene", 
                   type = "start_site")
#TSS = getBioRegion(TxDb=Tx_Pantro4, upstream=3000, downstream=3000, by = "gene", 
#                   type = "start_site")

#tagging of peaks by transcription start site. input is Epi_list
Epi_list_tagMatrix_k27 = lapply(Epi_list_k27, getTagMatrix, windows = TSS)
Epi_list_tagMatrix_k9 = lapply(Epi_list_k9, getTagMatrix, windows = TSS)

#plotting the transcription Start Site density meta plot.
plotAvgProf(Epi_list_tagMatrix_k27, xlim=c(-3000, 3000), ylab = "Count Frequency")
plotAvgProf(Epi_list_tagMatrix_k9, xlim=c(-3000, 3000), ylab = "Count Frequency")




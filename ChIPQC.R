rm(list = ls())
setwd("C:/Users/raoda/Desktop/R Stuff/CUT&Tag Processing")

# Load Files ------------------------------------------------------------------------
library("ChIPQC")
library("BiocParallel")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("GenomeInfoDb")


dataName_k27me3 <- c("ENCODE LV",'1_78-1_Accutase','2_78-1_Trypsin','4_78-1_Scraped') 
dataName_k27ac <- c('ENCODE LV','5_H21792_Accutase','6_H21792_Trypsin','10_78-1_Accutase',
                    '11_78-1_Scraped')
bamFiles <- c()
peakFiles <- c()

loadData <- function(){
        bamFiles <- choose.files()
        peakFiles <- choose.files()
        testExample <- ChIPQCsample(bamFiles,peakFiles, annotation = "hg38")
        return(testExample)
}

for (x in dataName_k27ac){
        H3K27ac <- loadData()
        assign(x, H3K27ac)
}
# Metrics Overview ------------------------------------------------------------------

#this is how to load all with the individual data sets
sampleSheet <- read.csv(choose.files())
        data_list <- list("ENCODE_LV" = `ENCODE LV`,
                         "1_78-1_Accutase" = `1_78-1_Accutase`,
                         "2_78-1_Trypsin" = `2_78-1_Trypsin`,
                         "4_78-1_Scraped" = `4_78-1_Scraped`)
        
        data_list <- list("ENCODE_LV" = `ENCODE LV`,
                          "5_H21792_Accutase" = `5_H21792_Accutase`,
                          "6_H21792_Trypsin" = `6_H21792_Trypsin`,
                          "10_78-1_Accutase" = `10_78-1_Accutase`,
                          "11_78-1_Scraped" = `11_78-1_Scraped`)
        
pooledSamples = ChIPQC(sampleSheet,annotation = 'hg38', samples = data_list)
plotCorHeatmap(pooledSamples, colScheme="Reds", margin=15)
#plotPrincomp(pooledSamples)

#run samples indiv, put into a list, and then run them as an experiemnt to ensure it goes 
#through all chromosomes, since i am not sure they do if i run them all together
#i can run QCmetrics() and indiv metrics,plot: corHeatmap(), coverageHist(), 
#Regi() (not sure this is helpful tho), plotSSD(), plotFrip(), plotRap(), plotPrincomp()

# Metrics and Plots -----------------------------------------------------------------

metrics <- QCmetrics(pooledSamples)
setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing")
metrics_fileName <- readline(prompt = "Enter the metrics excel filename: ")
metrics_fileName <- paste0(metrics_fileName,".csv")
write.csv(metrics, file = metrics_fileName)

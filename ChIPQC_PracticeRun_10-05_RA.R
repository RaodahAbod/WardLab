rm(list = ls())
setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing")

# Load Files ------------------------------------------------------------------------
library("ChIPQC")
library("BiocParallel")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("GenomeInfoDb")

num_of_files <- readline(prompt = "Enter the number of Files: ")
peak_caller <- readline(prompt = "Enter first letter of the Peak Caller (b n or s): ")
testfile <- rep(0,num_of_files)
testPeaks <- rep(0, num_of_files)

dataName <- c("Day0","Day0h","Day1","Day2","Day3","Day4","Day5","Day7","Day10","Day15","Day30") 

for (x in 1:num_of_files){
        # setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing/Bam Files")
        testfile[x] <- choose.files()
        if(peak_caller == "b"){
                setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing/MACS2 Broad Peaks")
        }else if(peak_caller == "n"){
                setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing/MACS2 Narrow Peaks")
                }else if(peak_caller == "s"){
                        setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing/SEACR Files")
                }
        testPeaks[x] <- choose.files()

}
# Metrics Overview ------------------------------------------------------------------

for(n in 1:num_of_files){
testExample <- ChIPQCsample(testfile[n],testPeaks[n], annotation = "hg38")

if(peak_caller == "b"){
        nameOfData <- paste0(dataName[n],"_H21792_Broad")
        }else if(peak_caller == 'n'){
        nameOfData <- paste0(dataName[n],"_H21792_Narrow")
                }else if(peak_caller == 's'){
                nameOfData <- paste0(dataName[n],"_H21792_SEACR")
                        #else{ nameOfData <- paste0(dataName[n+4])
                        #}
                }
assign(nameOfData, testExample)
}

#this is how to load all with the individual data sets
sampleSheet <- read.csv(choose.files())
        if(peak_caller == "b"){
                data_list <- list("1_D0_H3K27ac" = Day0_H21792_Broad,"2_D0h_H3K27ac" = Day0h_H21792_Broad,
                                  "3_D1_H3K27ac" = Day1_H21792_Broad,"4_D2_H3K27ac" = Day2_H21792_Broad,
                                  "5_D3_H3K27ac" = Day3_H21792_Broad,"6_D4_H3K27ac" = Day4_H21792_Broad,
                                  "7_D5_H3K27ac" = Day5_H21792_Broad,"8_D7_H3K27ac" = Day7_H21792_Broad,
                                  "9_D10_H3K27ac" = Day10_H21792_Broad,"10_D15_H3K27ac" = Day15_H21792_Broad,
                                  "11_D30_H3K27ac" = Day30_H21792_Broad)
        }else if(peak_caller == 'n'){
                data_list <- list("1_D0_H3K27ac" = Day0_H21792_Narrow,"2_D0h_H3K27ac" = Day0h_H21792_Narrow,
                                  "3_D1_H3K27ac" = Day1_H21792_Narrow,"4_D2_H3K27ac" = Day2_H21792_Narrow,
                                  "5_D3_H3K27ac" = Day3_H21792_Narrow,"6_D4_H3K27ac" = Day4_H21792_Narrow,
                                  "7_D5_H3K27ac" = Day5_H21792_Narrow,"8_D7_H3K27ac" = Day7_H21792_Narrow,
                                  "9_D10_H3K27ac" = Day10_H21792_Narrow,"10_D15_H3K27ac" = Day15_H21792_Narrow,
                                  "11_D30_H3K27ac" = Day30_H21792_Narrow)
        }else if(peak_caller == 's'){
                data_list <- list("1_D0_H3K27ac" = Day0_H21792_SEACR,"2_D0h_H3K27ac" = Day0h_H21792_SEACR,
                                  "3_D1_H3K27ac" = Day1_H21792_SEACR,"4_D2_H3K27ac" = Day2_H21792_SEACR,
                                  "5_D3_H3K27ac" = Day3_H21792_SEACR,"6_D4_H3K27ac" = Day4_H21792_SEACR,
                                  "7_D5_H3K27ac" = Day5_H21792_SEACR,"8_D7_H3K27ac" = Day7_H21792_SEACR,
                                  "9_D10_H3K27ac" = Day10_H21792_SEACR,"10_D15_H3K27ac" = Day15_H21792_SEACR,
                                  "11_D30_H3K27ac" = Day30_H21792_SEACR)
        }
pooledSamples = ChIPQC(sampleSheet,annotation = 'hg38', samples = data_list)


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
if(peak_caller == "b"){
        plotCorHeatmap(pooledSamples, colScheme="Blues")
        }else if(peak_caller == 'n'){
                plotCorHeatmap(pooledSamples, colScheme="Reds")
        }else if(peak_caller == 's'){
                plotCorHeatmap(pooledSamples, colScheme="Greens")
}


# Ill do it myself ------------------------------------------------------------------

extractDataIntoList <- function(ChIPQC_output){
 AvgPeakSig1 <- ChIPQC_output@AveragePeakSignal[1]
 AvgPeakSig2 <- ChIPQC_output@AveragePeakSignal[2]
 crossCov <- ChIPQC_output@CrossCoverage
 covHist <- ChIPQC_output@CoverageHistogram
 QC_metrics <- QCmetrics(ChIPQC_output)
 return(list('AvgPeakSig1' = AvgPeakSig1, 'AvgPeakSig2' = AvgPeakSig2, 
             'crossCov' = crossCov, 'covHist' = covHist, "ChIPQC Metrics" = QC_metrics))
}

data_list <- list("1_D0_H3K27ac" = Day0_H21792,"2_D0h_H3K27ac" = Day0h_H21792)

day0_H21792 <- extractDataIntoList(testExample)
day0h_H21792 <- extractDataIntoList(testExample2)



# lol i tired -----------------------------------------------------------------------
register(SerialParam())
exampleExp <- ChIPQC(file, annotation="hg38", consensus=TRUE)
maybe <- ChIPQCreport(exampleExp)

diffBindDBA <- dba(sampleSheet = choose.files())
anotherTry <- ChIPQC(diffBindDBA, annotation = 'hg38')

sampleSheet = read.csv(file.choose())

#View(ChIPQC)
     
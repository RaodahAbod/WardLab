rm(list = ls())
setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing")

# Load Data ---------------------------------------------------------------
loadFile <- function(){
  file <- read.table(choose.files())
  return(file)
}

day0_H <- loadFile()
day0h_H <- loadFile()
day1_H <- loadFile()
day2_H <- loadFile()
day3_H <- loadFile()
day4_H <- loadFile()
day5_H <- loadFile()
day7_H <- loadFile()
day10_H <- loadFile()
day15_H <- loadFile()
day30_H <- loadFile()

dataName <- readline(prompt = "Please enter the full Data Name: ")
autoDataName <- paste("autosomal",dataName)

# Removing Non Autosomal Data ---------------------------------------------
nonautosomal <- function(histonemod_day) {
  test <- subset(histonemod_day,histonemod_day$V1 == "chr1" | histonemod_day$V1 == "chr2" |
                   histonemod_day$V1 == "chr3" | histonemod_day$V1 == "chr4" |
                   histonemod_day$V1 == "chr5" | histonemod_day$V1 == "chr6"  |
                   histonemod_day$V1 == "chr7" | histonemod_day$V1 == "chr8" |
                   histonemod_day$V1 == "chr9" | histonemod_day$V1 == "chr10"  |
                   histonemod_day$V1 == "chr11" | histonemod_day$V1 == "chr12" |
                   histonemod_day$V1 == "chr13" | histonemod_day$V1 == "chr14"  |
                   histonemod_day$V1 == "chr15" | histonemod_day$V1 == "chr16" |
                   histonemod_day$V1 == "chr17" | histonemod_day$V1 == "chr18"  |
                   histonemod_day$V1 == "chr19" | histonemod_day$V1 == "chr20" |
                   histonemod_day$V1 == "chr21" | histonemod_day$V1 == "chr22")
  return(test)
}

nonautosomalNew <- function(histonemod_day) {
  test <- histonemod_day %>%
    select(-contains("chrUn"), -contains("chrX"), -contains("chrY"), -contains("chrM"))
  return(test)
}

orig <- nonautosomal(workingFile)
new <- nonautosomalNew(workingFile)

autos_workingFile <- nonautosomal(workingFile)

# Number of Peaks in Samples ----------------------------------------------
numPeaks <- nrow(workingFile)
autoNumPeaks <- nrow(autos_workingFile)

barplot(c(numPeaks, autoNumPeaks), ylab = "Number of Peaks", xlab = dataName,
        main = paste("Total Number of Peaks in", dataName), col = "pink",
        names = c("Non Autosomal (all)","Autosomal (filtered)"))

# Length of Peaks and Total Signal Value   ------------------------------------------------
ObtainPeakLength <- function(histonemod){
  PeakLength <- histonemod[,3] - histonemod[,2]
  histonemod <- cbind(histonemod,PeakLength)
  return(histonemod)
}
workingFile <- ObtainPeakLength(workingFile)
autos_workingFile <- ObtainPeakLength(autos_workingFile)

boxplot(workingFile$PeakLength, autos_workingFile$PeakLength, col = "pink",
        main = paste("Length of Peaks with Outliers in", dataName), ylab = "Peak Length",
        xlab = dataName, names = c("Non Autosomal (all)","Autosomal (filtered)"))

boxplot(workingFile$V7, autos_workingFile$V7, col = "pink",
        main = paste("Total Signal Values with Outliers in",dataName), ylab = "Signal Intensity",
        xlab = dataName, names= c("Non Autosomal (all)","Autosomal (filtered)"))

boxplot(workingFile$PeakLength, autos_workingFile$PeakLength, outline = FALSE,
        col = "pink", main = paste("Length of Peaks in",dataName), ylab = "Peak Length",
        xlab = dataName, names = c("Non Autosomal (all)","Autosomal (filtered)"))

boxplot(workingFile$V7, autos_workingFile$V7, outline = FALSE, col = "pink",
        main = paste("Total Signal Values in", dataName), ylab = "Signal Intensity",
        xlab = dataName, names= c("Non Autosomal (all)","Autosomal (filtered)"))

# Distribution of Peak Length (Autosomes) ---------------------------------------------

hist(autos_workingFile$PeakLength, main = paste("Distribution of Peak Length in", dataName),
     xlim = c(0,10000), breaks = 500, xlab = "Peak Length")


# Distribution of Peak Intensity (Autosomes) ------------------------------------------

hist(autos_workingFile$V7, main = paste("Distribution of Peak Intensity in", dataName),
     xlim = c(0,100), breaks = 50, xlab = "Peak Intensity")

# Save FIltered Autosomals ------------------------------------------------

newFileName <- readline(prompt = "Enter the name for new file: ")
newFileName <- paste0(newFileName,".bed")
write.table(autos_workingFile, newFileName ,row.names = F ,col.names = F ,
            sep="\t" , quote=FALSE)








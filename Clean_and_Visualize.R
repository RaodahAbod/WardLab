rm(list = ls())
setwd("C:/Users/raoda/Desktop/R Stuff/CUT&Tag Processing")

# Load Data ---------------------------------------------------------------
loadFile <- function(){
  file <- read.table(choose.files())
  return(file)
}

workingFile <- loadFile()

dataName <- readline(prompt = "Please enter the full Data Name: ")
autoDataName <- paste("autosomal",dataName)

# Removing Non Autosomal Data ---------------------------------------------
nonautosomalNew <- function(histonemod_day) {
  test <- histonemod_day %>%
    select(-contains("chrUn"), -contains("chrX"), -contains("chrY"), -contains("chrM"))
  return(test)
}

orig <- nonautosomal(workingFile)
new <- nonautosomalNew(workingFile)

autos_workingFile <- nonautosomal(workingFile)
autos_NewWorkingFile <- nonautosomalNew(workingFile)

# Save FIltered Autosomals ------------------------------------------------

newFileName <- readline(prompt = "Enter the name for new file: ")
newFileName <- paste0(newFileName,".bed")
write.table(autos_workingFile, newFileName ,row.names = F ,col.names = F ,
            sep="\t" , quote=FALSE)

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










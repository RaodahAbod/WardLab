rm(list = ls())
setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing")

#  Load File ------------------------------------------------------------------------
workingData <- read.table(choose.files())

# Split into Top groups -------------------------------------------------------------
#I need a top 1/4th, top 1/3rd, top 1/2, and top 2/3 :D
#get the criteria for the top whatever then use the select function (i think) to separate 
#based on equal to or greater than the critera

separatePeaks <- function(peak_file,topWhatever){
   sorted <- arrange(peak_file,-V7)
   index <- floor(length(sorted[,1])*topWhatever)
   #top_criteria <- peak_file$V5[index]
   topPeaks <- sorted[1:index,]
   #topPeaks <- sorted[index:length(sorted[,1]),]
   #topPeaks <- subset(sorted, sorted >= top_criteria)
   return(topPeaks)
}

topHalf_Peaks <- separatePeaks(workingData, 1/2) 
topThird_Peaks <- separatePeaks(workingData, 1/3) 
topFourth_Peaks <- separatePeaks(workingData, 1/4)
topTwoThird_Peaks <- separatePeaks(workingData, 2/3)

# Save Grouped Peaks ----------------------------------------------------------------

setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing/Quartiled Bed Files")


saveQuartileFile <- function(dataset){
   newFileName <- readline(prompt = "Enter the name for new file: ")
   newFileName <- paste0(newFileName,".bed")
   write.table(dataset, newFileName ,row.names = F ,col.names = F ,
               sep="\t" , quote=FALSE) 
}

saveQuartileFile(topHalf_Peaks)
saveQuartileFile(topThird_Peaks)
saveQuartileFile(topFourth_Peaks)
saveQuartileFile(topTwoThird_Peaks)


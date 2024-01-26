rm(list = ls())
setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing")

library("tidyverse")
library("readxl")

# Load Data -------------------------------------------------------------------------

workingData <- read.table(choose.files())

# Section into Quartiles.  ----------------------------------------------------------

quartilesRoad <- function(dataSet){
   firstQ <- quantile(dataSet$V5, 0.25)
   secondQ <- quantile(dataSet$V5, 0.5)
   thirdQ <- quantile(dataSet$V5, 0.75)
   
   Percent25 <- subset(dataSet, dataSet$V5 <= firstQ)
   Percent50 <- subset(dataSet, dataSet$V5 > firstQ & dataSet$V5 <= secondQ)
   Percent75 <- subset(dataSet, dataSet$V5 > secondQ & dataSet$V5 <= thirdQ)
   Percent100 <- subset(dataSet, dataSet$V5 > thirdQ)
   
   sortedQuartiles <- list(Percent25,Percent50,Percent75,Percent100)
   return(sortedQuartiles)
}
quartiles <- function(dataSet){
   firstQ <- quantile(dataSet$V4, 0.25)
   secondQ <- quantile(dataSet$V4, 0.5)
   thirdQ <- quantile(dataSet$V4, 0.75)
   
      Percent25 <- subset(dataSet, dataSet$V4 <= firstQ)
      Percent50 <- subset(dataSet, dataSet$V4 > firstQ & dataSet$V4 <= secondQ)
      Percent75 <- subset(dataSet, dataSet$V4 > secondQ & dataSet$V4 <= thirdQ)
      Percent100 <- subset(dataSet, dataSet$V4 > thirdQ)
      
   sortedQuartiles <- list(Percent25,Percent50,Percent75,Percent100)
   return(sortedQuartiles)
}

quartiled_data <- quartilesRoad(workingData)

quartiled_data <- quartiles(workingData)
Q1 <- quartiled_data[[1]]
Q2 <- quartiled_data[[2]]
Q3 <- quartiled_data[[3]]
Q4 <- quartiled_data[[4]]

# Save Quartile Information ---------------------------------------------------------
setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing/Quartiled Bed Files")

saveQuartileFile <- function(dataset){
newFileName <- readline(prompt = "Enter the name for new file: ")
newFileName <- paste0(newFileName,".bed")
write.table(dataset, newFileName ,row.names = F ,col.names = F ,
            sep="\t" , quote=FALSE) 
}

saveQuartileFile(Q1)
saveQuartileFile(Q2)
saveQuartileFile(Q3)
saveQuartileFile(Q4)
# Simple Split into Quartiles using new tricks -------------------------------------------------------

addQuartiles <- function(dataSet){
   firstQ <- quantile(dataSet$V4)
   secondQ <- quantile(dataSet$V4)
   thirdQ <- quantile(dataSet$V4)
   
   dataSet %>%
      dplyr::select(V1,V2, V3, V4) %>%
      mutate(Quartiles = if_else(V4 <= firstQ[2,1],0.25, 
                                 if_else(V4 > quantile(dataSet$V4,0.25)&& V4 <= quantile(dataSet$V4,0.5),0.5, 
                                         if_else(V4 > quantile(dataSet$V4,0.5)&& V4 <= quantile(dataSet$V4,0.75), 0.75),1.00))) -> QuartileData
   
   return(QuartileData)
}

data_two_quartiles <- addQuartiles(DataSet_two)


#doesnt work. 
# test <- function(dataSet){ 
#    firstQ <- quantile(dataSet$V4, 0.25)
#    secondQ <- quantile(dataSet$V4, 0.5)
#    thirdQ <- quantile(dataSet$V4, 0.75)
#    
#    ifelse(dataSet$V4 <= firstQ, dataSet$V5 <- 0.25, 
#           ifelse(dataSet$V4 > firstQ & dataSet$V4 <= secondQ, dataSet$V5 <- 0.5,
#                  ifelse(dataSet$V4 > secondQ & dataSet$V4 <= thirdQ, dataSet$V5 <- 0.75,
#                         NA)))
#    return(dataSet)
#    
#    }



#put these lines of code into the beginning of the chipseeker script
# Q1 <- quartile_dataOne[[1]]
# Q2 <- quartile_dataOne[[2]]
# Q3 <- quartile_dataOne[[3]]
# Q4 <- quartile_dataOne[[4]]


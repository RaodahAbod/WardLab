rm(list = ls())
setwd("C:/Users/raoda/Desktop/R Stuff/CUT&Tag Processing")
library(tidyverse)
library(dplyr)
library(ggplot2)

batch1 <- c("1_78-1_K27me3_Accutase", "2_78-1_K27me3_Trypsin",
            "3_78-1_K27me3_FF", "4_78-1_K27me3_Scraped")
batch2 <- c("5_H21792_K27ac_Accutase", "6_H21792_K27ac_Trypsin")
batch3 <- c("7_H21792_K9me3_Trypsin","8_H21792_K9me3_Accutase","9_H21792_K9me3_Scraped" )
batch4 <- c("10_78-1_K27ac_Accutase", "11_78-1_K27ac_Scraped")

batchIt <- function(namesOfData, batchNum){
  file <- c()
for (dataset in namesOfData){
  file = read.table(choose.files(), header = FALSE) %>% 
    mutate(`Data Name` = dataset, Batch = batchNum) %>%
    rbind(file, .)
}
  batch <- file
  return(batch)
}
batchOne <- batchIt(batch1, 1)
batchTwo <- batchIt(batch2, 2)
batchThree <- batchIt(batch3, 3)
batchFour <- batchIt(batch4, 4)

fullSample <- rbind(batchOne, batchTwo, batchThree, batchFour)

autoData <- subset(fullSample, !grepl("chrX|chrY|chrM|chrUn", fullSample$V1))

#only to run and save autosomal files
for (autos in namesOfData){
  origData <- read.table(choose.files(), header = FALSE)
  newFileName <- paste0(autos,"_Human_SEACR_H3K27ac.bed")
  autoData <- subset(origData, !grepl("chrX|chrY|chrM|chrUn", origData$V1))
  write.table(autoData, newFileName ,row.names = F ,col.names = F ,
            sep="\t" , quote=FALSE) 
}

# Plot Num of Peaks -------------------------------------------------------
#the golden format :)
ggplot(fullSample, aes(x = `Data Name`, fill = `Data Name`)) + geom_bar() +
  facet_wrap(~Batch, scales = "free_x", ncol = 4) +
  labs(title = "Lengths of Peaks Across M2Broad Samples Grouped by Batch") +
  geom_text(stat = "count", aes(label = after_stat(count)),
            position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Data Sample") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Boxplot of Peaks --------------------------------------------------------

ObtainPeakLength <- function(histonemod){
  histonemod <-  histonemod %>%
      mutate(histonemod,peakLength = histonemod[,3] - histonemod[,2]) %>%
  return()
}

autoData <- ObtainPeakLength(autoData)

ggplot(autoData, aes(x = `Data Name`, y = peakLength)) + geom_boxplot(outlier.shape = NA) +
  labs(title = "Distribution of Peak Lengths Across SEACR Samples") +
  ylim(0,2000) + ylab("Peak Length") + xlab("Data Sample") +
  theme_minimal()

# Boxplot of Signal Intensity ---------------------------------------------

#signal value for broad peaks is in column 7
ggplot(autoData, aes(x = `Data Name`, y = V4)) + geom_boxplot(outlier.shape = NA) +
  labs(title = "Distribution of Peak Signal Intensity Across SEACR Samples")+
  ylim(0,10000) + ylab("Peak Signal Intensity") + xlab("Data Sample") +
  theme_minimal()


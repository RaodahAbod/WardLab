rm(list = ls())
setwd("C:/Users/raoda/Desktop/R Stuff/CUT&Tag Processing")
library(tidyverse)
library(dplyr)
library(ggplot2)

batchIt <- function(namesOfData){
  file <- c()
for (dataset in namesOfData){
  file = read.table(choose.files(), header = FALSE) %>% 
    mutate(`Data Name` = dataset) %>%
    rbind(file, .)
}
  batch <- file
  return(batch)
}
chipseq <- c('Roadmap_LV_H3K27ac','Roadmap_LV_H3K9me3','Roadmap_iPSC20b_H3K27ac',
             'Roadmap_iPSC20b_H3K9me3')
encode <- c('Encode_LV_H3K27ac','Encode_iPSC20b_H3K27ac')
cuttag <- c('Day0 H3K27ac', 'Day0 H3K9me3','Day15 H3K27ac','Day15 H3K9me3','Day30 H3K27ac')

outsourced <- batchIt(chipseq)
encodestuff <- batchIt(encode)
wardSamples <- batchIt(cuttag)


lengthandauto <- function(sample){
  sample <-  sample %>%
    mutate(peakLength = .[,3] - .[,2]) %>%
      subset(!grepl("chrX|chrY|chrM|chrUn",V1)) %>%
return()
}

outsourced1 <- lengthandauto(outsourced)
encode1 <- lengthandauto(encodestuff)
ward1 <- lengthandauto(wardSamples)

encode1 <- select(encode1, V1,V2, V3, V4 = V7, `Data Name`, `peakLength` )
ward1 <- select(ward1, V1, V2, V3, V4,`Data Name`, `peakLength` )
fullSample <- rbind(outsourced1, ward1, encode1)
chipSamples <- rbind(outsourced1, encode1)

#only to run and save autosomal files
saveAutos <- function(autosName){
  ObtainPeakLength <- function(histonemod){
    histonemod <-  histonemod %>%
      mutate(histonemod,peakLength = histonemod[,3] - histonemod[,2]) %>%
      return()
  }
for (autos in autosName){
  origData <- read.table(choose.files(), header = FALSE) %>%
    ObtainPeakLength()
  newFileName <- paste0(autos,"_autosomal_H3K27ac.bed")
  autoData <- subset(origData, !grepl("chrX|chrY|chrM|chrUn", origData$V1))
  write.table(autoData, newFileName ,row.names = F ,col.names = F ,
            sep="\t" , quote=FALSE) 
}
}
saveAutos("EncodeiPS20b")


# Plot Num of Peaks -------------------------------------------------------
#the golden format :)
ggplot(fullSample, aes(x = `Data Name`, fill = `Data Name`)) + geom_bar() +
  #facet_wrap(~Batch, scales = "free_x", ncol = 4) +
  labs(title = "Number of Peaks") +
  geom_text(stat = "count", aes(label = after_stat(count)),
            position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Data Sample") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Boxplot of Peaks --------------------------------------------------------

# ObtainPeakLength <- function(histonemod){
#   histonemod <-  histonemod %>%
#       mutate(histonemod,peakLength = histonemod[,3] - histonemod[,2]) %>%
#   return()
# }
# 
# autoData <- ObtainPeakLength(autoData)

ggplot(fullSample, aes(x = `Data Name`, y = peakLength)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title = "Distribution of Peak Lengths") +
  ylim(0,2000) + ylab("Peak Length") + xlab("Data Sample") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# Boxplot of Signal Intensity ---------------------------------------------

#signal value for broad peaks is in column 7
ggplot(ward1, aes(x = `Data Name`, y = V4)) + geom_boxplot(outlier.shape = NA) +
  labs(title = "Distribution of Peak Signal Intensity") +
  ylab("Peak Signal Intensity") + xlab("Data Sample") + coord_cartesian(ylim = c(0,5000)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  
  

ggplot(chipSamples, aes(x = `Data Name`, y = V4)) + geom_boxplot(outlier.shape = NA) +
  labs(title = "Distribution of Peak Signal Intensity") +
  ylab("Peak Signal Intensity") + xlab("Data Sample") + coord_cartesian(ylim = c(0,30)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

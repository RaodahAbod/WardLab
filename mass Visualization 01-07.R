rm(list = ls())
setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing")

namesOfData <- c("Day0","Day0h","Day1","Day2","Day3","Day4","Day5","Day7","Day10", "Day15","Day30")

file <- c()

for (dataset in namesOfData){
  file = read.table(choose.files(), header = FALSE) %>% 
    mutate(`Data Name` = dataset) %>%
    rbind(file, .)
}

autoData <- subset(file, !grepl("chrX|chrY|chrM|chrUn", file$V1))

#only to run and save autosomal files
for (autos in namesOfData){
  origData <- read.table(choose.files(), header = FALSE)
  newFileName <- paste0(autos,"_Human_SEACR_H3K27ac.bed")
  autoData <- subset(origData, !grepl("chrX|chrY|chrM|chrUn", origData$V1))
  write.table(autoData, newFileName ,row.names = F ,col.names = F ,
            sep="\t" , quote=FALSE) 
}

# Plot Num of Peaks -------------------------------------------------------

ggplot(file, aes(x = `Data Name`, fill = `Data Name`)) + geom_bar() +
  labs(title = "Lengths of Peaks Across SEACR Samples") +
  geom_text(stat = "count", aes(label = after_stat(count)),
            position = position_dodge(width = 0.9), vjust = -0.5) + xlab("Data Sample")

ggplot(autoData, aes(x = `Data Name`, fill = `Data Name`)) + geom_bar() +
  labs(title = "Lengths of Peaks Across SEACR Samples After Autosomal Filtering") +
  geom_text(stat = "count", aes(label = after_stat(count)),
            position = position_dodge(width = 0.9), vjust = -0.5) + xlab("Data Sample")

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


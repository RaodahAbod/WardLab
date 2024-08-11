rm(list = ls())
setwd("C:/Users/raoda/Desktop/R Stuff/CUT&Tag Processing/Cardiotox Sequencing")
library("Rsubread")
library(ComplexHeatmap)
library(dplyr)
library(tibble)
library(tidyverse)
library(edgeR)
library(factoextra)
library(ggrepel)
library(ggfortify)
library(circlize)

#  Rearrange GFF file? --------------------------------------------------------------
loc <- choose.files()
file <- read.table(loc)
name <- basename(loc)
newFile <- select(file,V3, V1,V4,V5,V7,V9)
#fileName <- readline(prompt = "Enter the filename: ")
fileName <- paste0(name,".saf")
write.table(newFile, fileName ,row.names = F ,col.names = F ,
            sep="\t" , quote=FALSE) 

# FeatureCounts ---------------------------------------------------------------------

numOfFiles <- readline(prompt = "Enter the number of bam files you are using: ")
bam_files <- c()
for(c in 1:numOfFiles){
bam_files[c] <- choose.files()
}
# at this point, you can save your R workspace/R environment with the variable
# bam_files to avoid having to reload all files one by one if working with a large
# dataset. 

dataName <- c("87-1 DNR3", "87-1 DNR24",
              "87-1 DOX3", "87-1 DOX24", 
              "87-1 MTX3", "87-1 MTX24",
              "87-1 VEH3", "87-1 VEH24",
              
              "77-1 DNR3", "77-1 DNR24",
              "77-1 DOX3",
              "77-1 EPI3", "77-1 EPI24", 
              "77-1 MTX24",
              "77-1 VEH3", "77-1 VEH24",
              
              "79-1 DNR3", "79-1 DNR24",
              "79-1 DOX3", "79-1 DOX24",
              "79-1 EPI3", "79-1 EPI24", 
              "79-1 MTX3", "79-1 MTX24",
              "79-1 VEH3", "79-1 VEH24",
              
              "78-1 DNR3", "78-1 DNR24",
              "78-1 DOX3", "78-1 DOX24",
              "78-1 EPI3", "78-1 EPI24", 
              "78-1 MTX3", "78-1 MTX24",
              "78-1 VEH3", "78-1 VEH24",
              
              "71-1 DNR3", "71-1 DNR24",
              "71-1 DOX24",
              "71-1 EPI3", "71-1 EPI24", 
              "71-1 MTX3", "71-1 MTX24",
              "71-1 VEH3", "71-1 VEH24") 

characteristics <- data.frame(dataName)
characteristics <- mutate(characteristics, anthracycline =NA, time =NA, 
                          trt=NA, indiv=NA)

# Update the anthracycline column based on the conditions
characteristics <- characteristics %>%
 mutate(anthracycline = ifelse(grepl("VEH", dataName) | grepl("MTX", dataName), 'NO', 'YES'))

characteristics <- characteristics %>%
 mutate(time = ifelse(grepl("3", dataName), '3H', '24H'))

# Update the trt column based on the conditions
characteristics <- characteristics %>%
 mutate(trt = ifelse(grepl("DOX", dataName), 'DOX',
              ifelse(grepl("DNR", dataName), 'DNR',
              ifelse(grepl("EPI", dataName), 'EPI',
              ifelse(grepl("MTX", dataName), 'MTX', 'VEH')))))

# Update the trt column based on the conditions
characteristics <- characteristics %>%
 mutate(indiv = ifelse(grepl("87", dataName), '87-1',
                       ifelse(grepl("77", dataName), '77-1',
                              ifelse(grepl("79", dataName), '79-1',
                                     ifelse(grepl("78", dataName), '78-1', '71-1')))))

annot_file <- choose.files() #load in saf file
for (x in 1:numOfFiles){
featureCountsTrial <- featureCounts(files = bam_files[x], annot.ext = annot_file,
                                    isGTFAnnotationFile = FALSE, isPairedEnd = TRUE,
                                    nthreads = 1, verbose = TRUE)

assign(dataName[x],featureCountsTrial)
}

counts <- data.frame(`87-1 DNR3`[[1]],`87-1 DNR24`[[1]], 
                     `87-1 DOX3`[[1]],`87-1 DOX24`[[1]],
                     `87-1 MTX3`[[1]],`87-1 MTX24`[[1]],
                     `87-1 VEH3`[[1]],`87-1 VEH24`[[1]],
                     
                     `77-1 DNR3`[[1]],`77-1 DNR24`[[1]], 
                     `77-1 DOX3`[[1]],
                     `77-1 EPI3`[[1]],`77-1 EPI24`[[1]], 
                     `77-1 MTX24`[[1]],
                     `77-1 VEH3`[[1]],`77-1 VEH24`[[1]],
                     
                     `79-1 DNR3`[[1]],`79-1 DNR24`[[1]], 
                     `79-1 DOX3`[[1]],`79-1 DOX24`[[1]],
                     `79-1 EPI3`[[1]],`79-1 EPI24`[[1]], 
                     `79-1 MTX3`[[1]],`79-1 MTX24`[[1]],
                     `79-1 VEH3`[[1]],`79-1 VEH24`[[1]],
                     
                     `78-1 DNR3`[[1]],`78-1 DNR24`[[1]], 
                     `78-1 DOX3`[[1]],`78-1 DOX24`[[1]],
                     `78-1 EPI3`[[1]],`78-1 EPI24`[[1]], 
                     `78-1 MTX3`[[1]],`78-1 MTX24`[[1]],
                     `78-1 VEH3`[[1]],`78-1 VEH24`[[1]],
                     
                     `71-1 DNR3`[[1]],`71-1 DNR24`[[1]], 
                     `71-1 DOX24`[[1]],
                     `71-1 EPI3`[[1]],`71-1 EPI24`[[1]], 
                     `71-1 MTX3`[[1]],`71-1 MTX24`[[1]],
                     `71-1 VEH3`[[1]],`71-1 VEH24`[[1]])

colnames(counts) <- dataName
rownames(counts) <- featureLoc$chrLoc

row_means <- rowMeans(counts)
counts_filtered <- counts[row_means > 10,]

#-----
# Define color mappings with named levels
trt_colors <- c(DNR = "#F1B72B", DOX = "#8B006D", EPI = "#DF707E", MTX = "#3386DD", VEH = "#41B333")
indiv_colors <- c(`87-1` = "#1B9E77", `77-1` = "#D95F02", `79-1` = "#7570B3", `78-1` = "#E7298A", `71-1` = "#66A61E", indiv6 = "#E6AB02")
time_colors <- c(`3H` = "pink", `24H` = "chocolate4")
AC_colors <- c(`YES` = "yellow1", `NO` = "darkorange1")

# Define column annotations using named colors
heatChar <- HeatmapAnnotation(
 trt = anno_simple(characteristics$trt, col = trt_colors),
 indiv = anno_simple(characteristics$indiv, col = indiv_colors),
 time = anno_simple(characteristics$time, col = time_colors),
 AC = anno_simple(characteristics$anthracycline, col = AC_colors)
)
 
#these two heatmaps are identical
cpm_highConf <- counts %>% 
 cpm(., log = TRUE) %>% 
 cor(method = "spearman")
htmp3 <- Heatmap(as.matrix(cpm_highConf), width = unit(12, "cm"),column_title = 
                 'CPM Spearman Correlation of Cardiotox Treatment Groups
                   High Confidence Peak Set >= 2 Peaks', 
                top_annotation = heatChar)
htmp3@column_names_param[["gp"]][["fontsize"]] <- 8
htmp3@row_names_param[["gp"]][["fontsize"]] <- 8
htmp3

# Excluding Low peak Count Samples --------------------------------------------------

#this excludes samples with fewer than 200 peaks. those samples are
#77-1 VEH24, 79-1 DOX3, 79-1 EPI3, 79-1 MTX3, 79-1 VEH3
#79-1 VEH24, 78-1 MTX24, 78-1 VEH3, 78-1 VEH24

dataName2 <- c("87-1 DNR3", "87-1 DNR24",
               "87-1 DOX3", "87-1 DOX24", 
               "87-1 MTX3", "87-1 MTX24",
               "87-1 VEH3", "87-1 VEH24",
               
               "77-1 DNR3", "77-1 DNR24",
               "77-1 DOX3",
               "77-1 EPI3", "77-1 EPI24", 
               "77-1 MTX24",
               "77-1 VEH3",
               
               "79-1 DNR3", "79-1 DNR24",
               "79-1 DOX24",
               "79-1 EPI24", 
               "79-1 MTX24",
               
               "78-1 DNR3", "78-1 DNR24",
               "78-1 DOX3", "78-1 DOX24",
               "78-1 EPI3", "78-1 EPI24", 
               "78-1 MTX3",
               
               "71-1 DNR3", "71-1 DNR24",
               "71-1 DOX24",
               "71-1 EPI3", "71-1 EPI24", 
               "71-1 MTX3", "71-1 MTX24",
               "71-1 VEH3", "71-1 VEH24")

highConfChar <- data.frame(dataName2)
highConfChar <- mutate(highConfChar, anthracycline =NA, time =NA, 
                          trt=NA, indiv=NA)

# Update the anthracycline column based on the conditions
highConfChar <- highConfChar %>%
 mutate(anthracycline = ifelse(grepl("VEH", dataName2) | grepl("MTX", dataName2), 'NO', 'YES'))

highConfChar <- highConfChar %>%
 mutate(time = ifelse(grepl("3", dataName2), '3H', '24H'))

# Update the trt column based on the conditions
highConfChar <- highConfChar %>%
 mutate(trt = ifelse(grepl("DOX", dataName2), 'DOX',
                     ifelse(grepl("DNR", dataName2), 'DNR',
                            ifelse(grepl("EPI", dataName2), 'EPI',
                                   ifelse(grepl("MTX", dataName2), 'MTX', 'VEH')))))

# Update the trt column based on the conditions
highConfChar <- highConfChar %>%
 mutate(indiv = ifelse(grepl("87", dataName2), '87-1',
                       ifelse(grepl("77", dataName2), '77-1',
                              ifelse(grepl("79", dataName2), '79-1',
                                     ifelse(grepl("78", dataName2), '78-1', '71-1')))))


counts2 <- data.frame(`87-1 DNR3`[[1]],`87-1 DNR24`[[1]], 
                     `87-1 DOX3`[[1]],`87-1 DOX24`[[1]],
                     `87-1 MTX3`[[1]],`87-1 MTX24`[[1]],
                     `87-1 VEH3`[[1]],`87-1 VEH24`[[1]],
                     
                     `77-1 DNR3`[[1]],`77-1 DNR24`[[1]], 
                     `77-1 DOX3`[[1]],
                     `77-1 EPI3`[[1]],`77-1 EPI24`[[1]], 
                     `77-1 MTX24`[[1]],
                     `77-1 VEH3`[[1]],
                     
                     `79-1 DNR3`[[1]],`79-1 DNR24`[[1]], 
                     `79-1 DOX24`[[1]],
                     `79-1 EPI24`[[1]], 
                     `79-1 MTX24`[[1]],
                     
                     `78-1 DNR3`[[1]],`78-1 DNR24`[[1]], 
                     `78-1 DOX3`[[1]],`78-1 DOX24`[[1]],
                     `78-1 EPI3`[[1]],`78-1 EPI24`[[1]], 
                     `78-1 MTX3`[[1]],
                     
                     `71-1 DNR3`[[1]],`71-1 DNR24`[[1]], 
                     `71-1 DOX24`[[1]],
                     `71-1 EPI3`[[1]],`71-1 EPI24`[[1]], 
                     `71-1 MTX3`[[1]],`71-1 MTX24`[[1]],
                     `71-1 VEH3`[[1]],`71-1 VEH24`[[1]])

colnames(counts2) <- dataName2
rownames(counts2) <- featureLoc$chrLoc

cor_noLowConf <- cor(counts2, method = "spearman") # rounded to 2 decimals

# Define column annotations using named colors
heatChar_highconf <- HeatmapAnnotation(
  trt = anno_simple(highConfChar$trt, col = trt_colors),
  indiv = anno_simple(highConfChar$indiv, col = indiv_colors),
  time = anno_simple(highConfChar$time, col = time_colors),
  AC = anno_simple(highConfChar$anthracycline, col = AC_colors)
)

#-------

htmp2 <- Heatmap(as.matrix(cor_noLowConf), width = unit(12, "cm"),column_title = 
                   'Spearman Correlation of Cardiotox Treatment Groups 
                 Excluding Samples with Low Peak Count
                 High Conf Peak Set >= 3 Peaks', top_annotation = heatChar_highconf,
)
htmp2@column_names_param[["gp"]][["fontsize"]] <- 8
htmp2@row_names_param[["gp"]][["fontsize"]] <- 8
htmp2

# Extracting total 14 samples -------------------------------------------------------
#this excludes samples with fewer than 1.5K peaks. those samples are
#77-1 VEH24, 
#79-1 DOX3, 79-1 DOX24, 79-1 EPI3, 79-1 MTX3, 79-1 VEH3, 79-1 VEH24
#78-1 DNR24,78-1 DOX3, 78-1 EPI24, 78-1 MTX24, 78-1 VEH3, 78-1 VEH24
#71-1 VEH3

dataName3 <- c("87-1 DNR3", "87-1 DNR24",
               "87-1 DOX3", "87-1 DOX24", 
               "87-1 MTX3", "87-1 MTX24",
               "87-1 VEH3", "87-1 VEH24",
               
               "77-1 DNR3", "77-1 DNR24",
               "77-1 DOX3",
               "77-1 EPI3", "77-1 EPI24", 
               "77-1 MTX24",
               "77-1 VEH3",
               
               "79-1 DNR3", "79-1 DNR24",
               "79-1 EPI24", 
               "79-1 MTX24",
               
               "78-1 DNR3", 
               "78-1 DOX24",
               "78-1 EPI3",  
               "78-1 MTX3",
               
               "71-1 DNR3", "71-1 DNR24",
               "71-1 DOX24",
               "71-1 EPI3", "71-1 EPI24", 
               "71-1 MTX3", "71-1 MTX24",
               "71-1 VEH24")

highConfChar2 <- data.frame(dataName3)
highConfChar2 <- mutate(highConfChar2, anthracycline =NA, time =NA, 
                       trt=NA, indiv=NA)

# Update the anthracycline column based on the conditions
highConfChar2 <- highConfChar2 %>%
  mutate(anthracycline = ifelse(grepl("VEH", dataName3) | grepl("MTX", dataName3), 'NO', 'YES'))

highConfChar2 <- highConfChar2 %>%
  mutate(time = ifelse(grepl("3", dataName3), '3H', '24H'))

# Update the trt column based on the conditions
highConfChar2 <- highConfChar2 %>%
  mutate(trt = ifelse(grepl("DOX", dataName3), 'DOX',
                      ifelse(grepl("DNR", dataName3), 'DNR',
                             ifelse(grepl("EPI", dataName3), 'EPI',
                                    ifelse(grepl("MTX", dataName3), 'MTX', 'VEH')))))

# Update the trt column based on the conditions
highConfChar2 <- highConfChar2 %>%
  mutate(indiv = ifelse(grepl("87", dataName3), '87-1',
                        ifelse(grepl("77", dataName3), '77-1',
                               ifelse(grepl("79", dataName3), '79-1',
                                      ifelse(grepl("78", dataName3), '78-1', '71-1')))))


counts3 <- data.frame(`87-1 DNR3`[[1]],`87-1 DNR24`[[1]], 
                      `87-1 DOX3`[[1]],`87-1 DOX24`[[1]],
                      `87-1 MTX3`[[1]],`87-1 MTX24`[[1]],
                      `87-1 VEH3`[[1]],`87-1 VEH24`[[1]],
                      
                      `77-1 DNR3`[[1]],`77-1 DNR24`[[1]], 
                      `77-1 DOX3`[[1]],
                      `77-1 EPI3`[[1]],`77-1 EPI24`[[1]], 
                      `77-1 MTX24`[[1]],
                      `77-1 VEH3`[[1]],
                      
                      `79-1 DNR3`[[1]],`79-1 DNR24`[[1]],
                      `79-1 EPI24`[[1]], 
                      `79-1 MTX24`[[1]],
                      
                      `78-1 DNR3`[[1]],
                      `78-1 DOX24`[[1]],
                      `78-1 EPI3`[[1]], 
                      `78-1 MTX3`[[1]],
                      
                      `71-1 DNR3`[[1]],`71-1 DNR24`[[1]], 
                      `71-1 DOX24`[[1]],
                      `71-1 EPI3`[[1]],`71-1 EPI24`[[1]], 
                      `71-1 MTX3`[[1]],`71-1 MTX24`[[1]],
                      `71-1 VEH24`[[1]])

colnames(counts3) <- dataName3

cor_noLowConf2 <- cor(counts3, method = "spearman") # rounded to 2 decimals


heatChar_minus12 <- HeatmapAnnotation(
  trt = anno_simple(highConfChar2$trt, col = trt_colors),
  indiv = anno_simple(highConfChar2$indiv, col = indiv_colors),
  time = anno_simple(highConfChar2$time, col = time_colors),
  AC = anno_simple(highConfChar2$anthracycline, col = AC_colors)
)
htmp_minus12 <- Heatmap(as.matrix(cor_noLowConf2), width = unit(12, "cm"),column_title = 
                   'Spearman Correlation of Cardiotox Treatment Groups 
                 Excluding 12 Samples with <1.5K Peaks
                 High Conf Peak Set >= 3 Peaks', top_annotation = heatChar_minus12,
)
htmp_minus12@column_names_param[["gp"]][["fontsize"]] <- 8
htmp_minus12@row_names_param[["gp"]][["fontsize"]] <- 8
htmp_minus12

row_means2 <- rowMeans(counts3)
counts_filtered2 <- counts3[row_means2 > 10,]

cor_noLowConf2_filt <- cor(counts_filtered2, method = "spearman") # rounded to 2 decimals

htmp_minus12_filt <- Heatmap(as.matrix(cor_noLowConf2_filt), width = unit(12, "cm"),column_title = 
                          'Spearman Correlation of Cardiotox Treatment Groups 
                 Excluding 12 Samples with <1.5K Peaks
                 Filtered for Mean Peak Counts of >10
                 High Conf Peak Set >= 3 Peaks', top_annotation = heatChar_minus12,
)
htmp_minus12_filt@column_names_param[["gp"]][["fontsize"]] <- 8
htmp_minus12_filt@row_names_param[["gp"]][["fontsize"]] <- 8
htmp_minus12_filt

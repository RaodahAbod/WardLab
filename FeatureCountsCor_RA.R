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

#  Rearrange 'GFF' and change to an SAF file --------------------------------------------------------------
loc <- choose.files()
file <- read.table(loc)
name <- basename(loc)
newFile <- select(file,V3, V1,V4,V5,V7,V9)
#fileName <- readline(prompt = "Enter the filename: ")
fileName <- paste0(name,".saf")
write.table(newFile, fileName ,row.names = F ,col.names = F ,
            sep="\t" , quote=FALSE) 

# FeatureCounts ---------------------------------------------------------------------

# this will prompt you to enter a number
numOfFiles <- readline(prompt = "Enter the number of bam files you are using: ")
bam_files <- c()

# will iterate through the number you specified above and store all your bam file locations
# from your computer. 
for(c in 1:numOfFiles){
bam_files[c] <- choose.files()
}
# at this point, you can save your R workspace/R environment with the variable
# bam_files to avoid having to reload all files one by one if working with a large
# dataset. 


# a master variable containing all the sample names of your data set. 
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

# creates a base variable to store sample metadata/characteristics
# the following few lines will assign values to the columns based on the sample name. 
# i.e: 71-1 DOX 24 anthra col will be YES, time col will be 24H, trt will be DOX and indiv
# will be 5. 

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

# feature counts will run through all samples you have and store everything in its own
# renamed variable using the large dataname variable created above. 
annot_file <- choose.files() #load in saf file
for (x in 1:numOfFiles){
featureCountsTrial <- featureCounts(files = bam_files[x], annot.ext = annot_file,
                                    isGTFAnnotationFile = FALSE, isPairedEnd = TRUE,
                                    nthreads = 1, verbose = TRUE)

assign(dataName[x],featureCountsTrial)
}

# manually assigning all the count information into one big count matrix. 
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

# renames the columns to be the name of each of the samples. for easy identification. 
colnames(counts) <- dataName

# this is an optimization step that can allow us to filter out any lowly expressed regions
# across all samples. 
row_means <- rowMeans(counts)
counts_filtered <- counts[row_means > 10,]
# this filtered count matrix can be passed through the heat map function and be visualized, 
# similar to the other count matrices in this script. 

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
 
# generating the correlation heatmap. 
cpm_highConf <- counts %>% 
 cpm(., log = TRUE) %>% 
 cor(method = "spearman")
corHtmp <- Heatmap(as.matrix(cpm_highConf), width = unit(12, "cm"),column_title = 
                 'CPM Spearman Correlation of Cardiotox Treatment Groups
                   High Confidence Peak Set >= 3 Peaks', 
                top_annotation = heatChar)
corHtmp@column_names_param[["gp"]][["fontsize"]] <- 8
corHtmp@row_names_param[["gp"]][["fontsize"]] <- 8
corHtmp

# Excluding Low peak Count Samples --------------------------------------------------

#this excludes samples with fewer than 200 peaks. those samples are
#77-1 VEH24, 79-1 DOX3, 79-1 EPI3, 79-1 MTX3, 79-1 VEH3
#79-1 VEH24, 78-1 MTX24, 78-1 VEH3, 78-1 VEH24

dataName_lowPeaks <- c("87-1 DNR3", "87-1 DNR24",
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

char_lowPeaks <- data.frame(dataName_lowPeaks)
char_lowPeaks <- mutate(char_lowPeaks, anthracycline =NA, time =NA, 
                          trt=NA, indiv=NA)

# Update the anthracycline column based on the conditions
char_lowPeaks <- char_lowPeaks %>%
 mutate(anthracycline = ifelse(grepl("VEH", dataName_lowPeaks) | 
                                 grepl("MTX", dataName_lowPeaks), 'NO', 'YES'))

char_lowPeaks <- char_lowPeaks %>%
 mutate(time = ifelse(grepl("3", dataName_lowPeaks), '3H', '24H'))

# Update the trt column based on the conditions
char_lowPeaks <- char_lowPeaks %>%
 mutate(trt = ifelse(grepl("DOX", dataName_lowPeaks), 'DOX',
                     ifelse(grepl("DNR", dataName_lowPeaks), 'DNR',
                            ifelse(grepl("EPI", dataName_lowPeaks), 'EPI',
                                   ifelse(grepl("MTX", dataName_lowPeaks), 'MTX', 'VEH')))))

# Update the trt column based on the conditions
char_lowPeaks <- char_lowPeaks %>%
 mutate(indiv = ifelse(grepl("87", dataName_lowPeaks), '87-1',
                       ifelse(grepl("77", dataName_lowPeaks), '77-1',
                              ifelse(grepl("79", dataName_lowPeaks), '79-1',
                                     ifelse(grepl("78", dataName_lowPeaks), '78-1', '71-1')))))


counts_lowPeaks <- data.frame(`87-1 DNR3`[[1]],`87-1 DNR24`[[1]], 
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

colnames(counts_lowPeaks) <- dataName_lowPeaks

cor_lowPeaks <- cor(counts_lowPeaks, method = "spearman") # rounded to 2 decimals

# Define column annotations using named colors
heatChar_lowPeaks <- HeatmapAnnotation(
  trt = anno_simple(char_lowPeaks$trt, col = trt_colors),
  indiv = anno_simple(char_lowPeaks$indiv, col = indiv_colors),
  time = anno_simple(char_lowPeaks$time, col = time_colors),
  AC = anno_simple(char_lowPeaks$anthracycline, col = AC_colors)
)

#-------

htmp_lowPeaks <- Heatmap(as.matrix(cor_lowPeaks), width = unit(12, "cm"),column_title = 
                   'Spearman Correlation of Cardiotox Treatment Groups 
                 Excluding Samples with Low Peak Count
                 High Conf Peak Set >= 3 Peaks', top_annotation = heatChar_lowPeaks,
)
htmp_lowPeaks@column_names_param[["gp"]][["fontsize"]] <- 8
htmp_lowPeaks@row_names_param[["gp"]][["fontsize"]] <- 8
htmp_lowPeaks

# Extracting total 14 samples -------------------------------------------------------
#this excludes samples with fewer than 1.5K peaks. those samples are
#77-1 VEH24, 
#79-1 DOX3, 79-1 DOX24, 79-1 EPI3, 79-1 MTX3, 79-1 VEH3, 79-1 VEH24
#78-1 DNR24,78-1 DOX3, 78-1 EPI24, 78-1 MTX24, 78-1 VEH3, 78-1 VEH24
#71-1 VEH3

dataName_strict <- c("87-1 DNR3", "87-1 DNR24",
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

char_strict <- data.frame(dataName_strict)
char_strict <- mutate(char_strict, anthracycline =NA, time =NA, 
                       trt=NA, indiv=NA)

# Update the anthracycline column based on the conditions
char_strict <- char_strict %>%
  mutate(anthracycline = ifelse(grepl("VEH", dataName_strict) | 
                                  grepl("MTX", dataName_strict), 'NO', 'YES'))

char_strict <- char_strict %>%
  mutate(time = ifelse(grepl("3", dataName_strict), '3H', '24H'))

# Update the trt column based on the conditions
char_strict <- char_strict %>%
  mutate(trt = ifelse(grepl("DOX", dataName_strict), 'DOX',
                      ifelse(grepl("DNR", dataName_strict), 'DNR',
                             ifelse(grepl("EPI", dataName_strict), 'EPI',
                                    ifelse(grepl("MTX", dataName_strict), 'MTX', 'VEH')))))

# Update the trt column based on the conditions
char_strict <- char_strict %>%
  mutate(indiv = ifelse(grepl("87", dataName_strict), '87-1',
                        ifelse(grepl("77", dataName_strict), '77-1',
                               ifelse(grepl("79", dataName_strict), '79-1',
                                      ifelse(grepl("78", dataName_strict), '78-1', '71-1')))))


counts_strict <- data.frame(`87-1 DNR3`[[1]],`87-1 DNR24`[[1]], 
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

colnames(counts_strict) <- dataName_strict

cor_strict <- cor(counts_strict, method = "spearman") # rounded to 2 decimals


heatChar_strict <- HeatmapAnnotation(
  trt = anno_simple(char_strict$trt, col = trt_colors),
  indiv = anno_simple(char_strict$indiv, col = indiv_colors),
  time = anno_simple(char_strict$time, col = time_colors),
  AC = anno_simple(char_strict$anthracycline, col = AC_colors)
)
htmp_strict <- Heatmap(as.matrix(cor_strict), width = unit(12, "cm"),column_title = 
                   'Spearman Correlation of Cardiotox Treatment Groups 
                 Excluding 12 Samples with <1.5K Peaks
                 High Conf Peak Set >= 3 Peaks', top_annotation = heatChar_strict)

htmp_strict@column_names_param[["gp"]][["fontsize"]] <- 8
htmp_strict@row_names_param[["gp"]][["fontsize"]] <- 8
htmp_strict

row_means_strict <- rowMeans(counts_strict)
counts_filtered_strict <- counts_strict[row_means_strict > 10,]

cor_strict_filt <- cor(counts_filtered_strict, method = "spearman") # rounded to 2 decimals

htmp_strict_filt <- Heatmap(as.matrix(cor_strict_filt), width = unit(12, "cm"),column_title = 
                          'Spearman Correlation of Cardiotox Treatment Groups 
                 Excluding 12 Samples with <1.5K Peaks
                 Filtered for Mean Peak Counts of >10
                 High Conf Peak Set >= 3 Peaks', top_annotation = heatChar_strict)

htmp_strict_filt@column_names_param[["gp"]][["fontsize"]] <- 8
htmp_strict_filt@row_names_param[["gp"]][["fontsize"]] <- 8
htmp_strict_filt

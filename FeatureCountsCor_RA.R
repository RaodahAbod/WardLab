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
# annot_file <- readGFFFeatures(choose.files())

numOfFiles <- readline(prompt = "Enter the number of bam files you are using: ")
bam_files <- c()
for(c in 1:numOfFiles){
bam_files[c] <- choose.files()
}


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

chrPos <- data.frame(featureCountsTrial[1:2])
featureLoc <- chrPos %>% select(2:5) %>% unite("chrLoc",2:4,sep = ".")

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
#dim(x_filtered)

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
 
#-------

testit <- cor(counts, method = "spearman") # rounded to 2 decimals
htmp <- Heatmap(as.matrix(testit), width = unit(12, "cm"),column_title = 
                 'Spearman Correlation of Cardiotox Treatment Groups
                High Confidence Peak Set >= 3 Peaks', 
                top_annotation = heatChar)
htmp@column_names_param[["gp"]][["fontsize"]] <- 8
htmp@row_names_param[["gp"]][["fontsize"]] <- 8
htmp


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

#removing rowmeans = 0
cpmFilt_highConf <- counts_filtered %>% 
  cpm(., log = TRUE) %>% 
  cor(method = "spearman")
htmp_filt <- Heatmap(as.matrix(cpmFilt_highConf), width = unit(12, "cm"),column_title = 
                   'CPM Spearman Correlation of Cardiotox Treatment Groups
                   High Confidence Peak Set >= 2 Peaks', 
                 top_annotation = heatChar)
htmp_filt@column_names_param[["gp"]][["fontsize"]] <- 8
htmp_filt@row_names_param[["gp"]][["fontsize"]] <- 8
htmp_filt


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

#-----
# Define color mappings with named levels
trt_colors <- c(DNR = "#F1B72B", DOX = "#8B006D", EPI = "#DF707E", MTX = "#3386DD", VEH = "#41B333")
indiv_colors <- c(`87-1` = "#1B9E77", `77-1` = "#D95F02", `79-1` = "#7570B3", `78-1` = "#E7298A", `71-1` = "#66A61E", indiv6 = "#E6AB02")
time_colors <- c(`3H` = "pink", `24H` = "chocolate4")
AC_colors <- c(`YES` = "yellow1", `NO` = "darkorange1")

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


#####
cpm_highConf <- counts %>% 
                cpm(., log = TRUE) %>% 
                cor(method = "spearman")

highConfCounts <- data.frame(sample_list = colnames(cpm_highConf))

# PCA Plot -----------------------------------------------------------------------------------
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

pca_counts <- cpm(counts, log = TRUE) %>% as.matrix() 
pca_counts <- prcomp(t(counts), scale. = TRUE)

groups <- as.factor(characteristics$indiv)

drugPalette <- c("#F1B72B","#8B006D","#DF707E","#3386DD","#41B333")

pca_counts  %>%
 ggplot(.,aes(x = PC1, y = PC2, col=characteristics$trt,
              shape=characteristics$time, group=groups))+
 geom_point(size= 5)+
 scale_color_manual(values=drug_pal)+
 ggrepel::geom_text_repel(aes(label = groups))+
 labs(col = "Color Group", shape = "Shape Group") +
 #scale_shape_manual(name = "Time",values= c("3h"=0,"24h"=1))+
 ggtitle(expression("log2 cpm PCA plot"))+
 theme_bw()+
 guides(col="none", size =4)+
 #labs(y = "PC 2 (13.81%)", x ="PC 1 (16.58%)")+
 theme(plot.title=element_text(size= 14,hjust = 0.5),
       axis.title = element_text(size = 12, color = "black")+
        geom_label_repel())

pca_counts  %>%
 ggplot(.,aes(x = PC3, y = PC4, col=characteristics$trt,
              shape=characteristics$time, group=groups))+
 geom_point(size= 5)+
 scale_color_manual(values=drug_pal)+
 ggrepel::geom_text_repel(aes(label = groups))+
 #scale_shape_manual(name = "Time",values= c("3h"=0,"24h"=1))+
 ggtitle(expression("log2 cpm PCA plot"))+
 theme_bw()+
 guides(col="none", size =4)+
 #labs(y = "PC 2 (13.81%)", x ="PC 1 (16.58%)")+
 theme(plot.title=element_text(size= 14,hjust = 0.5),
       axis.title = element_text(size = 12, color = "black")+
        geom_label_repel())


groups2 <- as.factor(highConfChar$indiv)
pca_counts2 <- cpm(counts2, log = TRUE) %>% as.matrix() 
pca_counts2 <- prcomp(t(counts2), scale. = TRUE)
pca_counts2  %>%
 ggplot(.,aes(x = PC1, y = PC2, col=highConfChar$trt,
              shape=highConfChar$time, group=groups2))+
 geom_point(size= 5)+
 scale_color_manual(values=drug_pal)+
 ggrepel::geom_text_repel(aes(label = groups2))+
 #scale_shape_manual(name = "Time",values= c("3h"=0,"24h"=1))+
 ggtitle(expression("log2 PCA excluding low confidence peaksets"))+
 theme_bw()+
 guides(col="none", size =4)+
 #labs(y = "PC 2 (13.81%)", x ="PC 1 (16.58%)")+
 theme(plot.title=element_text(size= 14,hjust = 0.5),
       axis.title = element_text(size = 12, color = "black")+
        geom_label_repel())

pca_counts2  %>%
 ggplot(.,aes(x = PC3, y = PC4, col=highConfChar$trt,
              shape=highConfChar$time, group=groups2))+
 geom_point(size= 5)+
 scale_color_manual(values=drug_pal)+
 ggrepel::geom_text_repel(aes(label = groups2))+
 #scale_shape_manual(name = "Time",values= c("3h"=0,"24h"=1))+
 ggtitle(expression("log2 PCA excluding low confidence peaksets"))+
 theme_bw()+
 guides(col="none", size =4)+
 #labs(y = "PC 2 (13.81%)", x ="PC 1 (16.58%)")+
 theme(plot.title=element_text(size= 14,hjust = 0.5),
       axis.title = element_text(size = 12, color = "black")+
        geom_label_repel())



# boxplot ---------------------------------------------------------------------------
conditions <- c("chr8.11696723.11705488", "chr5.173227675.173237356",
                "chr7.5424533.5430372", "chr8.11679200.11684553",
                "chr11.65416268.65424443", "chr1.201376777.201377946",
                "chr1.236689729.236690453", "chr1.236699546.236700750", 
                "chr6.36675879.36685632") #actn2

conditionnames <- c("GATA4-1of2", "NKX2-5", "LINC03073", "GATA4-2of2",
                    "FRMD8","TNNT2","ACTN2-1of2","ACTN2-2of2", "CDKN1A")

cntr <- 1
for(a in conditions){
 cond <- grepl(a, rownames(counts))
 tmp <- counts %>% .[cond, ] %>% 
  mutate(Peak = row.names(.)) %>% 
   mutate(PeakName = conditionnames[cntr]) %>%
  pivot_longer(col = !Peak & !PeakName, names_to = "sample", values_to = "counts") %>% 
  cbind(.,characteristics)
 
 tmp$counts <- cpm(tmp$counts, log = TRUE)
 tmp$time <- factor(tmp$time, levels = c("3H","24H"))
 
 assign(paste0("box",cntr),tmp) 
 
 cntr <- cntr + 1
}

masterbox <- rbind(box1,box2,box3,box4,box5,box6,box7, box8, box9)

ggplot(masterbox, aes(x = time, y=counts))+
  geom_boxplot(aes(fill=trt))+
  geom_jitter(aes(color = indiv),
              position = position_jitter(width = 0.36),
              size = 2) +  # Jittered points with individual colors
  facet_wrap(~PeakName, ncol = 4) +
  labs(y = "log2 CPM", x = 'Time') +
  #ggtitle("top 3 DAR in 3 hour DOX")+
  scale_fill_manual(values = drugPalette) +
  scale_color_manual(values = indiv_colors) # Custom colors for individuals
  #theme_bw() 



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
rownames(counts3) <- featureLoc$chrLoc

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




#####

groups3 <- as.factor(highConfChar2$indiv)
pca_counts3 <- cpm(counts3, log = TRUE) %>% as.matrix() 
pca_counts3 <- prcomp(t(counts3), scale. = TRUE)
summary(pca_counts3)

pca_counts3  %>%
  ggplot(.,aes(x = PC1, y = PC2, col=highConfChar2$trt,
               shape=highConfChar2$time, group=groups3))+
  geom_point(size= 5)+
  scale_color_manual(values=drugPalette)+
  ggrepel::geom_text_repel(aes(label = groups3))+
  #scale_shape_manual(name = "Time",values= c("3h"=0,"24h"=1))+
  ggtitle(expression("log2 PCA excluding low confidence peaksets"))+
  theme_bw()+
  guides(col="none", size =4)+
  labs(y = "PC 2 (8.048%)", x ="PC 1 (68.92%)")+
  theme(plot.title=element_text(size= 14,hjust = 0.5),
        axis.title = element_text(size = 12, color = "black")+
          geom_label_repel())


pca_counts3  %>%
  ggplot(.,aes(x = PC2, y = PC3, col=highConfChar2$trt,
               shape=highConfChar2$time, group=groups3))+
  geom_point(size= 5)+
  scale_color_manual(values=drugPalette)+
  ggrepel::geom_text_repel(aes(label = groups3))+
  #scale_shape_manual(name = "Time",values= c("3h"=0,"24h"=1))+
  ggtitle(expression("log2 PCA excluding low confidence peaksets"))+
  theme_bw()+
  guides(col="none", size =4)+
  labs(x = "PC 2 (8.048%)", y ="PC 3 (4.305%)")+
  theme(plot.title=element_text(size= 14,hjust = 0.5),
        axis.title = element_text(size = 12, color = "black")+
          geom_label_repel())

#####



#cpmFilt_highConf
pca_counts4 <- prcomp(t(cor_noLowConf2_filt), scale. = TRUE)
summary(pca_counts4)

pca_counts4  %>%
  ggplot(.,aes(x = PC1, y = PC2, col=highConfChar2$trt,
               shape=highConfChar2$time, group=groups3))+
  geom_point(size= 5)+
  scale_color_manual(values=drugPalette)+
  ggrepel::geom_text_repel(aes(label = groups3))+
  #scale_shape_manual(name = "Time",values= c("3h"=0,"24h"=1))+
  ggtitle(expression("log2 PCA excluding low confidence peaksets"))+
  theme_bw()+
  guides(col="none", size =4)+
  labs(y = "PC 2 (23.11%)", x ="PC 1 (33.91%)")+
  theme(plot.title=element_text(size= 14,hjust = 0.5),
        axis.title = element_text(size = 12, color = "black")+
          geom_label_repel())


pca_counts4  %>%
  ggplot(.,aes(x = PC2, y = PC3, col=highConfChar2$trt,
               shape=highConfChar2$time, group=groups3))+
  geom_point(size= 5)+
  scale_color_manual(values=drugPalette)+
  ggrepel::geom_text_repel(aes(label = groups3))+
  #scale_shape_manual(name = "Time",values= c("3h"=0,"24h"=1))+
  ggtitle(expression("log2 PCA excluding low confidence peaksets"))+
  theme_bw()+
  guides(col="none", size =4)+
  labs(x = "PC 2 (23.11%)", y ="PC 3 (11.76%)")+
  theme(plot.title=element_text(size= 14,hjust = 0.5),
        axis.title = element_text(size = 12, color = "black")+
          geom_label_repel())



# DE genes for renees rna -----------------------------------------------------------

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


# cdkn1a, sort1, znf740, slc28a3, 
# rmi1, frs2, ppip5k2, celsr1, 
# gramd4, stxbp3, cnih1, Fermt2, yeats4 

DElocs_1of3 <- c("chr6.36675879.36685632","chr1.109392819.109393385",
                 "chr12.53180110.53181129","chr9.84342981.84343462")

DElocs_2of3 <- c("chr9.83979656.83980586", "chr12.69470347.69470926",
                 "chr5.103119908.103120803", "chr22.46377548.46378224")

DElocs_3of3 <- c("chr22.46577318.46578137", "chr1.108747097.108747410",
                 "chr14.54441009.54441400", "chr14.52950398.52951422",
                 "chr12.69404924.69405766")

DE_1of3 <- c("cdkn1a, sort1, znf740, slc28a3")

DE_2of3 <- c("rmi1, frs2, ppip5k2, celsr1")

DE_3of3 <- c("gramd4, stxbp3, cnih1, Fermt2, yeats4")

listDElocs <- list(DElocs_1of3,DElocs_2of3,DElocs_3of3)
listDEs <- list(DE_1of3, DE_2of3, DE_3of3)

cntr <- 1
for(a in DElocs_3of3){
  cond <- grepl(a, rownames(counts))
  tmp <- counts %>% .[cond, ] %>% 
    mutate(Peak = row.names(.)) %>% 
    mutate(PeakName = DE_3of3[cntr]) %>%
    pivot_longer(col = !Peak & !PeakName, names_to = "sample", values_to = "counts") %>% 
    cbind(.,characteristics)
  
  tmp$counts <- cpm(tmp$counts, log = TRUE)
  tmp$time <- factor(tmp$time, levels = c("3H","24H"))
  
  assign(paste0("set3_",cntr),tmp) 
  cntr <- cntr + 1
}

masterbox1 <- rbind(set1_1,set1_2, set1_3, set1_4)
masterbox2 <- rbind(set2_1,set2_2,set2_3,set2_4)
masterbox3 <- rbind(set3_1,set3_2,set3_3,set3_4,set3_5)

indiv_colors <- c(`87-1` = "#1B9E77", `77-1` = "#D95F02", 
                  `79-1` = "#7570B3", `78-1` = "#E7298A", 
                  `71-1` = "#66A61E")
drugPalette <- c("#F1B72B","#8B006D","#DF707E","#3386DD","#41B333")

ggplot(masterbox1, aes(x = time, y=counts))+
  geom_boxplot(aes(fill=trt))+
  geom_jitter(aes(color = indiv),
              position = position_jitter(width = 0.36),
              size = 2) +  # Jittered points with individual colors
  facet_wrap(~PeakName, ncol = 4) +
  labs(y = "log2 CPM", x = 'Time') +
  #ggtitle("top 3 DAR in 3 hour DOX")+
  scale_fill_manual(values = drugPalette) +
  scale_color_manual(values = indiv_colors) # Custom colors for individuals
#theme_bw() 

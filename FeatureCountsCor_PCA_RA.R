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

# extracts the chromosome position into a variable. This will be helpful
# in performing the counts boxplot later in the script. 
# it will aid in locating the nearest gene to the peak chromosome location
chrPos <- data.frame(featureCountsTrial[1:2])
featureLoc <- chrPos %>% select(2:5) %>% unite("chrLoc",2:4,sep = ".")

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
#renames the rows as the chromosome location rather than 'feature', which is nonspecific
rownames(counts) <- featureLoc$chrLoc


# this is an optimization step that can allow us to filter out any lowly expressed regions
# across all samples. row_means <- rowMeans(counts)
counts_filtered <- counts[row_means > 10,]
# this filtered count matrix can be passed through the heat map function and be visualized, 
# similar to the other count matrices in this script. 

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
cpm_cor_3pks <- counts %>% 
  cpm(., log = TRUE) %>% 
  cor(method = "spearman")
htmp_3pks <- Heatmap(as.matrix(cpm_cor_3pks), width = unit(12, "cm"),column_title = 
                 'Spearman Correlation of Cardiotox Treatment Groups
                High Confidence Peak Set >= 3 Peaks', 
                top_annotation = heatChar)
htmp_3pks@column_names_param[["gp"]][["fontsize"]] <- 8
htmp_3pks@row_names_param[["gp"]][["fontsize"]] <- 8
htmp_3pks

# PCA Plot -----------------------------------------------------------------------------------
pca_counts <- cpm(counts, log = TRUE) %>% as.matrix() 
pca_counts <- prcomp(t(counts), scale. = TRUE)

groups <- as.factor(characteristics$indiv)

drug_pal <- c("#F1B72B","#8B006D","#DF707E","#3386DD","#41B333")

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

# boxplot ---------------------------------------------------------------------------
conditions <- c("chr8.11696723.11705488", "chr5.173227675.173237356",
                "chr7.5424533.5430372", "chr8.11679200.11684553",
                "chr11.65416268.65424443", "chr1.201376777.201377946",
                "chr1.236689729.236690453", "chr1.236699546.236700750", 
                "chr6.36675879.36685632")

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

masterbox <- rbind(box1,box2,box3,box4,box5,box6,box7,box8,box9)

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


# DE genes from renees rna -----------------------------------------------------------

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

DE_1of3 <- c("CDKN1A", "SORT1", "ZNF740", "SLC28A3")

DE_2of3 <- c("RMI1", "FRS2", "PPIP5K2", "CELSR1")

DE_3of3 <- c("GRAMD4", "STXBP3", "CNIH1", "FERMT2", "YEATS4")

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

masterbox <- rbind(masterbox1, masterbox2, masterbox3)

indiv_colors <- c(`87-1` = "#1B9E77", `77-1` = "#D95F02", 
                  `79-1` = "#7570B3", `78-1` = "#E7298A", 
                  `71-1` = "#66A61E")
drugPalette <- c("#F1B72B","#8B006D","#DF707E","#3386DD","#41B333")

ggplot(masterbox1, aes(x = time, y=counts))+
  geom_boxplot(aes(fill=trt))+
  geom_jitter(aes(color = indiv),
              position = position_jitter(width = 0.36),
              size = 2) +  # Jittered points with individual colors
  facet_wrap(~PeakName, ncol = 5) +
  labs(y = "log2 CPM", x = 'Time') +
  #ggtitle("top 3 DAR in 3 hour DOX")+
  scale_fill_manual(values = drugPalette) +
  scale_color_manual(values = indiv_colors) # Custom colors for individuals
#theme_bw() 


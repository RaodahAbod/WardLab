rm(list = ls())
setwd("C:/Users/raoda/Desktop/R Stuff/CUT&Tag Processing")
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(readr)

# load in go enrich tibbles ---------------------------------------------------------

H3K27ac_Accutase_H21792 <- read_csv(choose.files())
H3K27ac_Trypsin_H21792 <- read_csv(choose.files())
H3K27ac_Accutase_78_1 <- read_csv(choose.files())
H3K27ac_Scraped_78_1 <- read_csv(choose.files())
H3K27ac_ENCODE_LV <- read_csv(choose.files())

# Building off the Cardiac Only Terms. Top 15 for Each. Except iPSC ---------------------
# this works. this is what i use. 

sortAndExtractHeartOnly <- function(dataSet){ 
 dataSet <- dataSet %>%
  filter(grepl("cardi|heart", Description, ignore.case = TRUE,)) %>%
  arrange(p.adjust) %>%
  slice(1:15) %>%
  select(Description, p.adjust)
 return(dataSet)
}

newAccutase_H21792 <- dplyr::select(H3K27ac_Accutase_H21792, Description, `Sample 1 p.adjust` = `p.adjust`)
newTrypsin_H21792 <- dplyr::select(H3K27ac_Trypsin_H21792, Description,`Sample 2 p.adjust` = `p.adjust`)
newAccutase_781 <- dplyr::select(H3K27ac_Accutase_78_1, Description, `Sample 3 p.adjust` = `p.adjust`)
newScraped_781 <- dplyr::select(H3K27ac_Scraped_78_1, Description, `Sample 4 p.adjust` = `p.adjust`)
newEncode <- dplyr::select(H3K27ac_ENCODE_LV, Description, `Sample 5 p.adjust` = `p.adjust`)

topTermsOnly <- rbind(sortAndExtractHeartOnly(H3K27ac_Accutase_H21792),
                      sortAndExtractHeartOnly(H3K27ac_Trypsin_H21792),
                      sortAndExtractHeartOnly(H3K27ac_Accutase_78_1),
                      sortAndExtractHeartOnly(H3K27ac_Scraped_78_1),
                      sortAndExtractHeartOnly(H3K27ac_ENCODE_LV)) %>%
  dplyr::select("Description")

masterTopHeartGOs <- topTermsOnly %>% 
  left_join(newAccutase_H21792, join_by(Description)) %>% 
  left_join(newTrypsin_H21792, join_by(Description)) %>% 
  left_join(newAccutase_781, join_by(Description)) %>% 
  left_join(newScraped_781, join_by(Description)) %>%
  left_join(newEncode, join_by(Description)) %>%
  unique()

colnames(masterTopHeartGOs) <- c("Description","H21792 H3K27ac Accutase", 
                                 "H21792 H3K27ac Trypsin", 
                                 "78-1 H3K27ac Accutase",
                                 "78-1 H3K27ac Scraped",
                                 "LV H3K27ac ENCODE")                           

mth_terms <- masterTopHeartGOs$Description
masterTopHeartGOs_noID <- subset(masterTopHeartGOs, select = -c(Description))
masterTopHeartGOs_noID <- as.data.frame(masterTopHeartGOs_noID)
rownames(masterTopHeartGOs_noID) <- mth_terms
logTopHeartGOs <- -log10(masterTopHeartGOs_noID)

value_colors <- c("NA" = "gray")
htmp <- Heatmap(as.matrix(logTopHeartGOs), na_col = value_colors)
htmp@row_names_param[["gp"]][["fontsize"]] <- 8
draw(htmp,column_title = "Heatmap of Top 15 Terms Containing `Heart` or `Cardiac`")

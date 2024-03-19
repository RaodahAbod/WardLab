rm(list = ls())
setwd("C:/Users/raoda/Desktop/R Stuff/CUT&Tag Processing")
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(readr)

# load in go enrich tibbles ---------------------------------------------------------

ipsc_ENCODE_enrichGO <- read_csv(choose.files())
Day0_enrichGO <- read_csv(choose.files())
Day15_enrichGO <- read_csv(choose.files())
Day30_enrichGO <- read_csv(choose.files())
lv_ENCODE_enrichGO <- read_csv(choose.files())
H21792_acc_enrichGO <- read_csv(choose.files())
ROA_scr_enrichGO <- read.csv(choose.files())
ROA_acc_enrichGO <- read_csv(choose.files())

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

newDay0 <- dplyr::select(Day0_enrichGO, Description, `Day 0 p.adjust` = `p.adjust`)
newDay15 <- dplyr::select(Day15_enrichGO, Description,`Day 15 p.adjust` = `p.adjust`)
newDay30 <- dplyr::select(Day30_enrichGO, Description, `Day 30 p.adjust` = `p.adjust`)
newENCODELV <- dplyr::select(lv_ENCODE_enrichGO, Description, `ENCODE p.adjust` = `p.adjust`)
newENCODEiPS <- dplyr::select(ipsc_ENCODE_enrichGO, Description, `ENCODE p.adjust` = `p.adjust`)
newROAscr <- dplyr::select(ROA_scr_enrichGO, Description, `ENCODE p.adjust` = `p.adjust`)
newROAacc <- dplyr::select(ROA_acc_enrichGO, Description, `ENCODE p.adjust` = `p.adjust`)
newH21792acc <- dplyr::select(H21792_acc_enrichGO, Description, `ENCODE p.adjust` = `p.adjust`)


topTermsOnly <- rbind(sortAndExtractHeartOnly(Day0_enrichGO),
                      sortAndExtractHeartOnly(Day15_enrichGO),
                      sortAndExtractHeartOnly(Day30_enrichGO),
                      sortAndExtractHeartOnly(lv_ENCODE_enrichGO),
                      sortAndExtractHeartOnly(ipsc_ENCODE_enrichGO),
                      sortAndExtractHeartOnly(ROA_scr_enrichGO),
                      sortAndExtractHeartOnly(ROA_acc_enrichGO),
                      sortAndExtractHeartOnly(H21792_acc_enrichGO)) %>%
  dplyr::select("Description")

masterTopHeartGOs <- topTermsOnly %>% 
  left_join(newDay0, join_by(Description)) %>% 
  left_join(newDay15, join_by(Description)) %>% 
  left_join(newDay30, join_by(Description)) %>% 
  left_join(newENCODELV, join_by(Description)) %>%
  left_join(newENCODEiPS, join_by(Description)) %>%
  left_join(newROAscr, join_by(Description)) %>%
  left_join(newROAacc, join_by(Description)) %>%
  left_join(newH21792acc, join_by(Description)) %>%
  unique()

colnames(masterTopHeartGOs) <- c("Description","Day 0 H3K27ac", "Day 15 H3K27ac", 
                                 "Day 30 H3K27ac", "ENCODE LV H3K27ac", 
                                 "ENCODE iPS H3K27ac",
                                 "ROA 78-1 Scr", 'ROA 78-1 Acc', 'H21792 Acc')                           

mth_terms <- masterTopHeartGOs$Description
masterTopHeartGOs_noID <- subset(masterTopHeartGOs, select = -c(Description))
masterTopHeartGOs_noID <- as.data.frame(masterTopHeartGOs_noID)
rownames(masterTopHeartGOs_noID) <- mth_terms
logTopHeartGOs <- -log10(masterTopHeartGOs_noID)

value_colors <- c("NA" = "gray")
htmp <- Heatmap(as.matrix(logTopHeartGOs), na_col = value_colors)
htmp@row_names_param[["gp"]][["fontsize"]] <- 8
draw(htmp,column_title = "Heatmap of Top 15 Terms Containing `Heart` or `Cardiac`")

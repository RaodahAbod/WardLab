rm(list = ls())
setwd("C:/Users/raoda/Desktop/R Stuff/CUT&Tag Processing")
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(gplots)

# load in go enrich tibbles ---------------------------------------------------------

ipsc_enrichGO <- read_csv(choose.files())
ipsc_cm_enrichGO <- read_csv(choose.files())
cm_enrichGO <- read_csv(choose.files())
lv_ENCODE_enrichGO <- read_csv(choose.files())

# Building off the Cardiac Only Terms. Top 12 for Each. Except iPSC ---------------------
# this works well enough, but the top 12 for ipsc skews the heatmap a bit. 

sortAndExtractHeartOnly <- function(dataSet){ 
 dataSet <- dataSet %>%
  filter(grepl("cardi|heart", Description, ignore.case = TRUE)) %>%
  arrange(p.adjust) %>%
  slice(1:15) %>%
  select(Description, p.adjust)
 return(dataSet)
}

newDay0 <- select(ipsc_enrichGO, Description, `Day 0 p.adjust` = `p.adjust`)
newDay15 <- select(ipsc_cm_enrichGO, Description,`Day 15 p.adjust` = `p.adjust`)
newDay30 <- select(cm_enrichGO, Description, `Day 30 p.adjust` = `p.adjust`)
newENCODE <- select(lv_ENCODE_enrichGO, Description, `ENCODE p.adjust` = `p.adjust`)

topTermsOnly <- rbind(sortAndExtractHeartOnly(ipsc_enrichGO),
                      sortAndExtractHeartOnly(ipsc_cm_enrichGO),
                      sortAndExtractHeartOnly(cm_enrichGO),
                      sortAndExtractHeartOnly(lv_ENCODE_enrichGO)) %>%
                select("Description")

masterTopHeartGOs <- topTermsOnly %>% 
  left_join(newDay0, join_by(Description)) %>% 
  left_join(newDay15, join_by(Description)) %>% 
  left_join(newDay30, join_by(Description)) %>% 
  left_join(newENCODE, join_by(Description)) %>%
  unique()

colnames(masterTopHeartGOs) <- c("Description","Day 0 H3K27ac", "Day 15 H3K27ac", 
                                 "Day 30 H3K27ac", "ENCODE LV H3K27ac")                           

mth_terms <- masterTopHeartGOs$Description
masterTopHeartGOs_noID <- subset(masterTopHeartGOs, select = -c(Description))
masterTopHeartGOs_noID <- as.data.frame(masterTopHeartGOs_noID)
rownames(masterTopHeartGOs_noID) <- mth_terms
logTopHeartGOs <- -log10(masterTopHeartGOs_noID)

value_colors <- c("NA" = "gray")
htmp <- Heatmap(as.matrix(logTopHeartGOs), na_col = value_colors)
htmp@row_names_param[["gp"]][["fontsize"]] <- 10
draw(htmp,column_title = "Heatmap of Top 15 Terms Containing `Heart` or `Cardiac`")

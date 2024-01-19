rm(list = ls())
setwd("C:/Users/raoda/Desktop/R Stuff/CUT&Tag Processing")
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

# load in go enrich tibbles ---------------------------------------------------------

ipsc_enrichGO <- read_csv(choose.files())
ipsc_cm_enrichGO <- read_csv(choose.files())
cm_enrichGO <- read_csv(choose.files())
lv_ENCODE_enrichGO <- read_csv(choose.files())

# Prepping Tibble Data and Merging to create master GO list --------------------------

newDay0 <- select(ipsc_enrichGO, Description, `Day 0 p.adjust` = `p.adjust`)
newDay15 <- select(ipsc_cm_enrichGO, Description,`Day 15 p.adjust` = `p.adjust`)
newDay30 <- select(cm_enrichGO, Description, `Day 30 p.adjust` = `p.adjust`)
newENCODE <- select(lv_ENCODE_enrichGO, Description, `ENCODE p.adjust` = `p.adjust`)

#idk if i use this tbh
masterCombinedGOs <- newDay0 %>% 
 full_join(newDay15) %>% 
 full_join(newDay30) %>% 
 full_join(newENCODE)


# Just Cardiac and Heart Terms ------------------------------------------------------

heartOnlyGOs <- masterCombinedGOs %>%
 filter(grepl("cardi|heart|ventricle", Description, ignore.case = TRUE))
#noNAs <- na.omit(heartOnlyGOs) this takes out alot of data...

heartOnlyTerms <- heartOnlyGOs$Description
heartOnlyGOs_noID <- subset(heartOnlyGOs, select = -c(Description))
colnames(heartOnlyGOs_noID) <- c("iPSC H3K27ac", "iPSC-CM H3K27ac", 
                                 "CM H3K27ac", "ENCODE LV H3K27ac")
heartOnlyGOs_noID <- as.data.frame(heartOnlyGOs_noID)
rownames(heartOnlyGOs_noID) <- heartOnlyTerms

logHeartGOs <- -log10(heartOnlyGOs_noID)
logHeartGOs[is.na(logHeartGOs)] <- 1

logHeartGOs_1 <- logHeartGOs %>% slice(1:33)
logHeartGOs_2 <- logHeartGOs %>% slice(34:66)
logHeartGOs_3 <- logHeartGOs %>% slice(67:98)

#produces useful information, however, there is still a significant amount 
#of NAs in the visualization

#heat_col <- colorRamp2(c(min(logHeartGOs_1), max(logHeartGOs_1)), c("black", "white"))
#htmp <- Heatmap(as.matrix(logHeartGOs_1), col = heat_col)

htmp <- Heatmap(as.matrix(logHeartGOs_1))
htmp@row_names_param[["gp"]][["fontsize"]] <- 10
draw(htmp,column_title = "Heatmap of Top Heart Terms (1/3) Containing `Heart`, `Cardiac`, or `Ventricle`")

# heat_col <- colorRamp2(c(min(logHeartGOs_2), max(logHeartGOs_2)), c("white", "red"))
# htmp2 <- Heatmap(as.matrix(logHeartGOs_2), col = heat_col)

htmp2 <- Heatmap(as.matrix(logHeartGOs_2))
htmp2@row_names_param[["gp"]][["fontsize"]] <- 10
draw(htmp2,column_title = "Heatmap of Top Heart Terms (2/3) Containing `Heart`, `Cardiac`, or `Ventricle`")

# heat_col <- colorRamp2(c(min(logHeartGOs_3), max(logHeartGOs_3)), c("black", "white"))
# htmp3 <- Heatmap(as.matrix(logHeartGOs_3), col = heat_col)

htmp3 <- Heatmap(as.matrix(logHeartGOs_3), na_col = "grey")
htmp3@row_names_param[["gp"]][["fontsize"]] <- 10
draw(htmp3,column_title = "Heatmap of Top Heart Terms (3/3) Containing `Heart`, `Cardiac`, or `Ventricle`")

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

sortAndExtract <- function(dataSet) {
 dataSet <- dataSet %>%
  arrange(p.adjust) %>%
  slice(1:15) %>%
  select(Description, p.adjust)
 return(dataSet)
}

newDay0 <- select(ipsc_enrichGO, Description, `Day 0 p.adjust` = `p.adjust`)
newDay15 <- select(ipsc_cm_enrichGO, Description,`Day 15 p.adjust` = `p.adjust`)
newDay30 <- select(cm_enrichGO, Description, `Day 30 p.adjust` = `p.adjust`)
newENCODE <- select(lv_ENCODE_enrichGO, Description, `ENCODE p.adjust` = `p.adjust`)

topTermsOnly <- rbind(sortAndExtract(ipsc_enrichGO),
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

colnames(masterTopHeartGOs) <- c("Description","iPSC H3K27ac", "iPSC-CM H3K27ac", 
                                 "CM H3K27ac", "ENCODE LV H3K27ac")                           

mth_terms <- masterTopHeartGOs$Description
masterTopHeartGOs_noID <- subset(masterTopHeartGOs, select = -c(Description))
masterTopHeartGOs_noID <- as.data.frame(masterTopHeartGOs_noID)
rownames(masterTopHeartGOs_noID) <- mth_terms

logTopHeartGOs <- -log10(masterTopHeartGOs_noID)
#logTopHeartGOs[is.na(logTopHeartGOs)] <- 1

# heat_col <- colorRamp2(c(min(logTopHeartGOs), max(logTopHeartGOs)), c("black", "white"))
# htmp <- Heatmap(as.matrix(logTopHeartGOs), col = heat_col)

htmp <- Heatmap(as.matrix(logTopHeartGOs), na_col = "black")
htmp@row_names_param[["gp"]][["fontsize"]] <- 10
draw(htmp,column_title = "Heatmap of Top 15 Terms Containing `Heart` or `Cardiac`
     (Includes iPSC Top General Terms)")
htmp

# Heatmap! -----------------------------------------
#not being used in this analysis
extractedTerms <- masterCombinedGOs$Terms
masterCombinedGOs_noID <- subset(masterCombinedGOs, select = -c(Terms,ID))
colnames(masterCombinedGOs_noID) <- c("Manual GSEA", "EnrichGO GSEA", "iPSC-CM H3K27ac", 
                                      "CM H3K27ac", "ENCODE LV")
masterCombinedGOs_noID <- as.data.frame(masterCombinedGOs_noID)
rownames(masterCombinedGOs_noID) <- extractedTerms

logMasterGOs <- -log10(masterCombinedGOs_noID)
logMasterGOs[is.na(logMasterGOs)] <- 1
pheatmap(as.matrix(logMasterGOs))
htmp <- Heatmap(as.matrix(logMasterGOs))
htmp@row_names_param[["gp"]][["fontsize"]] <- 4
htmp
draw(htmp,column_title = "Heatmap of All Extracted GSEA Terms")


# Nother take at some cardiac terms -------------------------------------------------
#not being used in this current analysis

# Create a new data frame including only rows where 'Description' contains 'Cardi', 
# 'development', 'contraction', 'heart', 'chamber', 'ventricle','band','morphogenesis',
# 'circulatory', 
selectedGOs <- masterCombinedGOs %>%
 filter(grepl("cardi|development|contraction|heart|chamber|ventricle|band|morphogenesis
               |circulatory", Description, ignore.case = TRUE))


# Selecting top 10 terms from Each  -------------------------------------------------
# this is actually terrible. not the greatest data visualization

sortAndExtract <- function(dataSet) {
 dataSet <- dataSet %>%
  arrange(p.adjust) %>%
  slice(1:12) %>%
  select(Description, p.adjust)
 return(dataSet)
}
sortedDay0 <- sortAndExtract(ipsc_enrichGO)
sortedDay15 <- sortAndExtract(ipsc_cm_enrichGO)
sortedDay30 <- sortAndExtract(cm_enrichGO)
sortedLV <- sortAndExtract(lv_ENCODE_enrichGO)

masterTopGOs <- sortedDay0 %>% 
 full_join(sortedDay15, join_by(Description)) %>% 
 full_join(sortedDay30, join_by(Description)) %>% 
 full_join(sortedLV, join_by(Description))




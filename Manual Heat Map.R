rm(list = ls())
setwd("C:/Users/raoda/Desktop/R Stuff/CUT&Tag Processing")
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(readr)

# load in go enrich tibbles ---------------------------------------------------------

H3K27me3_Accutase <- read_csv(choose.files())
H3K27me3_Trypsin <- read_csv(choose.files())
H3K27me3_FlashFroz <- read_csv(choose.files())
H3K27me3_Scraped <- read_csv(choose.files())
H3K27me3_ENCODE <- read_csv(choose.files())

# Building off the Cardiac Only Terms. Top 15 for Each. Except iPSC ---------------------
# this works. this is what i use. 

sortAndExtractHeartOnly <- function(dataSet){ 
  dataSet <- dataSet %>%
    filter(grepl("development|differentiation|morphogenesis",
                 Description, ignore.case = TRUE,)) %>%
    arrange(p.adjust) %>%
    slice(1:15) %>%
    select(Description, p.adjust)
  return(dataSet)
}

newAccutase <- dplyr::select(H3K27me3_Accutase, Description, `Sample 1 p.adjust` = `p.adjust`)
newTrypsin <- dplyr::select(H3K27me3_Trypsin, Description,`Sample 2 p.adjust` = `p.adjust`)
newFlashFroz <- dplyr::select(H3K27me3_FlashFroz, Description, `Sample 3 p.adjust` = `p.adjust`)
newScraped <- dplyr::select(H3K27me3_Scraped, Description, `Sample 4 p.adjust` = `p.adjust`)
newENCODE <- dplyr::select(H3K27me3_ENCODE, Description, `Sample 5 p.adjust` = `p.adjust`)

topTermsOnly <- rbind(sortAndExtractHeartOnly(H3K27me3_Accutase),
                      sortAndExtractHeartOnly(H3K27me3_Trypsin),
                      sortAndExtractHeartOnly(H3K27me3_FlashFroz),
                      sortAndExtractHeartOnly(H3K27me3_Scraped),
                      sortAndExtractHeartOnly(H3K27me3_ENCODE)) %>%
  dplyr::select("Description")

masterTopHeartGOs <- topTermsOnly %>% 
  left_join(newAccutase, join_by(Description)) %>% 
  left_join(newTrypsin, join_by(Description)) %>% 
  left_join(newFlashFroz, join_by(Description)) %>% 
  left_join(newScraped, join_by(Description)) %>%
  left_join(newENCODE, join_by(Description)) %>%
  unique()

colnames(masterTopHeartGOs) <- c("Description","78-1 H3K27me3 Accutase", 
                                 "78-1 H3K27me3 Trypsin", 
                                 "78-1 H3K27me3 Flash Frozen",
                                 "78-1 H3K27me3 Scraped",
                                 "LV ENCODE H3K27me3")                           

mth_terms <- masterTopHeartGOs$Description
masterTopHeartGOs_noID <- subset(masterTopHeartGOs, select = -c(Description))
masterTopHeartGOs_noID <- as.data.frame(masterTopHeartGOs_noID)
rownames(masterTopHeartGOs_noID) <- mth_terms
logTopHeartGOs <- -log10(masterTopHeartGOs_noID)

value_colors <- c("NA" = "gray")
htmp <- Heatmap(as.matrix(logTopHeartGOs), na_col = value_colors)
htmp@row_names_param[["gp"]][["fontsize"]] <- 8
draw(htmp,column_title = "Heatmap of Top Terms Containing 
`Development`,`Differentiation`, or `Morphogenesis`")

# For top General Terms Only ---------------------
# this works. this is what i use. 

sortAndExtract <- function(dataSet){ 
  dataSet <- dataSet %>%
    #filter(grepl("development|differentiation|morphogenesis",
    #            Description, ignore.case = TRUE,)) %>%
    arrange(p.adjust) %>%
    slice(1:20) %>%
    select(Description, p.adjust)
  return(dataSet)
}

# newAccutase <- dplyr::select(H3K27me3_Accutase, Description, `Sample 1 p.adjust` = `p.adjust`)
# newTrypsin <- dplyr::select(H3K27me3_Trypsin, Description,`Sample 2 p.adjust` = `p.adjust`)
# newFlashFroz <- dplyr::select(H3K27me3_FlashFroz, Description, `Sample 3 p.adjust` = `p.adjust`)
# newScraped <- dplyr::select(H3K27me3_Scraped, Description, `Sample 4 p.adjust` = `p.adjust`)
# newENCODE <- dplyr::select(H3K27me3_ENCODE, Description, `Sample 5 p.adjust` = `p.adjust`)

topGenTermsOnly <- rbind(sortAndExtract(H3K27me3_Accutase),
                      sortAndExtract(H3K27me3_Trypsin),
                      sortAndExtract(H3K27me3_FlashFroz),
                      sortAndExtract(H3K27me3_Scraped),
                      sortAndExtract(H3K27me3_ENCODE)) %>%
  dplyr::select("Description")

masterTopGOs <- topGenTermsOnly %>% 
  left_join(newAccutase, join_by(Description)) %>% 
  left_join(newTrypsin, join_by(Description)) %>% 
  left_join(newFlashFroz, join_by(Description)) %>% 
  left_join(newScraped, join_by(Description)) %>%
  left_join(newENCODE, join_by(Description)) %>%
  unique()

colnames(masterTopGOs) <- c("Description","78-1 H3K27me3 Accutase", 
                                 "78-1 H3K27me3 Trypsin", 
                                 "78-1 H3K27me3 Flash Frozen",
                                 "78-1 H3K27me3 Scraped",
                                 "LV ENCODE H3K27me3")                           

masterGen_terms <- masterTopGOs$Description
masterTopGOs_noID <- subset(masterTopGOs, select = -c(Description))
masterTopGOs_noID <- as.data.frame(masterTopGOs_noID)
rownames(masterTopGOs_noID) <- masterGen_terms
logTopGOs <- -log10(masterTopGOs_noID)

value_colors <- c("NA" = "gray")
htmp_gen <- Heatmap(as.matrix(logTopGOs), na_col = value_colors)
htmp_gen@row_names_param[["gp"]][["fontsize"]] <- 10
draw(htmp_gen,column_title = "Heatmap of Top Terms Containing 
`Development`,`Differentiation`, or `Morphogenesis`")

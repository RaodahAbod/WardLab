rm(list = ls())
setwd("C:/Users/raoda/Desktop/R Stuff/CUT&Tag Processing/Original Bed Files")

library(tidyverse)
library(dplyr)
library(readr)

fileName <- basename(choose.files()) %>% paste0("UCSC_",.)

read.table(choose.files()) %>% select(1:3) %>% write_tsv(file = fileName, col_names = FALSE, 
                                                         quote = "none")


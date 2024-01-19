rm(list = ls())
setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing")

library(tidyverse)
library(dplyr)
library(viridis)

sampleList = c("Day0 K27ac","Day15 K27ac", "Day0 K9me3", "Day15 K9me3")

#Reorder to IgG, Human0, Human2, Human5, Human15, Chimp0, Chimp2, Chimp5, Chimp15
#Using the color scheme below
CUTTag_Colors <- c("#883268","#3E7274","#94C470","#C03830")

histList = c("H3K27ac","H3K9me3")

fragLen = c()
for(hist in sampleList){
        histInfo = histList
        fragLen = read.table(choose.files(), header = FALSE) %>% 
                mutate(fragLen = V1 %>% as.numeric, 
                       fragCount = V2 %>% as.numeric, 
                       Weight = as.numeric(V2)/sum(as.numeric(V2)), 
                       Histone = histInfo[1], sampleInfo = hist) %>% rbind(fragLen, .) 
}

fragLen$sampleInfo = factor(fragLen$sampleInfo, levels = sampleList)
fragLen$Histone = factor(fragLen$Histone, levels = histList)

fragViolinPlot <- fragLen %>%
        ggplot(aes(x = sampleInfo, y = fragLen, weight = Weight, fill = sampleInfo)) +
        geom_violin(bw = 5, alpha = 0.8) +
        scale_y_continuous(breaks = seq(0, 800, 50)) +
        scale_fill_manual(values = CUTTag_Colors) +  # Use scale_fill_manual with your custom colors
        theme_bw(base_size = 20) +
        ylim(0, 1000) +
        ylab("Fragment Length") +
        xlab("")

fragViolinPlot

fragHistDist = fragLen %>% ggplot(aes(x = fragLen, y = fragCount, color = sampleInfo, 
                               group = sampleInfo, linetype = Histone)) +
        geom_line(size = 1) +
        scale_fill_manual(values = CUTTag_Colors) +  # Use scale_fill_manual with your custom colors        theme_bw(base_size = 20) +
        xlab("Fragment Length") +
        ylab("Count") +
        ylim(0,45000) +
        coord_cartesian(xlim = c(0, 500))
fragHistDist

fragHistDist <- fragLen %>% 
        ggplot(aes(x = fragLen, y = fragCount, color = sampleInfo, group = sampleInfo, linetype = Histone)) +
        geom_line(size = 1) +
        scale_color_manual(values = CUTTag_Colors) +  # Use scale_color_manual for line colors
        scale_linetype_manual(values = c("solid", "dashed")) +  # Use scale_linetype_manual for line types
        theme_bw(base_size = 20) +
        xlab("Fragment Length") +
        ylab("Count") +
        ylim(0, 20000) +
        coord_cartesian(xlim = c(0, 500))

fragHistDist




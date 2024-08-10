rm(list = ls())
setwd("C:/Users/raoda/Desktop/R Stuff/CUT&Tag Processing")
library(readr)
library(tidyverse)
library(dplyr)
library(ggplot2)


# Load Data -------------------------------------------------------------------------

stats <- data.frame(read.csv(choose.files()))
stats$Duplication.Percent <-stats$Duplication.Percent*100

# Add in columns for ggplot characteristics -----------------------------------------
stats <- stats %>%
 mutate(time = ifelse(grepl("3", stats$Sample), '3H', '24H'))

# Update the trt column based on the conditions
stats <- stats %>%
 mutate(trt = ifelse(grepl("DOX", stats$Sample), 'DOX',
                     ifelse(grepl("DNR", stats$Sample), 'DNR',
                            ifelse(grepl("EPI", stats$Sample), 'EPI',
                                   ifelse(grepl("MTX", stats$Sample), 'MTX', 'VEH')))))

# Update the trt column based on the conditions
stats <- stats %>%
 mutate(indiv = ifelse(grepl("87", stats$Sample), '87-1',
                       ifelse(grepl("77", stats$Sample), '77-1',
                              ifelse(grepl("79", stats$Sample), '79-1',
                                     ifelse(grepl("78", stats$Sample), '78-1', '71-1')))))

# Determining other annotation features such as order or color ----------------------
stats$indiv <- factor(stats$indiv, levels = c("87-1","77-1","79-1","78-1","71-1"))
stats$time <- factor(stats$time, levels = c("3H","24H"))
trt_colors <- c(DNR = "#F1B72B", DOX = "#8B006D", EPI = "#DF707E", MTX = "#3386DD", VEH = "#41B333")


# Plotting read info --------------------------------------------------------------------------

stats_reads_long <- stats %>%
 select(Demultiplexed.Reads, Concordantly.Aligned,Estimated.Lib.Size.after.Dedupl, Sample, indiv, trt) %>%
 pivot_longer(cols = c(Demultiplexed.Reads, Concordantly.Aligned, Estimated.Lib.Size.after.Dedupl), 
              names_to = "Event",
              values_to = "Reads")
desired_read_order <- c("Demultiplexed.Reads", "Concordantly.Aligned", "Estimated.Lib.Size.after.Dedupl")
stats_reads_long$Event <- factor(stats_reads_long$Event, levels = desired_read_order)

stats_percent_long <- stats %>%
 select(Percent.Alignment, Duplication.Percent,FRiP.Score, Sample, indiv,trt) %>%
 pivot_longer(cols = c(Percent.Alignment, Duplication.Percent, FRiP.Score), 
              names_to = "Event",
              values_to = "Percent")
desired_percent_order <- c("Percent.Alignment", "Duplication.Percent","FRiP.Score")
stats_percent_long$Event <- factor(stats_percent_long$Event, levels = desired_percent_order)


ggplot(stats_reads_long, aes(x = Sample, y = Reads, fill = trt , color = Event)) +
 geom_bar(stat = "identity", position = "dodge", size=1) +
 facet_wrap(~indiv,scales = "free", ncol=3) +
 scale_fill_manual(values=trt_colors) +
 scale_color_manual(values=c("darkgreen","cyan","black")) +
 #geom_text(aes(label = `Reads`), position = position_dodge(width = 1),
 #          vjust = -0.5, hjust = 0.5, size = 2.5) +
 labs(title = "Mapping Statistics Per Sample Grouped By Individual",
      x = "Sample", y = "Read Count") + 
 theme(axis.text.x = element_blank())

ggplot(stats_percent_long, aes(x = Sample, y = Percent, fill = trt, color = Event)) +
 geom_bar(stat = "identity",position = position_dodge2(width = 0.9, preserve = "single"), 
          size=1) +
 facet_wrap(~indiv, scales = "free_x", ncol=2) +
 scale_fill_manual(values=trt_colors) +
 scale_color_manual(values = c("darkred","violet","darkorange")) +
 labs(title = "Alignment and Duplication Percentage, and FRiP Score grouped by Individual",
      x = "Sample", y = "Percentage") + 
 coord_cartesian(ylim=c(0,105))+
 geom_text(aes(label = floor(`Percent`)), position = position_dodge2(width = 0.9, preserve = "single"),
           vjust = -0.5, hjust = 0.5, size = 3.25, color="black") +
 theme(axis.text.x = element_blank())












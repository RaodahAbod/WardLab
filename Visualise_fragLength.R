rm(list = ls())
setwd("C:/Users/raoda/Desktop/R Stuff/CUT&Tag Processing")
library(tidyverse)
library(dplyr)
library(viridis)


# Create a variable that stores the names for each individual sample.
# I have separated mine by individual to use as a means of grouping during
# final plotting

sampleList_87 <- c("87-1 DNR3", "87-1 DNR24",
                   "87-1 DOX3", "87-1 DOX24", 
                   "87-1 MTX3", "87-1 MTX24",
                   "87-1 VEH3", "87-1 VEH24")

sampleList_77 <- c("77-1 DNR3", "77-1 DNR24",
                    "77-1 DOX3",
                    "77-1 EPI3", "77-1 EPI24", 
                    "77-1 MTX24",
                    "77-1 VEH3", "77-1 VEH24")

sampleList_79 <- c("79-1 DNR3", "79-1 DNR24",
                    "79-1 DOX3", "79-1 DOX24",
                    "79-1 EPI3", "79-1 EPI24", 
                    "79-1 MTX3", "79-1 MTX24",
                    "79-1 VEH3", "79-1 VEH24")

sampleList_78 <- c("78-1 DNR3", "78-1 DNR24",
                   "78-1 DOX3", "78-1 DOX24",
                   "78-1 EPI3", "78-1 EPI24", 
                   "78-1 MTX3", "78-1 MTX24",
                   "78-1 VEH3", "78-1 VEH24")

sampleList_71 <- c("71-1 DNR3", "71-1 DNR24",
                   "71-1 DOX24",
                   "71-1 EPI3", "71-1 EPI24", 
                   "71-1 MTX3", "71-1 MTX24",
                   "71-1 VEH3", "71-1 VEH24")

# the following loops have the same setup
# create an empty variable where you will store your loaded metadata
fragLen_87 <- c()

# loop will iterate the variable 'sample' for each element in your sample name variable
for(sample in sampleList_87){
        histInfo = "H3K27ac"
        # read.table() will open a window every time for each new file. 
        # this will continue for the entire length of your variable.
        fragLen_87 = read.table(choose.files(), header = FALSE) %>% 
                mutate(fragLen = V1 %>% as.numeric, 
                       fragCount = V2 %>% as.numeric, 
                       Weight = as.numeric(V2)/sum(as.numeric(V2)), 
                       Histone = histInfo[1], sampleInfo = sample, Indv = '87-1') %>% 
                rbind(fragLen_87, .) 
}

fragLen_77 <- c()
for(sample in sampleList_77){
        histInfo = "H3K27ac"
        fragLen_77 = read.table(choose.files(), header = FALSE) %>% 
                mutate(fragLen = V1 %>% as.numeric, 
                       fragCount = V2 %>% as.numeric, 
                       Weight = as.numeric(V2)/sum(as.numeric(V2)), 
                       Histone = histInfo[1], sampleInfo = sample, Indv = '77-1') %>% 
                rbind(fragLen_77, .) 
}

fragLen_79 <- c()
for(sample in sampleList_79){
        histInfo = "H3K27ac"
        fragLen_79 = read.table(choose.files(), header = FALSE) %>% 
                mutate(fragLen = V1 %>% as.numeric, 
                       fragCount = V2 %>% as.numeric, 
                       Weight = as.numeric(V2)/sum(as.numeric(V2)), 
                       Histone = histInfo[1], sampleInfo = sample, Indv = '79-1') %>% 
                rbind(fragLen_79, .) 
}

fragLen_78 <- c()
for(sample in sampleList_78){
        histInfo = "H3K27ac"
        fragLen_78 = read.table(choose.files(), header = FALSE) %>% 
                mutate(fragLen = V1 %>% as.numeric, 
                       fragCount = V2 %>% as.numeric, 
                       Weight = as.numeric(V2)/sum(as.numeric(V2)), 
                       Histone = histInfo[1], sampleInfo = sample, Indv = '78-1') %>% 
                rbind(fragLen_78, .) 
}

fragLen_71 <- c()
for(sample in sampleList_71){
 histInfo = "H3K27ac"
 fragLen_71 = read.table(choose.files(), header = FALSE) %>% 
  mutate(fragLen = V1 %>% as.numeric, 
         fragCount = V2 %>% as.numeric, 
         Weight = as.numeric(V2)/sum(as.numeric(V2)), 
         Histone = histInfo[1], sampleInfo = sample, Indv = '71-1') %>% 
  rbind(fragLen_71, .) 
}

# combines all individuals together into one master variable. 
fragLen <- rbind(fragLen_87,fragLen_77,fragLen_79,fragLen_78,fragLen_71)

# transforms the master variable in a form more friendly to manipulation. 
holder <- as.data.frame(fragLen$sampleInfo)

# creating a separate variable to hold explicit sample metadata such as time and treatment
# for aid in final plotting. 
test <- holder %>% separate(col = `fragLen$sampleInfo`,
                            into = c(NA, "treatment"), sep = " ")
final <- test %>% mutate(treat = ifelse(grepl("DOX", test$treatment), 'DOX',
                       ifelse(grepl("DNR", test$treatment), 'DNR',
                              ifelse(grepl("EPI", test$treatment), 'EPI',
                                     ifelse(grepl("MTX", test$treatment), 'MTX', 'VEH'))))) %>%
                  mutate(time = ifelse(grepl("3", test$treatment), '3H', '24H'))

# combines variable holding metadata with the master variable holding fragment
# length information
fragLen <- cbind(fragLen, final)

# factoring will put metadata in the order I desire. It will order characteristics based
# on sequence. i.e: it will be ordered by individual first, then treatment, then time. 
fragLen$Indv <- factor(fragLen$Indv, levels = c("87-1","77-1","79-1","78-1","71-1"))
fragLen$treatment <- factor(fragLen$treatment, 
                            levels = c("DNR3","DOX3","EPI3","MTX3","VEH3",
                                       "DNR24","DOX24","EPI24","MTX24","VEH24"))
fragLen$time <- factor(fragLen$time, levels = c("3H","24H"))

# specified colors correlating to my metadata
trt_colors <- c(DNR = "#F1B72B", DOX = "#8B006D", EPI = "#DF707E", MTX = "#3386DD", VEH = "#41B333")


#ggplot plotting fragment length distribution. Tunable variable is fill, facet_wrap, and 
# scale fill manual. tuning other variables may alter the outcome of the ggplot figure. 
fragViolinPlot <- fragLen %>%
        ggplot(aes(x = time, y = fragLen, weight = Weight, fill = treat)) +
        geom_violin(bw = 5, alpha = 0.8) +
        facet_wrap(~Indv, scales = "free_x", ncol = 5) +
        scale_y_continuous(breaks = seq(0, 800, 50)) +
        theme_bw(base_size = 20) +
        scale_fill_manual(values=trt_colors)+
        ylim(0, 1000) +
        ylab("Fragment Length [bp]") +
        xlab("Treatment and Time") + theme_grey()

# this angles the text labels on the x axis
fragViolinPlot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
 #theme(legend.position = "none") +
 ggtitle("Fragment Length of Carditoxicity Libraries Grouped by Individual") 



fragHistDist = fragLen %>% ggplot(aes(x = fragLen, y = fragCount, color = sampleInfo, 
                               group = sampleInfo, linetype = Histone)) +
        geom_line(size = 1) +
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




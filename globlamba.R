library(ggplot2)
library(dplyr)

stats <- read.csv(choose.files(), header = TRUE) %>%
 pivot_longer(cols = c(`Original.M2Narrow`,`Original.M2Broad`,`Original.SEACR`, `New`),
              names_to = "Peak Caller", values_to = "Total Peaks")

stats[stats == 0] <- NA

h3k27ac <- subset(stats, Histone == 'H3K27ac')
h3k27me3 <- subset(stats, Histone == 'H3K27me3')
h3k9me3 <- subset(stats, Histone == 'H3K9me3')

ggplot(h3k27ac, aes(x = Sample.Name, y = `Total Peaks`, fill = `Peak Caller`)) +
 geom_bar(stat = "identity", position = position_dodge(width = 1)) +
 facet_wrap(~Sample.Name, scales = "free", ncol = 4) +
 geom_text(aes(label = `Total Peaks`), position = position_dodge(width = 1),
           vjust = 0.1, hjust = 0.5, size = 3.5) + 
 labs(title = "Comparison of Peak Count with Different Peak Callers - H3K27ac",
      x = "Sample", y = "Read Count") + 
 theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(h3k27me3, aes(x = Sample.Name, y = `Total Peaks`, fill = `Peak Caller`)) +
 geom_bar(stat = "identity", position = position_dodge(width = 1)) +
 facet_wrap(~Sample.Name, scales = "free", ncol = 2) +
 geom_text(aes(label = `Total Peaks`), position = position_dodge(width = 1),
           vjust = 0.1, hjust = 0.5) + 
 labs(title = "Comparison of Peak Count with Different Peak Callers - H3K27me3",
      x = "Sample", y = "Read Count") + 
 theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(h3k9me3, aes(x = Sample.Name, y = `Total Peaks`, fill = `Peak Caller`)) +
 geom_bar(stat = "identity", position = position_dodge(width = 1)) +
 facet_wrap(~Sample.Name, scales = "free", ncol = 3) +
 geom_text(aes(label = `Total Peaks`), position = position_dodge(width = 1),
           vjust = 0.1, hjust = 0.5, size = 3.5) + 
 labs(title = "Comparison of Peak Count with Different Peak Callers - H3K9me3",
      x = "Sample", y = "Read Count") + 
 theme(axis.text.x = element_text(angle = 45, hjust = 1))


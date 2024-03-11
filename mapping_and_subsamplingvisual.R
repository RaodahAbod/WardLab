rm(list = ls())
setwd("C:/Users/raoda/Desktop/R Stuff/CUT&Tag Processing")
library(tidyverse)
library(dplyr)

subsampling <- read.csv(choose.files())
filename <- choose.files() #renee
subsampling$Subsampling <- as.character(subsampling$Subsampling)

data_long <- subsampling %>%
             pivot_longer(cols = c(Bowtie2.Alignment.Reads, 
                                   Unique.Reads, 
                                   Peak.Count), 
                          names_to = "Metric", values_to = "Quantity")

# Reorder the levels of the Sample
# Define the desired order of levels
desired_order <- c("Bowtie2.Alignment.Reads", 
                   "Unique.Reads",
                   "Peak.Count")
subsample_order <- c("3.3mil","3mil","2.5mil","2mil","1.5mil","1mil","500k")

# Convert 'Sample' to a factor with the desired order
data_long$Metric <- factor(data_long$Metric, levels = desired_order)
data_long$Subsampling.Depth <- factor(data_long$Subsampling.Depth, levels = subsample_order)

data_long$Subsampling <- as.character(data_long$Subsampling)

dedupl_dataLong <- subset(data_long,data_long$Deduplicated. == "Yes")
noDedupl_dataLong <- subset(data_long, data_long$Deduplicated. == 'No')


ggplot(subsampling, aes(x = Sample, y = Paired.Reads, fill = Subsampling)) +
 geom_bar(stat = "identity", position = position_dodge(width = 1)) +
 facet_wrap(~Deduplicated., scales = "free_x") +
 geom_text(aes(label = Paired.Reads), position = position_dodge(width = 1),
           vjust = -0.5, hjust = 0.5) + 
 labs(title = "Subsampling Statistics Grouped by Metric Type", x = "Sample") + 
 scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
 theme(axis.text.x = element_text(angle = 45, hjust = 1))
print("Raoda is awesome!!")
#this one
ggplot(data_long, aes(x = Sample, y = Quantity, fill = Subsampling.Depth)) +
 geom_bar(stat = "identity", position = position_dodge(width = 1)) +
 facet_grid(Metric~Sample, scales = "free") +
 geom_text(aes(label = Quantity), position = position_dodge(width = 1),
           vjust = -0.5, hjust = 0.5, size = 5) + 
 labs(title = "Subsampling Statistics Grouped by Metric Type and Subsampling Depth", 
      x = "Sample") +
 scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
 theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dedupl_dataLong, aes(x = Sample, y = Quantity, fill = Subsampling)) +
 geom_bar(stat = "identity", position = position_dodge(width = 1)) +
 facet_grid(Metric~Deduplicated., scales = "free") +
 geom_text(aes(label = Quantity), position = position_dodge(width = 1),
           vjust = -0.5, hjust = 0.5, size = 3) + 
 labs(title = "Subsampling Paired Read Statistics", x = "Sample", y = "Read Count") +
 scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
 theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(noDedupl_dataLong, aes(x = Sample, y = Quantity, fill = Subsampling)) +
 geom_bar(stat = "identity", position = position_dodge(width = 1)) +
 facet_grid(Metric~Deduplicated., scales = "free") +
 geom_text(aes(label = Quantity), position = position_dodge(width = 1),
           vjust = -0.5, hjust = 0.5, size = 3) + 
 labs(title = "Subsampling Paired Read Statistics", x = "Sample", y = "Read Count") +
 scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
 theme(axis.text.x = element_text(angle = 45, hjust = 1))


# mapping stats ---------------------------------------------------------------------

mapping <- read.csv(choose.files())
mapping <- select(mapping,-Mapped)
mapping <- mutate(mapping, 'Concordantly Aligned %' = Concordantly.Aligned/Before.Alignment*100,
                 'After Filter %' = After.Filter/Before.Alignment*100,
                 'After Dedupl %' = After.Dedupl/Before.Alignment*100)
mapping_percentages <- select(mapping,Sample, Batch,Percent.Alignment, Duplication.Percent,
                              `Concordantly Aligned %`, `After Filter %`, `After Dedupl %`)
mapping_numbers <- select(mapping, Sample, Batch, Concordantly.Aligned, After.Filter, After.Dedupl)


long_mapping_percent <- pivot_longer(mapping_percentages,
                                     cols = c(Percent.Alignment, Duplication.Percent,
                                              `Concordantly Aligned %`, `After Filter %`, 
                                              `After Dedupl %`), 
                                     names_to = "Event", values_to = "Percentage")
long_mapping_percent$Percentage <- signif(long_mapping_percent$Percentage, 4) 

long_mapping_numbers <- pivot_longer(mapping_numbers,
                                     cols = c(Concordantly.Aligned, After.Filter, After.Dedupl),
                                     names_to = "Event", values_to = "Value")

# Define the desired order of levels
desired_order <- c("Duplication.Percent","Percent.Alignment", 'Concordantly Aligned %',
                   'Mapping Quality %', 'Remaining After Dedupl %')
# Convert 'Sample' to a factor with the desired order
long_mapping_percent$Event <- factor(long_mapping_percent$Event, levels = desired_order)

# Define the desired order of levels
desired_order <- c("Concordantly.Aligned", "After.Filter", "After.Dedupl")
# Convert 'Sample' to a factor with the desired order
long_mapping_numbers$Event <- factor(long_mapping_numbers$Event, levels = desired_order)

ggplot(long_mapping_percent, aes(x = Sample, y = Percentage, fill = Event)) +
 geom_bar(stat = "identity", position = position_dodge(width = 1)) +
 facet_grid(~Batch, scales = "free") +
 geom_text(aes(label = Percentage), position = position_dodge(width = 1),
           vjust = -0.5, hjust = 0.5) + 
 labs(title = "Subsampling Statistics Grouped by Metric Type and Deduplication Status", 
      x = "Sample") +
 scale_y_continuous(expand = expansion(mult = c(0, 0.3))) + 
 theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(long_mapping_numbers, aes(x = Sample, y = Value, fill = Event)) +
 geom_bar(stat = "identity", position = position_dodge(width = 1)) +
 facet_grid(~Batch, scales = "free") +
 geom_text(aes(label = Value), position = position_dodge(width = 1),
           vjust = -0.25, hjust = 0.5, size = 3.25) + 
 labs(title = "Subsampling Statistics Grouped by Metric Type and Deduplication Status", 
      x = "Sample") +
 scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
 theme(axis.text.x = element_text(angle = 45, hjust = 1))

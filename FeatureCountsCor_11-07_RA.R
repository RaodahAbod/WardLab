rm(list = ls())
setwd("C:/Users/raoda/OneDrive/Desktop/R Stuff/CUT&Tag Processing")
library("Rsubread")
library(ComplexHeatmap)

#  Rearrange GFF file? --------------------------------------------------------------
file <- read.table(choose.files())
newFile <- select(file,V3, V1,V4,V5,V7,V9)
fileName <- readline(prompt = "Enter the filename: ")
fileName <- paste0(fileName,".saf")
write.table(newFile, fileName ,row.names = F ,col.names = F ,
            sep="\t" , quote=FALSE) 

# FeatureCounts ---------------------------------------------------------------------
# annot_file <- readGFFFeatures(choose.files())
numOfFiles <- readline(prompt = "Enter the number of bam files you are using: ")
B_or_S <- readline(prompt = "Broad or SEACR?: ")
dataName <- c("Day0","Day0h","Day1","Day2","Day3","Day4","Day5","Day7","Day10","Day15","Day30") 
annot_file <- choose.files() #load in saf file

bam_files <- c()
for(c in 1:numOfFiles){
bam_files[c] <- choose.files()
}

for (x in 1:11){
featureCountsTrial <- featureCounts(files = bam_files[x], annot.ext = annot_file,
                                    isGTFAnnotationFile = FALSE, isPairedEnd = TRUE,
                                    nthreads = 1, verbose = TRUE)
newName <- paste(dataName[x],B_or_S)
assign(newName,featureCountsTrial)

}
counts <- data.frame(`Day0 S`[[1]],`Day0h S`[[1]], `Day1 S`[[1]],
                     `Day2 S`[[1]], `Day3 S`[[1]], `Day4 S`[[1]],
                     `Day5 S`[[1]],`Day7 S`[[1]],`Day10 S`[[1]],
                     `Day15 S`[[1]],`Day30 S`[[1]])
colnames(counts) <- c("Day0","Day0h","Day1","Day2","Day3","Day4","Day5","Day7","Day10","Day15","Day30")


testit <- cor(counts, method = "spearman") # rounded to 2 decimals
htmp <- Heatmap(as.matrix(testit), width = unit(12, "cm"))
htmp

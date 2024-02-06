rm(list = ls())
setwd("C:/Users/raoda/Desktop/R Stuff/CUT&Tag Processing")

library(ChIPseeker)
library(tibble)
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")
library(clusterProfiler)
library(DOSE)

# Load ChIPSeek File ----------------------------------------------------------------
loadFile_peakCall <- function(){
 file <- file.choose()
 file <- readPeakFile(file, header = FALSE)
 return(file)
}

peakFile <- loadFile_peakCall()
working_file <- readline(prompt = "Enter the name of the file you just selected: ")

# Protocol 1: ChIPseeker and epigenomic dataset prep ------------------------------------------
prepGRangeObj <- function(seek_object){
 seek_object$Peaks = seek_object$V4
 seek_object$level = seek_object$V5
 seek_object$V4 = seek_object$V5 = NULL
 return(seek_object)
}

peakFile <- prepGRangeObj(peakFile)

# Protocol 2: Annotation of Epigenomic datasets----------------------------------------------------------
TxDb_hg38 = TxDb.Hsapiens.UCSC.hg38.knownGene

Epi_list <- GRangesList(tibble_PeakFile = peakFile)

peakAnnoList <- lapply(Epi_list, annotatePeak, tssRegion = c(-2000,2000), TxDb = TxDb_hg38)
anno_genes_peaks <- as.data.frame(peakAnnoList$tibble_PeakFile)$geneId

# Protocol 4: Visualiztion of annotated results ----------------------------------------------
tibble_anno_enrichGO3 <- enrichGO(gene = anno_genes_peaks, OrgDb = "org.Hs.eg.db",
                                  ont = "ALL" , pvalueCutoff = 1, qvalueCutoff = 1,
                                  pAdjustMethod = 'none', minGSSize = 5,
                                  maxGSSize = 500)

# Save As Tibble --------------------------------------------------------------------
# View the top enriched GO terms as a tibble
GO_Tibble <- tibble_anno_enrichGO3 %>% as_tibble()
GO_Tibble
# Arrange by adjusted p value
GO_Tibble <- GO_Tibble %>%
 arrange(desc(p.adjust))
# Factor by P.adjusted value
GO_Tibble$Description <- factor(GO_Tibble$Description, 
                                levels = GO_Tibble$Description[order(GO_Tibble$p.adjust)])

setwd("C:/Users/raoda/Desktop/R Stuff/CUT&Tag Processing/Gene Ontology Enrichment Tibbles")
#enrichment_fileName <- readline(prompt = "Enter the EnrichGO Tibble filename: ")
enrichment_fileName <- paste0(working_file,"_EnrichmentTibble_p1.csv")
write.csv(GO_Tibble, file = enrichment_fileName)







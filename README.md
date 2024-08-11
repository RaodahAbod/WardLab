# WardLab

This folder (hopes to) consist of all files I have used in my CUT&Tag computation workflow. This majorly includes optimization and QC scripts. 

The most important/beneficial scripts are:

- FeatureCountsCor: generates featurecounts heat map and pca plot using a count matrix
- Manual Heat Map: uses a csv of extracted enriched GO terms (from the ChIPseeker_TibbleScript) to generate a list of enriched terms which is then portrayed in a heatmap
- Visualize_fragLength: uses a file generated from the 9th column of a .sam file to visualize the fragment length distribution in sequenced CUT&Tag Libraries
- _viusalization files: these are used to portray many sequencing statistics or peak calling metrics such as peak intensity/length, etc. 

I hope these scripts will be of benefit to you and serve as a ground zero to inspire more QC and optimization workflows!

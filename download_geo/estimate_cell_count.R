library(minfi)
library(GEOquery)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(FlowSorted.Blood.450k)
library(data.table)

# Sample
sample <- "GSE51032"

# Select idat
idat_file <- "GSE51032/idat/normal_mama"

# Open idat files as array
rgSet <- read.metharray.exp(idat_file)

# Quality control
## calculate the detection p-values
detP <- detectionP(rgSet)

# remove poor quality samples
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]

# Estimate cell count
## "Bcell", "CD4T", "CD8T", "Eos", "Gran", "Mono", "Neu", and "NK" 
cell_estimation <- estimateCellCounts(rgSet, compositeCellType = "Blood",
                   processMethod = "auto", probeSelect = "auto",
                   cellTypes = c("CD8T","CD4T","NK","Bcell","Mono","Gran"),
                   referencePlatform = c("IlluminaHumanMethylation450k",
                                         "IlluminaHumanMethylationEPIC",
                                         "IlluminaHumanMethylation27k"),
                   returnAll = FALSE, meanPlot = FALSE, verbose = TRUE)

cell_estimation_path <- paste(sample, "/", sample, "_bvalues.csv", sep="")
write.csv(file=cell_estimation_path, x=cell_estimation)
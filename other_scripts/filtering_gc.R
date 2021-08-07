# https://kasperdanielhansen.github.io/genbioconductor/html/minfi.html
# examples from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72308

# Packages needed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("minfi")
# BiocManager::install("GEOquery")
# BiocManager::install("IlluminaHumanMethylation450kmanifest")
# BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

library(minfi)
library(GEOquery)

library(knitr)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)

# Breast cancer datasets
# samples <- c("GSE72245", "GSE72251", "GSE72254")

sample <- "GSE72245"

raw <- paste(sample, "/", sample, "_RAW.tar", sep="")
idat_file <- paste(sample, "/", "idat", sep="")

rgSet <- read.metharray.exp(idat_file)
gc()
grSet <- preprocessQuantile(rgSet)

# Differential methylation analysis
# calculate the detection p-values
detP <- detectionP(rgSet)

# remove poor quality samples
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]

# remove poor quality samples from detection p-value table
detP <- detP[,keep]

# normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessQuantile(rgSet) 
# preprocessQuantile() is more suited for datasets where you do not expect 
# global differences between your samples, for example a single tissue

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
# table(keep) # 21821 removed
mSetSqFlt <- mSetSq[keep,]

# if your data includes males and females, remove probes on the sex chromosomes
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% 
                                                        c("chrX","chrY")])
# table(keep) # 11027 removed
mSetSqFlt <- mSetSqFlt[keep,]

# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt) # 15185 removed

# exclude cross reactive probes 
xReactiveProbes <- read.csv(file="GSE72245/48639-non-specific-probes-Illumina450k_part1.csv"
, stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
mSetSqFlt <- mSetSqFlt[keep,] # 26021 removed

# export result as csv file
bVals <- getBeta(mSetSqFlt)
filtered_output <- paste(sample, "/", sample, "_filtered_bvalues.csv", sep="")
write.csv(file=filtered_output, x=bVals)

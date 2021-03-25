library(minfi)
library(GEOquery)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(data.table)

# Sample
sample <- "GSE51032"

# Select idat
idat_file <- paste("GSE51032/idat/normal_mama", sep="")

# Open idat files as array
rgSet <- read.metharray.exp(idat_file)

# Quality control
## calculate the detection p-values
detP <- detectionP(rgSet)

# remove poor quality samples
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]

# Normalization
## normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessFunnorm(rgSet)
rm(idat_file)
rm(rgSet)
gc()

# Filtering
## ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

## remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
mSetSqFlt <- mSetSq[keep,]
rm(mSetSq)
rm(detP)
gc()

## remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)

## exclude cross reactive probes 
xReactiveProbes <- read.csv(file=paste("48639-non-specific-probes-Illumina450k_part1.csv", sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
mSetSqFlt <- mSetSqFlt[keep,]

## calculate M-values for statistical analysis and B-values for interpretation
rm(keep)
rm(xReactiveProbes)
bVals <- getBeta(mSetSqFlt)
# 394937    659
bvalues <- paste(sample, "/", sample, "_bvalues.csv", sep="")
write.csv(file=bvalues, x=bVals)
rm(bVals)
mVals <- getM(mSetSqFlt)
rda_mvals <- paste(sample, "/", sample, "_mVals.Rda", sep="")
save(mVals, file=rda_mvals)
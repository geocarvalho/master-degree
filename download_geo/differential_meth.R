library(minfi)
library(GEOquery)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(data.table)

sample <- "GSE51032"

# Select idat
idat_file <- paste("GSE51032/idat/normal_colon", sep="")

# Open idat files as array
rgSet <- read.metharray.exp(idat_file)
gc()

# Quality control
## calculate the detection p-values
detP <- detectionP(rgSet)

# remove poor quality samples
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]

# Predict sex
MSet <- preprocessRaw(rgSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)$predictedSex
fwrite(list(predictedSex), "sex_vector.csv")

# Normalization
## normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessFunnorm(rgSet)
rm(rgSet)
gc()

# Filtering
## ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

## remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
mSetSqFlt <- mSetSq[keep,]
rm(mSetSq)
gc()

#___________________ analyse if is necessary remove sex chr depending of the biological question
## if your data includes males and females, remove probes on the sex chromosomes
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
mSetSqFlt <- mSetSqFlt[keep,]

## remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)

## exclude cross reactive probes 
xReactiveProbes <- read.csv(file=paste(dataDirectory, "48639-non-specific-probes-Illumina450k_part1.csv", sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
mSetSqFlt <- mSetSqFlt[keep,] 

## calculate M-values for statistical analysis and B-values for interpretation
# 385750    659
mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)
bvalues <- paste(sample, "/", sample, "_bvalues.csv", sep="")
write.csv(file=bvalues, x=bVals)


# Probe-wise differential methylation analysis

# Differential methylation analysis

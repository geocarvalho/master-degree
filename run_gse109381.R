library(minfi)
library(GEOquery)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)

sample <- "GSE109381"

idat_file <- paste(sample, "/", "idat/major", sep="")

rgSet <- read.metharray.exp(idat_file)
gc()

# Predict sex
# MSet <- preprocessRaw(rgSet)
# RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
# GRset <- mapToGenome(RSet)
# predictedSex <- getSex(GRset, cutoff = -2)$predictedSex

# grSet <- preprocessQuantile(rgSet)

# Differential methylation analysis
# calculate the detection p-values
detP <- detectionP(rgSet)
# remove poor quality samples
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
# remove poor quality samples from detection p-value table
detP <- detP[,keep]
# normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessQuantile(rgSet) # 485512
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
# table(keep) 
mSetSqFlt <- mSetSqFlt[keep,]
# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt) # 402974

# exclude cross reactive probes 
csv_file <- "/home/watson/george/master-degree/48639-non-specific-probes-Illumina450k_part1.csv"
xReactiveProbes <- read.csv(csv_file, stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
mSetSqFlt <- mSetSqFlt[keep,] # 377799
bVals <- getBeta(mSetSqFlt)

# Download beta-values
bvalues <- paste(sample, "/", sample, "_bvalues_gbm.csv", sep="")
write.csv(file=bvalues, x=bVals)

# Download phenotype data
pheno <- paste(sample, "/", sample, "_all_phenotype.csv", sep="")
geoMat <- getGEO(sample)
pD.all <- pData(geoMat[[1]])
write.csv(file=pheno, x=pD.all)
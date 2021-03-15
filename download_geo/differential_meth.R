library(minfi)
library(limma)
library(GEOquery)
library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(data.table)

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

#___________________ analyse if is necessary remove sex chr depending of the biological question
## if your data includes males and females, remove probes on the sex chromosomes
###  I DON'T THINK IT IS A GREAT IDEA, BECAUSE WE WANT TO SEPARATE SUBGROUPS BY SEX IN THE FUTURE
# ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
# mSetSqFlt <- mSetSqFlt[keep,]

## remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)

## exclude cross reactive probes 
xReactiveProbes <- read.csv(file=paste("48639-non-specific-probes-Illumina450k_part1.csv", sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
mSetSqFlt <- mSetSqFlt[keep,] 

# Predict sex
# MSet <- preprocessRaw(mSetSqFlt)
# RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
# GRset <- mapToGenome(RSet)
# predictedSex <- getSex(GRset, cutoff = -2)$predictedSex
# fwrite(list(predictedSex), "sex_vector.csv")

## calculate M-values for statistical analysis and B-values for interpretation
# 385750    659
mVals <- getM(mSetSqFlt)
rm(keep)
rm(xReactiveProbes)
rm(mSetSqFlt)
save(mVals, file="mVals.Rda")
# bVals <- getBeta(mSetSqFlt)
# bvalues <- paste(sample, "/", sample, "_bvalues.csv", sep="")
# write.csv(file=bvalues, x=bVals)

# Show objects consuming memory
# sort( sapply(ls(),function(x){object.size(get(x))}))

# Probe-wise differential methylation analysis

# Differential methylation analysis
# geoMat <- getGEO(filename="./GSE51032/GSE51032_series_matrix.txt.gz", GSEMatrix=TRUE)
# pD.all <- pData(geoMat)

# Rename cancer type column
# names(pD.all)[names(pD.all) == "cancer type (icd-10):ch1"] <- "cancer_type"

# Replace NA with normal in cancer type column
# pD.all["cancer_type"][is.na(pD.all["cancer_type"])] <- "normal"

# Select just breast cancer and normal samples in 'cancer type (icd-10):ch1' column
# filtered <- pD.all[pD.all$cancer_type %in% c("normal", "C50"), ]

# Create the desing matrix
# cancerType <- factor(filtered$cancer_type)
# design <- model.matrix(~cancerType -1)
# colnames(design) <- c("C50", "normal")
# rm(cancerType)
# rm(filtered)
# rm(pD.all)
# rm(geoMat)
# gc()
# save(foo,file="design.Rda")

# Load design file
load("design.Rda")
load("mVals.Rda")
# fit the linear model 
fit <- lmFit(mVals, design)
fit <- eBayes(fit)

# The numbers of hyper-methylated (1) and hypo-methylated (-1)
summary(decideTests(fit))
save(fit, file="fit.Rda")
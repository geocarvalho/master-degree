suppressPackageStartupMessages({
    library(minfi)
    library(limma)
    library(missMethyl)
    library(GEOquery)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    library(IlluminaHumanMethylation450kmanifest)
    library(data.table)
})

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

## remove probes on the sex cromosomes
# ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
# table(keep)
# mSetSqFlt <- mSetSqFlt[keep,]
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
# 385750    659
bvalues <- paste(sample, "/", sample, "_bvalues.csv", sep="")
write.csv(file=bvalues, x=bVals)
rm(bVals)
mVals <- getM(mSetSqFlt)
# rda_mvals <- paste(sample, "/", sample, "_mVals.Rda", sep="")
# save(mVals, file=rda_mvals)

## create_design.R
filtered <- read.csv("../download_geo/GSE51032/GSE51032_classes_design.csv")

# Create the desing matrix
diagnosisTime <- factor(filtered$time_to_diagnosis)
cancerType <- factor(filtered$cancer_type)
design <- model.matrix(~0 + cancerType + diagnosisTime)

colnames(design) <- c("cancert_type_C50", "cancert_type_normal", "diagnosisTime_1", "diagnosisTime_6", "diagnosisTime_7", "diagnosisTime_8", "diagnosisTime_2", "diagnosisTime_3", "diagnosisTime_4", "diagnosisTime_5")

rm(cancerType)
rm(filtered)
rm(pD.all)
rm(geoMat)
gc()
# design_file = paste(sample, "/", sample, "_design.Rda", sep="")
# save(design,file=design_file)

## Create model
sample <- "GSE51032"
# design <- paste(sample, "/", sample, "_design.Rda", sep="")
# mVals <- paste(sample, "/", sample, "_mVals.Rda", sep="")
# Load design file
# load(design)
# load(mVals)
# fit the linear model 
fit <- lmFit(mVals, design)
fit <- eBayes(fit)

# The numbers of hyper-methylated (1) and hypo-methylated (-1)
summary(decideTests(fit))
# fit_output <- paste(sample, "/", sample, "_fit.Rda", sep="")
# save(fit, file=fit_output)

## Differential meth
# create a contrast matrix for specific comparisons
# contMatrix <- makeContrasts("C50-normal", levels=design)
contMatrix <- makeContrasts("cancer_type_C50", levels=design)

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2)) # 175827
diff <- decideTests(fit2)
diff_out <- paste(sample, "/", sample, "_diff_probes.csv", sep="")
write.table(diff, file=diff_out, sep=",", row.names=TRUE)

# get the table of results for the first contrast (naive - rTreg)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
DMPs_out <- paste(sample, "/", sample, "_DMPs.csv", sep="")
write.table(DMPs, file=DMPs_out, sep=",", row.names=FALSE)

# Get the significant CpG sites at less than 5% FDR
sigCpGs <- DMPs$Name[DMPs$adj.P.Val<0.05]
# length(sigCpGs) 175827
# Write a txt file with the significant CpGs
sig_out <- paste(sample, "/", sample, "_sigCpGs.csv", sep="")
fileConn <- file(sig_out)    
writeLines(sigCpGs, fileConn)    
close(fileConn)

# Get all the CpG sites used in the analysis to form the background
all <- DMPs$Name
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=TRUE)

# Top 10 GO categories
topGSA(gst, number=10)
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

# Age example datasets
samples <- c("GSE51032", "GSE42861", "GSE87571")

for (sample in samples) {
    # Download idat
    getGEOSuppFiles(sample)

    raw <- paste(sample, "/", sample, "_RAW.tar", sep="")
    idat_file <- paste(sample, "/", "idat", sep="")

    untar(raw, exdir = idat_file)
    head(list.files(idat_file, pattern = "idat"))

    idatFiles <- list.files(idat_file, pattern = "idat.gz$", full = TRUE)
    sapply(idatFiles, gunzip, overwrite = TRUE)

    rgSet <- read.metharray.exp(idat_file)
    gc()
    grSet <- preprocessQuantile(rgSet)

    # Download beta-values
    bvalues <- paste(sample, "/", sample, "_bvalues.csv", sep="")
    beta <- getBeta(grSet)
    write.csv(file=bvalues, x=beta)

    # Download phenotype data
    pheno <- paste(sample, "/", sample, "_all_phenotype.csv", sep="")
    geoMat <- getGEO(sample)
    pD.all <- pData(geoMat[[1]])
    write.csv(file=pheno, x=pD.all)
    
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


    # calculate M-values for statistical analysis
    mVals <- getM(mSetSqFlt)
    bVals <- getBeta(mSetSqFlt)

    # Probe-wise differencial methylation analysis
    # change columns name
    # names(pD.all)[names(pD.all) == "her2_status:ch1"] <- "her2_status"
    names(pD.all)[names(pD.all) == "subtype_ihc:ch1"] <- "subtype_ihc"
    # names(pD.all)[names(pD.all) == "age_bin:ch1"] <- "age_bin"
    # this is the factor of interest
    cellType <- factor(pD.all$subtype_ihc)
    # her2 = factor(pD.all$her2_status) # try age
    # age = factor(pD.all$age_bin)

    # use the above to create a design matrix
    design <- model.matrix(~cellType -1)
    colnames(design) <- c("Basal","HER2","LumA","LumB")

    # fit the linear model 
    fit <- lmFit(mVals, design)
    fit <- eBayes(fit)

    # The numbers of hyper-methylated (1) and hypo-methylated (-1)
    summary(decideTests(fit))

    # create a contrast matrix for specific comparisons
    # contMatrix <- makeContrasts("Basal-HER2", "Basal-LumA", "Basal-LumB", "HER2-LumA", "HER2-LumB", "LumA-LumB", levels=design)

    # fit the contrasts
    # fit2 <- contrasts.fit(fit, contMatrix)
    # fit2 <- eBayes(fit2)
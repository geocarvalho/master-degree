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

# "GSE72245", "GSE72251", "GSE72254"
sample <- "GSE72245"
idat_file <- "/home/genomika/george/master-degree/GSE72245/idat"

# Download idat
# getGEOSuppFiles(sample)

# raw <- paste(sample, "/", sample, "_RAW.tar", sep="")
# idat_file <- paste(sample, "/", "idat", sep="")

# untar(raw, exdir = idat_file)
# head(list.files(idat_file, pattern = "idat"))

# idatFiles <- list.files(idat_file, pattern = "idat.gz$", full = TRUE)
# sapply(idatFiles, gunzip, overwrite = TRUE)

# Reading data 
rgSet <- read.metharray.exp(idat_file)

# Phenotype data
geoMat <- getGEO(sample)
phenoData <- pData(geoMat[[1]])

# Create object with methylated/unmethylated signals for QC check
MSet <- preprocessRaw(rgSet)

# Create RatioSet to store b-values matrix
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
# beta <- getBeta(RSet)

# Add genomic coordinates to each probe with additional information
GRset <- mapToGenome(RSet)

# Extract and plot the QC information
# qc <- getQC(MSet)
# plotQC(qc)

# Preprocessing and normalization
# If we aren't comparing cancer/normal samples, there is large-scale differences
GRset.quantile <- preprocessQuantile(rgSet, fixOutliers = TRUE,
                                     removeBadSamples = TRUE, 
                                     badSampleCutoff = 10.5,
                                     quantileNormalize = TRUE, 
                                     stratified = TRUE, 
                                     mergeManifest = FALSE, 
                                     sex = NULL)
# But if we are
GRset.funnorm <- preprocessFunnorm(rgSet)


# To drop SNPs present in probes
GRset <- dropLociWithSnps(GRset.funnorm, snps=c("SBE","CpG"), maf=0)
beta <- getBeta(GRset)
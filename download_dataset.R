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

# Download idat
getGEOSuppFiles("GSE72245")
untar("GSE72245/GSE72245_RAW.tar", exdir = "GSE72245/idat")
head(list.files("GSE72245/idat", pattern = "idat"))

idatFiles <- list.files("GSE72245/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

rgSet <- read.metharray.exp("GSE72245/idat")
gc()
grSet <- preprocessQuantile(rgSet)

# Download beta-values
beta <- getBeta(grSet)
write.csv(file="GSE72245/GSE72245_bvalues.csv", x=beta)

# Download phenotype data
geoMat <- getGEO("GSE72245")
pD.all <- pData(geoMat[[1]])
write.csv(file="GSE72245/GSE72245_all_phenotype.csv", x=pD.all)
# https://kasperdanielhansen.github.io/genbioconductor/html/minfi.html
# examples from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72308

library(minfi)
library(GEOquery)

# Download idat
# getGEOSuppFiles("/home/george/Documents/Mestrado/projeto/GEO/GSE72251")
# untar("/home/george/Documents/Mestrado/projeto/GEO/GSE72251/GSE72251_RAW.tar", exdir = "/home/george/Documents/Mestrado/projeto/GEO/GSE72251/idat")
# head(list.files("/home/george/Documents/Mestrado/projeto/GEO/GSE72251/idat", pattern = "idat"))

# idatFiles <- list.files("/home/george/Documents/Mestrado/projeto/GEO/GSE72251/idat", pattern = "idat.gz$", full = TRUE)
# sapply(idatFiles, gunzip, overwrite = TRUE)

# Download phenotype data
# geoMat <- getGEO("GSE72251")
# pD.all <- pData(geoMat[[1]])
# write.csv(file="/home/george/Documents/Mestrado/projeto/GEO/GSE72251/GSE72251_all_phenotype.csv", x=pD.all, row.names=FALSE)
# pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1", "characteristics_ch1.2")]
# head(pD)

# Open dataset
rgSet <- read.metharray.exp("/home/george/Documents/Mestrado/projeto/GEO/GSE72251/idat")

# Preprocess
gc()
grSet <- preprocessQuantile(rgSet)
getBeta(grSet)


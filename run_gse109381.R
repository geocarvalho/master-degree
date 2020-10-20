library(minfi)
library(GEOquery)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)

sample <- "GSE109381"

idat_file <- paste(sample, "/", "idat/part4", sep="")

rgSet <- read.metharray.exp(idat_file)
gc()
grSet <- preprocessQuantile(rgSet)

# Download beta-values
bvalues <- paste(sample, "/", sample, "_bvalues_part4.csv", sep="")
beta <- getBeta(grSet)
write.csv(file=bvalues, x=beta)

# Download phenotype data
pheno <- paste(sample, "/", sample, "_all_phenotype.csv", sep="")
geoMat <- getGEO(sample)
pD.all <- pData(geoMat[[1]])
write.csv(file=pheno, x=pD.all)
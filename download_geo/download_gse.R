library(minfi)
library(GEOquery)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(data.table)

sample <- "GSE51032"
dir.create(sample)
setwd(sample)

# Download idat via wget 
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51032/suppl//GSE51032_RAW.tar
# or 
# getGEOSuppFiles(sample)

raw <- paste(sample, "_RAW.tar", sep="")
idat_file <- paste("idat", sep="")
untar(raw, exdir = idat_file)

head(list.files(idat_file, pattern = "idat"))
idatFiles <- list.files(idat_file, pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

rgSet <- read.metharray.exp(idat_file)
gc()

# Predict sex
MSet <- preprocessRaw(rgSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
predictedSex <- getSex(GRset, cutoff = -2)$predictedSex
fwrite(list(predictedSex), "sex_vector.csv")

# Download phenotype data
pheno <- paste(sample, "_all_phenotype.csv", sep="")
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51032/matrix/GSE51032_series_matrix.txt.gz
# geoMat <- getGEO(sample)
geoMat <- getGEO(filename="./GSE51032/GSE51032_series_matrix.txt.gz", GSEMatrix=TRUE)
pD.all <- pData(geoMat)
write.csv(file=pheno, x=pD.all)

# Differential methylation analysis

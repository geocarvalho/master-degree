library(minfi)
library(GEOquery)

# Download phenotype data
pheno <- paste(sample, "_all_phenotype.csv", sep="")
# Download Matrix
# geoMat <- getGEO(sample) # timeout, use download_GSE51032.sh
geoMat <- getGEO(filename="./GSE51032/GSE51032_series_matrix.txt.gz", GSEMatrix=TRUE)
pD.all <- pData(geoMat)
write.csv(file=pheno, x=pD.all)
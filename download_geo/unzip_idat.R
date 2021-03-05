library(minfi)
library(GEOquery)

sample <- "GSE51032"
# dir.create(sample)
# setwd(sample)

# Download idat via wget or:
# getGEOSuppFiles(sample) # timeout, use download_GSE51032.sh

# Unzip RAW file and idat files
raw <- paste(sample, "_RAW.tar", sep="")
idat_file <- paste("idat", sep="")
untar(raw, exdir = idat_file)

head(list.files(idat_file, pattern = "idat"))
idatFiles <- list.files(idat_file, pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)
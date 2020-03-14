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

samples <- c("GSE72245", "GSE72251", "GSE72254")
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
}

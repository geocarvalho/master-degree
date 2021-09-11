if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ENmix")

suppressPackageStartupMessages({
    library(minfi)
    library(ENmix)
})

# Sample
sample <- "GSE51032"

# Select idat
idat_file <- "GSE51032/idat/normal_mama"

# Open idat files as array
rgSet <- read.metharray.exp(idat_file)

# Calculate age
meth=getmeth(rgSet)
beta=getB(meth)
mage=methyAge(beta)

write.csv(mage, file = "./GSE51032/age_prediction_GSE109381.csv", row.names = FALSE)
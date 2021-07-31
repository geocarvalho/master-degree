library(GEOquery)

# Probe-wise differential methylation analysis
sample <- "GSE51032"
# Differential methylation analysis
geoMat <- getGEO(filename="./GSE51032/GSE51032_series_matrix.txt", GSEMatrix=TRUE)
pD.all <- pData(geoMat)

# Rename cancer type column
names(pD.all)[names(pD.all) == "cancer type (icd-10):ch1"] <- "cancer_type"
names(pD.all)[names(pD.all) == "gender:ch1"] <- "gender"

# Replace NA with normal in cancer type column
pD.all["cancer_type"][is.na(pD.all["cancer_type"])] <- "normal"

# Select just breast cancer and normal samples in 'cancer type (icd-10):ch1' column
filtered <- pD.all[pD.all$cancer_type %in% c("normal", "C50"), ]
filtered <- filtered[filtered$gender %in% c("F"), ]

# Create the desing matrix
cancerType <- factor(filtered$cancer_type)
design <- model.matrix(~cancerType -1)
colnames(design) <- c("C50", "normal")
rm(cancerType)
rm(filtered)
rm(pD.all)
rm(geoMat)
gc()
design_file = paste(sample, "/", sample, "_design.Rda", sep="")
save(design,file=design_file)
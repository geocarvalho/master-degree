library(GEOquery)

# Probe-wise differential methylation analysis
# sample <- "GSE51032"
# Differential methylation analysis
# geoMat <- getGEO(filename="../data/GSE51032_classes_design.csv", GSEMatrix=TRUE)
# pD.all <- pData(geoMat)

# Rename cancer type column
# names(pD.all)[names(pD.all) == "cancer type (icd-10):ch1"] <- "cancer_type"
# names(pD.all)[names(pD.all) == "gender:ch1"] <- "gender"
# names(pD.all)[names(pD.all) == "time to diagnosis:ch1"] <- "time_to_diagnosis"

# Replace NA with normal in cancer type column
# pD.all["cancer_type"][is.na(pD.all["cancer_type"])] <- "normal"

# Select just breast cancer and normal samples in 'cancer type (icd-10):ch1' column
# filtered <- pD.all[pD.all$cancer_type %in% c("normal", "C50"), ]
# filtered <- filtered[filtered$gender %in% c("F"), ]

# Read csv classes design
filtered <- read.csv("../download_geo/GSE51032/GSE51032_classes_design.csv")

# Create the desing matrix
diagnosisTime <- factor(filtered$time_to_diagnosis_classes)
cancerType <- factor(filtered$cancer_type)
design <- model.matrix(~0 + cancerType + diagnosisTime, data=filtered)
# design <- model.matrix(~0+ cancerType)
# colnames(design) <- c("C50", "normal")
# colnames(design) <- c("cancer_type_C50", "cancer_type_normal", "diagnosisTime>=0.04,<2.04", 
# "diagnosisTime>=10.04,<12.04", "diagnosisTime>=12.04,<14.04" "diagnosisTime>=14.04,<16.04", 
# "diagnosisTime>=2.04,<4.04", "diagnosisTime>=4.04,<6.04", "diagnosisTime>=6.04,<8.04", 
# "diagnosisTime>=8.04,<10.04")

colnames(design) <- c("cancer_type_C50", "cancer_type_normal", "diagnosisTime_1", "diagnosisTime_6", "diagnosisTime_7", "diagnosisTime_8", "diagnosisTime_2", "diagnosisTime_3", "diagnosisTime_4", "diagnosisTime_5")

rm(cancerType)
rm(filtered)
rm(pD.all)
rm(geoMat)
gc()
design_file = paste(sample, "/", sample, "_design.Rda", sep="")
save(design,file=design_file)
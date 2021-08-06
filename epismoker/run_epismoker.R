# http://htmlpreview.github.io/?https://github.com/sailalithabollepalli/EpiSmokEr/blob/master/vignettes/epismoker.html

suppressPackageStartupMessages({
library(EpiSmokEr)  
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(htmlTable)
library(rmarkdown)
})

# Careful with the path because it should be inside the docker image
idatpath <- "./GSE51032/idat/normal_mama"
rawdata <- loadData(idatpath)
samplesheet <- read.csv("./GSE51032/idat/normal_mama/samplesheet_GSE109381_simple.csv", header=TRUE, sep=",")
samplesheet2 <- samplesheet[,-1]
rownames(samplesheet2) <- samplesheet[,1]
# knitr::kable(head(samplesheet, 5)
dataset_QN <- normalizeData(RGset=rawdata, normMethod = "QN")

# dataset_SQN <- normalizeData(RGset=rawdata, normMethod = "SQN")
# Error in normalize.quantiles(mat[Index2, ]) : 
#   ERROR; return code from pthread_create() is 22

# dataset_ILM <- normalizeData(RGset=rawdata, normMethod = "ILM")

# dataset_ALL <- normalizeData(RGset=rawdata, normMethod = "ALL")
# str(dataset_ALL)
# dataset_QN <- dataset_ALL$dataset_QN
# dataset_ILM <- dataset_ALL$dataset_ILM
# dataset_SQN <- dataset_ALL$dataset_SQN

# Smoking Status (SSt) Prediction
result_SSt <- epismoker(dataset=dataset_QN, samplesheet = samplesheet2, method = "SSt")

# Smoking Score (SSc)
# result_SSc <- epismoker(dataset = dataset_SQN, method = "SSc")

# Methylation Score (MS)
# result_MS <- epismoker(dataset = dataset_ILM,  method = "MS")

# Compreshensive approach (all)
# result_All <- epismoker(dataset_QN=dataset_QN, dataset_ILM=dataset_ILM, dataset_SQN=dataset_SQN, samplesheet = samplesheet, method = "all")

write.csv(result_SSt, file = "./GSE51032/epismoker_SSt_GSE109381.csv", row.names = FALSE)
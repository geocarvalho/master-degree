library(limma)

sample <- "GSE51032"
design <- paste(sample, "/", sample, "_design.Rda", sep="")
mVals <- paste(sample, "/", sample, "_mVals.Rda", sep="")
# Load design file
load(design)
load(mVals)
# fit the linear model
fit <- lmFit(mVals, design)

contMatrix <- makeContrasts(cancerTypeC50-cancerTypenormal, levels=design) # 70107 + 82787
# contMatrix <- makeContrasts(cancer_type_C50-cancer_type_normal, # 6106 + 7982
#                             cancer_type_C50-diagnosisTime_1,
#                             cancer_type_C50-diagnosisTime_2,
#                             cancer_type_C50-diagnosisTime_3,
#                             cancer_type_C50-diagnosisTime_4,
#                             cancer_type_C50-diagnosisTime_6,
#                             cancer_type_C50-diagnosisTime_7,
#                             cancer_type_C50-diagnosisTime_8,
#                            levels=design)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
# fit <- eBayes(fit)

# The numbers of hyper-methylated (1) and hypo-methylated (-1)
summary(decideTests(fit2))
fit_output <- paste(sample, "/", sample, "_fit.Rda", sep="")
save(fit2, file=fit_output)

# Create output for the summary for the big matrix

diff <- summary(decideTests(fit2))
diff_out <- paste(sample, "/", sample, "_summary_big_matrix.csv", sep="")
write.table(diff, file=diff_out, sep="\t", row.names=TRUE)

# DMRs
library(doParallel)
registerDoParallel(cores = 30)

dmrs <- bumphunter(mSetSqFlt, design = design, cutoff = 0.95, B=0, type="Beta")

# if > 30k increase cutoff

dmrs <- bumphunter(mSetSqFlt, design = design, cutoff = 0.95, B=500, type="Beta")

# GO and KEGG terms

library(limma)

sample <- "GSE51032"
design <- paste(sample, "/", sample, "_design.Rda", sep="")
mVals <- paste(sample, "/", sample, "_mVals.Rda", sep="")
# Load design file
load(design)
load(mVals)
# fit the linear model 
fit <- lmFit(mVals, design)

contMatrix <- makeContrasts(cancerTypeC50-cancerTypenormal, levels=design)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
# fit <- eBayes(fit)

# The numbers of hyper-methylated (1) and hypo-methylated (-1)
summary(decideTests(fit2))
fit_output <- paste(sample, "/", sample, "_fit.Rda", sep="")
save(fit2, file=fit_output)

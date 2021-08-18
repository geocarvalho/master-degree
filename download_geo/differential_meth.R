suppressPackageStartupMessages({
    library(limma)
    library(missMethyl)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    library(data.table)
})

sample <- "GSE51032"

fit <- paste(sample, "/", sample, "_fit.Rda", sep="")
design <- paste(sample, "/", sample, "_design.Rda", sep="")
mVals <- paste(sample, "/", sample, "_mVals.Rda", sep="")
load(fit)
load(design)
load(mVals)


# create a contrast matrix for specific comparisons
# contMatrix <- makeContrasts("C50-normal", levels=design)
contMatrix <- makeContrasts("cancer_type_C50", levels=design)

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2)) # 152894 down+up
diff <- decideTests(fit2)
diff_out <- paste(sample, "/", sample, "_diff_probes_cn.csv", sep="")
write.table(diff, file=diff_out, sep=",", row.names=TRUE)

# get the table of results for the first contrast (naive - rTreg)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
DMPs_out <- paste(sample, "/", sample, "_DMPs_cn.csv", sep="")
write.table(DMPs, file=DMPs_out, sep=",", row.names=FALSE)

# Get the significant CpG sites at less than 5% FDR
sigCpGs <- DMPs$Name[DMPs$adj.P.Val<0.05]
# length(sigCpGs) 152894
# Write a txt file with the significant CpGs
sig_out <- paste(sample, "/", sample, "_sigCpGs_cn.csv", sep="")
fileConn <- file(sig_out)    
writeLines(sigCpGs, fileConn)    
close(fileConn)

# Get all the CpG sites used in the analysis to form the background
all <- DMPs$Name
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=TRUE)

# Top 10 GO categories
topGSA(gst, number=10)
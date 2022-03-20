suppressPackageStartupMessages({
    library(minfi)
    library(limma)
    library(missMethyl)
    library(GEOquery)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    library(IlluminaHumanMethylation450kmanifest)
    library(data.table)
    library(hash)
})

# Sample
sample <- "GSE51032"

# Select idat
idat_file <- "GSE51032/idat/normal_mama"

# Open idat files as array
rgSet <- read.metharray.exp(idat_file)

# Quality control
## calculate the detection p-values
detP <- detectionP(rgSet)

# remove poor quality samples
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]

# Normalization
## normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessFunnorm(rgSet)
rm(idat_file)
rm(rgSet)
gc()

# Filtering
## ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

## remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
mSetSqFlt <- mSetSq[keep,]

## remove probes on the sex cromosomes
# ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
# table(keep)
# mSetSqFlt <- mSetSqFlt[keep,]
rm(mSetSq)
rm(detP)
gc()

## remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)

## exclude cross reactive probes 
xReactiveProbes <- read.csv(file=paste("48639-non-specific-probes-Illumina450k_part1.csv", sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
mSetSqFlt <- mSetSqFlt[keep,]

## calculate M-values for statistical analysis and B-values for interpretation
rm(keep)
rm(xReactiveProbes)
bVals <- getBeta(mSetSqFlt)
# 385750    659
bvalues <- paste(sample, "/", sample, "_bvalues.csv", sep="")
write.csv(file=bvalues, x=bVals)
# rm(bVals)
mVals <- getM(mSetSqFlt)
# rda_mvals <- paste(sample, "/", sample, "_mVals.Rda", sep="")
# save(mVals, file=rda_mvals)

## create_design.R
filtered <- read.csv("../download_geo/GSE51032/GSE51032_classes_design.csv")

## https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
### dmpFinder: to find differentially methylated positions (DMPs)
diagnosisTime <- factor(filtered$time_to_diagnosis_classes)
cancerType <- factor(filtered$cancer_type)
# nk <- factor(filtered$NK)
# dmp <- dmpFinder(bVals, pheno=cancerType, type="categorical")
# dmp <- dmpFinder(bVals, pheno=nk, type="continuous")
# keep <- dmp$pval < 0.01
# Select most significant DMPs
# dmp[keep,]

# Create the desing matrix
design <- model.matrix(~0 + cancerType)
# design <- model.matrix(~0 + cancerType + diagnosisTime, data=filtered)

# colnames(design) <- c("cancer_type_C50", "cancer_type_normal", "diagnosisTime_1", "diagnosisTime_6", "diagnosisTime_7", "diagnosisTime_8", "diagnosisTime_2", "diagnosisTime_3", "diagnosisTime_4", "diagnosisTime_5")

rm(cancerType)
rm(filtered)
rm(pD.all)
rm(geoMat)
gc()
# design_file = paste(sample, "/", sample, "_design.Rda", sep="")
# save(design,file=design_file)

## Create model
# sample <- "GSE51032"
# design <- paste(sample, "/", sample, "_design.Rda", sep="")
# mVals <- paste(sample, "/", sample, "_mVals.Rda", sep="")
# Load design file
# load(design)
# load(mVals)
# fit the linear model 
fit <- lmFit(mVals, design)

# The numbers of hyper-methylated (1) and hypo-methylated (-1)
# summary(decideTests(fit))
# fit_output <- paste(sample, "/", sample, "_fit.Rda", sep="")
# save(fit, file=fit_output)

## Differential meth
# create a contrast matrix for specific comparisons
# contMatrix <- makeContrasts("C50-normal", levels=design)
contMatrix <- makeContrasts("cancerTypeC50-cancerTypenormal", levels=design)

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2)) # 70107 + 82787 = 152894
diff <- decideTests(fit2)
diff_out <- paste(sample, "/", sample, "_diff_probes.csv", sep="")
write.table(diff, file=diff_out, sep=",", row.names=TRUE)

# get the table of results for the first contrast (naive - rTreg)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
#head(DMPs)
DMPs_out <- paste(sample, "/", sample, "_DMPs.csv", sep="")
write.table(DMPs, file=DMPs_out, sep=",", row.names=FALSE)

# Get the significant CpG sites at less than 5% FDR
sigCpGs <- DMPs$Name[DMPs$adj.P.Val<0.05] #14088
# length(sigCpGs) 152894
# Write a txt file with the significant CpGs
sig_out <- paste(sample, "/", sample, "_sigCpGs.csv", sep="")
fileConn <- file(sig_out)    
writeLines(sigCpGs, fileConn)    
close(fileConn)

# Get all the CpG sites used in the analysis to form the background
## Pysubgroup
pysg <- hash()
pysg[["subgroup1"]] <- c('cg00015930', 'cg00262415', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01962937', 'cg02260340', 'cg02527528', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03397332', 'cg03801179', 'cg04104256', 'cg04285418', 'cg04551440', 'cg05006947', 'cg05046026', 'cg05664296', 'cg06097659', 'cg06668922', 'cg06770554', 'cg06775930', 'cg07062711', 'cg07610406', 'cg09449449', 'cg10318121', 'cg10437571', 'cg12071806', 'cg12510286', 'cg14651446', 'cg15717853', 'cg15859496', 'cg20972117', 'cg21279459', 'cg21881074', 'cg24154336', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')
pysg[["subgroup2"]] <- c('cg00015930', 'cg00262415', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01962937', 'cg02260340', 'cg02527528', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03397332', 'cg03801179', 'cg04104256', 'cg04285418', 'cg04551440', 'cg05006947', 'cg05046026', 'cg05664296', 'cg06097659', 'cg06668922', 'cg06770554', 'cg06775930', 'cg07062711', 'cg07610406', 'cg09449449', 'cg10318121', 'cg10437571', 'cg12071806', 'cg13613388', 'cg14651446', 'cg15717853', 'cg15859496', 'cg17386812', 'cg20972117', 'cg21279459', 'cg24154336', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')
pysg[["subgroup3"]] <- c('cg00015930', 'cg00262415', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01962937', 'cg02260340', 'cg02527528', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03397332', 'cg03801179', 'cg04104256', 'cg04285418', 'cg04551440', 'cg05006947', 'cg05046026', 'cg05664296', 'cg06097659', 'cg06668922', 'cg06770554', 'cg06775930', 'cg07062711', 'cg07610406', 'cg09449449', 'cg10318121', 'cg10437571', 'cg12071806', 'cg14651446', 'cg15717853', 'cg15859496', 'cg17386812', 'cg20972117', 'cg21279459', 'cg21327609', 'cg24154336', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')
pysg[["subgroup4"]] <- c('cg00015930', 'cg00262415', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01962937', 'cg02260340', 'cg02527528', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03397332', 'cg03801179', 'cg04104256', 'cg04285418', 'cg04551440', 'cg05006947', 'cg05046026', 'cg05664296', 'cg06097659', 'cg06668922', 'cg06770554', 'cg06775930', 'cg07062711', 'cg07610406', 'cg09449449', 'cg10318121', 'cg10437571', 'cg12071806', 'cg12510286', 'cg14651446', 'cg15717853', 'cg15859496', 'cg17386812', 'cg20972117', 'cg21279459', 'cg24154336', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')
pysg[["subgroup5"]] <- c('cg00015930', 'cg00262415', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01962937', 'cg02260340', 'cg02527528', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03397332', 'cg03801179', 'cg04104256', 'cg04285418', 'cg04551440', 'cg05006947', 'cg05046026', 'cg05664296', 'cg06097659', 'cg06668922', 'cg06770554', 'cg06775930', 'cg07062711', 'cg07610406', 'cg09449449', 'cg10318121', 'cg10437571', 'cg12071806', 'cg13613388', 'cg14651446', 'cg15717853', 'cg15859496', 'cg20972117', 'cg21279459', 'cg21881074', 'cg24154336', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')
pysg[["subgroup6"]] <- c('cg00015930', 'cg00262415', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01962937', 'cg02260340', 'cg02527528', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03397332', 'cg03801179', 'cg04104256', 'cg04285418', 'cg04551440', 'cg05006947', 'cg05046026', 'cg05664296', 'cg06097659', 'cg06668922', 'cg06770554', 'cg06775930', 'cg07062711', 'cg07610406', 'cg09449449', 'cg10318121', 'cg10437571', 'cg12071806', 'cg14651446', 'cg15717853', 'cg15859496', 'cg20972117', 'cg21279459', 'cg21327609', 'cg21881074', 'cg24154336', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')
pysg[["subgroup7"]] <- c('cg00015930', 'cg00262415', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01962937', 'cg02260340', 'cg02527528', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03397332', 'cg03801179', 'cg04104256', 'cg04285418', 'cg04551440', 'cg05006947', 'cg05046026', 'cg05664296', 'cg06097659', 'cg06668922', 'cg06770554', 'cg06775930', 'cg07062711', 'cg07610406', 'cg09449449', 'cg10318121', 'cg10858568', 'cg12071806', 'cg12510286', 'cg14651446', 'cg15717853', 'cg15859496', 'cg20972117', 'cg21279459', 'cg21881074', 'cg24154336', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')
pysg[["subgroup8"]] <- c('cg00015930', 'cg00262415', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01962937', 'cg02260340', 'cg02527528', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03397332', 'cg03801179', 'cg04104256', 'cg04285418', 'cg04551440', 'cg05006947', 'cg05046026', 'cg05664296', 'cg06097659', 'cg06668922', 'cg06770554', 'cg06775930', 'cg07062711', 'cg07610406', 'cg09449449', 'cg10318121', 'cg10858568', 'cg12071806', 'cg12510286', 'cg14651446', 'cg15717853', 'cg15859496', 'cg17386812', 'cg20972117', 'cg21279459', 'cg24154336', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')
pysg[["subgroup9"]] <- c('cg00015930', 'cg00262415', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01962937', 'cg02260340', 'cg02527528', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03397332', 'cg03801179', 'cg04104256', 'cg04285418', 'cg04551440', 'cg05006947', 'cg05046026', 'cg05664296', 'cg06097659', 'cg06668922', 'cg06770554', 'cg06775930', 'cg07062711', 'cg07610406', 'cg09449449', 'cg10318121', 'cg10858568', 'cg12071806', 'cg14651446', 'cg15717853', 'cg15859496', 'cg17386812', 'cg20972117', 'cg21279459', 'cg21327609', 'cg24154336', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')
pysg[["subgroup10"]] <- c('cg00015930', 'cg00262415', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01962937', 'cg02260340', 'cg02527528', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03397332', 'cg03801179', 'cg04104256', 'cg04285418', 'cg04551440', 'cg05006947', 'cg05046026', 'cg05664296', 'cg06097659', 'cg06668922', 'cg06770554', 'cg06775930', 'cg07062711', 'cg07610406', 'cg09449449', 'cg10318121', 'cg10858568', 'cg12071806', 'cg13613388', 'cg14651446', 'cg15717853', 'cg15859496', 'cg17386812', 'cg20972117', 'cg21279459', 'cg24154336', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')

# For each subgroup create kte and gte
for (v in ls(pysg)) {
    print(pysg[[v]])
    kte <- topGSA(gometh(sig.cpg=pysg[[v]], all.cpg=sigCpGs, plot.bias=FALSE, collection='KEGG'), number=10)
    gte <- topGSA(gometh(sig.cpg=pysg[[v]], all.cpg=sigCpGs, plot.bias=FALSE), number=10)

    kte_out <- paste(sample, "/filtered_results/sg_anno/", sample, "_", v, "_beam_kte.csv", sep="")
    write.table(kte, file=kte_out, sep=",", row.names=TRUE)
    # subgroups_anno filtered_results
    gte_out <- paste(sample, "/filtered_results/sg_anno/", sample, "_", v, "_beam_gte.csv", sep="")
    write.table(gte, file=gte_out, sep=",", row.names=TRUE)
} 

## SSDP+ (qg)
qg <- hash()
qg[["subgroup1"]] <- c('cg03610970', 'cg22793136')
qg[["subgroup2"]] <- c('cg08879482', 'cg24434368', 'cg09706833')
qg[["subgroup3"]] <- c('cg14154819', 'cg10668569')
# qg[["subgroup4"]] <- c('cg22947748', 'cg03534326')
qg[["subgroup5"]] <- c('cg24367316', 'cg16995299', 'cg06532611')
qg[["subgroup6"]] <- c('cg02674789', 'cg26684511', 'cg24434368')
qg[["subgroup7"]] <- c('cg12785643', 'cg10110581', 'cg03610970')
qg[["subgroup8"]] <- c('cg17820247', 'cg10872521')
qg[["subgroup9"]] <- c('cg24463471', 'cg14154819', 'cg25741023', 'cg09706833')
qg[["subgroup10"]] <- c('cg09444358', 'cg13551227')

# For each subgroup create kte and gte
for (v in ls(qg)) {
    print(qg[[v]])
    kte <- topGSA(gometh(sig.cpg=qg[[v]], all.cpg=sigCpGs, plot.bias=FALSE, collection='KEGG'), number=10)
    gte <- topGSA(gometh(sig.cpg=qg[[v]], all.cpg=sigCpGs, plot.bias=FALSE), number=10)

    kte_out <- paste(sample, "/filtered_results/sg_anno/", sample, "_", v, "_qg_kte.csv", sep="")
    write.table(kte, file=kte_out, sep=",", row.names=TRUE)
    # subgroups_anno filtered_results
    gte_out <- paste(sample, "/filtered_results/sg_anno/", sample, "_", v, "_qg_gte.csv", sep="")
    write.table(gte, file=gte_out, sep=",", row.names=TRUE)
} 

# Top 10 GO categories
topGSA(gst, number=10)

suppressPackageStartupMessages({
    library(minfi)
    library(limma)
    library(missMethyl)
    library(GEOquery)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    library(IlluminaHumanMethylation450kmanifest)
    library(data.table)
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
# design <- model.matrix(~0 + cancerType)
design <- model.matrix(~0 + cancerType + diagnosisTime, data=filtered)

colnames(design) <- c("cancer_type_C50", "cancer_type_normal", "diagnosisTime_1", "diagnosisTime_6", "diagnosisTime_7", "diagnosisTime_8", "diagnosisTime_2", "diagnosisTime_3", "diagnosisTime_4", "diagnosisTime_5")

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
summary(decideTests(fit2)) # 6106 + 7982
diff <- decideTests(fit2)
diff_out <- paste(sample, "/", sample, "_diff_probes.csv", sep="")
write.table(diff, file=diff_out, sep=",", row.names=TRUE)

# get the table of results for the first contrast (naive - rTreg)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
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
# subgroup1 <- c('cg00015930', 'cg00262415', 'cg00475161', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01946760', 'cg01962937', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03801179', 'cg04104256', 'cg04551440', 'cg05006947', 'cg05046026', 'cg05664296', 'cg06721712', 'cg06770554', 'cg06932756', 'cg07062711', 'cg07610406', 'cg07859799', 'cg08905629', 'cg09449449', 'cg10318121', 'cg10437571', 'cg12071806', 'cg12510286', 'cg12589307', 'cg14651446', 'cg15717853', 'cg15859496', 'cg17386812', 'cg20972117', 'cg21279459', 'cg24154336', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')
# subgroup2 <- c('cg00015930', 'cg00262415', 'cg00475161', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01946760', 'cg01962937', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03801179', 'cg04104256', 'cg04551440', 'cg05006947', 'cg05046026', 'cg06721712', 'cg06770554', 'cg06932756', 'cg07062711', 'cg07610406', 'cg07859799', 'cg08905629', 'cg09449449', 'cg10318121', 'cg10437571', 'cg12071806', 'cg12510286', 'cg12589307', 'cg14651446', 'cg15717853', 'cg15859496', 'cg16661579', 'cg17386812', 'cg20972117', 'cg21279459', 'cg24154336', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')
# subgroup3 <- c('cg00015930', 'cg00262415', 'cg00475161', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01946760', 'cg01962937', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03801179', 'cg04104256', 'cg04551440', 'cg05006947', 'cg05046026', 'cg05664296', 'cg06721712', 'cg06770554', 'cg06932756', 'cg07062711', 'cg07610406', 'cg07859799', 'cg08905629', 'cg09449449', 'cg10318121', 'cg10437571', 'cg12071806', 'cg12510286', 'cg12589307', 'cg14651446', 'cg15717853', 'cg15859496', 'cg20972117', 'cg21279459', 'cg21881074', 'cg24154336', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')
# subgroup4 <- c('cg00015930', 'cg00262415', 'cg00475161', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01946760', 'cg01962937', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03801179', 'cg04104256', 'cg04551440', 'cg05006947', 'cg05046026', 'cg06721712', 'cg06770554', 'cg06932756', 'cg07062711', 'cg07610406', 'cg07859799', 'cg08905629', 'cg09449449', 'cg10318121', 'cg10437571', 'cg12071806', 'cg12510286', 'cg12589307', 'cg14651446', 'cg15717853', 'cg15859496', 'cg16661579', 'cg20972117', 'cg21279459', 'cg21881074', 'cg24154336', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')
# subgroup5 <- c('cg00015930', 'cg00262415', 'cg00475161', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01946760', 'cg01962937', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03801179', 'cg04104256', 'cg04551440', 'cg05006947', 'cg05046026', 'cg06721712', 'cg06770554', 'cg06932756', 'cg07062711', 'cg07610406', 'cg07859799', 'cg08905629', 'cg09449449', 'cg10318121', 'cg10437571', 'cg12071806', 'cg12510286', 'cg12589307', 'cg14651446', 'cg15717853', 'cg15859496', 'cg20972117', 'cg21279459', 'cg21881074', 'cg24154336', 'cg24197303', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')
# subgroup6 <- c('cg00015930', 'cg00262415', 'cg00475161', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01946760', 'cg01962937', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03801179', 'cg04104256', 'cg04551440', 'cg05006947', 'cg05046026', 'cg06721712', 'cg06770554', 'cg06932756', 'cg07062711', 'cg07610406', 'cg07859799', 'cg08905629', 'cg09449449', 'cg10318121', 'cg10437571', 'cg12071806', 'cg12510286', 'cg12589307', 'cg14651446', 'cg15717853', 'cg15859496', 'cg20972117', 'cg21279459', 'cg21881074', 'cg24154336', 'cg25045746', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')
# subgroup7 <- c('cg00015930', 'cg00262415', 'cg00475161', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01946760', 'cg01962937', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03801179', 'cg04104256', 'cg04551440', 'cg05006947', 'cg05046026', 'cg06721712', 'cg06770554', 'cg06932756', 'cg07062711', 'cg07610406', 'cg07859799', 'cg08905629', 'cg09449449', 'cg10318121', 'cg10437571', 'cg12071806', 'cg12510286', 'cg12589307', 'cg14651446', 'cg15717853', 'cg15859496', 'cg17386812', 'cg20972117', 'cg21279459', 'cg24154336', 'cg25045746', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')
# subgroup8 <- c('cg00015930', 'cg00262415', 'cg00475161', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01946760', 'cg01962937', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03801179', 'cg04104256', 'cg04551440', 'cg05006947', 'cg05046026', 'cg06721712', 'cg06770554', 'cg06932756', 'cg07062711', 'cg07610406', 'cg07859799', 'cg08905629', 'cg09449449', 'cg10318121', 'cg10437571', 'cg12071806', 'cg12510286', 'cg12589307', 'cg14651446', 'cg15717853', 'cg15859496', 'cg17386812', 'cg20972117', 'cg21279459', 'cg24154336', 'cg24197303', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')
# subgroup9 <- c('cg00015930', 'cg00262415', 'cg00475161', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01946760', 'cg01962937', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03801179', 'cg04551440', 'cg05006947', 'cg05046026', 'cg05664296', 'cg06721712', 'cg06770554', 'cg06932756', 'cg07062711', 'cg07610406', 'cg07859799', 'cg08905629', 'cg09449449', 'cg10318121', 'cg10437571', 'cg12071806', 'cg12510286', 'cg12589307', 'cg14651446', 'cg15717853', 'cg15859496', 'cg16733855', 'cg17386812', 'cg20972117', 'cg21279459', 'cg24154336', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')
# subgroup10 <- c('cg00015930', 'cg00262415', 'cg00475161', 'cg00505073', 'cg01016662', 'cg01031032', 'cg01622006', 'cg01892567', 'cg01946760', 'cg01962937', 'cg02647825', 'cg03161767', 'cg03255749', 'cg03801179', 'cg04551440', 'cg05006947', 'cg05046026', 'cg05664296', 'cg06721712', 'cg06770554', 'cg06932756', 'cg07062711', 'cg07610406', 'cg07859799', 'cg08905629', 'cg09449449', 'cg10318121', 'cg10437571', 'cg12071806', 'cg12510286', 'cg12589307', 'cg14651446', 'cg15717853', 'cg15859496', 'cg16733855', 'cg20972117', 'cg21279459', 'cg21881074', 'cg24154336', 'cg25373297', 'cg26821539', 'cg27544799', 'ch.9.2042397F')

## SSDP+
# subgroup1 <- c('cg15515265', 'cg09470905', 'cg24302688', 'cg05027081', 'cg10318121', 'cg25045746', 'cg04285418', 'cg07506744', 'cg16090392', 'cg16325394', 'cg15131808', 'cg20300129', 'cg23634318', 'cg21961890', 'cg25323005', 'cg20705358', 'cg20079899', 'cg05019854', 'cg27628340', 'ch.9.2042397F', 'cg07577957', 'cg18591136', 'cg19693626', 'cg18262591', 'cg20935363', 'cg10166888', 'cg24828603', 'cg01103590', 'cg12404465', 'cg17158761', 'cg06948630', 'cg20601096', 'cg15479073', 'cg06603997', 'cg02214595', 'cg03051411', 'cg00647830', 'cg18413367', 'cg21033655', 'cg12090003', 'ch.19.1787100F', 'cg13407422', 'cg12947626', 'cg14179908')
# subgroup2 <- c('cg09470905', 'cg21906728', 'cg26695837', 'cg23634318', 'cg13831540', 'cg01653422', 'ch.9.2042397F', 'cg17751435', 'cg27459530', 'cg09920427', 'cg21033655', 'cg26846726', 'cg19069346', 'cg14179908')
# subgroup3 <- c('cg12576688', 'cg10318121', 'cg15637234', 'cg07506744', 'cg16090392', 'cg04654261', 'cg12722569', 'cg13314057', 'cg13831540', 'cg20079899', 'ch.9.2042397F', 'cg17751435', 'cg06603997', 'cg21033655', 'ch.19.1787100F', 'cg15593721', 'cg13407422')
# subgroup4 <- c('cg09470905', 'cg24302688', 'cg26029682', 'cg05300996', 'cg21906728', 'cg17469292', 'cg16090392', 'cg12722569', 'cg13314057', 'cg20705358', 'ch.9.2042397F', 'cg07577957', 'cg12404465', 'cg02214595', 'cg27459530', 'cg21033655', 'cg12090003', 'cg12947626')
# subgroup5 <- c('cg00542041', 'cg26029682', 'cg10318121', 'cg15637234', 'cg04285418', 'cg12837046', 'cg12722569', 'ch.9.2042397F', 'cg02214595', 'cg01633184', 'cg08138366', 'cg09920427', 'cg21033655', 'cg15593721')
# subgroup6 <- c('cg12576688', 'cg09470905', 'cg10318121', 'cg04285418', 'cg04654261', 'cg23634318', 'cg12722569', 'cg13314057', 'cg20705358', 'cg05019854', 'cg17751435', 'cg03503634', 'cg24828603', 'cg04968132', 'cg27459530', 'cg01633184', 'cg09920427', 'cg26846726')
# subgroup7 <- c('cg15515265', 'cg25045746', 'cg05300996', 'cg17469292', 'cg26695837', 'cg20300129', 'cg23634318', 'cg13314057', 'cg18591136', 'cg17751435', 'cg27459530', 'cg08138366', 'cg21033655')
# subgroup8 <- c('cg26029682', 'cg07712892', 'cg21906728', 'cg23634318', 'cg21961890', 'cg13831540', 'cg17509872', 'ch.9.2042397F', 'cg26846726')
# subgroup9 <- c('cg00542041', 'cg24302688', 'cg26029682', 'cg10318121', 'cg15637234', 'cg15131808', 'cg10166888', 'cg01103590', 'cg20601096', 'cg01633184', 'cg08138366', 'cg15593721')
# subgroup <- c('cg00542041', 'cg09470905', 'cg24596729', 'cg15637234', 'cg20822448', 'cg16090392', 'cg12722569', 'cg13831540', 'cg20705358', 'ch.9.2042397F', 'cg17751435', 'cg19693626', 'cg15859496', 'cg12404465', 'cg12090003', 'cg15593721')

# beam search 2 - first one is the common
# subgroup <- c('cg17934367', 'cg26885578', 'cg18497550', 'cg01745873', 'cg06529439', 'cg10456069', 'cg15730116', 'cg07704144', 'cg10587741', 'cg01096478', 'cg07003587', 'cg13554071', 'cg07063325', 'cg25117123', 'cg26045220', 'cg07839336', 'cg20694619', 'cg19456996==3 & cg02964474', 'cg26680520', 'cg22466771', 'cg00126034', 'cg00847892', 'cg02105261', 'cg07782002', 'cg02582925', 'cg22106401', 'cg09909351', 'cg04420932', 'cg05667817', 'cg12597983', 'cg11422861', 'cg10474712')
# subgroup <- c('cg00126034', 'cg00577361', 'cg00847892', 'cg01096478', 'cg01745873', 'cg02105261', 'cg02582925', 'cg02964474', 'cg04420932', 'cg05667817', 'cg06529439', 'cg07003587', 'cg07063325', 'cg07704144', 'cg07782002', 'cg07839336', 'cg09909351', 'cg10456069', 'cg10474712', 'cg10587741', 'cg11422861', 'cg12597983', 'cg13554071', 'cg15730116', 'cg17934367', 'cg18497550', 'cg19456996', 'cg20694619', 'cg22106401', 'cg22466771', 'cg25117123', 'cg26045220', 'cg26680520', 'cg26885578')
# subgroup <- c('cg00126034', 'cg00847892', 'cg01096478', 'cg01745873', 'cg02068989', 'cg02105261', 'cg02582925', 'cg02964474', 'cg04420932', 'cg05667817', 'cg06529439', 'cg07003587', 'cg07063325', 'cg07704144', 'cg07782002', 'cg07839336', 'cg09909351', 'cg10456069', 'cg10474712', 'cg10587741', 'cg11422861', 'cg12597983', 'cg13554071', 'cg15730116', 'cg17934367', 'cg18497550', 'cg19456996', 'cg20694619', 'cg22106401', 'cg22466771', 'cg25117123', 'cg26045220', 'cg26680520', 'cg26885578')
subgroup <- c('cg00126034', 'cg00776782', 'cg00847892', 'cg01096478', 'cg01745873', 'cg02105261', 'cg02582925', 'cg02964474', 'cg04420932', 'cg05667817', 'cg06529439', 'cg07003587', 'cg07063325', 'cg07704144', 'cg07782002', 'cg07839336', 'cg09909351', 'cg10456069', 'cg10474712', 'cg10587741', 'cg11422861', 'cg12597983', 'cg13554071', 'cg15730116', 'cg17934367', 'cg18497550', 'cg19456996', 'cg20694619', 'cg22106401', 'cg22466771', 'cg25117123', 'cg26045220', 'cg26680520', 'cg26885578')

# ssdp 2
# subgroup <- c('cg10892587', 'cg00396625', 'cg18581221', 'cg25117123', 'cg13711457', 'cg26457665', 'cg27249304', 'cg16821694', 'cg05667817', 'cg12666166', 'cg09889850', 'cg10024437', 'cg12384499', 'cg02105261', 'cg18497550', 'cg10590400', 'cg25814383', 'cg21775899', 'cg06080893', 'cg11422861')
# subgroup <- c('cg10892587', 'cg13711457', 'cg25503086', 'cg26457665', 'cg16068620', 'cg08303612', 'cg17200920', 'cg09500200', 'cg27249304', 'cg05667817', 'cg09934894', 'cg10024437', 'cg16228087', 'cg00106744', 'cg05339472')
# subgroup <- c('cg10892587', 'cg06529439', 'cg25117123', 'cg26457665', 'cg19456996', 'cg17200920', 'cg09909351', 'cg27249304', 'cg18497550', 'cg21775899', 'cg11422861')
# subgroup <- c('cg19844724', 'cg06918887', 'cg22466620', 'cg08303612', 'cg05184917', 'cg10590400', 'cg21775899')
# subgroup <- c('cg26521784', 'cg06529439', 'cg06918887', 'cg09909351', 'cg17270843', 'cg07899956')

kte <- topGSA(gometh(sig.cpg=subgroup, all.cpg=sigCpGs, plot.bias=FALSE, collection='KEGG'), number=10)
gte <- topGSA(gometh(sig.cpg=subgroup, all.cpg=sigCpGs, plot.bias=FALSE), number=10)

kte_out <- paste(sample, "/filtered_results/sg_anno/", sample, "_bs3_kte.csv", sep="")
write.table(kte, file=kte_out, sep=",", row.names=TRUE)
# subgroups_anno filtered_results
gte_out <- paste(sample, "/filtered_results/sg_anno/", sample, "_bs3_gte.csv", sep="")
write.table(gte, file=gte_out, sep=",", row.names=TRUE)


# Top 10 GO categories
topGSA(gst, number=10)

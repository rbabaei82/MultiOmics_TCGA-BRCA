library(curatedTCGAData)
library(MultiAssayExperiment)
library(GenomicRanges)
library(RaggedExperiment)
library(GenomicDataCommons)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(TCGAutils)
library(UpSetR)
library(mirbase.db)
library(AnnotationFilter)
library(EnsDb.Hsapiens.v86)
library(survival)
library(survminer)
library(pheatmap)
library(ggplot2)

curatedTCGAData(diseaseCode = "BRCA")
BRCAmae <- curatedTCGAData("BRCA", c("GISTIC_Peaks","Mutation", "RNASeq2GeneNorm", "RPPAArray"),
                           dry.run = FALSE)


BRCAmae
dim(colData(BRCAmae))
head(colnames(colData(BRCAmae)))
head(metadata(colData(BRCAmae))[["subtypes"]])
colData(BRCAmae)[1:4,1:4]
table(BRCAmae$histological_type)

summary(complete.cases(BRCAmae)) 
listAssays <- assays(BRCAmae)

Mutation <- listAssays[[1]]
colnames(Mutation) <- gsub("(\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[0-9]+)$", "", colnames(Mutation))

RNASeq <- listAssays[[2]]
colnames(RNASeq) <- gsub("(\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[0-9]+)$", "", colnames(RNASeq))

RPPA <- listAssays[[3]]
dfRPPA <- data.frame(matrix(nrow = 226))
for(i in 1:937){dfRPPA <- cbind(dfRPPA, RPPA@listData[i])}
dfRPPA <- dfRPPA[,-1]
rownames(dfRPPA) <- RPPA@rownames
colnames(dfRPPA) <- gsub("(\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[0-9]+)$", "", colnames(dfRPPA))
rownames(dfRPPA) <- gsub("\\-|\\_", ".", rownames(dfRPPA))


copyNumber <- listAssays[[4]]
colnames(copyNumber) <- gsub("(\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[0-9]+)$", "", colnames(copyNumber))



brcaMatched <- intersectColumns(BRCAmae) # select complete cases and rearrange each ExperimentList element so its columns correspond exactly to rows of colData in the same order
colnames(brcaMatched)

# common row names across assays
brcaMatched <- intersectRows(brcaMatched) # keeps only rows that are common to each assay, and aligns them in identical order
rownames(brcamatched)


upsetSamples(BRCAmae)# 780 samples have all assays and three are missing RNA seq

# kaplan-meier plot for overall survival stratified by nodal stage
Surv(BRCAmae$days_to_death, BRCAmae$vital_status)# The colData provides clinical data for things like a Kaplan-Meier plot for overall survival stratified by nodal stage.

## remove patients missing overall survival information
BRCAmaesurv <- BRCAmae[, complete.cases(BRCAmae$days_to_death, BRCAmae$vital_status), ]


fit <- survfit(Surv(days_to_death, vital_status)~pathology_N_stage, data = colData(BRCAmaesurv))
ggsurvplot(fit, data = colData(BRCAmaesurv), risk.table = T)


###### correlation between copy number and RNA-Seq
subAssay1 <- BRCAmae[,, c("BRCA_RNASeq2GeneNorm-20160128" ,"BRCA_GISTIC_Peaks-20160128")]

subAssay1 <- intersectColumns(subAssay1)# keep samples with both assays available
subAssay1 <- intersectRows(subAssay1)

subAssay1_list <- assays(subAssay1) # list of matrices

subAssay1_list[[1]] <- log2(subAssay1_list[[1]] +1) # log2 transform of RNA-Seq assay

# transpose genes in column
subAssay1_list <- lapply(subAssay1_list, t)

## correlation between assay
corres <- cor(subAssay1_list[[1]], subAssay1_list[[2]])## incompatible dimenssions
hist(diag(corres))

hist(corres[upper.tri(corres)])


which.max(diag(corres))## input for next line
df <- wideFormat(subAssay1["from which.max",,])
boxplot(RNASeq2GeneNorm_EIF4E ~ gistict_EIF4E,
        data=df, varwidth=TRUE,
        xlab="GISTIC Relative Copy Number Call",
        ylab="RNA-seq counts")


########## Correlated Principal Components
# custom function for PCA
doPCA <- function(x, ncomp = 10, dolog =FALSE, center = TRUE, scale. = TRUE){
  if(dolog){
    x <- log2(x + 1)
  }
  pc = prcomp(x, center = center, scale. = scale.)
  return(t(pc$rotation[, 1:10]))
}

pcaBRCA <- intersectColumns(BRCAmae)
pcaBRCA <- c(pcaBRCA, gisticPCA = doPCA(listAssays[[1]], center = FALSE, scale. = FALSE), mapFrom = 1L)




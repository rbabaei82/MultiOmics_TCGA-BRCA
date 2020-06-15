library(curatedTCGAData)
#library(MultiAssayExperiment)
#library(GenomicRanges)
#library(RaggedExperiment)
#library(GenomicDataCommons)
#library(SummarizedExperiment)
#library(SingleCellExperiment)
#library(TCGAutils)
#library(UpSetR)
#library(mirbase.db)
#library(AnnotationFilter)
#library(EnsDb.Hsapiens.v86)
#library(survival)
#library(survminer)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(tictoc)


curatedTCGAData(diseaseCode = "BRCA")
BRCAmae <- curatedTCGAData("BRCA", c("GISTIC_ThresholdedByGene","Mutation", "RNASeq2GeneNorm", "RPPAArray"),
                           dry.run = FALSE)

#####
BRCAmae
dim(colData(BRCAmae))
head(colnames(colData(BRCAmae)))
head(metadata(colData(BRCAmae))[["subtypes"]])
colData(BRCAmae)[1:4,1:4]
table(BRCAmae$histological_type)# infiltrating ductal and lobular carcinomas are the most prevalent type, we select only these for analysis

## col_data for graph annotation
dfcolData <- data.frame(matrix(nrow = nrow(colData(BRCAmae))))
for(i in 1:ncol(colData(BRCAmae))){dfcolData <- cbind(dfcolData, BRCAmae@colData@listData[i])}
dfcolData <- dfcolData[,-1]
rownames(dfcolData) <- BRCAmae@colData@rownames
colnames(dfcolData) <- gsub("\\-|\\_", ".", colnames(dfcolData))
rownames(dfcolData) <- gsub("\\-|\\_", ".", rownames(dfcolData))

# select patients with only ductal or lobular carcinoma
dfcolData <- dfcolData %>% filter(histological.type == "infiltrating ductal carcinoma" | histological.type == "infiltrating lobular carcinoma")


######
summary(complete.cases(BRCAmae)) 
listAssays <- assays(BRCAmae)

#######
copyNumber <- listAssays[[1]]
colnames(copyNumber) <- gsub("(\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[0-9]+)$", "", colnames(copyNumber))
copyNumberdf <- as.data.frame(copyNumber)
colnames(copyNumberdf) <- gsub("\\-|\\_", ".", colnames(copyNumberdf))

# select patients with only ductal or lobular carcinoma
copyNumberdf <- copyNumberdf[,intersect(colnames(copyNumberdf), rownames(dfcolData))]
dim(copyNumberdf)
length(unique(colnames(copyNumberdf))) ## there are unique variables
sum(is.na(copyNumberdf)) # there is no missing value





Mutation <- listAssays[[2]]
colnames(Mutation) <- gsub("(\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[0-9]+)$", "", colnames(Mutation))

# reshape the sparse matrix and set binary values for mutations
Mutationdf <- as.data.frame(Mutation)
tmp <- data.frame(GeneSymbol = rep(NA, nrow(Mutationdf)))
for(i in 1:ncol(Mutationdf)){
  tmp$GeneSymbol <- coalesce(Mutationdf[,i], tmp$GeneSymbol)
  Mutationdf$GeneSymbol <- tmp$GeneSymbol
}

sum(is.na(Mutationdf$GeneSymbol))# no missing value


tmp <- data.frame(GeneSymbol = unique(Mutationdf$GeneSymbol))# select unique gene symbols
sum(is.na(tmp$GeneSymbol))# no NA

library(foreach)
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

tmp_data <- foreach(i = c(1:993), j = c(1:993))%dopar%{
  tmp[,paste0("V", i)] <- ifelse(tmp$GeneSymbol%in%Mutationdf[,j], 1, 0)
}
stopCluster(cl)

tmp_df <- data.frame(matrix(unlist(tmp_data), ncol = length(tmp_data)))# convert to data frame
colnames(tmp_df) <- colnames(Mutationdf[,1:993])
rownames(tmp_df) <- tmp$GeneSymbol
Mutationdf <- tmp_df
colnames(Mutationdf) <- gsub("\\-|\\_", ".", colnames(Mutationdf))

# select patients with only ductal or lobular carcinoma
Mutationdf <- Mutationdf[,intersect(colnames(Mutationdf), rownames(dfcolData))]
dim(Mutationdf)

length(unique(colnames(Mutationdf)))# unique variables
rm(tmp, tmp_data, tmp_df)




RNASeq <- listAssays[[3]]
colnames(RNASeq) <- gsub("(\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[0-9]+)$", "", colnames(RNASeq))
# log transformation
RNASeq <- log2(RNASeq +1)

RNASeqdf <- as.data.frame(RNASeq)


dim(RNASeqdf)
length(unique(colnames(RNASeqdf)))# there are duplicated variables

RNASeqdf <- RNASeqdf[,unique(colnames(RNASeqdf))] # remove duplicated variables
colnames(RNASeqdf) <- gsub("\\-|\\_", ".", colnames(RNASeqdf))

# select patients with only ductal or lobular carcinoma
RNASeqdf <- RNASeqdf[,intersect(colnames(RNASeqdf), rownames(dfcolData))]
dim(RNASeqdf)



RPPA <- listAssays[[4]]
dfRPPA <- data.frame(matrix(nrow = 226))
for(i in 1:937){dfRPPA <- cbind(dfRPPA, RPPA@listData[i])}
dfRPPA <- dfRPPA[,-1]
rownames(dfRPPA) <- RPPA@rownames
colnames(dfRPPA) <- gsub("(\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[0-9]+)$", "", colnames(dfRPPA))
rownames(dfRPPA) <- gsub("\\-|\\_", ".", rownames(dfRPPA))

length(unique(colnames(dfRPPA)))# there are duplicated variables
dfRPPA <- dfRPPA[,unique(colnames(dfRPPA))]
colnames(dfRPPA) <- gsub("\\-|\\_", ".", colnames(dfRPPA))

# select patients with only ductal or lobular carcinoma
dfRPPA <- dfRPPA[,intersect(colnames(dfRPPA), rownames(dfcolData))]
dim(dfRPPA)


## list of data frames with same variables
dataList <- list(CN =copyNumberdf,MU= Mutationdf,RS= RNASeqdf,PP= dfRPPA)
intersectID <- Reduce(intersect, lapply(dataList, colnames))
dataList <- lapply(dataList, function(x) x[intersectID])

## filter colData by shared variables
dfcolData <- dfcolData[intersect(rownames(dfcolData), colnames(dataList$CN)),]

## list of data frames with same genes (excluding RPPA)
omicsList <- dataList[1:3]
intersectID <- Reduce(intersect, lapply(omicsList, rownames))
omicsList <- lapply(omicsList, function(x) x[intersectID,])

dfRPPA <- dataList$PP

save(omicsList, dfRPPA, dfcolData, file = "Data.RData")

###
# clustering in relation to histological types
anno_col <- data.frame(type = as.factor(dfcolData$histological.type))
rownames(anno_col) <- rownames(dfcolData)

### heatmap of gene expression
tic("heat map GE")
pheatmap(as.matrix(omicsList$RS), annotation_col = anno_col,
         show_rownames = F, show_colnames =  F, main = "Gene Expression") # clustered by genes , not by subtypes

toc() # 34 minutes


# pheatmap(as.matrix(omicsList$MU), annotation_col = anno_col,
#         show_rownames = F, show_colnames =  F, main = "Mutation")

pheatmap(as.matrix(omicsList$CN), annotation_col = anno_col,
         show_rownames = F, show_colnames =  F, main = "Copy Number")



# correlation between gene expression and copy number
corr <- cor(t(omicsList$CN), t(omicsList$RS))
which.max(diag(corr))

# top ten correlated genes
top10 <- order(diag(corr), decreasing = T)[1:10]
k <- arrayInd(top10, dim(corr))
topGenes <- rownames(corr)[k[,1]]
topGenes
# boxplot of the top gene
dfTop <- data.frame(RNASeq = t(omicsList$RS[topGenes[1],]), CopyNumber = t(omicsList$CN[topGenes[1],]))
colnames(dfTop) <- c("RNASeq", "CopyNumber")
boxplot(RNASeq~CopyNumber, data = dfTop, xlab = "GISTIC Relative Copy Number Call", ylab = "RNA-Seq log FC", main = "KIAA1967")


################################

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





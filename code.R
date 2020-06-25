required_packages = c("curatedTCGAData", "iClusterPlus", "GenomicRanges", "foreach", "doParallel",
                      "pheatmap", "ggplot2", "tidyverse", "tictoc", "missForest", "NMF", "enrichR")
suppressPackageStartupMessages(lapply(required_packages, require, character.only =T))

library(curatedTCGAData)
library(iClusterPlus)
#library(MultiAssayExperiment)
library(GenomicRanges)
library(foreach)
library(doParallel)
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


#########
upsetSamples(BRCAmae)# 780 samples have all assays and three are missing RNA seq

# kaplan-meier plot for overall survival stratified by nodal stage
Surv(BRCAmae$days_to_death, BRCAmae$vital_status)# The colData provides clinical data for things like a Kaplan-Meier plot for overall survival stratified by nodal stage.

## remove patients missing overall survival information
BRCAmaesurv <- BRCAmae[, complete.cases(BRCAmae$days_to_death, BRCAmae$vital_status), ]


fit <- survfit(Surv(days_to_death, vital_status)~pathology_N_stage, data = colData(BRCAmaesurv))
ggsurvplot(fit, data = colData(BRCAmaesurv), risk.table = T)



## col_data for graph annotation
dfcolData <- data.frame(matrix(nrow = nrow(colData(BRCAmae))))
for(i in 1:ncol(colData(BRCAmae))){dfcolData <- cbind(dfcolData, BRCAmae@colData@listData[i])}
dfcolData <- dfcolData[,-1]
rownames(dfcolData) <- BRCAmae@colData@rownames
colnames(dfcolData) <- gsub("\\-|\\_", ".", colnames(dfcolData))
rownames(dfcolData) <- gsub("\\-|\\_", ".", rownames(dfcolData))

toBeRemoved <- which(data.frame(map(dfcolData, ~sum(is.na(.))))> nrow(dfcolData)*0.3)# just to reduce the dimension, not for imputing
dfcolData <- dfcolData[,-(toBeRemoved)]

# select patients with only ductal or lobular carcinoma
 dfcolData <- dfcolData %>% filter(histological.type == "infiltrating ductal carcinoma" | histological.type == "infiltrating lobular carcinoma")


######
summary(complete.cases(BRCAmae)) 
listAssays <- assays(BRCAmae)

#######
## RNA-Seq
RNASeq <- listAssays[[3]]
colnames(RNASeq) <- gsub("(\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[0-9]+)$", "", colnames(RNASeq))
# log transformation
RNASeq <- log2(RNASeq +1)

dim(RNASeq)
length(unique(colnames(RNASeq)))# there are duplicated variables

RNASeq <- RNASeq[,unique(colnames(RNASeq))] # remove duplicated variables
colnames(RNASeq) <- gsub("\\-|\\_", ".", colnames(RNASeq))

# select patients with only ductal or lobular carcinoma
 RNASeq <- RNASeq[,intersect(colnames(RNASeq), rownames(dfcolData))]
# dim(RNASeq)

RNASeq <- t(RNASeq)


## remove low variable genes
SDs <- apply(RNASeq,2,sd)
orderSD <- order(SDs, decreasing = T)[1:2000]
RNASeq <- RNASeq[,orderSD]
dim(RNASeq)




##copy number
copyNumber <- listAssays[[1]]
colnames(copyNumber) <- gsub("(\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[0-9]+)$", "", colnames(copyNumber))
colnames(copyNumber) <- gsub("\\-|\\_", ".", colnames(copyNumber))

# select patients with only ductal or lobular carcinoma
 copyNumber <- copyNumber[,intersect(colnames(copyNumber), rownames(dfcolData))]
 dim(copyNumber)
length(unique(colnames(copyNumber))) ## there are unique variables
sum(is.na(copyNumber)) # there is no missing value
copyNumber <- t(copyNumber)
copyNumber <- copyNumber[,intersect(colnames(copyNumber), colnames(RNASeq))] # relative to RNASeq
dim(copyNumber)

## mutation
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
Mutation <- t(as.matrix(Mutationdf))

# filter out genes with too few mutations (less than 5%)
mut.rate <- apply(Mutation, 2, mean)
Mutation <- Mutation[,which(mut.rate>0.02)]
dim(Mutation)


## PP
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
mRPPA <- t(as.matrix(dfRPPA))
dim(mRPPA)
sum(is.na(mRPPA))# there are missing values
toBeRemoved <- which(data.frame(map(as.data.frame(mRPPA), ~sum(is.na(.))))> nrow(mRPPA)/2) # three variables with more than 50% missing values, to be removed

mRPPA <- mRPPA[,-c(toBeRemoved)]
# impute the rest missing values
set.seed(4578)
registerDoParallel(cores = 3)
impute_forest <- missForest(data.frame(mRPPA), maxiter = 10, ntree = 100, parallelize = "forests")
stopImplicitCluster()

impute_forest$OOBerror # 90% accuracy
mRPPA <- as.matrix(impute_forest$ximp)
sum(is.na(mRPPA))# no missing value



rm(BRCAmae, cl, dfRPPA, listAssays, Mutationdf, RPPA, mut.rate, SDs)

#### list of omics with same patients
dataList <- list(CN =copyNumber,MU= Mutation,RS= RNASeq,PP= mRPPA)
intersectID <- Reduce(intersect, lapply(dataList, rownames))
dataList <- lapply(dataList, function(x) x[intersectID,]) # data is selected for 705 patients

## filter colData by shared rows
dfcolData <- dfcolData[intersect(rownames(dfcolData), rownames(dataList$CN)),]


## list of data frames with same genes (excluding RPPA)
#omicsList <- dataList[1:3]
#intersectID <- Reduce(intersect, lapply(omicsList, rownames))
#omicsList <- lapply(omicsList, function(x) x[intersectID,])
#dfRPPA <- dataList$PP

#save(omicsList, dfRPPA, dfcolData, file = "Data.RData")

########
# clustering in relation to clinical data
anno_col <- data.frame(type = as.factor(dfcolData$histological.type))
rownames(anno_col) <- rownames(dfcolData)

### heatmap of gene expression
# tic("heatmap GE")
# pheatmap(t(dataList$RS), annotation_col = anno_col,
#          show_rownames = F, show_colnames =  F, main = "Gene Expression") # clustered by genes , not by subtypes
# 
# toc()  
# 
# 
# tic("heatmap mutation")
# pheatmap(t(dataList$MU), annotation_col = anno_col,
#         show_rownames = F, show_colnames =  F, main = "Mutation")
# 
# toc()
# 
# 
# tic("heatmap copy number")
# pheatmap(t(dataList$CN), annotation_col = anno_col,
#          show_rownames = F, show_colnames =  F, main = "Copy Number")
# toc()
# 
# 
# tic("heatmap Phosphopeptides")
# pheatmap(t(dataList$PP), annotation_col = anno_col,
#          show_rownames = F, show_colnames =  F, main = "Phosphopeptides")
#  toc()

# correlation between gene expression and copy number
tic("Correlation")
corr <- cor(dataList$CN, dataList$RS)
toc()
which.max(diag(corr))

# top ten correlated genes
top10 <- order(diag(corr), decreasing = T)[1:10]
k <- arrayInd(top10, dim(corr))
topGenes <- rownames(corr)[k[,1]]
topGenes
# boxplot of the top gene
dfTop <- data.frame(RNASeq = dataList$RS[,topGenes[1]], CopyNumber = dataList$CN[,topGenes[1]])
colnames(dfTop) <- c("RNASeq", "CopyNumber")
boxplot(RNASeq~CopyNumber, data = dfTop, xlab = "GISTIC Relative Copy Number Call", ylab = "RNA-Seq log FC", main = topGenes[1])



####### Matrix factorization


# check if all samples are in the same order for the four datasets
all(rownames(dataList$CN) == rownames(dataList$MU)) # TRUE
all(rownames(dataList$CN) == rownames(dataList$RS)) # TRUE
all(rownames(dataList$CN) == rownames(dataList$PP)) # TRUE


#################
# set.seed(5214)
# data()
# cl <- makeCluster(3)
# registerDoParallel(cl)
# foreach(k = c(1:5))%dopar%{
#   require(iClusterPlus)
# i.cluster <- tune.iClusterPlus(cpus = 1, dt1 = dataList$RS, dt2 = dataList$MU, dt3 = dataList$CN, dt4 = dataList$PP,
#                           type = c("gaussian", "binomial", "multinomial", "gaussian"), K = k, n.lambda = 307,
#                           scale.lambda = c(1, 0.05, 1, 1), maxiter = 20)
# save(i.cluster, file = paste("i.cluster.k", k, ".RData", sep = ""))
# }
# stopCluster(cl)
# date()   
###############
# i.cluster <- iClusterPlus(dataList$RS, dataList$MU, dataList$CN, dataList$PP,
#                           type = c("gaussian", "binomial", "multinomial", "gaussian"), K=1,
#                           lambda =c(0.9, 0.03, 0.1, 0.6), maxiter = 10)
# 
# 
# icluster.z <- i.cluster$meanZ
# icluster.ws <- i.cluster$beta
# icluster.clusters <- i.cluster$clusters
# 
# icp_df <- as.data.frame(icluster.z)
# colnames(icp_df) <- c("dim1", "dim2")
# rownames(icp_df) <- rownames(dataList$CN)
# 
# icp_df$subtype <- factor(anno_col$type)
# ggplot(icp_df, aes(x=dim1, y = dim2, color = subtype)) + geom_point() + ggtitle("Scatter plot of iCluster+ factors") # ignore


# K-means clustering
# rownames(icluster.z) <- rownames(anno_col)
# clusters <- kmeans(icluster.z, 2)$cluster
# names(clusters) <- rownames(icluster.z)
# anno_cluster <- data.frame(Cluster = factor(clusters), subtype = factor(anno_col$type))
# pheatmap(t(icluster.z[order(clusters),]),
#          cluster_cols = F, cluster_rows = F, show_colnames = F, annotation_col = anno_cluster,
#          main = "iCluster factors with clusters and Histolofical types")

# #### Feature selection
# features <- alist()
# features[[1]] <- colnames(dataList$RS)
# features[[2]] <- colnames(dataList$MU)
# features[[3]] <- colnames(dataList$CN)
# features[[4]] <- colnames(dataList$PP)
# 
# sig.features <- alist()
# for(i in 1:4){
#   rowSum = apply(abs(i.cluster$beta[[i]]),1,sum)
#   upper = quantile(rowSum, probs = 0.75)
#   sig.features[[i]] = (features[[i]])[which(rowSum>upper)]
# }
# names(sig.features) = c("expression", "mutation", "copy number", "Phospho peptides")


############ NMF ( non-negative multi factor)
# Feature-normalize the data
CN.featnorm <- t(dataList$CN)/rowSums(t(abs(dataList$CN)))# take abs to avoid rowsums of zero
MU.featnorm <- t(dataList$MU)/rowSums(t(dataList$MU))
RS.featnorm <- t(dataList$RS)/rowSums(t(dataList$RS))
PP.featnorm <- t(dataList$PP)/rowSums(t(abs(dataList$PP)))# take abs to avoid rowsums of zero

# Normalize by each omics type's frobenius norm
CN.featnorm.frobnorm <- CN.featnorm/norm(as.matrix(CN.featnorm), type = "F")
MU.featnorm.frobnorm <- MU.featnorm/norm(as.matrix(MU.featnorm), type = "F")
RS.featnorm.frobnorm <- RS.featnorm/norm(as.matrix(RS.featnorm), type = "F")
PP.featnorm.frobnorm <- PP.featnorm/norm(as.matrix(PP.featnorm), type = "F")

# Split the features of the CN and PP matrix into two non-negative features each
split_neg_columns <- function(x) {
  new_cols <- list()
  for(i in seq_len(dim(x)[2])) {
    new_cols[[paste0(colnames(x)[i],'+')]] <- sapply(X = x[,i], 
                                                     function(x) max(0,x))
    new_cols[[paste0(colnames(x)[i],'-')]] <- sapply(X = -x[,i],
                                                     function(x) max(0,x))
  }
  new_cols
  return(do.call(cbind, new_cols))
}

CN.nonneg <- t(split_neg_columns(t(CN.featnorm.frobnorm)))
PP.nonneg <- t(split_neg_columns(t(PP.featnorm.frobnorm)))


r.nmf <- nmf(t(rbind(CN.nonneg, MU.featnorm.frobnorm, RS.featnorm.frobnorm, PP.nonneg)), 2, method = "Frobenius")
nmf.h <- basis(r.nmf)
nmf.w <- NMF::coef(r.nmf)


nmf_df <- as.data.frame(nmf.h)
colnames(nmf_df) <- c("dim1", "dim2")
nmf_df$subtypes <- factor(dfcolData[rownames(nmf_df),]$histological.type)
table(nmf_df$subtypes)

ggplot(nmf_df, aes(x = dim1, y = dim2, color = subtypes)) + geom_point() + ggtitle("2-component NMF for multi-omics integratoin")
 
pheatmap(t(nmf_df[,1:2]),annotation_col = anno_col, show_colnames=FALSE, main="Heatmap of 2-component NMF")

# One-hot clustering

nmf.cluster <- max.col(nmf.h)
names(nmf.cluster) <- rownames(nmf.h)
anno_nmf_cl <- data.frame(nmf.cluster = factor(nmf.cluster), subtypes = factor(dfcolData[rownames(nmf.h),]$histological.type))

pheatmap(t(nmf.h[order(nmf.cluster),]), cluster_cols = F, cluster_rows = F, annotation_col = anno_nmf_cl,
         show_colnames = F, main = "Joint NMF factors with clusters and subtypes")

### biological interpretation of latent factors
nmfw <- t(nmf.w)
data_anno <- data.frame(omics=c(rep('expression',dim(dataList$RS)[2]), rep('mut',dim(dataList$MU)[2]),
                                rep('cnv',dim(CN.nonneg)[1]), rep("phosphopeptides", dim(PP.nonneg)[1])))
rownames(data_anno) <- c(rownames(t(dataList$RS)), paste0("mut:", rownames(t(dataList$MU))),paste0("CN:", rownames(CN.nonneg)), rownames(PP.nonneg))
rownames(nmfw) <- rownames(data_anno)

pheatmap::pheatmap(nmfw,
                   cluster_cols = FALSE,
                   annotation_row = data_anno,
                   main="NMF coefficients",
                   clustering_distance_rows = "manhattan",
                   fontsize_row = 1)


## feature inspection
topExp <- data.frame(nmfw[colnames(dataList$RS),])
colnames(topExp) <- c("F1","F2")
knitr::kable(topExp[order(topExp$F1, decreasing = T),][1:3,])
knitr::kable(topExp[order(topExp$F2, decreasing = T),][1:3,])



topCN <- data.frame(nmfw[paste0("CN:", rownames(CN.nonneg)),])
colnames(topCN) <- c("F1","F2")
knitr::kable(topCN[order(topCN$F1, decreasing = T),][1:3,])
knitr::kable(topCN[order(topCN$F2, decreasing = T),][1:3,])



topMU <- data.frame(nmfw[paste0("mut:", colnames(dataList$MU)),])
colnames(topMU) <- c("F1","F2")
knitr::kable(topMU[order(topMU$F1, decreasing = T),][1:3,])
knitr::kable(topMU[order(topMU$F2, decreasing = T),][1:3,])


topPP <- data.frame(nmfw[rownames(PP.nonneg),])
colnames(topPP) <- c("F1","F2")
knitr::kable(topPP[order(topPP$F1, decreasing = T),][1:3,])
knitr::kable(topPP[order(topPP$F2, decreasing = T),][1:3,])


#### Enrichment analysis
# gene expression affect only the first factor

genes.factor <- rownames(topExp[order(topExp$F1, decreasing = T),][1:100,])# analysis for top 100 genes
go.factor <- enrichr(genes.factor, databases = c("GO_Biological_Process_2018"))$GO_Biological_Process_2018

ggplot(head(go.factor), aes(y=-log10(P.value), x = Term))+
  geom_bar(stat = "identity") + coord_flip() + ggtitle("Top GO-Terms\nassociated with\nNMF factor")

################################







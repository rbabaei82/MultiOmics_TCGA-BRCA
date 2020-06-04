source("Module_A.R")
source("Module_B.R")
clinBiosp <- DownloadBiospecimenClinicalData(cancerType = "BRCA", saveFolderName = "D:/DataScience/Profile/Multiomics_BRCA/MultiOmics_TCGA-BRCA/BiospecimenClinicalData",
                                             outputFileName = "test")
CNA <- DownloadCNAData(cancerType = "BRCA", assayPlatform = "cna_cnv.hg19",
                       saveFolderName ="D:/DataScience/Profile/Multiomics_BRCA/MultiOmics_TCGA-BRCA/CNA")
RNASeq <- DownloadRNASeqData(cancerType = "BRCA", assayPlatform = "gene.normalized_RNAseq",
                             saveFolderName = "D:/DataScience/Profile/Multiomics_BRCA/MultiOmics_TCGA-BRCA/RNAseq")
somaticMutation <- DownloadSomaticMutationData(cancerType = "BRCA", assayPlatform = "somaticMutation_DNAseq",
                                               saveFolderName ="D:/DataScience/Profile/Multiomics_BRCA/MultiOmics_TCGA-BRCA/somaticMutation")
phosphoProteom <- DownloadCPTACData(cancerType = "BRCA", assayPlatform = "phosphoproteome_iTRAQ",
                                    saveFolderName ="D:/DataScience/Profile/Multiomics_BRCA/MultiOmics_TCGA-BRCA/iTRAQ")

## process
GeneLevel.CNA <- ProcessCNAData(inputFilePath = CNA[1], outputFileName = "GeneLevelCNA", outputFileFolder = "D:/DataScience/Profile/Multiomics_BRCA/MultiOmics_TCGA-BRCA",
                                refGenomeFile = "SupportingFiles/Hg19GenePosition.txt")
processedRNAseq <- ProcessRNASeqData(inputFilePath = RNASeq[1], outputFileName = "RNASeqProcessed",
                                     outputFileFolder = "D:/DataScience/Profile/Multiomics_BRCA/MultiOmics_TCGA-BRCA", dataType = "geneExp", verType = "RNASeqV2")
procesedMutation <- ProcessSomaticMutationData(inputFilePath = somaticMutation[1], outputFileName = "MutationProcessed",
                                               outputFileFolder = "D:/DataScience/Profile/Multiomics_BRCA/MultiOmics_TCGA-BRCA")
processedProteome <- ProcessCPTACData(inputFilePath = phosphoProteom[1], outputFileName = "ProteomeProcessed",
                                      outputFileFolder = "D:/DataScience/Profile/Multiomics_BRCA/MultiOmics_TCGA-BRCA")

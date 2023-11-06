# ---------- Part 0: Preparations ----------
rm(list=ls()); gc()
setwd('D:/Projects/CRC_expData_LDA_Carcinogenesis')
set.seed(42)
random.seeds <- sample(1:1000000, 10000)


# ---------- Part I: Select colon samples from XENA data ----------
processExpAndphtData <- function(phenotypeFile, expDataFile, # pht = phenotype
                                 primarySite=NULL, sampleType=NULL,skip=0,
                                 sampleID.ncol=1, sampleType.ncol=2){
  library(data.table)
  library(stringr)
  phenotype <- data.frame(fread(phenotypeFile), check.names = F, stringsAsFactors = F)
  phenotypeSelected <- na.omit(phenotype)
  if (!is.null(primarySite)) {
    phenotypeSelected <-
      phenotypeSelected[phenotypeSelected$primarySite %in% primarySite,]
  }
  if (!is.null(sampleType)) {
    phenotypeSelected <- subset(phenotypeSelected, sampleType %in% sampleType)
  }
  phenotypeSelected <- phenotypeSelected[, c(sampleID.ncol, sampleType.ncol)]
  colnames(phenotypeSelected) <- c('sampleID', 'sampleType')
  phenotypeSelected$sampleID <- str_replace_all(phenotypeSelected$sampleID, '-', '\\.')
  
  expData <- data.frame(fread(expDataFile, stringsAsFactors=F, header=T, sep='\t',
                              nThread = 12,skip = skip),stringsAsFactors = F)
  #fread() get data.table class, data frame transformation is needed
  #for downstream analysis.
  expData <- na.omit(expData)
  sampleIDsSelected <- intersect(phenotypeSelected$sampleID, colnames(expData))
  phenotypeSelected <- phenotypeSelected[phenotypeSelected$sampleID %in% sampleIDsSelected,]
  expDataSelected <- expData[, c(1, which(colnames(expData) %in% sampleIDsSelected))]
  rownames(expDataSelected) <- expDataSelected[, 1]
  expDataSelected <- expDataSelected[, -1]
  expDataSelected <- t(expDataSelected)
  phenotypeSelected <- phenotypeSelected[match(rownames(expDataSelected),
                                               phenotypeSelected$sampleID), ]
  return(list(phenotype=phenotypeSelected, expData=expDataSelected))
}

getSampleNumbers <- function(sampleInfo){
  library(stringr)
  sampleInfo[,2] <- as.character(sampleInfo[,2])
  nNor <- paste0('Nor (n=', sum(str_detect(sampleInfo[,2], 'Nor')), ')')
  nAdj <- paste0('Adj (n=', sum(str_detect(sampleInfo[,2], 'Adj')), ')')
  nCan <- paste0('Can (n=', sum(str_detect(sampleInfo[,2], 'Can')), ')')
  sampleInfo[str_detect(sampleInfo[,2], 'Nor'),2] <- nNor
  sampleInfo[str_detect(sampleInfo[,2], 'Adj'),2] <- nAdj
  sampleInfo[str_detect(sampleInfo[,2], 'Can'),2] <- nCan
  sampleInfo[,2] <- factor(sampleInfo[,2], levels = c(nNor, nAdj, nCan))
  return(list(sampleInfo, c(nNor, nAdj, nCan)))
}

# Process XENA
expDataAndSampleInfo.XENA.CRC <-
  processExpAndphtData(
    'TcgaTargetGTEX_phenotype.txt',
    'E:/PublicData/DataFromXena/TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2',
    c('Colon', 'Rectum'), c('Normal Tissue', 'Primary Tumor', 'Solid Tissue Normal'),
    sampleID.ncol = 1, sampleType.ncol = 4)
sampleInfo.XENA.CRC <- expDataAndSampleInfo.XENA.CRC$phenotype
sampleInfo.XENA.CRC$sampleType <-
  ifelse(sampleInfo.XENA.CRC$sampleType=='Normal Tissue', 'Nor',
         ifelse(sampleInfo.XENA.CRC$sampleType=='Primary Tumor', 'Can', 'Adj'))
sampleInfo.XENA.CRC <- getSampleNumbers(sampleInfo.XENA.CRC)[[1]]
sampleInfo.XENA.CRC$sampleType <-
  factor(sampleInfo.XENA.CRC$sampleType,
         levels = getSampleNumbers(sampleInfo.XENA.CRC)[[2]])
sampleInfo.XENA.CRC <-
  sampleInfo.XENA.CRC[order(sampleInfo.XENA.CRC$sampleType),]
expData.XENA.CRC <- expDataAndSampleInfo.XENA.CRC$expData
expData.XENA.CRC <- expData.XENA.CRC[match(sampleInfo.XENA.CRC$sampleID,
                                           rownames(expData.XENA.CRC)),]
colnames(expData.XENA.CRC) <- str_sub(colnames(expData.XENA.CRC), 1, 15)
rm(expDataAndSampleInfo.XENA.CRC); gc()


# ---------- Part II: Clean the data from XENA ----------
cleanData <- function(expData, sampleInfo){
  # Remove the genes whose expression level are zero in more than 90% samples.
  # Remove the samples whose expression levels are zero in more than 90% genes.
  library(stringr)
  tumorIndex <- str_detect(sampleInfo$sampleType, 'Can')
  adjacentIndex <- str_detect(sampleInfo$sampleType, 'Adj')
  normalIndex <- str_detect(sampleInfo$sampleType, 'Nor')
  genesNotExp <-
    colSums(expData[tumorIndex,]!=0)<=0.1*sum(tumorIndex)|
    colSums(expData[adjacentIndex,]!=0)<=0.1*sum(adjacentIndex)|
    colSums(expData[normalIndex,]!=0)<=0.1*sum(normalIndex)
  samplesNotExp <- rowSums(expData==0) >= 0.9*dim(expData)[2]
  filteredExpData <- expData[!samplesNotExp, !genesNotExp]

  # Update the sample information.
  commonSampleIDs <- intersect(sampleInfo$sampleID, rownames(filteredExpData))
  filteredsampleInfo <- sampleInfo[sampleInfo$sampleID %in% commonSampleIDs,]
  filteredExpData <- filteredExpData[rownames(filteredExpData) %in% commonSampleIDs, ]
  return(list(sampleInfo=filteredsampleInfo, expData=filteredExpData))
}

# Clean XENA data
data.XENA.CRC.clean <- cleanData(expData.XENA.CRC, sampleInfo.XENA.CRC)
sampleInfo.XENA.CRC <- data.XENA.CRC.clean[[1]]
sampleInfo.XENA.CRC$sampleType <- as.character(sampleInfo.XENA.CRC$sampleType)
sampleInfo.XENA.CRC <- getSampleNumbers(sampleInfo.XENA.CRC)[[1]]
sampleInfo.XENA.CRC$sampleType <-
  factor(sampleInfo.XENA.CRC$sampleType,
         levels = getSampleNumbers(sampleInfo.XENA.CRC)[[2]])
expData.XENA.CRC <- data.XENA.CRC.clean[[2]]
rm(data.XENA.CRC.clean)


# ---------- Part III: LDA based on XENA data ----------
plotLDAmodel <- function(ldaModel, trainSet, trainGroup,
                         title, legend.position = c(0.75,0.05)){
  library(stringr)
  trainGroup[,2] <- as.character(trainGroup[,2])
  trainGroup <- getSampleNumbers(trainGroup)[[1]]
  trainGroup[,2] <-
    factor(trainGroup[,2],
           levels = getSampleNumbers(trainGroup)[[2]])
  z <- data.frame(trainSet %*% ldaModel$scaling,
                  'Sample Type'=trainGroup[,2],
                  Samples=1:nrow(trainGroup), check.names = F)
  library(ggplot2)
  LDAmodelPlot <-
    ggplot(data=z, mapping=aes(x=LD1, y=LD2, colour=`Sample Type`, shape=`Sample Type`)) +
    geom_point(size=0.8, alpha=0.7) + ggtitle(title) +
    theme(text=element_text(family = NULL),
          legend.position = legend.position, legend.justification = c(0, 0.1),
          legend.key.size = unit(3, 'pt'), legend.title = element_text(size = 14),
          legend.text = element_text(size=12), axis.title=element_text(size=14, face="bold"),
          axis.text=element_text(size=12, face="bold"),
          plot.title = element_text(face = "bold", size = 14, vjust = 1))
  return(LDAmodelPlot)
}

plotROC <- function(LDAmodel, testSet, testGroup, title){
  #windowsFonts(A=windowsFont("Times New Roman"))
  library(pROC)
  library(stringr)
  testGroup[,2] <- as.character(testGroup[,2])
  testGroup <- getSampleNumbers(testGroup)[[1]]
  testGroup[,2] <- factor(testGroup[,2], 
                          levels = getSampleNumbers(testGroup)[[2]])
  pred <- predict(LDAmodel, testSet)$posterior
  groups <- levels(testGroup[,2])
  ROCs <- list()
  for (i in 1:(length(groups)-1)) {
    for (j in (i+1):length(groups)) {
      roc.object <- roc(testGroup[,2], pred[,i],
                        levels = c(groups[i], groups[j]))
      AUC <- as.numeric(auc(roc.object))
      ROC <-
        ggroc(roc.object, legacy.axes=T) +
        xlab('1 - specificity') +
        ylab('Sensitivity') +
        ggtitle(paste0(title, '\n', groups[i], ' VS ', groups[j])) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        #    scale_y_continuous(expand=c(0,0)) +
        #    scale_x_continuous(expand=c(0,0)) +
        theme(axis.title=element_text(size=14, face="bold")) +
        theme(axis.text=element_text(size=12, face="bold")) +
        theme(plot.title = element_text(face = "bold",
                                        size = 14, vjust = 1)) +
        annotate('text', x=0.5, y=0.5, size=5, fontface="bold",
                 label=paste0('AUC=', round(AUC,digits = 4)))
      ROCs[[length(ROCs)+1]] <- ROC
    }
  }
  return(ROCs)
}

getTrainIndex <- function(random.seed, source){
  set.seed(random.seed)
  trainIndex <- sample(nrow(source), floor(nrow(source)*0.7))
  return(trainIndex)
}

trainIndex.XENA.CRC <- getTrainIndex(random.seeds[1], expData.XENA.CRC)
library(MASS)
ldaModel.XENA.CRC <-
  lda(expData.XENA.CRC[trainIndex.XENA.CRC,],
      grouping=sampleInfo.XENA.CRC$sampleType[trainIndex.XENA.CRC])
save(ldaModel.XENA.CRC, file = 'ldaModel.XENA.CRC.RData')
LDAplot.XENA.CRC <-
  plotLDAmodel(ldaModel = ldaModel.XENA.CRC, title = 'XENA.CRC',
               trainSet = expData.XENA.CRC[trainIndex.XENA.CRC,],
               trainGroup = sampleInfo.XENA.CRC[trainIndex.XENA.CRC,],
               legend.position = c(0.05, 0.05))
ROCplot.XENA.CRC <-
  plotROC(LDAmodel=ldaModel.XENA.CRC, title='XENA.CRC',
          testSet=expData.XENA.CRC[-trainIndex.XENA.CRC,],
          testGroup=sampleInfo.XENA.CRC[-trainIndex.XENA.CRC,])


# ---------- Part IV: Recursive feature elimination (iRFE) based on XENA data ----------
# Recursive features elimination and get optional genes.
multiRfe <- function(dataSet, group, sizeRatios, clusterNo=3){
  library(caret)
  library(doSNOW)
  cl <- makeCluster(clusterNo, type='SOCK')
  # Memory space limits the number of threads.
  registerDoSNOW(cl)
  ldaModels <- list()
  count <- 1
  while (length(ldaModels) < 100 & count < 200) {
    skip2next=F
    ldaProfile <- tryCatch(
      {rfe(dataSet, group,
           # Grouping parameter must be factor, or R will throw 
           # "Error during wrapup: subscript out of bounds"
           sizes = ceiling(ncol(dataSet)*sizeRatios),
           rfeControl = rfeControl(functions = ldaFuncs,
                                   method = "cv", rerank=T))},
      error = function(e) {skip2next <- T})
    count <- count + 1
    if (isTRUE(ldaProfile)) {
      next
    } else {
      ldaModels[[length(ldaModels)+1]] <- ldaProfile
      print(paste('Model', length(ldaModels), 'completed!'))
    }
  }
  stopCluster(cl)
  return(ldaModels)
}

getConsistantGenes <- function(rfeModels, threshold){
  optEnsemblIDs <- data.frame(table(Reduce('c', sapply(rfeModels, `[[`, 6))))
  names(optEnsemblIDs) <- c('id', 'Freq')
  consistantEnsemblIDs <- optEnsemblIDs[optEnsemblIDs$Freq>threshold, 1]
  return(consistantEnsemblIDs)
}

# RFE of XENA data
ldaModels.RFE.XENA.CRC <-
  multiRfe(as.data.frame(expData.XENA.CRC),
           sampleInfo.XENA.CRC$sampleType,
           c(0.5^(1:6), 0.5^6*(0.7^(1:6)),
             0.5^6*0.7^6*0.9^(1:6)))
save(ldaModels.RFE.XENA.CRC, file='ldaModels.RFE.XENA.CRC.RData')

geneNo <- c()
for (i in 1:100) {
  geneNo <- c(geneNo,length(getConsistantGenes(rfeModels = ldaModels.RFE.XENA.CRC, i)))
}
plot(geneNo)

ensemblID.RFE.XENA.CRC <-
  getConsistantGenes(rfeModels = ldaModels.RFE.XENA.CRC, 65)
ensemblID.RFE.XENA.CRC <- str_sub(ensemblID.RFE.XENA.CRC, 1, 15)
expData.RFE.XENA.CRC <-
  expData.XENA.CRC[, colnames(expData.XENA.CRC) %in%
                     ensemblID.RFE.XENA.CRC]
trainIndex.XENA.CRC <- getTrainIndex(random.seeds[2],
                                     expData.XENA.CRC)
ldaModel.RFE.XENA.CRC <-
  lda(expData.RFE.XENA.CRC[trainIndex.XENA.CRC,],
      grouping=sampleInfo.XENA.CRC$sampleType[trainIndex.XENA.CRC])
LDAplot.RFE.XENA.CRC <-
  plotLDAmodel(
    ldaModel = ldaModel.RFE.XENA.CRC, title ='RFE.XENA.CRC',
    trainSet = expData.RFE.XENA.CRC[trainIndex.XENA.CRC,],
    trainGroup =  sampleInfo.XENA.CRC[trainIndex.XENA.CRC,],
    legend.position = c(0.05, 0.05))
ROCplot.RFE.XENA.CRC <-
  plotROC(LDAmodel=ldaModel.RFE.XENA.CRC, title='RFE.XENA.CRC',
          testSet=expData.RFE.XENA.CRC[-trainIndex.XENA.CRC,],
          testGroup=sampleInfo.XENA.CRC[-trainIndex.XENA.CRC,])


# ---------- Part X: iRFE analysis ----------
ensemblID.iRFE.XENA.CRC <- ensemblID.RFE.XENA.CRC
expData.iRFE.XENA.CRC <- expData.RFE.XENA.CRC
ensemblID.iRFE.No1 <- 0
ensemblID.iRFE.No2 <- length(ensemblID.iRFE.XENA.CRC)
iRFEModels.XENA.CRC <- list(); length(iRFEModels.XENA.CRC) <- 100
while ((ensemblID.iRFE.No1 != ensemblID.iRFE.No2) &
       (length(iRFEModels.XENA.CRC)==100) &
       (length(ensemblID.iRFE.XENA.CRC) != 0)) {
  iRFEModels.XENA.CRC <-
    multiRfe(dataSet = as.data.frame(expData.iRFE.XENA.CRC),
             group = sampleInfo.XENA.CRC$sampleType,
             sizeRatios = c(0.8^(1:7), 0.8^7*0.9^(1:7)),
             clusterNo = 11)
  
  ensemblID.iRFE.No1 <- ensemblID.iRFE.No2
  tmp <- c()
  for (i in 1:100) {
    tmp <- c(tmp, length(getConsistantGenes(iRFEModels.XENA.CRC, i)))
  }
  ensemblID.iRFE.XENA.CRC <-
    getConsistantGenes(iRFEModels.XENA.CRC,
                       sum(tmp==max(tmp))+1)
  ensemblID.iRFE.No2 <- length(ensemblID.iRFE.XENA.CRC)
  expData.iRFE.XENA.CRC <-
    subset(expData.iRFE.XENA.CRC,
           select=colnames(expData.iRFE.XENA.CRC) %in%
             ensemblID.iRFE.XENA.CRC)
}
save.image()

ensemblID.iRFE.XENA.CRC <- levels(ensemblID.iRFE.XENA.CRC)
expData.iRFE.XENA.CRC <-
  subset(expData.RFE.XENA.CRC,
         select=colnames(expData.RFE.XENA.CRC) %in%
           ensemblID.iRFE.XENA.CRC)

trainIndex.XENA.CRC <-
  getTrainIndex(random.seeds[3], expData.XENA.CRC)
ldaModel.iRFE.XENA.CRC <-
  lda(expData.iRFE.XENA.CRC[trainIndex.XENA.CRC,],
      grouping=sampleInfo.XENA.CRC$sampleType[trainIndex.XENA.CRC])
LDAplot.iRFE.XENA.CRC <-
  plotLDAmodel(ldaModel = ldaModel.iRFE.XENA.CRC,
               trainSet = expData.iRFE.XENA.CRC[trainIndex.XENA.CRC,],
               trainGroup = sampleInfo.XENA.CRC[trainIndex.XENA.CRC,],
               title = "iRFE genes, XENA CRC, train set",
               legend.position = c(0.05, 0.05))

ROCplot.iRFE.XENA.CRC <-
  plotROC(LDAmodel = ldaModel.iRFE.XENA.CRC,
          testSet = expData.iRFE.XENA.CRC[-trainIndex.XENA.CRC,],
          testGroup = sampleInfo.XENA.CRC[-trainIndex.XENA.CRC,],
          title = 'iRFE genes, XENA CRC, test set')


# ---------- Part XII: Craw gene summary from genecards ----------
getGeneSummary <- function(ensemblID, geneID){
  library(rvest)
  library(httr)
  url <- paste('https://www.genecards.org/cgi-bin/carddisp.pl?gene=',
               geneID, '&keywords=', ensemblID, sep='')
  webPage <- url %>%
    session(add_headers(`User-Agent`='Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/105.0.0.0 Safari/537.36 Edg/105.0.1343.50'))
  if ('try-error' %in% class(webPage)) {
    geneSummary <- 'Can not open url.'
    return(geneSummary)
  }
  location <- html_nodes(webPage, css='section#summaries')
  geneSummary <- try(paste(html_text(location), collapse=''), silent = T)
  geneSummary <- gsub(' +', ' ', geneSummary)
  geneSummary <- gsub('(\r\n( )*)+', '\n', geneSummary)
  geneSummary <- gsub('^\n', '', geneSummary)
  geneSummary <- gsub('Jump to sectionAliasesDisordersDomains & FamiliesDrugs & CompoundsExpressionFunctionGenomicsLocalizationOrthologsParalogsPathways & InteractionsProductsProteinsPublicationsSourcesTranscriptsVariants\n',
                      '', geneSummary)
  if('try-error' %in% class(geneSummary)){
    geneSummary <- 'try-error'
  }
  return(geneSummary)
}

crawGeneSummaries <- function(ensemblID_geneID_df){
  geneSummaries <- c()
  for (i in 1:nrow(ensemblID_geneID_df)) {
    ensemblID <- ensemblID_geneID_df[i, 1]
    geneID <- ensemblID_geneID_df[i, 2]
    geneSummary <- getGeneSummary(ensemblID, geneID)
    geneSummaries <- c(geneSummaries, geneSummary)
    Sys.sleep(1)
  }
  ensemblID_geneID_df[['geneSummary']] <- geneSummaries
  return(ensemblID_geneID_df)
}

library(clusterProfiler)
library(org.Hs.eg.db)
gene.enID.symbol.XENA.CRC <-
  bitr(ensemblID.iRFE.XENA.CRC, fromType = 'ENSEMBL',
       toType = 'SYMBOL',OrgDb = org.Hs.eg.db)
geneSummary.iRFE.XENA.CRC <- crawGeneSummaries(gene.enID.symbol.XENA.CRC)

library(openxlsx)
write.xlsx(geneSummary.iRFE.XENA.CRC, 'geneSummary.iRFE.XENA.CRC.xlsx')


# ---------- Part V: Enrichment analysis based on the optional genes from RFE ----------
# GO and KEGG enrichment analysis.
GO_KEGG_enrichment <- function(geneset, pvalueCutoff, qvalueCutoff){
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  ego_cc <- enrichGO(gene=geneset,
                     keyType='ENSEMBL',
                     OrgDb=org.Hs.eg.db,
                     ont='CC',
                     pAdjustMethod='BH',
                     pvalueCutoff=pvalueCutoff,
                     qvalueCutoff=qvalueCutoff,
                     readable=T)
  ego_bp <- enrichGO(gene=geneset,
                     keyType='ENSEMBL',
                     OrgDb=org.Hs.eg.db,
                     ont='BP',
                     pAdjustMethod='BH',
                     pvalueCutoff=pvalueCutoff,
                     qvalueCutoff=qvalueCutoff,
                     readable=T)
  ego_mf <- enrichGO(gene=geneset,
                     keyType='ENSEMBL',
                     OrgDb=org.Hs.eg.db,
                     ont='MF',
                     pAdjustMethod='BH',
                     pvalueCutoff=pvalueCutoff,
                     qvalueCutoff=qvalueCutoff,
                     readable=T)
  entrezID <- bitr(geneset,
                   fromType='ENSEMBL',
                   toType='ENTREZID',
                   OrgDb=org.Hs.eg.db)
  kk <- enrichKEGG(gene=entrezID$ENTREZID,
                   organism='human',
                   pvalueCutoff=pvalueCutoff,
                   qvalueCutoff=qvalueCutoff)
  results <- list()
  results$CC <- ego_cc
  results$BP <- ego_bp
  results$MF <- ego_mf
  results$KK <- kk
  return(results)
}

plotEnrichResults <- function(enrichmentResults, txtwidth=27, showNum=20,
                              pcutoff=0.05, qcutoff=0.05){
  library(Rmisc)
  library(ggplot2)
  library(stringr)
  if (!(is.null(enrichmentResults$CC)) &&
      ((min(enrichmentResults$CC@result$p.adjust)<=pcutoff &&
        min(enrichmentResults$CC@result$qvalue)<=qcutoff) ||
       (min(enrichmentResults$CC@result$p.adjust)<=pcutoff &&
        sum(is.na(enrichmentResults$CC@result$qvalue))>0))) {
    CCenrich <- clusterProfiler::dotplot(enrichmentResults$CC,
                                         title='CCenrichment',
                                         showCategory=showNum) +
      scale_y_discrete(labels=function(x) str_wrap(x, width=txtwidth))
  } else {
    CCenrich <- NULL
  }
  
  if (!(is.null(enrichmentResults$BP)) &&
      ((min(enrichmentResults$BP@result$p.adjust)<=pcutoff &&
        min(enrichmentResults$BP@result$qvalue)<=qcutoff) ||
       (min(enrichmentResults$BP@result$p.adjust)<=pcutoff &&
        sum(is.na(enrichmentResults$BP@result$qvalue))>0))) {
    BPenrich <- clusterProfiler::dotplot(enrichmentResults$BP,
                                         title='BPenrichment',
                                         showCategory=showNum) +
      scale_y_discrete(labels=function(x) str_wrap(x, width=txtwidth))
  } else {
    BPenrich <- NULL
  }
  
  if (!(is.null(enrichmentResults$MF)) &&
      ((min(enrichmentResults$MF@result$p.adjust)<=pcutoff &&
        min(enrichmentResults$MF@result$qvalue)<=qcutoff) ||
       (min(enrichmentResults$MF@result$p.adjust)<=pcutoff &&
        sum(is.na(enrichmentResults$MF@result$qvalue))>0))) {
    MFenrich <- clusterProfiler::dotplot(enrichmentResults$MF,
                                         title='MFenrichment',
                                         showCategory=showNum) +
      scale_y_discrete(labels=function(x) str_wrap(x, width=txtwidth))
  } else {
    MFenrich <- NULL
  }
  
  if (!(is.null(enrichmentResults$KK)) &&
      ((min(enrichmentResults$KK@result$p.adjust)<=pcutoff &&
        min(enrichmentResults$KK@result$qvalue)<=qcutoff) ||
       (min(enrichmentResults$KK@result$p.adjust)<=pcutoff &&
        sum(is.na(enrichmentResults$KK@result$qvalue))>0))) {
    KKenrich <- clusterProfiler::dotplot(enrichmentResults$KK,
                                         title='KKenrichment',
                                         showCategory=showNum) +
      scale_y_discrete(labels=function(x) str_wrap(x, width=txtwidth))
  } else {
    KKenrich <- NULL
  }
  return(list(CCenrich, BPenrich, MFenrich, KKenrich))
}

# GO and KEGG enrichment analysis of the genes from XENA data set.
enrich.RFE.XENA.CRC <- GO_KEGG_enrichment(ensemblID.RFE.XENA.CRC, 0.05, 0.05)
enrichplot.RFE.XENA.CRC <- plotEnrichResults(enrich.RFE.XENA.CRC, 27, 20)

enrich.iRFE.XENA.CRC <- GO_KEGG_enrichment(ensemblID.iRFE.XENA.CRC, 0.05, 0.05)
enrichplot.iRFE.XENA.CRC <- plotEnrichResults(enrich.iRFE.XENA.CRC, 27, 20)


# ---------- Part VI: Analysis of the genes from ribosome and proteasome ----------
library(clusterProfiler)
geneSets <-
  read.gmt('E:/PublicData/msigdb_v2022.1.Hs_GMTs/msigdb.v2022.1.Hs.entrez.gmt')
getEnsemblIDsOfGeneSet <- function(geneSetName){
  entrezIDs <- geneSets[geneSets[,1]==geneSetName,2]
  library(clusterProfiler)
  ensemblIDs <- bitr(entrezIDs, fromType = 'ENTREZID',
                     toType = 'ENSEMBL', OrgDb = 'org.Hs.eg.db')
  return(ensemblIDs$ENSEMBL)
}

getGeneSetCor <- function(geneSetNames, expData,
                          sampleInfo, geneSetLabels=NA){
  expData.selected.df <- data.frame(matrix(nrow = nrow(expData)))
  expData.selected.df[,1] <- NULL
  expData.selected.list <- list()
  for (i in 1:length(geneSetNames)) {
    assign(paste0('ensemblIDs.set', i),
           getEnsemblIDsOfGeneSet(geneSetNames[[i]]))
    if (i>1) {
      for (j in 1:(i-1)) {
        assign(paste0('ensemblIDs.set', i),
               setdiff(get(paste0('ensemblIDs.set', i)),
                       get(paste0('ensemblIDs.set', j))))
      }
    }
    assign(paste0('expData', i),
           expData[, colnames(expData) %in%
                     get(paste0('ensemblIDs.set', i)), drop=F])
    assign(paste0('expData', i),
           get(paste0('expData', i))[,apply(get(paste0('expData', i)),
                                            2,sd)!=0, drop=F])
    expData.selected.df <- cbind(expData.selected.df,
                                 get(paste0('expData', i)))
    expData.selected.list[[length(expData.selected.list)+1]] <-
      get(paste0('expData', i))
  }
  ifelse(all(is.na(geneSetLabels)),
         names(expData.selected.list) <- geneSetNames,
         names(expData.selected.list) <- geneSetLabels)
  CC.geneSets.matrix <- cor(expData.selected.df)

  CC.geneSets.vec <- list()
  for (i in 1:(length(expData.selected.list)-1)) {
    for (j in names(expData.selected.list)[(i+1):length(expData.selected.list)]) {
      assign(paste0(names(expData.selected.list[i]), ' VS ', j),
                    c(cor(expData.selected.list[[names(expData.selected.list)[i]]],
                          expData.selected.list[[get('j')]])))
      CC.geneSets.vec[[length(CC.geneSets.vec)+1]] <-
        get(paste0(names(expData.selected.list)[i], ' VS ', j))
      names(CC.geneSets.vec)[length(CC.geneSets.vec)] <-
        paste0(names(expData.selected.list)[i], ' VS ', j)}
  }
  if (length(CC.geneSets.vec)==1) {
    CC.geneSets.vec <- CC.geneSets.vec[[1]]
  }
  CC.geneSets.vec <- CC.geneSets.vec[!is.na(CC.geneSets.vec)]

  library(ComplexHeatmap)
  library(circlize)
  if (all(is.na(geneSetLabels))) {
    en.group <- factor(rep(geneSetNames,
                           times=unname(sapply(expData.selected.list, ncol))),
                       levels = geneSetNames)
  } else {
    en.group <- factor(rep(geneSetLabels,
                    times=unname(sapply(expData.selected.list, ncol))),
                    levels = geneSetLabels)
  }
  sampleType <- sampleInfo$sampleType
  col_fun.cor = colorRamp2(c(min(CC.geneSets.matrix, na.rm = T), 0,
                             max(CC.geneSets.matrix, na.rm = T)),
                           c('blue', "#ffffff", "#ff0000"))
  if (all(is.na(geneSetLabels))) {
    hm.labels.en <-
      HeatmapAnnotation(cluster=anno_block(labels = geneSetNames,
                                           gp = gpar(fill = c("#2fa1dd", "#f87669")),
                                           labels_gp = gpar(col = "white", fontsize = 10)))
    hm.labels.en.row <-
      rowAnnotation(cluster=anno_block(labels = geneSetNames,
                                       gp = gpar(fill = c("#2fa1dd", "#f87669")),
                                       labels_gp = gpar(col = "white", fontsize = 10)))
  } else {
    hm.labels.en <-
      HeatmapAnnotation(cluster=anno_block(labels = geneSetLabels,
                                           gp = gpar(fill = c("#2fa1dd", "#f87669")),
                                           labels_gp = gpar(col = "white", fontsize = 10)))
    hm.labels.en.row <-
      rowAnnotation(cluster=anno_block(labels = geneSetLabels,
                                       gp = gpar(fill = c("#2fa1dd", "#f87669")),
                                       labels_gp = gpar(col = "white", fontsize = 10)))
  }
  
  hm.labels.sampleType <-
    rowAnnotation(cluster=anno_block(
      labels = sapply(str_split(unique(sampleInfo$sampleType),
                                ' '),paste,collapse='\n'),
      gp = gpar(fill = c("#2fa1dd", "#f87669", '#40e0d0')),
      labels_gp = gpar(col = "white", fontsize = 10)))
  hm.CC <- Heatmap(CC.geneSets.matrix, name = " ",
                   top_annotation = hm.labels.en,
                   left_annotation = hm.labels.en.row,
                   column_split = en.group, row_split = en.group,
                   show_heatmap_legend = T, border = F,
                   show_column_names = F, show_row_names = F,
                   column_title = NULL, row_title = NULL, col=col_fun.cor,
                   cluster_row_slices = F, cluster_column_slices = F)
  hm.exp <- Heatmap(expData.selected.df, name = " ",
                    top_annotation = hm.labels.en,
                    left_annotation = hm.labels.sampleType,
                    column_split = en.group, row_split = sampleType,
                    show_heatmap_legend = T, border = F,
                    show_column_names = F, show_row_names = F,
                    column_title = NULL, row_title = NULL,
                    cluster_row_slices = F, cluster_column_slices = F)
  return(list(CC=CC.geneSets.vec, hm.cor=hm.CC, hm.exp=hm.exp))
}

cor.ribo.prot.XENA.CRC <-
  getGeneSetCor(geneSetNames = c('GOCC_RIBOSOME', 'GOCC_PROTEASOME_COMPLEX'),
                expData = expData.XENA.CRC, sampleInfo = sampleInfo.XENA.CRC,
                geneSetLabels = c('Ribosome', 'Proteasome'))
CC.ribo.prot.XENA.CRC <- cor.ribo.prot.XENA.CRC[[1]]
hm.CC.ribo.prot.XENA.CRC <- cor.ribo.prot.XENA.CRC[[2]]
hm.exp.ribo.prot.XENA.CRC <- cor.ribo.prot.XENA.CRC[[3]]
rm(cor.ribo.prot.XENA.CRC); gc()

# Correlation analysis between randomly selected gene sets and specified gene set.
geneSetNames <- unique(geneSets[,1])
geneSetNames <-
  str_subset(geneSetNames, '^GOBP_|^GOCC_|GOMF_|KEGG_|HALLMARK_')
get.largest.random.CC.df <-
  function(random.seeds, expdata, sampleInfo,
           largestNo=10, specified.geneSet.name,
           specified.label, specified.df,
           group.levels){
  CC.random.geneSets <- list()
  medians.CC <- c()
  random.geneSetNames <- specified.geneSet.name
  j <- max(random.seeds)
  for (i in random.seeds) {
    set.seed(i)
    randomlySelectedGeneSets <-
      as.character(geneSetNames[sample(1:length(geneSetNames), 1)])
    while (randomlySelectedGeneSets %in% random.geneSetNames) {
      j <- j + 1
      set.seed(j)
      randomlySelectedGeneSets <-
        as.character(geneSetNames[sample(1:length(geneSetNames), 1)])
    }
    CC.random.geneSet <-
      getGeneSetCor(geneSetNames = c(specified.geneSet.name,
                                     randomlySelectedGeneSets),
                    expData = expdata, sampleInfo = sampleInfo)
    CC.random.geneSets[[length(CC.random.geneSets)+1]] <- CC.random.geneSet$CC
    names(CC.random.geneSets)[length(CC.random.geneSets)] <- randomlySelectedGeneSets
    medians.CC <- c(medians.CC, median(CC.random.geneSet$CC))
    random.geneSetNames <- c(random.geneSetNames, randomlySelectedGeneSets)
  }
  CC.random.geneSets <-
    CC.random.geneSets[order(medians.CC, decreasing = T)[1:largestNo]]
  
  random.CC.group <- c()
  for (i in 1:length(CC.random.geneSets)) {
    random.CC.group <- c(random.CC.group,
                         rep(paste0(specified.label,' VS ',
                                    names(CC.random.geneSets)[i]),
                             length(CC.random.geneSets[[i]])))
  }
  
  random.CC.df <- data.frame(comPairs=random.CC.group,
                             CC=unlist(CC.random.geneSets),
                             Group=paste0(specified.label, ' VS Rand'))
  random.CC.df <- rbind(specified.df, random.CC.df)
  random.CC.df$comPairs <- factor(random.CC.df$comPairs,
                                  levels = unique(random.CC.df$comPairs))
  random.CC.df$Group <- factor(random.CC.df$Group, levels = group.levels)
  return(random.CC.df)
}

CC.ribo.prot.df.XENA.CRC <-
  data.frame(comPairs='Ribo VS Prot', CC=CC.ribo.prot.XENA.CRC,
             Group='Ribo VS Prot')
random.CC.ribo.df.XENA.CRC <-
  get.largest.random.CC.df(random.seeds = random.seeds[4:103],
                           specified.geneSet.name = 'GOCC_RIBOSOME',
                           specified.label = 'Ribo',
                           expdata = expData.XENA.CRC,
                           sampleInfo = sampleInfo.XENA.CRC,
                           specified.df = CC.ribo.prot.df.XENA.CRC,
                           group.levels = c('Ribo VS Prot', 'Ribo VS Rand'))

random.CC.prot.df.XENA.CRC <-
  get.largest.random.CC.df(random.seeds = random.seeds[104:203],
                           specified.geneSet.name = 'GOCC_PROTEASOME_COMPLEX',
                           specified.label = 'Prot',
                           expdata = expData.XENA.CRC,
                           sampleInfo = sampleInfo.XENA.CRC,
                           specified.df = CC.ribo.prot.df.XENA.CRC,
                           group.levels = c('Ribo VS Prot', 'Prot VS Rand'))

getComPairs <- function(geneSets, pair1=NULL, pair2=NULL, comPairs=NULL){
  comPairs <- list()
  if (is.null(pair1) & is.null(pair2) & is.null(pairs)) {
    for (i in 1:(length(geneSets)-1)) {
      for (j in (i+1):length(geneSets)) {
        comPairs[[length(comPairs)+1]] <-
          c(as.character(geneSets[i]),
            as.character(geneSets[j]))
      }
    }
    return(comPairs)
  } else if (is.null(pairs)) {
    for (i in pair1) {
      for (j in pair2) {
        comPairs[[length(comPairs)+1]] <-
          c(as.character(geneSets[i]),
            as.character(geneSets[j]))
      }
    }
    return(comPairs)
  } else if (!is.null(pairs)) {
    for (i in pairs) {
      comPairs[[length(comPairs)+1]] <- c(geneSets[i[1]], geneSets[i[2]])
    }
    return(comPairs)
  }
}

plotBox <- function(CC.df, show.legend=T, xlabel=NULL,
                    ylabel='Correlation coefficient',
                    title='XENA.CRC', comPair1=NULL,
                    comPair2=NULL, as.title=T, comPairIndex=NULL){
  library(ggplot2)
  library(ggpubr)
  library(stringr)
  colnames(CC.df) <- c('comPairs', 'CC', 'Group')
  if (is.null(comPairIndex)) {
    test.comPairs <-
      getComPairs(as.character(unique(CC.df$comPairs)),
                  comPair1, comPair2)
  } else {
    test.comPairs <-
      getComPairs(as.character(unique(CC.df$comPairs)),
                  pairs = comPairIndex)
  }
  if (as.title) {
    xlabels <- str_wrap(str_to_title(str_replace_all(unique(CC.df[,1]), '_', ' ')), width=12)
  } else {
    xlabels <- str_wrap(str_replace_all(unique(CC.df[,1]), '_', ' '), width=12)
  }
  boxplot.CC <-
    ggplot(data=CC.df, aes(y=CC, x=comPairs, fill = Group)) +
    geom_boxplot(notch = FALSE, outlier.alpha = 1) +
    #    scale_fill_brewer(palette = "Set2") +
    geom_signif(comparisons = test.comPairs, step_increase = 0.1,
                map_signif_level = F, test = t.test) +
    theme(axis.title = element_text(size=10, face = 'plain'),
          axis.text = element_text(size=8, face = 'plain'),
          plot.title = element_text(hjust = 0)) +
    labs(x=xlabel, y=str_wrap(str_replace_all(ylabel, '_', ' '),25),
         title=title) +
    geom_hline(
      yintercept=median(CC.df[,2][CC.df[,1]==unique(CC.df[,1])[1]]),
      colour='red', linetype='dashed') +
    scale_x_discrete(labels=xlabels)
  if (show.legend==T) {
    return(boxplot.CC)
  } else {
    return(boxplot.CC + guides(fill='none'))
  }
}

boxplot.CC.ribo.random.XENA.CRC <-
  plotBox(random.CC.ribo.df.XENA.CRC, comPair1 = 1, comPair2 = 2:11)
boxplot.CC.prot.random.XENA.CRC <-
  plotBox(random.CC.prot.df.XENA.CRC, comPair1 = 1, comPair2 = 2:11)

# The difference between protein synthesis and degradation
get_PCA_scores <- function(expData, geneSetName, sampleType, title){
  library(factoextra)
  ensemblIDs.set <- getEnsemblIDsOfGeneSet(geneSetName)
  expData <- expData[, colnames(expData) %in% ensemblIDs.set]
  PCAres <- prcomp(as.matrix(expData), scale. = T)
  PCAplot <- fviz_pca_ind(PCAres, col.ind=sampleType, mean.point=F,
                          addEllipses = T, repel=T, geom='point',
                          legend.title="Sample Type", title=title)
  PCAscores <- PCAres$x[,1]
  if (abs(sum(apply(expData, 2, cor, PCAscores)[apply(expData, 2, cor, PCAscores)<0])) >
      sum(apply(expData,2,cor,PCAscores)[apply(expData,2,cor,PCAscores)>0])) {
    PCAscores <- -PCAscores
  }
  return(list(PCAscores=PCAscores, PCAres=PCAres, PCAplot=PCAplot))
}

ensemblID.ribo.prot <-
  c(getEnsemblIDsOfGeneSet('GOCC_RIBOSOME'),
    getEnsemblIDsOfGeneSet('GOCC_PROTEASOME_COMPLEX'))
expData.ribo.prot.XENA.CRC <-
  expData.XENA.CRC[, colnames(expData.XENA.CRC) %in%
                     ensemblID.ribo.prot]
PCAscores.ribo.XENA.CRC <-
  get_PCA_scores(
    expData = expData.ribo.prot.XENA.CRC,
    sampleType = sampleInfo.XENA.CRC$sampleType,
    geneSetName = 'GOCC_RIBOSOME',
    title = 'Ribo.XENA.CRC')
PCAscores.prot.XENA.CRC <-
  get_PCA_scores(
    expData = expData.ribo.prot.XENA.CRC,
    sampleType = sampleInfo.XENA.CRC$sampleType,
    geneSetName = 'GOCC_PROTEASOME_COMPLEX',
    title = 'Prot.XENA.CRC')

PCAscores.ribo.prot.XENA.CRC <-
  data.frame(Ribosome=PCAscores.ribo.XENA.CRC$PCAscores,
             Proteasome=PCAscores.prot.XENA.CRC$PCAscores,
             sampleType=sampleInfo.XENA.CRC$sampleType)

plotCor.Classes <- function(twoDimData, title){
  library(ggpmisc)
  library(ggpubr)
  corPlot.Classes <-
    ggplot(data=twoDimData,
           aes_string(x=colnames(twoDimData)[1],
                      y=colnames(twoDimData)[2],
                      color=colnames(twoDimData)[3])) +
    geom_smooth(method = 'lm', se=F, formula = y~x, size=0.6) +
    stat_poly_eq(formula = y~x,
                 aes(label = paste(..eq.label.., ..rr.label..,
                                   ..p.value.label.., ..n.label..,
                                   sep='~~~')), 
                 parse = T, label.x.npc = 0.04,
                 label.y.npc = seq(0.98,0,-0.06)[1:length(levels(twoDimData[,3]))]) +
    geom_point(size=0.8, alpha=0.5) + geom_rug() +
    labs(title=title)
  return(corPlot.Classes)
}

corPlot.ribo.prot.XENA.CRC <-
  plotCor.Classes(PCAscores.ribo.prot.XENA.CRC, title = 'XENA.CRC')

PCAscores.diff.ribo.prot.XENA.CRC <-
  data.frame(
    sampleType=PCAscores.ribo.prot.XENA.CRC$sampleType,
    diff=PCAscores.ribo.XENA.CRC$PCAscores -
      PCAscores.prot.XENA.CRC$PCAscores, Group=3)
groups <- levels(PCAscores.diff.ribo.prot.XENA.CRC$sampleType)
boxplot.diff.ribo.prot.XENA.CRC <-
  plotBox(CC.df = PCAscores.diff.ribo.prot.XENA.CRC,
          comPairIndex = list(c(1,2), c(1,3), c(2,3)), show.legend = F,
          ylabel = 'PCAscores.Ribosome - PCAscores.Proteasome')


# ---------- Part IX: Analysis of all the gene sets about protein synthesis and degradation ----------
Pr.Fold <-
  c('GOBP_PROTEIN_FOLDING',
    'GOBP_DE_NOVO_PROTEIN_FOLDING',
    'GOBP_PROTEIN_QUALITY_CONTROL_FOR_MISFOLDED_OR_INCOMPLETELY_SYNTHESIZED_PROTEINS',
    'GOBP_PROTEIN_MATURATION_BY_PROTEIN_FOLDING',
    'GOBP_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_CELLULAR_RESPONSE_TO_UNFOLDED_PROTEIN',
    'GOBP_PROTEIN_FOLDING_IN_ENDOPLASMIC_RETICULUM',
    'GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_ATF6_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_PROTEIN_REFOLDING',
    'GOBP_CHAPERONE_COFACTOR_DEPENDENT_PROTEIN_REFOLDING',
    'GOBP_RESPONSE_TO_MISFOLDED_PROTEIN',
    'GOBP_CHAPERONE_MEDIATED_PROTEIN_FOLDING',
    'GOBP_ER_ASSOCIATED_MISFOLDED_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_REGULATION_OF_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_NEGATIVE_REGULATION_OF_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_POSITIVE_REGULATION_OF_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_REGULATION_OF_PROTEIN_FOLDING',
    'GOBP_REGULATION_OF_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_NEGATIVE_REGULATION_OF_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_POSITIVE_REGULATION_OF_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_REGULATION_OF_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_NEGATIVE_REGULATION_OF_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    'GOMF_PROTEIN_FOLDING_CHAPERONE',
    'GOMF_UNFOLDED_PROTEIN_BINDING',
    'GOMF_MISFOLDED_PROTEIN_BINDING',
    'GOMF_ATP_DEPENDENT_PROTEIN_FOLDING_CHAPERONE')

Pr.Cata <-
  c('GOBP_GLYCOPROTEIN_CATABOLIC_PROCESS',
    'GOBP_PROTEIN_CATABOLIC_PROCESS_IN_THE_VACUOLE',
    'GOBP_PROTEASOMAL_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_UBIQUITIN_DEPENDENT_SMAD_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_SCF_DEPENDENT_PROTEASOMAL_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_REGULATION_OF_PROTEASOMAL_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_NEGATIVE_REGULATION_OF_PROTEASOMAL_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_POSITIVE_REGULATION_OF_PROTEASOMAL_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_MITOCHONDRIAL_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_LIPOPROTEIN_CATABOLIC_PROCESS',
    'GOBP_REGULATION_OF_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_NEGATIVE_REGULATION_OF_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_AMYLOID_PRECURSOR_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS_VIA_THE_MULTIVESICULAR_BODY_SORTING_PATHWAY',
    'GOBP_PROTEIN_TRANSPORT_TO_VACUOLE_INVOLVED_IN_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS_VIA_THE_MULTIVESICULAR_BODY_SORTING_PATHWAY',
    'GOBP_POSITIVE_REGULATION_OF_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_PROTEOLYSIS_INVOLVED_IN_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_ER_ASSOCIATED_MISFOLDED_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_PROTEIN_DEUBIQUITINATION_INVOLVED_IN_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_UBIQUITIN_INDEPENDENT_PROTEIN_CATABOLIC_PROCESS_VIA_THE_MULTIVESICULAR_BODY_SORTING_PATHWAY',
    'GOBP_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS_VIA_THE_C_END_DEGRON_RULE_PATHWAY',
    'GOBP_NEGATIVE_REGULATION_OF_PROTEASOMAL_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_REGULATION_OF_ASPARTIC_TYPE_ENDOPEPTIDASE_ACTIVITY_INVOLVED_IN_AMYLOID_PRECURSOR_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_REGULATION_OF_AMYLOID_PRECURSOR_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_NEGATIVE_REGULATION_OF_AMYLOID_PRECURSOR_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_POSITIVE_REGULATION_OF_AMYLOID_PRECURSOR_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_REGULATION_OF_PROTEOLYSIS_INVOLVED_IN_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_NEGATIVE_REGULATION_OF_PROTEOLYSIS_INVOLVED_IN_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_POSITIVE_REGULATION_OF_PROTEOLYSIS_INVOLVED_IN_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_REGULATION_OF_ER_ASSOCIATED_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_POSITIVE_REGULATION_OF_ER_ASSOCIATED_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_REGULATION_OF_PROTEIN_CATABOLIC_PROCESS_IN_THE_VACUOLE',
    'GOBP_POSITIVE_REGULATION_OF_PROTEIN_CATABOLIC_PROCESS_IN_THE_VACUOLE',
    'GOBP_LYSOSOMAL_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_REGULATION_OF_LYSOSOMAL_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_RIBOSOME_ASSOCIATED_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_REGULATION_OF_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_NEGATIVE_REGULATION_OF_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_POSITIVE_REGULATION_OF_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS',
    'GOMF_PROTEIN_KINASE_A_CATALYTIC_SUBUNIT_BINDING',
    'GOMF_CATALYTIC_ACTIVITY_ACTING_ON_A_GLYCOPROTEIN',
    'GOBP_AUTOPHAGY_OF_MITOCHONDRION',
    'GOBP_AUTOPHAGY_OF_NUCLEUS',
    'GOBP_AUTOPHAGY_OF_PEROXISOME',
    'GOBP_CHAPERONE_MEDIATED_AUTOPHAGY',
    'GOBP_LYSOSOMAL_MICROAUTOPHAGY',
    'GOBP_MACROAUTOPHAGY',
    'GOBP_NEGATIVE_REGULATION_OF_AUTOPHAGY',
    'GOBP_NEGATIVE_REGULATION_OF_AUTOPHAGY_OF_MITOCHONDRION',
    'GOBP_NEGATIVE_REGULATION_OF_MACROAUTOPHAGY',
    'GOBP_PIECEMEAL_MICROAUTOPHAGY_OF_THE_NUCLEUS',
    'GOBP_POSITIVE_REGULATION_OF_AUTOPHAGY',
    'GOBP_POSITIVE_REGULATION_OF_AUTOPHAGY_OF_MITOCHONDRION',
    'GOBP_POSITIVE_REGULATION_OF_MACROAUTOPHAGY',
    'GOBP_REGULATION_OF_AUTOPHAGY',
    'GOBP_REGULATION_OF_AUTOPHAGY_OF_MITOCHONDRION',
    'GOBP_REGULATION_OF_AUTOPHAGY_OF_MITOCHONDRION_IN_RESPONSE_TO_MITOCHONDRIAL_DEPOLARIZATION',
    'GOBP_REGULATION_OF_CHAPERONE_MEDIATED_AUTOPHAGY',
    'GOBP_REGULATION_OF_MACROAUTOPHAGY',
    'GOBP_SELECTIVE_AUTOPHAGY',
    'KEGG_REGULATION_OF_AUTOPHAGY',
    'GOBP_LYSOSOMAL_MICROAUTOPHAGY',
    'GOBP_LYSOSOMAL_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_REGULATION_OF_LYSOSOMAL_PROTEIN_CATABOLIC_PROCESS',
    'GOCC_AUTOLYSOSOME',
    'GOCC_PHAGOLYSOSOME',
    'GOMF_LYSOZYME_ACTIVITY',
    'KEGG_LYSOSOME'
  )

ER <-
  c('GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I_VIA_ER_PATHWAY',
    'GOBP_PROTEIN_RETENTION_IN_ER_LUMEN',
    'GOBP_ER_OVERLOAD_RESPONSE',
    'GOBP_ER_NUCLEUS_SIGNALING_PATHWAY',
    'GOBP_PREASSEMBLY_OF_GPI_ANCHOR_IN_ER_MEMBRANE',
    'GOBP_PROTEIN_INSERTION_INTO_ER_MEMBRANE',
    'GOBP_PROTEIN_INSERTION_INTO_ER_MEMBRANE_BY_STOP_TRANSFER_MEMBRANE_ANCHOR_SEQUENCE',
    'GOBP_VESICLE_TARGETING_ROUGH_ER_TO_CIS_GOLGI',
    'GOBP_REGULATION_OF_ER_TO_GOLGI_VESICLE_MEDIATED_TRANSPORT',
    'GOBP_ER_ASSOCIATED_MISFOLDED_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_TAIL_ANCHORED_MEMBRANE_PROTEIN_INSERTION_INTO_ER_MEMBRANE',
    'GOBP_REGULATION_OF_ER_ASSOCIATED_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_POSITIVE_REGULATION_OF_ER_ASSOCIATED_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS',
    'GOBP_REGULATION_OF_RETROGRADE_PROTEIN_TRANSPORT_ER_TO_CYTOSOL',
    'GOMF_ER_RETENTION_SEQUENCE_BINDING',
    'GOCC_ER_UBIQUITIN_LIGASE_COMPLEX',
    'GOCC_ER_TO_GOLGI_TRANSPORT_VESICLE_MEMBRANE',
    'GOCC_COPII_COATED_ER_TO_GOLGI_TRANSPORT_VESICLE',
    'GOCC_ER_MEMBRANE_INSERTION_COMPLEX',
    'GOBP_POST_TRANSLATIONAL_PROTEIN_TARGETING_TO_ENDOPLASMIC_RETICULUM_MEMBRANE',
    'GOBP_ENDOPLASMIC_RETICULUM_TO_GOLGI_VESICLE_MEDIATED_TRANSPORT',
    'GOBP_RETROGRADE_VESICLE_MEDIATED_TRANSPORT_GOLGI_TO_ENDOPLASMIC_RETICULUM',
    'GOBP_ENDOPLASMIC_RETICULUM_ORGANIZATION',
    'GOBP_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_ENDOPLASMIC_RETICULUM_CALCIUM_ION_HOMEOSTASIS',
    'GOBP_NEGATIVE_REGULATION_OF_ENDOPLASMIC_RETICULUM_CALCIUM_ION_CONCENTRATION',
    'GOBP_PROTEIN_EXIT_FROM_ENDOPLASMIC_RETICULUM',
    'GOBP_PROTEIN_FOLDING_IN_ENDOPLASMIC_RETICULUM',
    'GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS',
    'GOBP_MAINTENANCE_OF_PROTEIN_LOCALIZATION_IN_ENDOPLASMIC_RETICULUM',
    'GOBP_REGULATION_OF_TRANSLATION_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS',
    'GOBP_REGULATION_OF_TRANSLATION_INITIATION_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS',
    'GOBP_POSITIVE_REGULATION_OF_TRANSLATION_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS',
    'GOBP_ENDOPLASMIC_RETICULUM_LOCALIZATION',
    'GOBP_ENDOPLASMIC_RETICULUM_PLASMA_MEMBRANE_TETHERING',
    'GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS',
    'GOBP_REGULATION_OF_PROTEIN_EXIT_FROM_ENDOPLASMIC_RETICULUM',
    'GOBP_NEGATIVE_REGULATION_OF_PROTEIN_EXIT_FROM_ENDOPLASMIC_RETICULUM',
    'GOBP_POSITIVE_REGULATION_OF_PROTEIN_EXIT_FROM_ENDOPLASMIC_RETICULUM',
    'GOBP_PROTEIN_LOCALIZATION_TO_ENDOPLASMIC_RETICULUM',
    'GOBP_PROTEIN_LOCALIZATION_TO_ENDOPLASMIC_RETICULUM_EXIT_SITE',
    'GOBP_ENDOPLASMIC_RETICULUM_TUBULAR_NETWORK_ORGANIZATION',
    'GOBP_ESTABLISHMENT_OF_PROTEIN_LOCALIZATION_TO_ENDOPLASMIC_RETICULUM',
    'GOBP_ENDOPLASMIC_RETICULUM_MEMBRANE_ORGANIZATION',
    'GOBP_REGULATION_OF_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_NEGATIVE_REGULATION_OF_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_POSITIVE_REGULATION_OF_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_REGULATION_OF_ENDOPLASMIC_RETICULUM_STRESS_INDUCED_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY',
    'GOBP_NEGATIVE_REGULATION_OF_ENDOPLASMIC_RETICULUM_STRESS_INDUCED_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY',
    'GOBP_POSITIVE_REGULATION_OF_ENDOPLASMIC_RETICULUM_STRESS_INDUCED_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY',
    'GOBP_REGULATION_OF_ENDOPLASMIC_RETICULUM_TUBULAR_NETWORK_ORGANIZATION',
    'GOBP_ENDOPLASMIC_RETICULUM_TO_CYTOSOL_TRANSPORT',
    'GOBP_RELEASE_OF_SEQUESTERED_CALCIUM_ION_INTO_CYTOSOL_BY_ENDOPLASMIC_RETICULUM',
    'GOBP_NEGATIVE_REGULATION_OF_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS',
    'GOBP_NEGATIVE_REGULATION_OF_ENDOPLASMIC_RETICULUM_STRESS_INDUCED_EIF2_ALPHA_PHOSPHORYLATION',
    'GOBP_ENDOPLASMIC_RETICULUM_MANNOSE_TRIMMING',
    'GOBP_REGULATION_OF_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS',
    'GOBP_POSITIVE_REGULATION_OF_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS',
    'GOBP_POSITIVE_REGULATION_OF_TRANSCRIPTION_FROM_RNA_POLYMERASE_II_PROMOTER_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS',
    'GOBP_MITOCHONDRION_ENDOPLASMIC_RETICULUM_MEMBRANE_TETHERING',
    'GOCC_ENDOPLASMIC_RETICULUM_LUMEN',
    'GOCC_SMOOTH_ENDOPLASMIC_RETICULUM',
    'GOCC_ROUGH_ENDOPLASMIC_RETICULUM',
    'GOCC_ENDOPLASMIC_RETICULUM_GOLGI_INTERMEDIATE_COMPARTMENT',
    'GOCC_ROUGH_ENDOPLASMIC_RETICULUM_MEMBRANE',
    'GOCC_SMOOTH_ENDOPLASMIC_RETICULUM_MEMBRANE',
    'GOCC_ENDOPLASMIC_RETICULUM_PALMITOYLTRANSFERASE_COMPLEX',
    'GOCC_INTRINSIC_COMPONENT_OF_ENDOPLASMIC_RETICULUM_MEMBRANE',
    'GOCC_CORTICAL_ENDOPLASMIC_RETICULUM',
    'GOCC_ENDOPLASMIC_RETICULUM_GOLGI_INTERMEDIATE_COMPARTMENT_MEMBRANE',
    'GOCC_ENDOPLASMIC_RETICULUM_CHAPERONE_COMPLEX',
    'GOCC_NUCLEAR_OUTER_MEMBRANE_ENDOPLASMIC_RETICULUM_MEMBRANE_NETWORK',
    'GOCC_EXTRINSIC_COMPONENT_OF_ENDOPLASMIC_RETICULUM_MEMBRANE',
    'GOCC_MITOCHONDRIA_ASSOCIATED_ENDOPLASMIC_RETICULUM_MEMBRANE',
    'GOCC_ENDOPLASMIC_RETICULUM_QUALITY_CONTROL_COMPARTMENT',
    'GOCC_ENDOPLASMIC_RETICULUM_EXIT_SITE',
    'GOCC_INTEGRAL_COMPONENT_OF_CYTOPLASMIC_SIDE_OF_ENDOPLASMIC_RETICULUM_MEMBRANE',
    'GOCC_ENDOPLASMIC_RETICULUM_TUBULAR_NETWORK',
    'GOCC_PERINUCLEAR_ENDOPLASMIC_RETICULUM',
    'GOCC_LUMENAL_SIDE_OF_ENDOPLASMIC_RETICULUM_MEMBRANE',
    'GOCC_CYTOPLASMIC_SIDE_OF_ENDOPLASMIC_RETICULUM_MEMBRANE',
    'GOCC_ENDOPLASMIC_RETICULUM_TUBULAR_NETWORK_MEMBRANE',
    'GOCC_ENDOPLASMIC_RETICULUM_PLASMA_MEMBRANE_CONTACT_SITE',
    'GOCC_ENDOPLASMIC_RETICULUM_PROTEIN_CONTAINING_COMPLEX')

UPR <-
  c('GOBP_ATF6_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_NEGATIVE_REGULATION_OF_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_NEGATIVE_REGULATION_OF_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_NEGATIVE_REGULATION_OF_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_POSITIVE_REGULATION_OF_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_POSITIVE_REGULATION_OF_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_REGULATION_OF_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_REGULATION_OF_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_REGULATION_OF_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE')

UPR.proliferation.apoptosis.geneSets <-
  c('GOCC_RIBOSOME',
    'GOBP_RIBOSOME_DISASSEMBLY',
    'GOBP_RIBOSOME_BIOGENESIS',
    'GOBP_RIBOSOME_ASSEMBLY',
    'GOBP_REGULATION_OF_RIBOSOME_BIOGENESIS',
    'GOBP_NEGATIVE_REGULATION_OF_RIBOSOME_BIOGENESIS',
    'GOBP_RIBOSOME_ASSOCIATED_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS',
    
    'GOCC_PROTEASOME_COMPLEX',
    'GOBP_PROTEASOME_ASSEMBLY',
    
    'GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_ATF6_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    'GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    
    'GOBP_POSITIVE_REGULATION_OF_CELL_POPULATION_PROLIFERATION',
    'GOBP_NEGATIVE_REGULATION_OF_CELL_POPULATION_PROLIFERATION',
    'GOBP_REGULATION_OF_CELL_POPULATION_PROLIFERATION',
    'GOBP_EPITHELIAL_CELL_PROLIFERATION',
    'GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION',
    'GOBP_NEGATIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION',
    
    'GOBP_EXECUTION_PHASE_OF_APOPTOSIS',
    'GOBP_REGULATION_OF_EXECUTION_PHASE_OF_APOPTOSIS',
    'GOBP_NEGATIVE_REGULATION_OF_EXECUTION_PHASE_OF_APOPTOSIS',
    'GOBP_POSITIVE_REGULATION_OF_EXECUTION_PHASE_OF_APOPTOSIS',
    'GOBP_REGULATION_OF_CYSTEINE_TYPE_ENDOPEPTIDASE_ACTIVITY_INVOLVED_IN_EXECUTION_PHASE_OF_APOPTOSIS',
    'KEGG_APOPTOSIS',
    
    'GOBP_CELL_GROWTH',
    'GOBP_POSITIVE_REGULATION_OF_CELL_GROWTH',
    'GOBP_NEGATIVE_REGULATION_OF_CELL_GROWTH',
    'GOBP_GROWTH',
    'GOBP_REGULATION_OF_GROWTH',
    'GOBP_REGULATION_OF_GROWTH_RATE',
    'GOBP_NEGATIVE_REGULATION_OF_GROWTH',
    'GOBP_POSITIVE_REGULATION_OF_GROWTH',
    'GOBP_REGULATION_OF_EXTENT_OF_CELL_GROWTH')


get.CC.df <- function(geneSets, expData, sampleInfo, group){
  library(stringr)
  cor.selected <<-
    getGeneSetCor(geneSetNames = geneSets,
                  expData = expData,
                  sampleInfo = sampleInfo)
  
  cc.selected <- cor.selected$CC
  CC.median <- median(cc.selected)
  i <- sapply(lapply(lapply(str_split(geneSets, '_'), `[`, (-1)),
                     str_to_title), paste, collapse=' ')
  df.selected <- data.frame(comPairs=paste(i, collapse = ' VS '),
                            CC=cc.selected,
                            Group=str_wrap(paste0(i[1], ' VS ', group),
                                           width = 12))
  return(list(median=CC.median, df=df.selected))
}

plotBoxOne2manySets <- function(specified.geneSet='GOCC_RIBOSOME',
                                geneSets, expData=expData.XENA.CRC,
                                sampleInfo=sampleInfo.XENA.CRC,
                                group='Protein Folding Gene Sets',
                                CC.ribo.prot.df=CC.ribo.prot.df.XENA.CRC,
                                title='XENA.CRC'){
  CC.df.median.list <- list()
  for (i in geneSets) {
    CC.df.median <- get.CC.df(c(specified.geneSet, i),
                              expData, sampleInfo,
                              group=group)
    CC.df.median.list[[length(CC.df.median.list)+1]] <- CC.df.median
  }
  
  CC.medians <- c()
  CC.df.list <- list()
  for (j in CC.df.median.list) {
    CC.medians <- c(CC.medians, j[[1]])
    CC.df.list[[length(CC.df.list)+1]] <- j[[2]]
  }
  CC.df.list <- CC.df.list[order(CC.medians, decreasing = T)][1:10]
  CC.df <- CC.ribo.prot.df
  for (k in CC.df.list) {
    CC.df <- rbind(CC.df, k)
  }
  CC.df[,1] <- factor(CC.df[,1], levels = unique(CC.df[,1]))
  CC.df[,3] <- factor(CC.df[,3], levels = unique(CC.df[,3]))
  return(plotBox(CC.df = CC.df, comPair1 = 1,
                 comPair2 = 2:length(unique(CC.df[,1])),
                 show.legend = T, title=title))
}

for (i in c('Pr.Fold', 'Pr.Cata',
            'ER', 'UPR')) {
  for (j in list(list('expData.XENA.CRC', 'sampleInfo.XENA.CRC'))) {
    for (k in c('GOCC_RIBOSOME', 'GOCC_PROTEASOME_COMPLEX')) {
      assign(paste('boxplot', str_split(k, '_')[[1]][2], i,
                   str_split(j[[1]], '\\.')[[1]][2],sep = '.'),
             plotBoxOne2manySets(specified.geneSet = k, geneSets = get(i),
                                 expData = get(j[[1]]), sampleInfo = get(j[[2]]),
                                 group = i))
    }
  }
}

cor.hm.UPR.ribo.prot.XENA.CRC <-
  getGeneSetCor(geneSetNames =
                  c('GOCC_RIBOSOME', 'GOCC_PROTEASOME_COMPLEX',
                    'GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
                    'GOBP_ATF6_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
                    'GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE'),
                expData = expData.XENA.CRC,
                sampleInfo = sampleInfo.XENA.CRC,
                geneSetLabels = c('Ribosome', 'Proteasome', 'IRE1',
                                  'ATF6', 'PERK'))

getUPRdf <- function(List, listNo) {
  df <- rbind(data.frame(comPairs=names(List$CC)[listNo],
                         CC=List$CC[[listNo]],
                         Group=names(List$CC)[listNo]))
  #df$Group <- factor(df$Group, levels = group.levels)
  
  return(df)
}

plotBox.RiboProtUPR <-
  function(UPR.ribo.prot, selectedCol, random.seeds,
           expData, sampleInfo, specified.geneSet.name,
           specified.label, title='XENA.CRC', group.levels){
  UPR.ribo.prot.df <- getUPRdf(UPR.ribo.prot, selectedCol)
  random.CC.ribo.prot.df <-
    get.largest.random.CC.df(
      random.seeds = random.seeds,
      expdata = expData, sampleInfo = sampleInfo,
      specified.geneSet.name = specified.geneSet.name,
      specified.label = specified.label,
      specified.df = UPR.ribo.prot.df)
  random.CC.ribo.prot.df$Group <- factor(random.CC.ribo.prot.df$Group,
                                         levels = group.levels)
  boxplot.CC.ribo.prot.UPR <-
    plotBox(CC.df = random.CC.ribo.prot.df,
            comPair1 = 1, comPair2 = 2:11,
            title = title)
  return(boxplot.CC.ribo.prot.UPR)
}

boxplot.CC.ribo.ATF6.XENA.CRC <-
  plotBox.RiboProtUPR(
    UPR.ribo.prot = cor.hm.UPR.ribo.prot.XENA.CRC,
    random.seeds = random.seeds[404:503],
    expData = expData.XENA.CRC, specified.label = 'ATF6',
    sampleInfo = sampleInfo.XENA.CRC, selectedCol = 3,
    specified.geneSet.name = 'GOBP_ATF6_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    group.levels = c('Ribosome VS ATF6', 'ATF6 VS Rand'))
boxplot.CC.ribo.IRE1.XENA.CRC <-
  plotBox.RiboProtUPR(
    UPR.ribo.prot = cor.hm.UPR.ribo.prot.XENA.CRC,
    random.seeds = random.seeds[504:603],
    expData = expData.XENA.CRC, specified.label = 'IRE1',
    sampleInfo = sampleInfo.XENA.CRC, selectedCol = 2,
    specified.geneSet.name = 'GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    group.levels = c('Ribosome VS IRE1', 'IRE1 VS Rand'))
boxplot.CC.ribo.PERK.XENA.CRC <-
  plotBox.RiboProtUPR(
    UPR.ribo.prot = cor.hm.UPR.ribo.prot.XENA.CRC,
    random.seeds = random.seeds[604:703],
    expData = expData.XENA.CRC, specified.label = 'PERK', 
    sampleInfo = sampleInfo.XENA.CRC, selectedCol = 4,
    specified.geneSet.name = 'GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    group.levels = c('Ribosome VS PERK', 'PERK VS Rand'))
boxplot.CC.prot.ATF6.XENA.CRC <-
  plotBox.RiboProtUPR(
    UPR.ribo.prot = cor.hm.UPR.ribo.prot.XENA.CRC,
    random.seeds = random.seeds[704:803],
    expData = expData.XENA.CRC, specified.label = 'ATF6',
    sampleInfo = sampleInfo.XENA.CRC, selectedCol = 6,
    specified.geneSet.name = 'GOBP_ATF6_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    group.levels = c('Proteasome VS ATF6', 'ATF6 VS Rand'))
boxplot.CC.prot.IRE1.XENA.CRC <-
  plotBox.RiboProtUPR(
    UPR.ribo.prot = cor.hm.UPR.ribo.prot.XENA.CRC,
    random.seeds = random.seeds[804:903],
    expData = expData.XENA.CRC, specified.label = 'IRE1',
    sampleInfo = sampleInfo.XENA.CRC, selectedCol = 5,
    specified.geneSet.name = 'GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    group.levels = c('Proteasome VS IRE1', 'IRE1 VS Rand'))
boxplot.CC.prot.PERK.XENA.CRC <-
  plotBox.RiboProtUPR(
    UPR.ribo.prot = cor.hm.UPR.ribo.prot.XENA.CRC,
    random.seeds = random.seeds[904:1003],
    expData = expData.XENA.CRC, specified.label = 'PERK',
    sampleInfo = sampleInfo.XENA.CRC, selectedCol = 7,
    specified.geneSet.name = 'GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
    group.levels = c('Proteasome VS PERK', 'PERK VS Rand'))


# ---------- Part VIII: Analysis of other cancer types in XENA ----------
library(data.table)
library(stringr)

ensemblID.related <-
  unique(c(ensemblID.RFE.XENA.CRC,
    ensemblID.ribo.prot,
    unlist(lapply(Pr.Fold, getEnsemblIDsOfGeneSet)),
    unlist(lapply(Pr.Cata, getEnsemblIDsOfGeneSet)),
    unlist(lapply(UPR, getEnsemblIDsOfGeneSet)),
    unlist(lapply(ER, getEnsemblIDsOfGeneSet)),
    unlist(lapply(UPR.proliferation.apoptosis.geneSets,
                  getEnsemblIDsOfGeneSet))))

get.theMostRelevant.geneSets <- function(specifiedGeneSetName, expData, sampleInfo){
  cor.with.speSet <- list()
  geneSetNames <- setdiff(geneSetNames, specifiedGeneSetName)
  for (i in geneSetNames) {
    cor.i <-
      getGeneSetCor(geneSetNames = c(specifiedGeneSetName, i),
                    expData = expData,
                    sampleInfo = sampleInfo)
    if (length(cor.i$CC)!=0) {
      if (median(abs(cor.i$CC))>0.8) {
        cor.with.speSet[[i]] <- cor.i
      }
    }
  }
  return(cor.with.speSet)
}

phenotype.GTExTCGA <- fread('TcgaTargetGTEX_phenotype.txt',
                            check.names = F, data.table = F)
expGTExTCGA <-
  fread(
    'E:/PublicData/DataFromXena/TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2',
    stringsAsFactors=F, header=T, sep='\t', nThread = 12, skip = 0, data.table = F)
rownames(expGTExTCGA) <- expGTExTCGA$V1
expGTExTCGA$V1 <- NULL
rownames(expGTExTCGA) <- str_sub(rownames(expGTExTCGA), 1, 15)
expGTExTCGA.related <- expGTExTCGA[rownames(expGTExTCGA) %in%
                                     ensemblID.related,]
rm(expGTExTCGA); gc()

processGTExTCGAall <- function(primaryDisease){
  sampleType=c('Normal Tissue', 'Solid Tissue Normal',
               'Primary Tumor')
  phenotype <- phenotype.GTExTCGA[phenotype.GTExTCGA$sampleType %in%
                                    sampleType, ]
  primarySite <-
    unique(phenotype$primarySite[phenotype$primaryDisease==primaryDisease])
  phenotype <- phenotype[phenotype$primaryDisease==primaryDisease|
                           (phenotype$primarySite==primarySite&
                              phenotype$sampleType == 'Normal Tissue'),]
  #fread() get data.table class, data frame transformation is needed
  #for downstream analysis.
  expData <- expGTExTCGA.related
  sampleIDsSelected <- intersect(phenotype$sampleID,
                                 colnames(expData))
  if (length(sampleIDsSelected)!=0) {
    phenotype <- phenotype[phenotype$sampleID %in% sampleIDsSelected,]
    expData <- expData[, colnames(expData) %in% sampleIDsSelected]
    expData <- t(expData)
    nNor <- paste0('Nor (n=', sum(phenotype$sampleType=='Normal Tissue'), ')')
    nAdj <- paste0('Adj (n=', sum(phenotype$sampleType=='Solid Tissue Normal'), ')')
    nCan <- paste0('Can (n=', sum(phenotype$sampleType=='Primary Tumor'), ')')
    phenotype$sampleType <-
      ifelse(phenotype$sampleType=='Normal Tissue', nNor,
             ifelse(phenotype$sampleType=='Primary Tumor',
                    nCan, nAdj))
    phenotype$sampleType <- factor(phenotype$sampleType,
                                   levels = c(nNor, nAdj, nCan))
    phenotype <- phenotype[order(phenotype$sampleType),]
    expData <- expData[match(phenotype$sampleID,
                             rownames(expData)),]
    library(preprocessCore)
    return(list(phenotype=phenotype, expData=expData))
  } else {
    return('Invalid data set')
  }
}

for (i in unique(phenotype.GTExTCGA$primaryDisease)) {
  tmp <- processGTExTCGAall(primaryDisease = i)
  if (!is.character(tmp)) {
    if (length(unique(tmp$phenotype$sampleType))==3) {
      if (min(table(tmp$phenotype$sampleType)) > 3) {
        assign(paste0('exp.pht.', i), tmp)
      }
    }
  }
}

library(factoextra)
random.seeds.start <- 1004
for (i in ls()[str_detect(ls(), '^exp.pht.+')]) {
  exp.pht <- get(i)
  i <- str_split(i, '\\.')[[1]][3]
  exp.pht$expData <-
    exp.pht$expData[, sapply(apply(exp.pht$expData,
                                   2,unique), length)!=1]
  
  # PCA and outliers remove.
  assign(paste0('PCAres.XENA.', i, '.allSamples'),
         prcomp(as.matrix(exp.pht$expData), scale. = T))
  assign(paste0('PCAplot.XENA.', i, '.allSamples'),
         fviz_pca_ind(get(paste0('PCAres.XENA.', i, '.allSamples')),
                      col.ind=exp.pht$phenotype$sampleType,
                      mean.point=F, addEllipses = T, repel=T,
                      geom='point', legend.title="Sample Type",
                      title=paste0('XENA.', i, '\nAll samples')))
  outliers.index <-
    which(abs(get(paste0('PCAres.XENA.', i, '.allSamples'))[[5]][,1])>250)
  if (length(outliers.index)!=0) {
    exp.pht$phenotype <- exp.pht$phenotype[-outliers.index,]
    exp.pht$expData <- exp.pht$expData[-outliers.index,]
  }
  #exp.pht <- cleanData(expData = exp.pht$expData,
  #                    sampleInfo = exp.pht$phenotype)
  #exp.pht$phenotype <- exp.pht$sampleInfo
  nNor <- paste0('Nor (n=',
                 sum(str_detect(exp.pht$phenotype$sampleType, 'Nor')), ')')
  nAdj <- paste0('Adj (n=',
                 sum(str_detect(exp.pht$phenotype$sampleType, 'Adj')), ')')
  nCan <- paste0('Can (n=',
                 sum(str_detect(exp.pht$phenotype$sampleType, 'Can')), ')')
  exp.pht$phenotype$sampleType <-
    ifelse(str_detect(exp.pht$phenotype$sampleType, 'Nor'), nNor,
           ifelse(str_detect(exp.pht$phenotype$sampleType, 'Can'),
                  nCan, nAdj))
  exp.pht$phenotype$sampleType <-
    factor(exp.pht$phenotype$sampleType, levels = c(nNor, nAdj, nCan))
  assign(paste0('PCAres.XENA.', i, '.outliers.removed'),
         prcomp(as.matrix(exp.pht$expData),scale. = T))
  assign(paste0('PCAplot.XENA.', i, '.outliers.removed'),
         fviz_pca_ind(get(paste0('PCAres.XENA.', i,
                                 '.outliers.removed')),
                      col.ind=exp.pht$phenotype$sampleType,
                      mean.point=F, addEllipses = T, repel=T,
                      geom='point', legend.title="Sample Type",
                      title=paste0('XENA.', i, '\nOutliers removed')))
  
  # Plot ROC and LDA model.
  expData.RFE <- exp.pht$expData[, colnames(exp.pht$expData) %in% ensemblID.RFE.XENA.CRC]
  sampleInfo.RFE <- exp.pht$phenotype[,c(1,4)]
  outliers <- c(names(boxplot(predict(ldaModel.RFE.XENA.CRC, expData.RFE)$x[,1])$out),
                names(boxplot(predict(ldaModel.RFE.XENA.CRC, expData.RFE)$x[,2])$out))
  expData.RFE.outliers.removed <- expData.RFE[!(rownames(expData.RFE) %in% outliers),]
  sampleInfo.RFE.outliers.removed <- sampleInfo.RFE[!(sampleInfo.RFE$sampleID %in% outliers),]
  sampleInfo.RFE.outliers.removed <- getSampleNumbers(sampleInfo.RFE.outliers.removed)[[1]]
  if (min(table(sampleInfo.RFE.outliers.removed$sampleType))==0) {
    assign(paste0('ROCplot.RFE.XENA.', i),
           plotROC(LDAmodel=ldaModel.RFE.XENA.CRC,
                   title=paste0('RFE.XENA.', i),
                   testSet=expData.RFE, testGroup=sampleInfo.RFE))
    assign(paste0('LDAplot.RFE.XENA.', i),
           plotLDAmodel(ldaModel = ldaModel.RFE.XENA.CRC,
                        trainSet = expData.RFE,
                        trainGroup = sampleInfo.RFE,
                        title = paste0('RFE.XENA.', i),
                        legend.position = c(0.05, 0.05)))
  } else {
    assign(paste0('ROCplot.RFE.XENA.', i),
           plotROC(LDAmodel=ldaModel.RFE.XENA.CRC,
                   title=paste0('RFE.XENA.', i),
                   testSet=expData.RFE.outliers.removed,
                   testGroup=sampleInfo.RFE.outliers.removed))
    assign(paste0('LDAplot.RFE.XENA.', i),
           plotLDAmodel(ldaModel = ldaModel.RFE.XENA.CRC,
                        trainSet = expData.RFE.outliers.removed,
                        trainGroup = sampleInfo.RFE.outliers.removed,
                        title = paste0('RFE.XENA.', i),
                        legend.position = c(0.05, 0.05)))
  }
  expData.iRFE <- exp.pht$expData[, colnames(exp.pht$expData) %in% ensemblID.iRFE.XENA.CRC]
  sampleInfo.iRFE <- exp.pht$phenotype[,c(1,4)]
  outliers <- c(names(boxplot(predict(ldaModel.iRFE.XENA.CRC, expData.iRFE)$x[,1])$out),
                names(boxplot(predict(ldaModel.iRFE.XENA.CRC, expData.iRFE)$x[,2])$out))
  expData.iRFE.outliers.removed <- expData.iRFE[!(rownames(expData.iRFE) %in% outliers),]
  sampleInfo.iRFE.outliers.removed <- sampleInfo.iRFE[!(sampleInfo.iRFE$sampleID %in% outliers),]
  sampleInfo.iRFE.outliers.removed <- getSampleNumbers(sampleInfo.iRFE.outliers.removed)[[1]]
  if (min(table(sampleInfo.iRFE.outliers.removed$sampleType))==0) {
    assign(paste0('ROCplot.iRFE.XENA.', i),
           plotROC(LDAmodel=ldaModel.iRFE.XENA.CRC,
                   title=paste0('iRFE.XENA.', i),
                   testSet=expData.iRFE, testGroup=sampleInfo.iRFE))
    assign(paste0('LDAplot.iRFE.XENA.', i),
           plotLDAmodel(ldaModel = ldaModel.iRFE.XENA.CRC,
                        trainSet = expData.iRFE,
                        trainGroup = sampleInfo.iRFE,
                        title = paste0('iRFE.XENA.', i),
                        legend.position = c(0.05, 0.05)))
  } else {
    assign(paste0('ROCplot.iRFE.XENA.', i),
           plotROC(LDAmodel=ldaModel.iRFE.XENA.CRC,
                   title=paste0('iRFE.XENA.', i),
                   testSet=expData.iRFE.outliers.removed,
                   testGroup=sampleInfo.iRFE.outliers.removed))
    assign(paste0('LDAplot.iRFE.XENA.', i),
           plotLDAmodel(ldaModel = ldaModel.iRFE.XENA.CRC,
                        trainSet = expData.iRFE.outliers.removed,
                        trainGroup = sampleInfo.iRFE.outliers.removed,
                        title = paste0('iRFE.XENA.', i),
                        legend.position = c(0.05, 0.05)))
  }
  
  if (F) {
    # Get the correlation coefficients between all gene sets and ribo/prot genes.
    assign(paste0('cor.with.ribosome.', i),
           get.theMostRelevant.geneSets(specifiedGeneSetName = 'GOCC_RIBOSOME',
                                        expData = exp.pht$expData,
                                        sampleInfo = exp.pht$phenotype))
    assign(paste0('cor.with.proteasome.', i),
           get.theMostRelevant.geneSets(specifiedGeneSetName = 'GOCC_PROTEASOME_COMPLEX',
                                        expData = exp.pht$expData,
                                        sampleInfo = exp.pht$phenotype)) 
    
  }
  
  # Get the CC between ribosome and proteasome.
  assign(paste0('ribo.prot.cor.hm.', i),
         getGeneSetCor(
           geneSetNames = c('GOCC_RIBOSOME', 'GOCC_PROTEASOME_COMPLEX'),
           geneSetLabels = c('Ribosome', 'Proteasome'),
           expData = exp.pht$expData, sampleInfo = exp.pht$phenotype))
  
  # Get the CC between random selected gene sets and ribo/prot genes.
  cor.hm <- get(paste0('ribo.prot.cor.hm.', i))
  CC.ribo.prot.df <- data.frame(comPairs='Ribo VS Prot',
                                CC=cor.hm$CC, Group='Ribo VS Prot')
  random.CC.ribo.df <-
    get.largest.random.CC.df(
      random.seeds = random.seeds[random.seeds.start:
                                    (random.seeds.start+99)],
      specified.geneSet.name = 'GOCC_RIBOSOME',
      specified.label = 'Ribo', expdata = exp.pht$expData,
      sampleInfo = exp.pht$phenotype, specified.df = CC.ribo.prot.df,
      group.levels = c('Ribo VS Prot', 'Ribo VS Rand'))
  assign(paste0('boxplot.CC.ribo.', i),
         plotBox(CC.df = random.CC.ribo.df,
                 comPair1 = 1, comPair2 = 2:11,
                 show.legend = T, title=i))
  
  random.CC.prot.df <-
    get.largest.random.CC.df(
      random.seeds = random.seeds[(random.seeds.start+100):
                                    (random.seeds.start+199)],
      specified.geneSet.name = 'GOCC_PROTEASOME_COMPLEX',
      specified.label = 'Prot', expdata = exp.pht$expData,
      sampleInfo = exp.pht$phenotype, specified.df = CC.ribo.prot.df,
      group.levels = c('Ribo VS Prot', 'Prot VS Rand'))
  assign(paste0('boxplot.CC.prot.', i),
         plotBox(CC.df = random.CC.prot.df,
                 comPair1 = 1, comPair2 = 2:11, show.legend = T,
                 title=i))
  
  # Plot corplot and boxplot of PCA scores.
  assign(paste0('PCAscores.ribo.', i),
         get_PCA_scores(
           expData = exp.pht$expData,
           sampleType = exp.pht$phenotype$sampleType,
           geneSetName = 'GOCC_RIBOSOME',
           title = paste0('Ribo.', i)))
  assign(paste0('PCAscores.prot.', i),
         get_PCA_scores(
           expData = exp.pht$expData,
           sampleType = exp.pht$phenotype$sampleType,
           geneSetName = 'GOCC_PROTEASOME_COMPLEX',
           title = paste0('Prot', i)))
  PCAscores.ribo.prot <-
    data.frame(
      PCAscores.Ribosome=get(paste0('PCAscores.ribo.', i))[[1]],
      PCAscores.Proteasome=get(paste0('PCAscores.prot.', i))[[1]],
      sampleType=exp.pht$phenotype$sampleType)
  library(ggpmisc)
  library(ggpubr)
  assign(paste0('corPlot.ribo.prot.', i),
         plotCor.Classes(PCAscores.ribo.prot, title = i))
  
  PCAscores.diff.ribo.prot <-
    data.frame(sampleType=exp.pht$phenotype$sampleType,
               diff=PCAscores.ribo.prot$PCAscores.Ribosome -
                 PCAscores.ribo.prot$PCAscores.Proteasome, Group=1)
  assign(paste0('boxplot.diff.ribo.prot.', i),
         plotBox(CC.df = PCAscores.diff.ribo.prot,
                 comPairIndex = list(c(1, 2), c(1, 3), c(2, 3)),
                 show.legend = F,
                 title=tail(unlist(str_split(i,'\\.')),1),
                 ylabel = 'PCAscores.Ribosome - PCAscores.Proteasome'))
  
  # Plot boxplot of CC of related gene sets.
  for (j in c('Pr.Fold', 'Pr.Cata', 'ER')) {
    for (k in c('GOCC_RIBOSOME', 'GOCC_PROTEASOME_COMPLEX')) {
      assign(
        paste('boxplot', str_split(k, '_')[[1]][2], j,
              i, sep = '.'),
        plotBoxOne2manySets(
          specified.geneSet = k, geneSets = get(j),
          expData = exp.pht$expData, sampleInfo = exp.pht$phenotype,
          group = j, CC.ribo.prot.df = CC.ribo.prot.df, title = i))
    }
  }

  for (j in c('GOBP_ATF6_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
              'GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE',
              'GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE')) {
    geneSetJ <- str_split(j, '_')[[1]][2]
    tmp <- get_PCA_scores(
      expData = exp.pht$expData, geneSetName = j,
      sampleType = exp.pht$phenotype$sampleType,
      title = paste0(geneSetJ, '.XENA\n', i))
    tmp <- data.frame(sampleID=names(tmp[[1]]),
                      PCAscores=tmp[[1]], group=1)
    tmp <- merge(exp.pht$phenotype[,c(1,4)], tmp, by='sampleID')
    tmp$sampleID <- NULL

    assign(paste0('boxplot.PCAscoresof', j, 'in.', i),
           plotBox(CC.df = tmp, show.legend = F,
                   ylabel = paste('PCAscores of Gene Exp Levels\n(',
                                  geneSetJ, 'Mediated UPR)'), title = i))
  }
  print(paste(i, 'completed!'))
  random.seeds.start <- random.seeds.start+200
}

reduce <- function(cor.with){
  for (i in names(cor.with)) {
    cor.with[[i]] <- cor.with[[i]][['CC']]
  }
  return(cor.with)
}
for (i in ls()[str_detect(ls(), 'cor.with')]) {
  assign(i, reduce(get(i)))
}

save.image()


# ---------- Part XI: Analyze validate data sets with 3 sample types ----------
processGEOdataset <-
  function(expDataFile, expDataSkip, GPLFile=NULL,
           GPLSkip=NULL, GPLSelectedCol=NULL,
           GPLobject=NULL, sampleInfoFile,
           sampleTypeLevels, geneIDtype='ENTREZID', takeLog=F){
    library(org.Hs.eg.db)
    library(data.table)
    library(clusterProfiler)
    expData <- fread(file = expDataFile, skip = expDataSkip,
                     header = T, nThread = 10, data.table = F)

    if (is.null(GPLFile) & (!is.null(GPLobject))) {
      GPL <- GPLobject
      commonIDs <- intersect(expData[,1], GPL[,1])
      expData <- expData[expData[,1] %in% commonIDs, ]
      GPL <- GPL[GPL[,1] %in% commonIDs,]
      GPL <- GPL[match(expData[,1], GPL[,1]), ]
      geneIDs.ensemblIDs <- bitr(GPL[,2],fromType = geneIDtype,
                                 toType = 'ENSEMBL',OrgDb = org.Hs.eg.db)
      expData[,1] <-
        sapply(sapply(GPL[,2],
                      function(x) geneIDs.ensemblIDs[geneIDs.ensemblIDs[,1]==x, 2]),
               `[`, 1)
    } else if ((!is.null(GPLFile)) & is.null(GPLobject)) {
      GPL <-
        fread(file = GPLFile, skip = GPLSkip, header = T, nThread = 10)
      GPL <- na.omit(data.frame(GPL, check.names = F)[, GPLSelectedCol])
      commonIDs <- intersect(expData[,1], GPL[,1])
      expData <- expData[expData[,1] %in% commonIDs, ]
      GPL <- GPL[GPL[,1] %in% commonIDs,]
      GPL <- GPL[match(expData[,1], GPL[,1]), ]
      geneIDs.ensemblIDs <- bitr(GPL[,2],fromType = geneIDtype,
                                 toType = 'ENSEMBL',OrgDb = org.Hs.eg.db)
      expData[,1] <-
        sapply(sapply(GPL[,2],
                      function(x) geneIDs.ensemblIDs[geneIDs.ensemblIDs[,1]==x, 2]),
               `[`, 1)
    } else {
      geneIDs.ensemblIDs <-
        bitr(expData[,1], fromType = geneIDtype,
             toType = 'ENSEMBL', OrgDb = org.Hs.eg.db)
      expData[,1] <-
        sapply(sapply(expData[,1],
                      function(x) geneIDs.ensemblIDs[geneIDs.ensemblIDs[,1]==x, 2]),
               `[`, 1)
    }
    expData <- expData[!is.na(expData[,1]),]
    expData <- aggregate(expData, FUN=mean,
                         na.rm=T, by=list(expData[,1]))
    rownames(expData) <- expData$Group.1
    expData <- expData[,-c(1,2)]
    library(preprocessCore)
    expData <- normalize.quantiles.robust(
      as.matrix(expData), remove.extreme = 'mean',
      use.log2 = takeLog, keep.names = T)
    expData <- t(expData)
    sampleInfo <-
      read.csv(sampleInfoFile, header = T)
    sampleInfo$sampleType <- factor(sampleInfo$sampleType,
                                    levels = c(sampleTypeLevels))
    sampleInfo <- sampleInfo[order(sampleInfo$sampleType),]
    expData <- expData[match(sampleInfo[,1], rownames(expData)),]
    if (takeLog==T) {
      expData <- log2(expData + 1e-10)
    }
    expData <- expData[, colnames(expData) %in% ensemblID.related]
    expData <- expData[, !(sapply(apply(expData,2,unique), length)==1)]
    return(list(expData=expData, sampleInfo=sampleInfo))
}

corAnalysis <-
  function(expData, sampleInfo, dataSetName,
           geneSet1='GOCC_RIBOSOME',
           geneSet2='GOCC_PROTEASOME_COMPLEX',
           setName1='ribosome', setName2='proteasome',
           comPair1=NULL, comPair2=NULL,
           comPairs=NULL, as.title=T){
    library(stringr)
    if ((!is.null(geneSet1)) & (!is.null(geneSet2))) {
      if (length(intersect(getEnsemblIDsOfGeneSet(geneSet1), colnames(expData)))>1 &
          length(intersect(getEnsemblIDsOfGeneSet(geneSet2), colnames(expData)))>1) {
        cor <-
          getGeneSetCor(expData = expData, sampleInfo = sampleInfo,
                        geneSetNames = c(geneSet1, geneSet2),
                        geneSetLabels = c(setName1, setName2))
        
        PCAscores.set1 <-
          get_PCA_scores(
            expData = expData,
            sampleType = sampleInfo$sampleType,
            geneSetName = geneSet1,
            title = setName1)
        PCAscores.set2 <-
          get_PCA_scores(
            expData = expData,
            sampleType = sampleInfo$sampleType,
            geneSetName = geneSet2,
            title = setName2)
        PCAscores.set1.set2 <-
          data.frame(set1=PCAscores.set1[[1]],
                     set2=PCAscores.set2[[1]],
                     sampleType=sampleInfo$sampleType)
        library(ggpmisc)
        library(ggpubr)
        colnames(PCAscores.set1.set2) <-
          c(paste0('PCAscores.', setName1),
            paste0('PCAscores.', setName2),
            'sampleType')
        corPlot.PCAscores <- plotCor.Classes(PCAscores.set1.set2,
                                             title = dataSetName)
        colnames(PCAscores.set1.set2) <- c('set1', 'set2', 'sampleType')
        PCAscores.diff <-
          data.frame(sampleType=PCAscores.set1.set2$sampleType,
                     diff=PCAscores.set1.set2$set1 -
                       PCAscores.set1.set2$set2, Group=1)
        if (is.null(comPair1) & is.null(comPair2) & is.null(comPairs)){
          comPairs <- list()
          for (i in 1:(length(levels(sampleInfo$sampleType))-1)) {
            for (j in (i+1):length(levels(sampleInfo$sampleType))) {
              comPairs[[length(comPairs)+1]] <- c(i, j)
            }
          }
          boxplot.PCAscores.diff <-
            plotBox(CC.df = PCAscores.diff,
                    comPairIndex = comPairs, show.legend = F,
                    ylabel = paste0('PCAscores.',setName1,
                                    '\n- PCAscores.',setName2),
                    title = dataSetName)
        } else if (is.null(comPairs)) {
          boxplot.PCAscores.diff <-
            plotBox(
              CC.df = PCAscores.diff, comPair1 = comPair1,
              comPair2 = comPair2, show.legend = F,
              ylabel = paste0('PCAscores.',setName1,
                              '\n- PCAscores.',setName2),
              title = dataSetName)
        } else if (!is.null(comPairs)) {
          boxplot.PCAscores.diff <-
            plotBox(
              CC.df = PCAscores.diff, comPairIndex = comPairs,
              show.legend = F,
              ylabel = paste0('PCAscores.',setName1,
                              '\n- PCAscores.',setName2),
              title = dataSetName)
        }
        return(list(cor, corPlot.PCAscores,
                    boxplot.PCAscores.diff))
      }
    } else if (is.null(geneSet2)) {
      if (length(intersect(getEnsemblIDsOfGeneSet(geneSet1), colnames(expData)))>1) {
        PCAscores.set1 <-
          get_PCA_scores(
            expData = expData,
            sampleType = sampleInfo$sampleType,
            geneSetName = geneSet1,
            title = setName1)
        PCAscores.set1.df <-
          data.frame(sampleType=sampleInfo$sampleType,
                     set1=PCAscores.set1[[1]], Group=1)
        colnames(PCAscores.set1.df)[1] <- setName1
        if (is.null(comPairs1) & is.null(comPairs2) & is.null(comPairs)) {
          comPairs <- list()
          for (i in 1:(length(levels(sampleInfo$sampleType))-1)) {
            for (j in (i+1):length(levels(sampleInfo$sampleType))) {
              comPairs[[length(comPairs)+1]] <- c(i, j)
            }
          }
          boxplot.PCAscores <-
            plotBox(
              CC.df = PCAscores.set1.df,
              comPairIndex = comPairs,
              show.legend = F,
              ylabel = paste0('PCAscores.',
                              str_wrap(str_replace_all(setName1,
                                                       '_', ' '), 20)),
              title = dataSetName, as.title = F)
        } else if (is.null(comPairs)){
          boxplot.PCAscores <-
            plotBox(CC.df = PCAscores.set1.df,
                    comPair1 = comPair1,
                    comPair2 = comPair2,
                    show.legend = F,
                    ylabel = paste0('PCAscores.',
                                    str_wrap(str_replace_all(setName1,
                                                             '_', ' '), 20)),
                    title = dataSetName, as.title = F)
        } else if (!is.null(comPairs)) {
          boxplot.PCAscores <-
            plotBox(CC.df = PCAscores.set1.df,
                    comPairIndex = comPairs,
                    show.legend = F,
                    ylabel = paste0('PCAscores.',
                                    str_wrap(str_replace_all(setName1,
                                                             '_', ' '), 20)),
                    title = dataSetName, as.title = F)
        }
        return(boxplot.PCAscores)
      }
    }
}

setwd('D:/Projects/CRC_expData_LDA_Carcinogenesis/validataData/')
# GSE10927
GSE10927.data <- 
  processGEOdataset(
    expDataFile = 'GSE10927/GSE10927_series_matrix.txt',
    expDataSkip = 86, GPLFile = 'GSE10927/GPL570-55999.txt', 
    GPLSkip = 16, GPLSelectedCol = c(1,12), geneIDtype = 'ENTREZID',
    sampleInfoFile = 'GSE10927/GSE10927_sampleInfo.csv',
    sampleTypeLevels = c('Nor (n=10)', 'Ade (n=22)', 'Can (n=33)'))
cor.GSE10927.ribo.prot <-
  corAnalysis(GSE10927.data$expData,
              GSE10927.data$sampleInfo,
              dataSetName = 'GSE10927\nAdrenocortical Carcinomas')

# GSE14323
GSE14323.GPL571.data <- 
  processGEOdataset(
    expDataFile = 'GSE14323/GSE14323-GPL571_series_matrix.txt',
    expDataSkip = 77, GPLFile = 'GSE14323/GPL571-17391.txt', 
    GPLSkip = 16, GPLSelectedCol = c(1,12), geneIDtype = 'ENTREZID',
    sampleInfoFile = 'GSE14323/GSE14323_GPL571_sampleInfo.csv',
    sampleTypeLevels = c('Nor (n=19)', 'crr (n=41)',
                         'crrHCC (n=17)', 'HCC (n=47)'))
GSE14323.GPL96.data <-
  processGEOdataset(
    expDataFile = 'GSE14323/GSE14323-GPL96_series_matrix.txt',
    expDataSkip = 77, GPLFile = 'GSE14323/GPL96-57554.txt',
    GPLSkip = 16, GPLSelectedCol = c(1,12), geneIDtype = 'ENTREZID',
    sampleInfoFile = 'GSE14323/GSE14323_GPL96_sampleInfo.csv',
    sampleTypeLevels = c('Nor (n=19)', 'crr (n=41)',
                         'crrHCC (n=17)', 'HCC (n=47)'))
GSE14323.data <- list()
GSE14323.data$expData <- rbind(GSE14323.GPL571.data$expData,
                               GSE14323.GPL96.data$expData)
GSE14323.data$sampleInfo <- rbind(GSE14323.GPL571.data$sampleInfo,
                                  GSE14323.GPL96.data$sampleInfo)
rm(list=c('GSE14323.GPL571.data', 'GSE14323.GPL96.data'))
cor.GSE14323.ribo.prot <-
  corAnalysis(GSE14323.data$expData,
              GSE14323.data$sampleInfo,
              dataSetName = 'GSE14323\nHepatocellular Carcinoma')

# GSE27678
GSE27678.GPL570.data <- 
  processGEOdataset(
    expDataFile = 'GSE27678/GSE27678-GPL570_series_matrix.txt',
    expDataSkip = 81, GPLFile = 'GSE27678/GPL570-55999.txt', 
    GPLSkip = 16, GPLSelectedCol = c(1,12), geneIDtype = 'ENTREZID',
    sampleInfoFile = 'GSE27678/GSE27678_GPL570_sampleInfo.csv',
    sampleTypeLevels = c('Nor (n=15)', 'LSIL (n=11)',
                         'HSIL (n=21)', 'Can (n=28)'))
GSE27678.GPL571.data <-
  processGEOdataset(
    expDataFile = 'GSE27678/GSE27678-GPL571_series_matrix.txt',
    expDataSkip = 81, GPLFile = 'GSE27678/GPL571-17391.txt',
    GPLSkip = 16, GPLSelectedCol = c(1,12), geneIDtype = 'ENTREZID',
    sampleInfoFile = 'GSE27678/GSE27678_GPL571_sampleInfo.csv',
    sampleTypeLevels = c('Nor (n=15)', 'LSIL (n=11)',
                         'HSIL (n=21)', 'Can (n=28)'))
commonEnsemblIDs <- intersect(colnames(GSE27678.GPL570.data$expData),
                              colnames(GSE27678.GPL571.data$expData))
GSE27678.GPL570.data$expData <-
  GSE27678.GPL570.data$expData[, colnames(GSE27678.GPL570.data$expData) %in%
                                 commonEnsemblIDs]
GSE27678.GPL571.data$expData <-
  GSE27678.GPL571.data$expData[, colnames(GSE27678.GPL571.data$expData) %in%
                                 commonEnsemblIDs]

GSE27678.data <- list()
GSE27678.data$expData <- rbind(GSE27678.GPL570.data$expData,
                               GSE27678.GPL571.data$expData)
GSE27678.data$sampleInfo <- rbind(GSE27678.GPL570.data$sampleInfo,
                                  GSE27678.GPL571.data$sampleInfo)
GSE27678.data$sampleInfo <-
  GSE27678.data$sampleInfo[order(GSE27678.data$sampleInfo$sampleType),]
GSE27678.data$expData <-
  GSE27678.data$expData[match(GSE27678.data$sampleInfo$sampleID,
                              rownames(GSE27678.data$expData)),]
rm(list=c('GSE27678.GPL570.data', 'GSE27678.GPL571.data'))
cor.GSE27678.ribo.prot <-
  corAnalysis(GSE27678.data$expData,
              GSE27678.data$sampleInfo,
              dataSetName = 'GSE27678\nCervical Carcinoma')

# GSE41657
GSE41657.data <- 
  processGEOdataset(
    expDataFile = 'GSE41657/GSE41657_series_matrix.txt',
    expDataSkip = 71, GPLFile = 'GSE41657/GPL6480-9577.txt', 
    GPLSkip = 17, GPLSelectedCol = c(1,6), geneIDtype = 'ENTREZID',
    sampleInfoFile = 'GSE41657/sampleInfo_GSE41657.csv',
    sampleTypeLevels = c('Nor (n=12)', 'AdeL (n=21)',
                         'AdeH (n=30)', 'Can (n=25)'))
cor.GSE41657.ribo.prot <-
  corAnalysis(expData = GSE41657.data$expData,
              sampleInfo = GSE41657.data$sampleInfo,
              dataSetName = 'GSE41657\nColorectal Cancer')

# GSE44076
GSE44076.data <-
  processGEOdataset(
    expDataFile = 'GSE44076/GSE44076_series_matrix.txt',
    expDataSkip = 80, GPLFile = 'GSE44076/GPL13667-15572.txt',
    GPLSkip = 43, GPLSelectedCol = c(1,21), geneIDtype = 'ENTREZID',
    sampleInfoFile = 'GSE44076/sampleInfo_GSE44076.csv',
    sampleTypeLevels = c('Nor (n=50)', 'Adj (n=98)', 'Can (n=98)'))
cor.GSE44076.ribo.prot <-
  corAnalysis(expData = GSE44076.data$expData,
              sampleInfo = GSE44076.data$sampleInfo,
              dataSetName = 'GSE44076\nColon Adenocarcinoma')

# GSE64041
GPL6244 <- fread('GSE64041/GPL6244-17930.txt',
                 skip = 12, data.table = F)[,c(1,10)]
GPL6244$gene_assignment <-
  sapply(str_split(GPL6244$gene_assignment,
                   ' // '), `[`, 2)
GSE64041 <- fread('GSE64041/GSE64041_series_matrix.txt',
                  skip = 57, data.table = F)
GSE64041 <- merge(GSE64041, GPL6244, by.x='ID_REF', by.y='ID')
GSE64041$ID_REF <- GSE64041$gene_assignment
names(GSE64041)[1] <- 'GeneSymbol'
GSE64041$gene_assignment <- NULL
GSE64041 <- na.omit(GSE64041)
write.table(GSE64041, 'GSE64041/GSE64041.txt',
            sep = '\t', row.names = F)
GSE64041.data <-
  processGEOdataset(expDataFile = 'GSE64041/GSE64041.txt',
                    expDataSkip = 0, geneIDtype = 'SYMBOL',
                    sampleInfoFile = 'GSE64041/GSE64041_sampleInfo.csv',
                    sampleTypeLevels = c('Nor (n=5)', 'Adj (n=60)',
                                         'Can (n=60)'))
cor.GSE64041.ribo.prot <-
  corAnalysis(expData = GSE64041.data$expData,
              sampleInfo = GSE64041.data$sampleInfo,
              dataSetName = 'GSE64041\nHepatocellular Carcinoma')

# GSE91035
GSE91035.data <-
  processGEOdataset(
    expDataFile = 'GSE91035/GSE91035_normalized_data_with_gene_symbol.txt',
    expDataSkip = 0, geneIDtype = 'SYMBOL',
    sampleInfoFile = 'GSE91035/GSE91035_sampleInfo.csv',
    sampleTypeLevels = c('Nor (n=8)', 'Ben (n=15)', 'Can (n=25)'))
cor.GSE91035.ribo.prot <-
  corAnalysis(expData = GSE91035.data$expDat,
              sampleInfo = GSE91035.data$sampleInfo,
              dataSetName = 'GSE91035\nPancreatic Ductal Adenocarcinoma')

# GSE95132
GSE95132.data <-
  processGEOdataset(
    expDataFile = 'GSE95132/GSE95132_RNAseq_read_counts.txt',
    expDataSkip = 0, geneIDtype = 'SYMBOL',
    sampleInfoFile = 'GSE95132/GSE95132_sampleInfo.csv',
    sampleTypeLevels = c('Nor (n=11)', 'Adj (n=10)', 'Can (n=10)'),
    takeLog = T)
cor.GSE95132.ribo.prot <-
  corAnalysis(expData = GSE95132.data$expData,
              sampleInfo = GSE95132.data$sampleInfo,
              dataSetName = 'GSE95132\nColorectal Cancer')

# GSE100179
GPL.GSE100179 <- read.csv('GSE100179/GPL17586-45144.csv')
GPL.GSE100179 <-
  cbind(
    GPL.GSE100179[,1],
    unname(sapply(GPL.GSE100179[,8],
                  function(x) str_remove_all(str_split(x, ' //')[[1]][5], ' '))))
GSE100179.data <-
  processGEOdataset(
    expDataFile = 'GSE100179/GSE100179_series_matrix.txt',
    expDataSkip = 65, GPLobject = GPL.GSE100179, geneIDtype = 'ENTREZID',
    sampleInfoFile = 'GSE100179/sampleInfo_GSE100179.csv',
    sampleTypeLevels = c('Nor (n=20)', 'Ade (n=20)', 'Can (n=20)'))
cor.GSE100179.ribo.prot <-
  corAnalysis(expData = GSE100179.data$expData,
              sampleInfo = GSE100179.data$sampleInfo,
              dataSetName = 'GSE100179\nColorectal Cancer')

# GSE147352
GSE147352.data <-
  t(fread('GSE147352/GSE147352_DESeq_normalized_counts.txt',
          data.table = F,header = T))
colnames(GSE147352.data) <- GSE147352.data[1,]
GSE147352.data <- GSE147352.data[-1,]
rownames.GSE147352 <- rownames(GSE147352.data)
GSE147352.data <- apply(GSE147352.data, 2, as.numeric)
rownames(GSE147352.data) <- rownames.GSE147352
GSE147352.data <- log2(GSE147352.data+1e-10)
sampleInfo.GSE147352 <- fread('GSE147352/GSE147352_sampleInfo.csv',
                              data.table = F)
sampleInfo.GSE147352$sampleType <-
  factor(sampleInfo.GSE147352$sampleType,
         levels=c('Nor (n=15)', 'GliL (n=18)', 'Glb (n=85)'))
sampleInfo.GSE147352 <-
  sampleInfo.GSE147352[match(rownames(GSE147352.data),
                             sampleInfo.GSE147352$sampleID),]
cor.GSE147352.ribo.prot <-
  corAnalysis(expData = GSE147352.data,
              sampleInfo = sampleInfo.GSE147352,
              dataSetName = 'GSE147352\nGlioblastomas')


CCs.validated <- matrix(ncol = 2, nrow = 0,
                        dimnames = list(NULL, c('Correlation Coefficients', 'group')))
for (i in ls()[str_detect(ls(), 'cor.GSE')]) {
  tmp <- get(i)
  dataSetName <- str_split(i, '\\.')[[1]][2]
  if (dataSetName %in% c('GSE27678', 'GSE91035', 'GSE10927', 'GSE64041')) {
    CCs.validated <- rbind(CCs.validated,
                           data.frame(group='Decrease', check.names = F,
                                      `Correlation Coefficients`=tmp[[1]]$CC,
                                      `Gene Set Name`='Decrease',stringsAsFactors = T))
  }
  else {
    CCs.validated <- rbind(CCs.validated,
                           data.frame(group='Increase', check.names = F,
                                      `Correlation Coefficients`=tmp[[1]]$CC,
                                      `Gene Set Name`='Increase',stringsAsFactors = T))
  }
}
boxplot.CC.validateData <- plotBox(CCs.validated,show.legend = F,title = 'Validate data sets')


# ---------- Part: protein synthesis inhibited data sets analysis ----------
setwd('D:/Projects/CRC_expData_LDA_Carcinogenesis/ribosome inhibited data/')
# GSE8597
GSE8597.data <- 
  processGEOdataset(
    expDataFile = 'GSE8597/GSE8597_series_matrix.txt',
    expDataSkip = 69, GPLFile = 'GSE8597/GPL570-55999.txt',
    GPLSkip = 16, GPLSelectedCol = c(1,12), geneIDtype = 'ENTREZID',
    sampleInfoFile = 'GSE8597/GSE8597_Info.csv',
    sampleTypeLevels = c('DMSO (n=8)', 'CHX 35.54mM,1h (n=8)'))
for (i in UPR.proliferation.apoptosis.geneSets) {
  assign(paste0('ribo.inhib.boxplot.GSE8597.', i),
         corAnalysis(expData = GSE8597.data$expData,
                     sampleInfo = GSE8597.data$sampleInfo,
                     dataSetName = 'GSE8597\nBreast Cancer(MCF7)',
                     geneSet1 = i, geneSet2 = NULL,
                     setName1 = i, setName2 = NULL,
                     comPair1 = 1, comPair2 = 2, as.title = F))
}

# GSE51068
GSE51068.data <- 
  processGEOdataset(
    expDataFile = 'GSE51068/GSE51068_series_matrix.txt',
    expDataSkip = 64, GPLFile = 'GSE51068/GPL13667-15572.txt',
    GPLSkip = 43, GPLSelectedCol = c(1,21), geneIDtype = 'ENTREZID',
    sampleInfoFile = 'GSE51068/GSE51068_Info.csv',
    sampleTypeLevels = c('DMSO 6h (n=8)', 'CHX 0.264μM,6h (n=3)', 'CHX 5μM,6h (n=3)',
                         'DMSO 12h (n=8)', 'CHX 0.264μM,12h (n=3)', 'CHX 5μM,12h (n=3)',
                         'DMSO 24h (n=8)', 'CHX 0.264μM,24h (n=3)', 'CHX 5μM,24h (n=3)'))
for (i in UPR.proliferation.apoptosis.geneSets) {
  assign(paste0('ribo.inhib.boxplot.GSE51068.', i),
         corAnalysis(expData = GSE51068.data$expData,
                     sampleInfo = GSE51068.data$sampleInfo,
                     dataSetName = 'GSE51068\nDiffuse Large B-Cell Lymphoma(OCI-Ly3)',
                     geneSet1 = i, geneSet2 = NULL,
                     setName1 = i, setName2 = NULL,
                     comPairs = list(c(1,2), c(1,3), c(2,3),
                                  c(4,5), c(4,6), c(5,6),
                                  c(7,8), c(7,9), c(8,9))))
}

# GSE124828
GSE124828.data <-
  processGEOdataset(
    expDataFile = 'GSE124828/GSE124828_series_matrix.txt',
    expDataSkip = 68, GPLFile = 'GSE124828/GPL23572-33493.txt',
    GPLSkip = 16, GPLSelectedCol = c(1,4), geneIDtype = 'ENTREZID', 
    sampleInfoFile = 'GSE124828/GSE124828_sampleInfo.csv',
    sampleTypeLevels = c('DMSO (n=11)', 'CHX 50μM (n=11)'))
for (i in UPR.proliferation.apoptosis.geneSets) {
  assign(paste0('ribo.inhib.boxplot.GSE124828.', i),
         corAnalysis(expData = GSE124828.data$expData,
                     sampleInfo = GSE124828.data$sampleInfo,
                     dataSetName = 'GSE124828\nCervical Cancer(HeLa)',
                     geneSet1 = i, geneSet2 = NULL,
                     setName1 = i, setName2 = NULL))
}

processGEOdataset2 <- function(expDataPath, sampleInfoPath,
                               levels, takeLog=F){
  expData <- fread(expDataPath, header = T, data.table = F)
  rownames(expData) <- expData[,1]
  expData <- expData[,-(1:6)]
  colnames(expData) <- sapply(str_split(colnames(expData), '_'), `[`, 1)
  expData <- t(expData)
  expData <- aggregate(expData, by=list(rownames(expData)), FUN=mean)
  rownames(expData) <- expData$Group.1
  expData$Group.1 <- NULL

  ensemblID <- bitr(rownames(expData),
                    fromType = 'SYMBOL',
                    toType = 'ENSEMBL',
                    OrgDb = org.Hs.eg.db)
  rownames <- data.frame(seqNo=1:nrow(expData),
                         SYMBOL=rownames(expData))
  ensemblID <- merge(rownames, ensemblID,
                     by='SYMBOL', all.x=T,
                     all.y=F, sort=F)
  ensemblID <- ensemblID[order(ensemblID$seqNo),]
  ensemblID <- ensemblID[!duplicated(ensemblID$seqNo),]
  expData <- expData[!(duplicated(ensemblID$ENSEMBL)|
                         is.na(ensemblID$ENSEMBL)),]
  ensemblID <- ensemblID[!(duplicated(ensemblID$ENSEMBL)|
                             is.na(ensemblID$ENSEMBL)),]
  rownames(expData) <- ensemblID$ENSEMBL
  library(preprocessCore)
  expData <- normalize.quantiles.robust(
    as.matrix(expData), remove.extreme = 'mean',
    use.log2 = takeLog, keep.names = T)
  expData <- data.frame(t(expData))
  
  sampleInfo <- fread(sampleInfoPath, header = T, data.table = F)
  commonSamples <- intersect(sampleInfo$sampleID,
                             rownames(expData))
  expData <- expData[rownames(expData) %in% commonSamples,]
  sampleInfo <- sampleInfo[sampleInfo$sampleID %in% commonSamples,]
  sampleInfo$sampleType <- factor(sampleInfo$sampleType,
                                  levels = levels)
  sampleInfo <- sampleInfo[order(sampleInfo$sampleType),]
  expData <- expData[match(sampleInfo$sampleID,
                           rownames(expData)),]
  expData <- expData[, sapply(apply(expData,2,unique), length)!=1]
  if (takeLog==T) {
    expData <- log(expData + 1e-10)
  }
  return(list(expData=expData, sampleInfo=sampleInfo))
}

# GSE162855
GSE162855.data <-
  processGEOdataset2(
    expDataPath = 'GSE162855/GSE162855_httr_mcf7_pilot_count_data.csv',
    sampleInfoPath = 'GSE162855/GSE162855_sampleInfo.csv',
    levels = c('DMSO (n=18)', 'CHX 0.03μM (n=3)',
               'CHX 0.1μM (n=3)', 'CHX 0.3μM (n=3)',
               'CHX 1μM (n=3)', 'CHX 3μM (n=3)',
               'CHX 10μM (n=3)', 'CHX 30μM (n=3)',
               'CHX 100μM (n=3)'), takeLog = T)
for (i in UPR.proliferation.apoptosis.geneSets) {
  assign(paste0('ribo.inhib.boxplot.GSE162855.', i),
         corAnalysis(expData = GSE162855.data$expData,
                     sampleInfo = GSE162855.data$sampleInfo,
                     dataSetName = 'GSE162855\nBreast Cancer(MCF7)',
                     geneSet1 = i, geneSet2 = NULL,
                     setName1 = i, setName2 = NULL,
                     comPairs = list(c(1,2), c(1,3), c(1,4), c(1,5),
                                     c(1,6), c(1,7), c(1,8), c(1,9))))
}

# GSE200845
GSE200845.data <-
  processGEOdataset2(
    expDataPath = 'GSE200845/GSE200845_httr_u2os_pilot_GEO_count_data.csv',
    sampleInfoPath = 'GSE200845/GSE200845_sampleInfo.csv',
    levels = c('DMSO (n=44)', 'CHX 0.01μΜ (n=4)', 'CHX 0.03μΜ (n=4)',
               'CHX 0.1μΜ (n=4)', 'CHX 0.3μΜ (n=4)', 'CHX 1μΜ (n=4)',
               'CHX 3μΜ (n=4)', 'CHX 10μΜ (n=4)'), takeLog = T)
for (i in UPR.proliferation.apoptosis.geneSets) {
  assign(paste0('ribo.inhib.boxplot.GSE200845.', i),
         corAnalysis(expData = GSE200845.data$expData,
                     sampleInfo = GSE200845.data$sampleInfo,
                     dataSetName = 'GSE200845\nOsteosarcoma(U-2 OS)',
                     geneSet1 = i, geneSet2 = NULL,
                     setName1 = i, setName2 = NULL,
                     comPairs = list(1:2, 2:3, 3:4, 4:5, 5:6, 6:7, 7:8)))
}


# ---------- Part: protein degradation inhibited data sets analysis ----------
setwd('D:/Projects/CRC_expData_LDA_Carcinogenesis/proteasome inhibited data/')
processGEOdataset3 <-
  function(expDataPath, sampleInfoPath,
           row.omitted=0, col.omitted=NULL,
           geneIDType, levels, takeLog=F){
  expData <- fread(expDataPath, header = T, data.table = F,
                   skip = row.omitted)
  if (!is.null(col.omitted)) {
    expData <- expData[,-col.omitted]
  }
  
  if (geneIDType!='ENSEMBL') {
    library(clusterProfiler)
    library(org.Hs.eg.db)
    tmp <- bitr(expData[,1],
                fromType = geneIDType,
                toType = 'ENSEMBL',
                OrgDb = org.Hs.eg.db)
    expData <- merge(tmp, expData, by.x=geneIDType,
                     by.y=colnames(expData)[1], all.y=T)
    expData <- expData[, -1]
  }
  expData$ENSEMBL <- sapply(str_split(expData$ENSEMBL,
                                      '\\.'), `[`, 1)
  expData <- aggregate(expData, FUN=mean,
                       by=list(expData$ENSEMBL))
  rownames(expData) <- expData$Group.1
  expData <- expData[,-(1:2)]
  library(preprocessCore)
  expData <- normalize.quantiles.robust(
    as.matrix(expData), remove.extreme = 'mean',
    use.log2 = takeLog, keep.names = T)
  expData <- t(expData)
  rownames(expData) <- str_replace(rownames(expData), '-', '_')

  sampleInfo <- fread(sampleInfoPath, header = T, data.table = F)
  commonSamples <- intersect(rownames(expData),
                             sampleInfo$sampleID)
  expData <- expData[rownames(expData) %in% commonSamples,]
  sampleInfo <- sampleInfo[sampleInfo$sampleID %in% commonSamples,]
  sampleInfo$sampleType <- factor(sampleInfo$sampleType,
                                  levels = levels)
  sampleInfo <- sampleInfo[order(sampleInfo$sampleType),]
  expData <- expData[match(sampleInfo$sampleID, rownames(expData)),]
  expData <- expData[, sapply(apply(expData,2,unique), length)!=1]
  if (takeLog==T) {
    expData <- log(expData + 1e-10)
  }
  return(list(expData=expData, sampleInfo=sampleInfo))
}

# GSE115973
GSE115973.data <-
  processGEOdataset3(
    expDataPath = 'GSE115973/GSE115973_NormalizedData_withAnnotations.xls',
    sampleInfoPath = 'GSE115973/GSE115973_sampleInfo.csv',
    col.omitted = c(1:3,5,6), geneIDType = 'ENTREZID',
    levels = c('DMSO 0h (n=4)', 'Bort 25nM,6h (n=4)', 'Bort 25nM,10h (n=4)'), takeLog=T)
for (i in UPR.proliferation.apoptosis.geneSets) {
  assign(paste0('prot.inhib.boxplot.GSE115973.', i),
         corAnalysis(expData = GSE115973.data$expData,
                     sampleInfo = GSE115973.data$sampleInfo,
                     dataSetName = 'GSE115973\nOsteosarcoma(U-2 OS)',
                     geneSet1 = i, geneSet2 = NULL,
                     setName1 = i, setName2 = NULL))
}

# GSE119841
GSE119841.data <-
  processGEOdataset3(
    expDataPath = 'GSE119841/GSE119841_L.Hoch_and_al_2018-normalized_gene_abundance_measurements_RPM.txt',
    sampleInfoPath = 'GSE119841/GSE119841_sampleInfo.csv',
    row.omitted = 1, geneIDType = 'SYMBOL',
    levels = c('DMSO (n=3)', 'Bort 30nM (n=3)', 'MG132 1mM (n=3)'), takeLog = T)
for (i in UPR.proliferation.apoptosis.geneSets) {
  assign(paste0('prot.inhib.boxplot.GSE119841.', i),
         corAnalysis(expData = GSE119841.data$expData,
                     sampleInfo = GSE119841.data$sampleInfo,
                     dataSetName = 'GSE119841\nFibroblast(R77C)',
                     geneSet1 = i, geneSet2 = NULL,
                     setName1 = i, setName2 = NULL,
                     comPairs = list(c(1,2), c(1,3))))
}

# GSE141858
GSE141858.data <-
  processGEOdataset(
    expDataFile = 'GSE141858/GSE141858_series_matrix.txt',
    expDataSkip = 59, GPLFile = 'GSE141858/GPL4133-12599.txt',
    GPLSkip = 22, GPLSelectedCol = c(1,9),
    sampleInfoFile = 'GSE141858/GSE141858_Info.csv',
    sampleTypeLevels = c('DMSO 0h (n=3)', 'MG132 1μM,4h (n=3)',
                         'MG132 1μM,24h (n=3)'),
    geneIDtype = 'ENTREZID', takeLog = T)
for (i in UPR.proliferation.apoptosis.geneSets) {
  assign(paste0('prot.inhib.boxplot.GSE141858.', i),
         corAnalysis(expData = GSE141858.data$expData,
                     sampleInfo = GSE141858.data$sampleInfo,
                     dataSetName = 'GSE141858\nBreast Cancer(MCF7)',
                     geneSet1 = i, geneSet2 = NULL,
                     setName1 = i, setName2 = NULL))
}

# GSE164788
GSE164788.data <-
  fread('GSE164788/GSE164788_normalized_counts.csv',
        header=T, data.table = F)
library(reshape2)
GSE164788.data <- dcast(data = GSE164788.data,
                        formula = gene_name~sample)
write.csv(GSE164788.data, row.names = F,
          file = 'GSE164788/GSE164788.data.csv')

GSE164788.data <-
  processGEOdataset3(expDataPath = 'GSE164788/GSE164788.data.csv',
                     sampleInfoPath = 'GSE164788/GSE164788_Info.csv',
                     geneIDType = 'SYMBOL',
                     levels = c('DMSO (n=23)', 'Bort 0.3μM (n=3)',
                                'Bort 1μM (n=3)', 'Bort 3μM (n=3)',
                                'Bort 10μM (n=3)', 'MG132 1μM (n=3)',
                                'MG132 10μM (n=3)'), takeLog = T)
for (i in UPR.proliferation.apoptosis.geneSets) {
  assign(paste0('prot.inhib.boxplot.GSE164788.', i),
         corAnalysis(expData = GSE164788.data$expData,
                     sampleInfo = GSE164788.data$sampleInfo,
                     dataSetName = 'GSE164788\nNeural and Glial Cells(ReNcell VM)',
                     geneSet1 = i, geneSet2 = NULL, setName1 = i, setName2 = NULL,
                     comPairs = list(c(1,2), c(1,3), c(1,4), c(1,5), c(1,6), c(1,7))))
}

# GSE165325
GSE165325.data <-
  processGEOdataset3(
    expDataPath = 'GSE165325/GSE165325_gene_te_counts.tsv',
    sampleInfoPath = 'GSE165325/GSE165325_Info.csv',
    geneIDType = 'ENSEMBL', takeLog = T,
    levels = c('Ctrl (n=4)','Bort 100nM,3h (n=4)','MG132 20μM,3h (n=4)'))
for (i in UPR.proliferation.apoptosis.geneSets) {
  assign(paste0('prot.inhib.boxplot.GSE165325.', i),
         corAnalysis(expData = GSE165325.data$expData,
                     sampleInfo = GSE165325.data$sampleInfo,
                     dataSetName = 'GSE165325\nColorectal Cancer(HCT-116)',
                     geneSet1 = i, geneSet2 = NULL,
                     setName1 = i, setName2 = NULL,
                     comPairs = list(c(1,2), c(1,3))))
}

# GSE184029
GSE184029.data <-
  processGEOdataset3(expDataPath = 'GSE184029/GSE184029_Raw_gene_counts_matrix.txt',
                     sampleInfoPath = 'GSE184029/GSE184029_sampleInfo.csv',
                     col.omitted = 1, geneIDType = 'ENSEMBL',
                     levels = c('H-929, DMSO (n=3)', 'H-929, Cart (n=3)',
                                'RPMI, DMSO (n=3)', 'RPMI, Cart (n=3)',
                                'DLD-1, DMSO (n=3)', 'DLD-1, Cart (n=3)',
                                'RKO, DMSO (n=3)', 'RKO, Cart (n=3)',
                                'CAPAN-2, DMSO (n=3)', 'CAPAN-2, Cart (n=3)',
                                'PANC-1, DMSO (n=3)', 'PANC-1, Cart (n=3)',
                                'H-1299, DMSO (n=3)', 'H-1299, Cart (n=3)',
                                'H-23, DMSO (n=3)', 'H-23, Cart (n=3)'), takeLog = T)
for (i in UPR.proliferation.apoptosis.geneSets) {
  assign(paste0('prot.inhib.boxplot.GSE184029.', i),
         corAnalysis(expData = GSE184029.data$expData,
                     sampleInfo = GSE184029.data$sampleInfo,
                     dataSetName = 'GSE184029\nMulti Cancers(Multi Cell Lines)',
                     geneSet1 = i, geneSet2 = NULL, setName1 = i, setName2 = NULL,
                     comPairs = list(1:2, 3:4, 5:6, 7:8, 9:10, 11:12, 13:14, 15:16)))
}
setwd('D:/Projects/CRC_expData_LDA_Carcinogenesis/')


# ---------- Part XII: Plot figures ----------
setwd('D:/Projects/CRC_expData_LDA_Carcinogenesis/figures and tables 20221214/')
library(cowplot)

# Figure S1: PCA plots of all cancer types in XENA.
figureS1 <- list()
for (i in ls()[str_detect(ls(), 'PCAplot')&str_detect(ls(),'allSamples')]) {
  figureS1[[length(figureS1)+1]] <- get(i)
}
png('Figure S1 PCAplots XENA all samples.png',
    width = 16000, height = 16000, res = 600)
plot_grid(plotlist = figureS1,
          nrow=4, labels=LETTERS[1:14], label_size = 16)
dev.off()

# Figure S2: PCA plots of all cancer types in XENA--outliers removed.
figureS2 <- list()
for (i in ls()[str_detect(ls(), 'PCAplot')&str_detect(ls(),'outliers.removed')]) {
  figureS2[[length(figureS2)+1]] <- get(i)
}
png('Figure S2 PCAplots XENA outliers removed.png',
    width = 16000, height = 16000, res = 600)
plot_grid(plotlist = figureS2,
          nrow=4, labels=LETTERS[1:14], label_size = 16)
dev.off()

# Figure 1: LDA and ROC plots.
png('Figure 1 LDAplots ROCplots XENA CRC.png',
    width = 12000, height = 9000, res = 600)
plot_grid(LDAplot.XENA.CRC, ROCplot.XENA.CRC[[1]],
          ROCplot.XENA.CRC[[2]], ROCplot.XENA.CRC[[3]],
          LDAplot.RFE.XENA.CRC, ROCplot.RFE.XENA.CRC[[1]],
          ROCplot.RFE.XENA.CRC[[2]], ROCplot.RFE.XENA.CRC[[3]],
          LDAplot.iRFE.XENA.CRC, ROCplot.iRFE.XENA.CRC[[1]],
          ROCplot.iRFE.XENA.CRC[[2]], ROCplot.iRFE.XENA.CRC[[3]],
          nrow = 3, labels = LETTERS[1:12], label_size = 16)
dev.off()

# Figure S3-S15: ROC plots of the other 13 cancer types in XENA.
figureIndex <- 3
library(ggplotify)
for (i in str_sub(ls()[str_detect(ls(), 'ribo.prot.cor.hm.')], 18)) {
  if (i != 'Colon Adenocarcinoma') {
    png(paste0('Figure S', figureIndex, ' ', i, ' ROC plots.png'),
        width = 12000, height = 6000, res = 600)
    top.row <- plot_grid(get(paste0('LDAplot.RFE.XENA.', i)),
                         as.grob(get(paste0('ROCplot.RFE.XENA.', i))[[1]]),
                         as.grob(get(paste0('ROCplot.RFE.XENA.', i))[[2]]),
                         as.grob(get(paste0('ROCplot.RFE.XENA.', i))[[3]]),
                         nrow = 1, labels = LETTERS[1:4], label_size = 16)
    bottom.row <- plot_grid(get(paste0('LDAplot.iRFE.XENA.', i)),
                            get(paste0('ROCplot.iRFE.XENA.', i))[[1]],
                            get(paste0('ROCplot.iRFE.XENA.', i))[[2]],
                            get(paste0('ROCplot.iRFE.XENA.', i))[[3]],
                            nrow = 1, labels = LETTERS[5:8], label_size = 16)
    print(plot_grid(top.row, bottom.row, nrow = 2,
                    labels = NULL))
    dev.off()
    figureIndex <- figureIndex + 1
  }
}

# Figure 2: Enrichment of RFE genes.
png('Figure 2 enrichment plots.png',
    width=14000, height=7000, res=600)
plot_grid(enrichplot.RFE.XENA.CRC[[1]],
          enrichplot.RFE.XENA.CRC[[2]],
          enrichplot.RFE.XENA.CRC[[3]],
          enrichplot.RFE.XENA.CRC[[4]],
          nrow = 1, labels = LETTERS[1:4],
          label_size = 16)
dev.off()

# Figure 3: The correlation between ribosome and proteasome -- CRC in XENA.
png('Figure 3 correlation plots XENA CRC.png',
    width = 10000, height = 12000, res = 600)
top.row <- plot_grid(as.grob(hm.exp.ribo.prot.XENA.CRC),
                     as.grob(hm.CC.ribo.prot.XENA.CRC),
                     nrow = 1, labels = LETTERS[1:2], label_size = 16)
bottom.row <- plot_grid(boxplot.CC.ribo.random.XENA.CRC,
                        boxplot.CC.prot.random.XENA.CRC,
                        nrow = 2, labels = c('C', 'D'), label_size = 16)
plot_grid(top.row, bottom.row, nrow = 2,
          labels = NULL, rel_heights = c(1,1.5))
dev.off()

# Figure S16 ~ figure S28: The correlation between ribosome and
# proteasome -- other 13 cancer types in XENA.
for (i in str_sub(ls()[str_detect(ls(), 'ribo.prot.cor.hm.')], 18)) {
  if (i != 'Colon Adenocarcinoma') {
    png(paste0('Figure S', figureIndex, ' ', i, ' correlation plots.png'),
        width = 12000, height = 12000, res = 600)
    top.row <- plot_grid(as.grob(get(paste0('ribo.prot.cor.hm.', i))[[3]]),
                         as.grob(get(paste0('ribo.prot.cor.hm.', i))[[2]]),
                         nrow = 1, labels = LETTERS[1:2], label_size = 16)
    bottom.row <- plot_grid(get(paste0('boxplot.CC.ribo.', i)),
                            get(paste0('boxplot.CC.prot.', i)),
                            nrow = 2, labels = c('C', 'D'), label_size = 16)
    print(plot_grid(top.row, bottom.row, nrow = 2,
                    labels = NULL, rel_heights = c(1,1.5)))
    dev.off()
    figureIndex <- figureIndex + 1
  }
}

# Figure S29 ~ figure S30: PCA of ribosome and proteasome genes -- XENA.
PCAplots.ribo <- list()
PCAplots.prot <- list()
for (i in str_sub(ls()[str_detect(ls(), 'ribo.prot.cor.hm.')], 18)) {
  PCAplots.ribo[[length(PCAplots.ribo)+1]] <-
    get(paste0('PCAscores.ribo.', i))
  PCAplots.prot[[length(PCAplots.prot)+1]] <-
    get(paste0('PCAscores.prot.', i))
}

png(paste0('Figure S', figureIndex , ' PCAplots of ribo gene sets.png'),
    width = 12000, height = 12000, res = 600)
plot_grid(plotlist = lapply(PCAplots.ribo, `[[`, 3), ncol = 4,
          label_size = 16, labels = LETTERS[1:length(PCAplots.ribo)])
dev.off()
figureIndex <- figureIndex + 1

png(paste0('Figure S', figureIndex , ' PCAplots of prot gene sets.png'),
    width = 12000, height = 12000, res = 600)
plot_grid(plotlist = lapply(PCAplots.prot, `[[`, 3), ncol = 4,
          label_size = 16, labels = LETTERS[1:length(PCAplots.prot)])
dev.off()
figureIndex <- figureIndex + 1

# Figure 4: The correlation between PC1s of ribosome and proteasome genes -- XENA.
PCAscores.corPlot <- list()
for (i in str_sub(ls()[str_detect(ls(), 'ribo.prot.cor.hm.')], 18)) {
  PCAscores.corPlot[[length(PCAscores.corPlot)+1]] <- get(paste0('corPlot.ribo.prot.', i))
}
png('Figure 4 PCAscores corplot of ribo and prot.png',
    width = 12000, height = 9000, res = 600)
plot_grid(plotlist = PCAscores.corPlot, nrow = 4, label_size = 16,
          labels = LETTERS[1:length(PCAscores.corPlot)])
dev.off()

# Figure 5: The box plots of the difference between the PC1s of
# ribosome and proteasome genes -- XENA.
boxplots.diff.list <- list()
for (i in str_sub(ls()[str_detect(ls(), 'ribo.prot.cor.hm.')], 18)) {
  boxplots.diff.list[[length(boxplots.diff.list)+1]] <-
    get(paste0('boxplot.diff.ribo.prot.', i))
}
png('Figure 5 PCAscores diff between ribo and prot.png',
    width = 8000, height = 8500, res = 600)
plot_grid(plotlist = boxplots.diff.list, ncol=4, label_size = 16,
          labels=LETTERS[1:length(boxplots.diff.list)])
dev.off()

# Figure S31: The box plots of the difference between the PC1s of
# ribosome and proteasome genes -- GEO data sets.
plotList <- list()
for (i in ls()[str_detect(ls(), 'cor.GSE')]) {
  tmp <- get(i)
  i <- str_split(i,'\\.')[[1]][2]
  plotList[[length(plotList)+1]] <- as.grob(tmp[[1]][[2]])
  plotList[[length(plotList)+1]] <- tmp[[3]]
}
plotList <- plotList[c(1:2, 5:8,11:14,19:20,15:16,3:4,9:10,17:18)]
plotList <- plotList[c(seq(1,9,2), seq(2,10,2), seq(11,19,2), seq(12,20,2))]
png('Figure S31 validate data.png',
    width = 12000, height = 9000, res = 600)
plot_grid(plotlist = plotList, ncol = 5, label_size = 16,
          labels = LETTERS[1:length(plotList)], byrow = T)
dev.off()
figureIndex <- figureIndex + 1


# Figure S32: Correlation coefficients between ribosome and proteasome genes
# in validate data sets.
png('Figure S32 Correlation coefficients of validate data.png',
    width = 1500, height = 2500, res = 600)
plot_grid(boxplot.CC.validateData)
dev.off()
figureIndex <- figureIndex + 1

# Figure 6: Mechanism exploration.
mechanism.list <- list()
for (j in c('Pr.Cata', 'Pr.Fold', 'ER')) {
  for (k in c('GOCC_RIBOSOME', 'GOCC_PROTEASOME_COMPLEX')) {
    mechanism.list[[length(mechanism.list)+1]] <- 
      get(paste('boxplot', str_split(k, '_')[[1]][2], j,
                'Colon Adenocarcinoma', sep = '.'))
  }
}
png('Figure 6 Mechanism exploration.png',
    width = 12000, height = 12000, res = 600)
print(plot_grid(plotlist = mechanism.list, nrow = 3,
                labels = LETTERS[1:6], label_size = 16))
dev.off()


# Figure S33 -- figure S45: Mechanism exploration.
for (i in ls()[str_detect(ls(), '^exp.pht.+')]) {
  if (i!='exp.pht.Colon Adenocarcinoma') {
    mechanism.list <- list()
    for (j in c('Pr.Cata', 'Pr.Fold', 'ER')) {
      for (k in c('GOCC_RIBOSOME', 'GOCC_PROTEASOME_COMPLEX')) {
        mechanism.list[[length(mechanism.list)+1]] <- 
          get(paste('boxplot', str_split(k, '_')[[1]][2], j,
                    str_split(i, '\\.')[[1]][3], sep = '.'))
      }
    }
    png(paste0('Figure S',figureIndex, ' ', i, '.png'),
        width = 12000, height = 12000, res = 600)
    print(plot_grid(plotlist = mechanism.list, nrow = 3,
                    labels = LETTERS[1:6], label_size = 16))
    dev.off()
    figureIndex <- figureIndex + 1
  }
}


# Figure 7:  The box plots of the PC1s of some gene sets
# after ribosome or proteasome inhibition.
png('Figure 7 PCAscores of UPR pathways after protein degradation or synthesis inhibited.png',
    width = 15000, height = 18000, res = 600)
PCAscores.UPR.pathways <- list()
for (i in c('GOBP_IRE1', 'GOBP_ATF6', 'GOBP_PERK')) {
  for (j in c('GSE115973', 'GSE119841', 'GSE141858', 'GSE164788')) {
      PCAscores.UPR.pathways[[length(PCAscores.UPR.pathways)+1]] <-
        get(paste0('prot.inhib.boxplot.', j, '.', i,
                   '_MEDIATED_UNFOLDED_PROTEIN_RESPONSE'))
  }
}
for (i in c('GOBP_IRE1', 'GOBP_ATF6', 'GOBP_PERK')) {
  for (j in c('GSE124828', 'GSE8597', 'GSE162855', 'GSE51068')) {
    PCAscores.UPR.pathways[[length(PCAscores.UPR.pathways)+1]] <-
      get(paste0('ribo.inhib.boxplot.', j, '.', i,
                 '_MEDIATED_UNFOLDED_PROTEIN_RESPONSE'))
  }
}
names(PCAscores.UPR.pathways) <- 1:24
narrow_cols <- plot_grid(plotlist = PCAscores.UPR.pathways[c(1,2,5,6,9,10,13,14,17,18,21,22)],
                         ncol = 2, label_size = 16,
                         labels = c('A', 'B', 'E', 'F', 'I', 'J', 'M', 'N', 'Q', 'R', 'U', 'V'))
wide_cols <- plot_grid(plotlist = PCAscores.UPR.pathways[c(3,4,7,8,11,12,15,16,19,20,23,24)],
                       ncol = 2, label_size = 16,
                       labels = c('C','D','G','H','K','L','O','P','S','T','W','X'))

plot_grid(narrow_cols, wide_cols, ncol = 2, label_size = 16, rel_widths = c(3,5))
dev.off()

# Figure 8: the apoptosis levels after inhibition of ribosome and proteasome.
png('Figure 8 the apoptosis levels after inhibition of ribosome and proteasome.png',
    width = 12000, height = 6000, res = 600)
narrow_cols <-
  plot_grid(prot.inhib.boxplot.GSE115973.KEGG_APOPTOSIS,
            prot.inhib.boxplot.GSE119841.KEGG_APOPTOSIS,
            
            ribo.inhib.boxplot.GSE8597.KEGG_APOPTOSIS,
            ribo.inhib.boxplot.GSE124828.KEGG_APOPTOSIS,
            nrow = 2, label_size = 16,
            labels = c('A', 'B', 'E', 'F'))
wide_cols <-
  plot_grid(prot.inhib.boxplot.GSE141858.KEGG_APOPTOSIS,
            prot.inhib.boxplot.GSE164788.KEGG_APOPTOSIS,
            
            ribo.inhib.boxplot.GSE51068.KEGG_APOPTOSIS,
            ribo.inhib.boxplot.GSE162855.KEGG_APOPTOSIS,
            nrow = 2, label_size = 16,
            labels = c('C', 'D', 'G', 'H'))
plot_grid(narrow_cols, wide_cols, ncol = 2, rel_widths = c(2,5))
dev.off()

# Figure S46 -- S47

# Figure S46
png('Figure S46 epithelial cell proliferation regulation and apoptosis after proteasome inhibition.png',
    width = 10000, height = 6000, res = 600)
narrow_cols <- plot_grid(prot.inhib.boxplot.GSE115973.GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION,
                         prot.inhib.boxplot.GSE115973.GOBP_NEGATIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION,

                         prot.inhib.boxplot.GSE119841.GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION,
                         prot.inhib.boxplot.GSE119841.GOBP_NEGATIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION,

                         prot.inhib.boxplot.GSE141858.GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION,
                         prot.inhib.boxplot.GSE141858.GOBP_NEGATIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION,

                         ncol = 3, label_size = 16, labels = LETTERS[c(1:3,5:7,9:11)], byrow = F)

wide_cols <- plot_grid(prot.inhib.boxplot.GSE164788.GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION,
                       prot.inhib.boxplot.GSE164788.GOBP_NEGATIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION,
                       ncol = 1, label_size = 16, labels = c('D', 'H', 'L'), byrow = F)
plot_grid(narrow_cols, wide_cols, ncol = 2, rel_widths = c(4,3))
dev.off()


# Figure S47
png('Figure S47 epithelial cell proliferation regulation and apoptosis after ribosome inhibition.png',
    width = 13000, height = 6000, res = 600)
narrow_cols <- plot_grid(
  ribo.inhib.boxplot.GSE8597.GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION,
  ribo.inhib.boxplot.GSE8597.GOBP_NEGATIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION,

  ribo.inhib.boxplot.GSE124828.GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION,
  ribo.inhib.boxplot.GSE124828.GOBP_NEGATIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION,

  ncol = 2, label_size = 16, labels = c('A', 'B', 'E', 'F', 'I', 'J'), byrow = F)
          
wide_cols <- plot_grid(
  ribo.inhib.boxplot.GSE51068.GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION,
  ribo.inhib.boxplot.GSE51068.GOBP_NEGATIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION,

  ribo.inhib.boxplot.GSE162855.GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION,
  ribo.inhib.boxplot.GSE162855.GOBP_NEGATIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION,

  ncol = 2, label_size = 16, labels = c('C','D', 'G', 'H', 'K', 'L'), byrow = F)
plot_grid(narrow_cols, wide_cols, ncol = 2, rel_widths = c(1,3))
dev.off()
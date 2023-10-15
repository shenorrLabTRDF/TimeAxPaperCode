#folder = "Add the source code folder here"
folder = "C:/Users/Shennor/Dropbox/research_documents/Postdoc/Shai/multipleAlignment/Figures/New/Revision/V2/sourceCode"
source(file.path(folder, "enrichmentTest.R"))
library(TimeAx)
library(ggplot2)

dataPlotter = function(currMetaData, currTraj, currMeta, additionalCol = NULL, box = T, yValues = NULL, xValues = NULL, regLine = F, DoF = 4){
  x = currMetaData[,currMeta]
  indsToUse = which(!is.na(x) & !is.na(currTraj))
  metaDataLevels = x[indsToUse]
  currTrajLevels = currTraj[indsToUse]
  currCorValue = cor(metaDataLevels, currTrajLevels)
  if(length(unique(metaDataLevels))<6 & box){
    if(length(unique(metaDataLevels)) == 2 & min(table(metaDataLevels))>3){
      currCorValue = t.test(currTrajLevels[metaDataLevels == unique(metaDataLevels)[1]],
                            currTrajLevels[metaDataLevels == unique(metaDataLevels)[2]])$p.value
    }
    ggplot(data = NULL, aes(x = factor(metaDataLevels), y = currTrajLevels))+geom_boxplot()+ggtitle(paste(currMeta, currCorValue,sep = " _ "))+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  }else{
    if(is.null(xValues)){
      xValues = c(min(currTrajLevels),max(currTrajLevels))
    }
    if(is.null(yValues)){
      yValues = c(min(metaDataLevels),max(metaDataLevels))
    }
    regLineFunc = NULL
    if(regLine){
      regLineFunc = geom_smooth(method = "lm", formula = y ~ poly(x, DoF), se = F)
    }
    if(!is.null(additionalCol)){
      additionalCol = factor(additionalCol[indsToUse])
      ggplot(data = NULL, aes(x = currTrajLevels, y = metaDataLevels))+geom_point(size = 3, aes(color = additionalCol))+ggtitle(paste(currMeta, currCorValue,sep = " _ "))+xlim(xValues)+ylim(yValues)+regLineFunc+
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    }else{
      ggplot(data = NULL, aes(x = currTrajLevels, y = metaDataLevels))+geom_point(size = 3)+ggtitle(paste(currMeta, currCorValue,sep = " _ "))+xlim(xValues)+ylim(yValues)+regLineFunc+
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    }
  }
}
proteinCodingGenes = as.character(as.matrix(read.table(file.path(folder,"proteinCodingGenes.txt"),sep = "\t", row.names = NULL, col.names = F)))

#### Load and norm data ####

# Main cohort #
UBCFullData = readRDS(file.path(folder,"UBCFullDataRaw.rds"))
annotations = UBCFullData$annotations
GEData = UBCFullData$GEData
annotationDataFull = read.table(file.path(folder, "annotationTableFull.txt"),sep = "\t", row.names = 1, header = T, check.names = F)[colnames(GEData),]
annotationDataFull = annotationDataFull[colnames(GEData),]
sampleNames = annotations$Patient_ID
trajectory = annotations$Tumor_Number
batchInfo = read.table(file.path(folder,"batchInformation.txt"),sep="\t", row.names = 1, header = T, check.names = F)
batchInfo = batchInfo[colnames(GEData),]
GEDataCombat = sva::ComBat(GEData, batchInfo)

trajectoryFixed = unlist(lapply(unique(sampleNames),function(currName){ 1:length(which(currName == sampleNames))}))
trajectoryTime = as.numeric(annotations$Months_from_Primary)/12

GEDataToUse = GEDataCombat

trainSamples = names(which(table(sampleNames)>3))
testSamples = unique(sampleNames)[!(unique(sampleNames) %in% trainSamples)]
trainIndexes = which(sampleNames %in% trainSamples)
testIndexes = which(sampleNames %in% testSamples)
sampleNamesTrain = sampleNames[trainIndexes]
sampleNamesTest = sampleNames[testSamples]
GEDataTrain = GEDataToUse[,trainIndexes]
GEDataTest = GEDataToUse[,testSamples]
trajectoryFixedTrain = trajectoryFixed[trainIndexes]
trajectoryFixedTest = trajectoryFixed[testSamples]
trajectoryTimeTrain = trajectoryTime[trainIndexes]
trajectoryTimeTest = trajectoryTime[testSamples]
annotationsTrain = annotations[trainIndexes,]
annotationsTest = annotations[testSamples,]

metaDataNumTrain = apply(annotationsTrain,2,function(x){
  xFixed = as.numeric(x)
  if(length(which(is.na(xFixed)))>30){
    xFixed = as.numeric(as.factor(x))
  }
  xFixed[x=="N/A"] = NA
  xFixed
})
metaDataNumTest = apply(annotationsTest,2,function(x){
  xFixed = as.numeric(x)
  if(length(which(is.na(xFixed)))>30){
    xFixed = as.numeric(as.factor(x))
  }
  xFixed[x=="N/A"] = NA
  xFixed
})

# Microarray cohort #
dataGSE83586Full = readRDS(file.path(folder,"UBCFullDataRawGSE83586.rds"))
dataGSE83586 = dataGSE83586Full$GEData
annotationsGSE83586 = dataGSE83586Full$annotations
batchInfoGSE83586 = read.table(file.path(folder,"GSE83586_batches.txt"),sep="\t", row.names = 1, header = T, check.names = F)
batchInfoGSE83586 = batchInfoGSE83586[colnames(dataGSE83586),]
dataGSE83586 = dataGSE83586[intersect(row.names(dataGSE83586),row.names(GEDataToUse)),]

dataGSE83586 = sva::ComBat(dataGSE83586, batchInfoGSE83586)

# TCGA cohort #
dataFullTCGA = readRDS(file.path(folder,"TCGAData.rds"))
dataTCGA = dataFullTCGA$data
colnames(dataTCGA) = sapply(colnames(dataTCGA),function(x){unlist(strsplit(x,"\\."))[1]})
dataTCGA = dataTCGA[which(!is.na(rowSums(dataTCGA))),]
annotationsTCGA = dataFullTCGA$clinical
dataTCGANorm = as.matrix(log(dataTCGA+1))
dataTCGANorm = dataTCGANorm[which(!is.nan(rowSums(dataTCGANorm))),]
indexesTCGAOneSample = match(names(which(table(annotationsTCGA$case_submitter_id)==1)),annotationsTCGA$case_submitter_id)
dataTCGANorm = dataTCGANorm[,indexesTCGAOneSample]
annotationsTCGA = annotationsTCGA[indexesTCGAOneSample,]

#### Calculate model and pseudotime ####
proteinCodingGenes = intersect(proteinCodingGenes, row.names(GEDataTrain))
multiAlignment = modelCreation(GEDataTrain[proteinCodingGenes,], sampleNamesTrain, numOfTopFeatures = 100)
seedGenes = multiAlignment$seed
metPoint = 0.7

trajToUseProfileTrain = predictByConsensus(multiAlignment,GEDataTrain)$predictions
trajToUseProfile = predictByConsensus(multiAlignment,GEDataTest)$predictions
annotationsToUseProfile = annotationsTest
trajToUseProfileGSE83586 = predictByConsensus(multiAlignment,dataGSE83586)$predictions
trajToUseProfileTCGA = predictByConsensus(multiAlignment,dataTCGANorm)$predictions


##### Figure 2 #####
#### 2E #####
trajectoryTimeTrainForVis = trajectoryTimeTrain
trajectoryTimeTrainForVis[which.max(trajectoryTimeTrainForVis)] = 8
chosenPatients = c("20","26","30","40","57","65","67","72")
tmpData = as.data.frame(cbind(trajectoryTimeTrainForVis,trajToUseProfileTrain,sampleNamesTrain))
colnames(tmpData) = c("Pre","Post","Sample")
tmpData = tmpData[tmpData$Sample %in% chosenPatients,]
ggplot(data = tmpData, aes(x = Pre,y = Post, colour = factor(Sample)))+geom_line(size = 1)+ylim(0,1)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### 2F ####
regPre = summary(lm(t(GEDataTrain)~poly(trajectoryTimeTrain,2)))
genePValuesPre = sapply(regPre, function(x){
  fstat <- x$fstatistic
  pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
  #x[[4]][2,4]
})

regAligned = summary(lm(t(GEDataTrain)~poly(trajToUseProfileTrain,2)))
genePValuesAligned = sapply(regAligned, function(x){
  fstat <- x$fstatistic
  pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
  #x[[4]][2,4]
})

geneRegPreLog = -log(genePValuesPre,base=10)
geneRegAlignedLog = -log(genePValuesAligned,base=10)
newCutOff = -log(0.01, base=10)
geneColourProVsAligned = rep(0, length(geneRegPreLog))
geneColourProVsAligned[geneRegPreLog>newCutOff & geneRegAlignedLog<newCutOff] = 1
geneColourProVsAligned[geneRegPreLog<newCutOff & geneRegAlignedLog>newCutOff] = 2
geneColourProVsAligned[geneRegPreLog>newCutOff & geneRegAlignedLog>newCutOff] = 3
ggplot(data = NULL, aes(x = geneRegPreLog, y = geneRegAlignedLog, colour = factor(geneColourProVsAligned)))+geom_point()+
  geom_hline(yintercept = newCutOff, linetype = "dashed")+geom_vline(xintercept = newCutOff, linetype = "dashed")+
  scale_color_manual(values = c("0" = "lightgray","1" = "darkgray","2" = "purple", "3" = "orange"))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### S3A ####
factorToPlot = "CCL2"
#factorToPlot = "IFITM2"
#factorToPlot = "SGPL1"
dataPlotter(t(GEDataTrain), trajToUseProfileTrain,factorToPlot,regLine = T, DoF=3)
dataPlotter(t(GEDataTrain), trajectoryTimeTrain,factorToPlot,regLine = T, DoF=3)

#### S3B ####
geneSets = readRDS(file.path(folder, "geneSets.rds"))
corGeneTrajTime = cor(t(GEDataTrain), trajectoryTimeTrain)[,1]
corGeneTrajPseudo = cor(t(GEDataTrain), trajToUseProfileTrain)[,1]
enrichmentTime = runEnrichment(sort(corGeneTrajTime),geneSets)
enrichmentPseudo = runEnrichment(sort(corGeneTrajPseudo),geneSets)
minPoint = min(enrichmentTime,enrichmentPseudo)
maxPoint = max(enrichmentTime,enrichmentPseudo)
ggplot(data = NULL, aes(x = enrichmentTime, y = enrichmentPseudo))+geom_point()+
  xlim(minPoint,maxPoint)+ylim(minPoint,maxPoint)+geom_hline(yintercept = 0, linetype = "dashed")+geom_vline(xintercept = 0, linetype = "dashed")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

##### Figure 3 #####
#### 3A ####
stagesFull = annotations$Stage
stagesFull[stagesFull == "N/A" | stagesFull == "cis" | stagesFull == "Met" | stagesFull == "T1M1"] = "Unknown"
stagesFull[stagesFull == "Ta"] = "T0"
stagesFull[grep("T3",stagesFull)] = "T3"
stages = stagesFull[trainIndexes]
stageCont = as.numeric(gsub("T","",stages[stages!="Unknown"]))
tStage = summary(lm(stageCont~trajToUseProfileTrain[stages!="Unknown"]))[[4]][2,4]
ggplot(data = NULL, aes(x = factor(stages[stages!="Unknown"]), y = trajToUseProfileTrain[stages!="Unknown"]))+geom_boxplot()+ggtitle(paste("Main",tStage,sep=" _ "))+geom_jitter(height = 0, width = 0.1)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### S5A ####
tTime = summary(lm(stageCont~trajectoryTimeTrain[stages!="Unknown"]))[[4]][2,4]
ggplot(data = NULL, aes(x = factor(stages[stages!="Unknown"]), y = trajectoryTimeTrain[stages!="Unknown"]))+geom_boxplot()+ggtitle(tTime)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### S5B ####
stagesTest = annotationsToUseProfile$Stage
stagesTest[stagesTest == "N/A" | stagesTest == "cis" | stagesTest == "Met" | stagesTest == "T1M1"] = "Unknown"
stagesTest[stagesTest == "Ta"] = "T0"
stagesTest[grep("T3",stagesTest)] = "T3"
stageContTest = as.numeric(gsub("T","",stagesTest[stagesTest!="Unknown"]))
tStageTest = summary(lm(stageContTest~trajToUseProfile[stagesTest!="Unknown"]))[[4]][2,4]
ggplot(data = NULL, aes(x = factor(stagesTest[stagesTest!="Unknown"]), y = trajToUseProfile[stagesTest!="Unknown"]))+geom_boxplot()+ggtitle(paste("Main",tStageTest,sep=" _ "))+geom_jitter(height = 0, width = 0.1)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

stagesTCGA = annotationsTCGA$ajcc_pathologic_t
stagesTCGA[stagesTCGA == "--"] = "Unknown"
stagesTCGA[stagesTCGA == ""] = "Unknown"
stagesTCGA[stagesTCGA == "TX"] = "Unknown"
stagesTCGA[grep("T1",stagesTCGA)] = "T1"
stagesTCGA[grep("T2",stagesTCGA)] = "T2"
stagesTCGA[grep("T3",stagesTCGA)] = "T3"
stagesTCGA[grep("T4",stagesTCGA)] = "T4"
stageTCGACont = as.numeric(gsub("T","",stagesTCGA[stagesTCGA!="Unknown"]))
tStageTCGA = summary(lm(stageTCGACont~trajToUseProfileTCGA[stagesTCGA!="Unknown"]))[[4]][2,4]
ggplot(data = NULL, aes(x = factor(stagesTCGA[stagesTCGA!="Unknown"]), y = trajToUseProfileTCGA[stagesTCGA!="Unknown"]))+geom_boxplot()+geom_jitter(height = 0, width = 0.1)+ggtitle(paste("TCGA",tStageTCGA,sep=" _ "))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


stagesGSE83586 = annotationsGSE83586$Pathologically_Reviewed_TURB_Sample_Stage
stagesGSE83586[stagesGSE83586 == "pTa" | stagesGSE83586 == "pTis"] = "Stage0"
stagesGSE83586[stagesGSE83586 == "pT1"] = "Stage1"
stagesGSE83586[stagesGSE83586 == "pT2"] = "Stage2"
stagesGSE83586[stagesGSE83586 == "pT3"] = "Stage3"
stagesGSE83586[stagesGSE83586 == "pT4"] = "Stage4"
stagesGSE83586[stagesGSE83586 == "pTx"] = "Unknown"
stageGSE83586Cont = as.numeric(gsub("Stage","",stagesGSE83586[stagesGSE83586!="Unknown"]))
tStageGSE83586 = summary(lm(stageGSE83586Cont~trajToUseProfileGSE83586[stagesGSE83586!="Unknown"]))[[4]][2,4]
ggplot(data = NULL, aes(x = factor(stagesGSE83586[stagesGSE83586!="Unknown"]), y = trajToUseProfileGSE83586[stagesGSE83586!="Unknown"]))+geom_boxplot()+geom_jitter(height = 0, width = 0.1)+ggtitle(paste("GSE83586",tStageGSE83586,sep=" _ "))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### S5C ####
clinicalFullTCGA = read.table(file.path(folder,"clinicalFull.tsv"),sep="\t",row.names = 1,header = T,check.names = F)
clinicalFullTCGA = clinicalFullTCGA[intersect(annotationsTCGA$case_id, row.names(clinicalFullTCGA)),]
clinicalTCGAIndexes = which(!is.na(match(annotationsTCGA$case_submitter_id,clinicalFullTCGA$case_submitter_id)))
trajForTCGAGroups = trajToUseProfileTCGA[clinicalTCGAIndexes]

cancerTypeTCGA = clinicalFullTCGA$primary_diagnosis
cancerTypeTCGA[grep("transitional",tolower(cancerTypeTCGA),invert = T)] = "Unknown"
tTypeTCGA = t.test(trajForTCGAGroups[cancerTypeTCGA=="Papillary transitional cell carcinoma"],
                   trajForTCGAGroups[cancerTypeTCGA=="Transitional cell carcinoma"],alternative = "less")$p.value
ggplot(data = NULL, aes(x = factor(cancerTypeTCGA[cancerTypeTCGA!="Unknown"]), y = trajForTCGAGroups[cancerTypeTCGA!="Unknown"]))+geom_boxplot()+geom_jitter(height = 0, width = 0.1)+ggtitle(paste("TCGA",tTypeTCGA,sep=" _ "))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### S5D ####
cancerPriorTCGA = clinicalFullTCGA$prior_malignancy
tPriorTCGA = t.test(trajForTCGAGroups[cancerPriorTCGA=="no"],trajForTCGAGroups[cancerPriorTCGA=="yes"],alternative = "less")$p.value
ggplot(data = NULL, aes(x = factor(cancerPriorTCGA), y = trajForTCGAGroups))+geom_boxplot()+geom_jitter(height = 0, width = 0.1)+ggtitle(paste("TCGA",tPriorTCGA,sep=" _ "))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### 3B ####
factorToPlot = "ESTIMATE_purity"
dataPlotter(metaDataNumTrain, trajToUseProfileTrain,factorToPlot,regLine = T, DoF=3)
dataPlotter(metaDataNumTrain, trajectoryTimeTrain,factorToPlot,regLine = T, DoF=3)

#### S5E-F ####
ageTCGA = annotationsTCGA$age_at_index
beforeAfterPointTCGA = rep(-1, length(trajToUseProfileTCGA))
beforeAfterPointTCGA[trajToUseProfileTCGA>metPoint] = 1
agePTCGA = t.test(ageTCGA[beforeAfterPointTCGA==-1],ageTCGA[beforeAfterPointTCGA==1])$p.value
ggplot(data = NULL, aes(x = factor(beforeAfterPointTCGA), y = ageTCGA))+geom_boxplot()+geom_jitter(height = 0, width = 0.1)+ggtitle(agePTCGA)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

sexTableTCGA = table(annotationsTCGA$gender,beforeAfterPointTCGA)
fisher.test(sexTableTCGA)

beforeAfterPoint = rep(-1, length(trajToUseProfileTrain))
beforeAfterPoint[trajToUseProfileTrain>metPoint] = 1
sexTableTrain = table(annotationsTrain$Sex,beforeAfterPoint)
fisher.test(sexTableTrain)

fullSexTable = rbind(t(t(sexTableTCGA)/colSums(sexTableTCGA)),t(t(sexTableTrain)/colSums(sexTableTrain)))
row.names(fullSexTable) = c("F-TCGA","M-TCGA","F-Long","M-Long")
sexTableDF = as.data.frame(reshape2::melt(fullSexTable))
colnames(sexTableDF) = c("General","PrePost","Ratio")
sexTableDF$Sex = sapply(as.character(sexTableDF$General),function(x){unlist(strsplit(x,"-"))[1]})
sexTableDF$Data = sapply(as.character(sexTableDF$General),function(x){unlist(strsplit(x,"-"))[2]})
ggplot(data = sexTableDF, aes(x = interaction(PrePost,Data), y = Ratio, fill = factor(Sex)))+geom_bar(stat="identity")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### 3C ####
source(file.path(folder,"CIBERSORT_LibLinear.R"))
LM22 = read.table(file.path(folder, "LM22Matrix.txt"),sep="\t", header = T, row.names = 1, check.names = F)
LM22 = LM22[intersect(row.names(LM22),row.names(GEDataTrain)),]
singleCellData = readRDS(file.path(folder, "GSM4307111_SingleCellData.rds"))
singleCellInfo = read.table(file.path(folder, "GSM4307111_CellInfo.txt"),sep="\t", header = T, row.names = 1, check.names = F)
cellTypes = singleCellInfo[colnames(singleCellData),"cell_type"]
selectedCells = grep("Unknown",cellTypes,invert = T)
singleCellData = singleCellData[,selectedCells]
cellTypes = cellTypes[selectedCells]

geneMarkers = unique(unlist(lapply(unique(cellTypes), function(currCellType){
  currCells = singleCellData[,cellTypes==currCellType]
  otherCells = singleCellData[,cellTypes!=currCellType]
  names(sort(rowMeans(currCells)-rowMeans(otherCells),decreasing = T)[1:20])
})))
signatureMatrix = do.call(cbind,lapply(unique(cellTypes), function(currCellType){
  rowMeans(singleCellData[geneMarkers,cellTypes==currCellType])
}))
row.names(signatureMatrix) = geneMarkers
colnames(signatureMatrix) = unique(cellTypes)

geneMarkersMain = intersect(geneMarkers, row.names(GEDataTrain))
cibersortRes = CIBERSORT(signatureMatrix[geneMarkersMain,],GEDataTrain[geneMarkersMain,])
cibersortResLM22 = CIBERSORT(LM22,GEDataTrain[geneMarkersMain,])
cibersortResLM22 = cbind(cibersortResLM22,cibersortResLM22[,"T cells CD4 memory resting"]+cibersortResLM22[,"T cells CD4 memory activated"])
colnames(cibersortResLM22)[dim(cibersortResLM22)[2]] = "T cells CD4 memory"

ciberCorWithTraj = cor(cibersortRes,trajToUseProfileTrain)[,1]
ciberCorWithTrajLM22 = cor(cibersortResLM22,trajToUseProfileTrain)[,1]

cellTypeOrder = names(sort(ciberCorWithTraj))
cellTypeOrderLM22 = names(sort(ciberCorWithTrajLM22))

cellTypeDeconPlots = lapply(cellTypeOrder, function(currCellType){
  currDeconLevels = cibersortRes[,currCellType]
  trajDisc = rep(1, length(currDeconLevels))
  trajDisc[trajToUseProfileTrain<metPoint] = 0
  currPlot = dataPlotter(as.data.frame(trajDisc), currDeconLevels,"trajDisc")
  currPlot$labels$title = gsub("trajDisc", currCellType, currPlot$labels$title)
  currPlot
})

cellTypeDeconPlotsLM22 = lapply(cellTypeOrderLM22, function(currCellType){
  currDeconLevels = cibersortResLM22[,currCellType]
  trajDisc = rep(1, length(currDeconLevels))
  trajDisc[trajToUseProfileTrain<metPoint] = 0
  currPlot = dataPlotter(as.data.frame(trajDisc), currDeconLevels,"trajDisc")
  currPlot$labels$title = gsub("trajDisc", currCellType, currPlot$labels$title)
  currPlot
})

#### S5G ####
cibersortResTest = CIBERSORT(signatureMatrix[geneMarkersMain,],GEDataTest[geneMarkersMain,])
cibersortResGSE83586 = CIBERSORT(signatureMatrix[geneMarkersMain,],dataGSE83586[geneMarkersMain,])
cibersortResTCGA = CIBERSORT(signatureMatrix[intersect(row.names(dataTCGANorm),geneMarkersMain),],dataTCGANorm[intersect(row.names(dataTCGANorm),geneMarkersMain),])
ciberCorWithTrajTest = cor(cibersortResTest,trajToUseProfile)[,1]
ciberCorWithTrajGSE83586 = cor(cibersortResGSE83586,trajToUseProfileGSE83586)[,1]
ciberCorWithTrajTCGA = cor(cibersortResTCGA,trajToUseProfileTCGA)[,1]

cibersortResTestLM22 = CIBERSORT(LM22,GEDataTest[row.names(LM22),])
cibersortResGSE83586LM22 = CIBERSORT(LM22,dataGSE83586[row.names(LM22),])
cibersortResTCGALM22 = CIBERSORT(LM22[intersect(row.names(dataTCGANorm),row.names(LM22)),],dataTCGANorm[intersect(row.names(dataTCGANorm),row.names(LM22)),])
ciberCorWithTrajTestLM22 = cor(cibersortResTestLM22,trajToUseProfile)[,1]
ciberCorWithTrajGSE83586LM22 = cor(cibersortResGSE83586LM22,trajToUseProfileGSE83586)[,1]
ciberCorWithTrajTCGALM22 = cor(cibersortResTCGALM22,trajToUseProfileTCGA)[,1]

fullCorMatrix = t(cbind(c(ciberCorWithTraj,ciberCorWithTrajLM22[1:(length(ciberCorWithTrajLM22)-1)]),c(ciberCorWithTrajTest,ciberCorWithTrajTestLM22),
                        c(ciberCorWithTrajGSE83586,ciberCorWithTrajGSE83586LM22),c(ciberCorWithTrajTCGA,ciberCorWithTrajTCGALM22)))
row.names(fullCorMatrix) = c("Train", "Test","GSE83586", "TCGA")

#### S5H ####
cellTypeDeconPlots[[which(cellTypeOrder=="Basal tumor cells")]]

#### 3D ####
library(survminer)
library(RTCGA.clinical)
library(survival)

deathMatrix = cbind(clinicalFullTCGA$days_to_last_follow_up,annotationsTCGA$days_to_death[clinicalTCGAIndexes])
survTime = apply(deathMatrix,1,function(x){max(x[!is.na(x)])})
survGroup = as.numeric(!is.na(deathMatrix[,2]))
pseudoGroups = as.numeric(trajForTCGAGroups>metPoint)
stageGroups = as.numeric(gsub("T","",stagesTCGA[clinicalTCGAIndexes]))
classGroups = gsub("-Inf","",TCGAClassifier[clinicalTCGAIndexes,]$LundTax_subclasses)
TCGAAgeAtIndex = clinicalFullTCGA$age_at_index
TCGAGender = clinicalFullTCGA$gender
survDF = data.frame(time = survTime, vital = survGroup, pseudo = trajForTCGAGroups, inflection = pseudoGroups, stage = stageGroups, age = TCGAAgeAtIndex, sex = TCGAGender)
survDF = survDF[survTime>-Inf,]
survFit <- survfit(Surv(time, vital) ~ inflection, data = survDF)
surv_pvalue(survFit)
ggsurvplot(survFit, data = survDF, risk.table = F)

#### 3E ####
TCGAClassifier = read.table(file.path(folder,"LundTax_and TCGA_subtypes_TCGA_data.txt"), sep="\t", row.names=1, header = T, check.names = F)
TCGAClassifier = TCGAClassifier[annotationsTCGA$case_submitter_id,]
uroASamplesTCGA = which(TCGAClassifier$LundTax_subclasses=="Basal_squamous")
trajUroASamplesTCGA = trajToUseProfileTCGA[uroASamplesTCGA]
annotationsTCGAUroA = annotationsTCGA[uroASamplesTCGA,]
highUro = as.numeric(trajUroASamplesTCGA>metPoint)
survTable = table(annotationsTCGAUroA$vital_status,highUro)
survTableP = fisher.test(survTable)$p.value

survDFFull = as.data.frame(reshape2::melt(t(t(survTable)/colSums(survTable))))
ggplot(data = survDFFull, aes(x = factor(highUro), y = value, fill = factor(Var1)))+geom_col()+ggtitle(survTableP)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### S5I-J ####
classForSurvival = "Urothelial-like"
#classForSurvival = "Basal_squamous"
specificIndexUro = which(TCGAClassifier[clinicalTCGAIndexes,]$LundTax_classes==classForSurvival)
survFitSpecific <- survfit(Surv(time, vital) ~ inflection, data = survDF[specificIndexUro,])
surv_pvalue(survFitSpecific)
ggsurvplot(survFitSpecific, data = survDF[specificIndexUro,], risk.table = F)

##### Figure 4 #####
#### S6A ####
tumorSubType = annotationsTrain$Subtype_LundTax_RNA
patientsNoBasal = unique(sampleNamesTrain)[sapply(unique(sampleNamesTrain),function(currPatient){
  length(which(c("Ba/Sq") %in% tumorSubType[sampleNamesTrain == currPatient]))==0
})]
GEDataTrainNoBasal = GEDataTrain[,sampleNamesTrain %in% patientsNoBasal]
sampleNamesTrainNoBasal = sampleNamesTrain[sampleNamesTrain %in% patientsNoBasal]

multiAlignmentNoBasal = modelCreation(GEDataTrainNoBasal[proteinCodingGenes,], sampleNamesTrainNoBasal, numOfTopFeatures = 100)
trajToUseProfileTrainNoBasal = predictByConsensus(multiAlignmentNoBasal,GEDataTrain)$predictions
noBasalOverlap = cor(trajToUseProfileTrain,trajToUseProfileTrainNoBasal)
ggplot(data = NULL, aes(x = trajToUseProfileTrain,y = trajToUseProfileTrainNoBasal))+geom_point()+ggtitle(noBasalOverlap)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### S6B ####
regAlignedNoBasal = summary(lm(t(GEDataTrain)~trajToUseProfileTrainNoBasal))
genePValuesAlignedNoBasal = sapply(regAlignedNoBasal, function(x){
  x[[4]][2,4]
})

geneRegAlignedNoBasalLog = -log(genePValuesAlignedNoBasal,base=10)
geneColourNoBasalVsAligned = rep(0, length(geneRegAlignedNoBasalLog))
geneColourNoBasalVsAligned[geneRegAlignedNoBasalLog>newCutOff & geneRegAlignedLog<newCutOff] = 1
geneColourNoBasalVsAligned[geneRegAlignedNoBasalLog<newCutOff & geneRegAlignedLog>newCutOff] = 2
geneColourNoBasalVsAligned[geneRegAlignedNoBasalLog>newCutOff & geneRegAlignedLog>newCutOff] = 3
ggplot(data = NULL, aes(x = geneRegAlignedNoBasalLog, y = geneRegAlignedLog, colour = factor(geneColourNoBasalVsAligned)))+geom_point()+
  geom_hline(yintercept = newCutOff, linetype = "dashed")+geom_vline(xintercept = newCutOff, linetype = "dashed")+
  scale_color_manual(values = c("0" = "lightgray","1" = "darkgray","2" = "purple", "3" = "orange"))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### 4A ####
lundTrain = annotationsTrain$Subtype_LundTax_RNA
tLund =aov(trajToUseProfileTrain~lundTrain)
summary(tLund)[[1]][["Pr(>F)"]][1]
ggplot(data = NULL, aes(x = reorder(factor(lundTrain),trajToUseProfileTrain), y = trajToUseProfileTrain))+geom_boxplot()+geom_jitter(height = 0, width = 0.1)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### 4B ####
tumorSubType = annotationsTrain$Subtype_LundTax_RNA

inflectionVector = rep("Before",length(tumorSubType))
inflectionVector[trajToUseProfileTrain>metPoint] = "After"
inflectionTypeTable = table(tumorSubType,inflectionVector)
inflectionTypeTable = t(t(inflectionTypeTable)/colSums(inflectionTypeTable))

inflectionDF = as.data.frame(reshape2::melt(inflectionTypeTable))
ggplot(data = inflectionDF, aes(x = reorder(factor(tumorSubType),value,na.rm = TRUE), y = value, fill = factor(inflectionVector)))+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### 4C ####
tumorSubTypeTrain = gsub(" ", "", annotationsTrain$Prim_Rec_Prog)
indexesForFigure4 = which(trajToUseProfileTrain<metPoint & tumorSubTypeTrain!="prog")
primRecT = t.test(trajToUseProfileTrain[tumorSubTypeTrain=="primary"],trajToUseProfileTrain[tumorSubTypeTrain=="rec"],alternative = "less")$p.value
ggplot(data = NULL, aes(x = reorder(factor(tumorSubTypeTrain[indexesForFigure4]),trajToUseProfileTrain[indexesForFigure4],na.rm = TRUE), y = trajToUseProfileTrain[indexesForFigure4]))+geom_boxplot()+geom_jitter(height = 0, width = 0.1)+ggtitle(primRecT)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### 4D ####
tumorSubType = annotationsTrain$Subtype_LundTax_RNA
ggplot(data = NULL, aes(x = reorder(factor(tumorSubType[indexesForFigure4]),trajToUseProfileTrain[indexesForFigure4],na.rm = TRUE), y = trajToUseProfileTrain[indexesForFigure4]))+geom_boxplot()+geom_jitter(height = 0, width = 0.1)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### 4E ####
classForSurvival = "Urothelial-like"
specificIndexUro = which(TCGAClassifier[clinicalTCGAIndexes,]$LundTax_classes==classForSurvival)
survDFUro = survDF[specificIndexUro,]
survDFUroPre = survDFUro[survDFUro$inflection==0,]
survDFForBins = survDFUroPre
survBinUroA = sapply(seq(0.1,metPoint,0.1),function(i){
  currSurTable = table(survDFForBins$vital[survDFForBins$pseudo>(i-0.1) & survDFForBins$pseudo<i])
  currSurTable[1]/sum(currSurTable)
})
uroASurvTrendP = summary(lm(survBinUroA~seq(0.1,metPoint,0.1)))[[4]][2,4]
ggplot(data = NULL, aes(x = seq(0.1,metPoint,0.1), y = survBinUroA-0.4))+geom_col()+ggtitle(uroASurvTrendP)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### S6C ####
ggplot(data = NULL, aes(x = lundTrain[stages!="Unknown"], y = stageCont))+geom_violin()+geom_jitter(position=position_jitter(height=0, width=0.3),size = 1)+ggtitle("Train")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### S6D ####
lundTCGAGlobal = TCGAClassifier$LundTax_subclasses
lundTCGAForComp = lundTCGAGlobal
lundTCGAForComp[grep("Basal_squamous",lundTCGAForComp)] = "Ba/Sq"
lundTCGAForComp[grep("GU",lundTCGAForComp)] = "GU"
lundTCGAForComp[grep("Mes",lundTCGAForComp)] = "Mes"
lundTCGAForComp[grep("ScNE",lundTCGAForComp)] = "ScNE"
stageContTCGA = as.numeric(gsub("T","",stagesTCGA[stagesTCGA!="Unknown"]))
lundTCGAForComp = lundTCGAForComp[stagesTCGA!="Unknown"]
stageContTCGA = stageContTCGA[which(lundTCGAForComp %in% unique(lundGlobal))]
lundTCGAForComp = lundTCGAForComp[which(lundTCGAForComp %in% unique(lundGlobal))]
ggplot(data = NULL, aes(x = lundTCGAForComp[!is.na(stageContTCGA)], y = stageContTCGA))+geom_violin()+geom_jitter(position=position_jitter(height=0, width=0.3),size = 1)+ggtitle("TCGA")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### S6E ####
tumorSubType = annotationsTrain$Consensus_Subtype_RNA
ggplot(data = NULL, aes(x = reorder(factor(tumorSubType[indexesForFigure4]),trajToUseProfileTrain[indexesForFigure4],na.rm = TRUE), y = trajToUseProfileTrain[indexesForFigure4]))+geom_boxplot()+geom_boxplot()+geom_jitter(height = 0, width = 0.1)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### S6F ####
classForSurvival = "Urothelial-like"
specificIndexUro = which(TCGAClassifier[clinicalTCGAIndexes,]$LundTax_classes==classForSurvival)
survDFUro = survDF[specificIndexUro,]
survDFUroPre = survDFUro[survDFUro$inflection==0,]
progSetScoresTCGA = trajToUseProfileTCGA[clinicalTCGAIndexes][specificIndexUro][survDFUro$inflection==0]
survDFUroPre$progSet = progSetScoresTCGA
survDFUroPre$progSetClass = as.numeric(survDFUroPre$progSet>0.25)

survTableUroEarlyLate = table(survDFUroPre$vital,survDFUroPre$progSetClass)
survDFEarlyLate = as.data.frame(reshape2::melt(t(t(survTableUroEarlyLate)/colSums(survTableUroEarlyLate))))
survTablePEarlyLate = fisher.test(survTableUroEarlyLate)$p.value
ggplot(data = survDFEarlyLate, aes(x = factor(Var2), y = value, fill = factor(Var1)))+geom_col()+ggtitle(survTablePEarlyLate)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### S6G ####
uroASamples = which(annotationsTrain$Subtype_LundTax_RNA=="UroA-Prog" & trajToUseProfileTrain<metPoint)
trajUroASamples = trajToUseProfileTrain[uroASamples]
uroSpecGenes = c("CCND1","FGFR3","FOXA1","RB1","CDKN2A","GATA3","ERBB2","PPARG","XBP1")
geneDataForUro = GEDataTrain[uroSpecGenes,uroASamples]
geneDataForUro = geneDataForUro[,order(trajUroASamples)]
geneDataForUroWithMean = cbind(geneDataForUro,rowMeans(geneDataForUro))
geneDataForUroWithMean

##### Figure 5 #####
#### 5A ####
GEDataUroAllGenes = GEDataTrain[,uroASamples]
GEDataUroAllGenes = GEDataUroAllGenes[,order(trajUroASamples)]
trajUroASorted = sort(trajUroASamples)
uroAGenePAll = apply(GEDataUroAllGenes,1,function(x){
  t.test(x[trajUroASorted<0.25],x[trajUroASorted>0.25])$p.value
})
uroAGenePAllQ = p.adjust(uroAGenePAll, method = "fdr")
#signGeneUroA = names(uroAGenePAllQ[uroAGenePAll<0.05])
signGeneUroA = names(uroAGenePAllQ[uroAGenePAllQ<0.05])
signGeneUroA = signGeneUroA[!(signGeneUroA %in% seedGenes)]
signGenesCorMatrix = cor(t(GEDataUroAllGenes[signGeneUroA,]))

createGeneHeatmap = function(coExp, prvStats = NULL){
  library(WGCNA)
  
  if(is.null(prvStats)){
    dissTOM=1-coExp
    hierADJ=hclust(as.dist(dissTOM), method="average")
    currOrder = hierADJ$order
  }else{
    genesToCheck = prvStats$genes
    currOrder = prvStats$order
    coExp = coExp[genesToCheck,genesToCheck]
  }
  collectGarbage()
  
  pd.m <- reshape2::melt(coExp[currOrder,currOrder], id.vars = "Gene1", variable.name = "Gene2")
  currPlot = ggplot( pd.m, aes(Var1, Var2) ) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(low="blue", mid = 'white',midpoint = 0, high="red")+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  list(plot = currPlot, order = currOrder)
}
signGenesHeatmap = createGeneHeatmap(signGenesCorMatrix)

#### S6H-I ####
firstGroup = 1:1055
secondGroup = (max(firstGroup)+1):dim(signGenesCorMatrix)[1]

smallGroup = row.names(signGenesCorMatrix)[signGenesHeatmap$order][firstGroup]
bigGroup = row.names(signGenesCorMatrix)[signGenesHeatmap$order][secondGroup]

#geneTotalScore = colMeans(GEDataTrain[bigGroup,])
geneTotalScore = colMeans(GEDataTrain[smallGroup,])
ggplot(data = NULL, aes(x = trajToUseProfileTrain[uroASamples], y = geneTotalScore[uroASamples]))+geom_point(size = 4)+geom_smooth()+geom_vline(xintercept = metPoint,linetype="dashed")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### 5B ####
lcRNATypes = c("^RNU","^MIR","^RN7SK","^RN7SL","^RNA5S")

pseudoEnrichment = function(groups,groupForEnrich, PW = F, groupCutOff = 0){
  sapply(groups,function(currGroup){
    if(PW){
      groupGenes = currGroup
    }else{
      groupGenes = row.names(GEDataTrain)[grep(currGroup,row.names(GEDataTrain))] 
    }
    groupSelected = length(intersect(groupForEnrich,groupGenes))
    if(groupSelected<groupCutOff){
      return(NA)
    }
    groupTotal = length(groupGenes)
    nonGroupTotal = dim(GEDataTrain)[1]-groupTotal
    phyper(groupSelected, groupTotal, nonGroupTotal, length(bigGroup), lower.tail = F, log.p = FALSE)
  })
}

set.seed(2)
smallGroupEnrich = -log(p.adjust(pseudoEnrichment(lcRNATypes,smallGroup),method = "fdr"),base=10)
bigGroupEnrich = -log(p.adjust(pseudoEnrichment(lcRNATypes,bigGroup),method = "fdr"),base = 10)
randomGroupEnrich = -log(p.adjust(pseudoEnrichment(lcRNATypes,sample(row.names(GEDataTrain),length(bigGroup))),method = "fdr"),base = 10)

groupEnrichDF = as.data.frame(reshape2::melt(cbind(bigGroupEnrich,smallGroupEnrich,randomGroupEnrich)))
ggplot(data = groupEnrichDF)+geom_bar(aes(x = Var1, y = value, fill = Var2), stat='identity',position="dodge")+ coord_flip()+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### 5C ####
geneSetsList = apply(geneSets,1,function(x){unlist(strsplit(as.character(as.matrix(x)),","))})
selectedGeneSetsLargeIndexes = c(grep("^HALLMARK",names(geneSetsList)),grep("^REACTOME",names(geneSetsList)),grep("^KEGG",names(geneSetsList)))
geneSetsListSmall = geneSetsList[selectedGeneSetsLargeIndexes]
bigGroupNoPseudo = bigGroup[!(bigGroup %in% bigGroup[unlist(sapply(lcRNATypes, function(currPseudo){grep(currPseudo,bigGroup)}))])]
bigGroupPWEnrich = pseudoEnrichment(geneSetsListSmall,bigGroupNoPseudo, PW = T, groupCutOff = 10)
bigGroupPWEnrichQ = p.adjust(bigGroupPWEnrich, method = "fdr")
bigGroupPWEnrichQSign = bigGroupPWEnrichQ[which(bigGroupPWEnrichQ<0.05)]

selectedBigGroupPWs = c("REACTOME_GPCR_LIGAND_BINDING","REACTOME_POTASSIUM_CHANNELS","REACTOME_OLFACTORY_SIGNALING_PATHWAY")
selectedBigGroupPWsDF = as.data.frame(as.matrix(sort(bigGroupPWEnrich[selectedBigGroupPWs])))
selectedBigGroupPWsDF$Name = row.names(selectedBigGroupPWsDF)
selectedBigGroupPWsDF$enrichLog = -log(selectedBigGroupPWsDF$V1,base=10)
ggplot(data = selectedBigGroupPWsDF, aes(x = reorder(Name,enrichLog), y = enrichLog))+geom_col()+coord_flip()+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### 5D ####
smallGroupPWEnrich = pseudoEnrichment(geneSetsListSmall,smallGroup, PW = T, groupCutOff = 10)
smallGroupPWEnrichQ = p.adjust(smallGroupPWEnrich, method = "fdr")
smallGroupPWEnrichQSign = smallGroupPWEnrichQ[which(smallGroupPWEnrichQ<0.05)]

selectedSmallGroupPWs = c("REACTOME_RNA_POLYMERASE_II_PRE_TRANSCRIPTION_EVENTS","REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_DNA_REPAIR_GENES","REACTOME_MITOTIC_SPINDLE_CHECKPOINT","REACTOME_RESOLUTION_OF_SISTER_CHROMATID_COHESION","REACTOME_INTRA_GOLGI_AND_RETROGRADE_GOLGI_TO_ER_TRAFFIC","REACTOME_MACROAUTOPHAGY","KEGG_OXIDATIVE_PHOSPHORYLATION","REACTOME_CELL_CYCLE","KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS","REACTOME_REGULATION_OF_TP53_ACTIVITY")
selectedSmallGroupPWsDF = as.data.frame(as.matrix(sort(smallGroupPWEnrichQSign[selectedSmallGroupPWs])))
selectedSmallGroupPWsDF$Name = row.names(selectedSmallGroupPWsDF)
selectedSmallGroupPWsDF$enrichLog = -log(selectedSmallGroupPWsDF$V1,base=10)
ggplot(data = selectedSmallGroupPWsDF, aes(x = reorder(Name,enrichLog), y = enrichLog))+geom_col()+coord_flip()+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### 5E ####
uroASamples = which(annotationsTrain$Subtype_LundTax_RNA=="UroA-Prog" & trajToUseProfileTrain<metPoint)
trajUroASamples = trajToUseProfileTrain[uroASamples]
kinetochoreGenes = c("NSL1","SKA2","HAUS2","CLASP2","RSF1","PPP2R5A","NUP37","NUP133","NUP214","RSF1")
geneDataForUro = GEDataTrain[kinetochoreGenes,uroASamples]
geneDataForUroNorm = t(apply(geneDataForUro,1,function(x){(x-mean(x))/sd(x)}))
colnames(geneDataForUroNorm) = as.numeric(sort(trajUroASamples)>0.25)
geneDataForUroDF = as.data.frame(reshape2::melt(geneDataForUroNorm[kinetochoreGenes,]))
uroGenesPs = apply(geneDataForUroNorm[kinetochoreGenes,],1,function(x){
  t.test(x[colnames(geneDataForUroNorm)==0],x[colnames(geneDataForUroNorm)==1])$p.value
})
uroGenesP = paste(formatC(uroGenesPs, digits = 4, format = "f"), collapse = " _ ")
ggplot(data = geneDataForUroDF, aes(x = factor(Var1), y = value, colour = factor(Var2)))+geom_boxplot()+geom_hline(yintercept = 0, linetype = "dashed")+ggtitle(uroGenesP)+
  theme_bw() + theme(title = element_text(size=9), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
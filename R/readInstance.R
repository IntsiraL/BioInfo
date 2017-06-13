library(limma)
#library(illuminaio)

idatfilesPath<-getwd()
#######################################################################################################
#######################################################################################################
############################## LECTURE DES FICHIERS IDATS
#######################################################################################################
#######################################################################################################
idatFiles44 <- list.files(paste(idatfilesPath,"/../Data/200729890044/",sep = ""),all.files=FALSE,pattern=".idat")
idatFiles64 <- list.files(paste(idatfilesPath,"/../Data/200796240064/",sep = ""),all.files=FALSE,pattern=".idat")
bgxfile=paste(idatfilesPath,"/../Data/HumanHT-12_V4_0_R2_15002873_B.bgx",sep = "")
controlProfil  <- read.table(paste(idatfilesPath,"/../Data/controlProfil.txt",sep = ""), sep = "\t", header = TRUE ,as.is = TRUE)
idatFiles = c()
for (i in idatFiles44){
  idatFiles <- c(idatFiles,paste(idatfilesPath,"/../Data/200729890044/",i,sep = ""))
}
for (i in idatFiles64){
  idatFiles <- c(idatFiles,paste(idatfilesPath,"/../Data/200796240064/",i,sep = ""))
}

#lectures des fichiers idat par ordre alphabétique des noms des fichiers avec read.idat de limma
obj<-read.idat(idatFiles, bgxfile, dateinfo=TRUE,annotation = "Symbol",tolerance=0L, verbose = TRUE)

#renomer les colonnes de la matrice des intensités
nameCol <- c("P1_G","P1_MAL","P1_GLPS","P1_GALLPS","P1_M","P1_GAL","P1_MLPS","P1_MALLPS")

  for (i in c("","AL","LPS","ALLPS")){
    nameCol <- c(nameCol, paste("P2_G",i,sep = ""))
  }

  for (i in c("","AL","LPS","ALLPS")){
    nameCol <- c(nameCol, paste("P2_M",i,sep = ""))
  }

  for (i in c("","AL","LPS","AL_LPS")){
    nameCol <- c(nameCol, paste("P3_G",i,sep = ""))
  }

  for (i in c("","AL","LPS","ALLPS")){
    nameCol <- c(nameCol, paste("P3_M",i,sep = ""))
  }
colnames(obj$E) <- nameCol
colnames(obj$other$NumBeads) <- nameCol
colnames(obj$other$STDEV) <- nameCol

#juste un petit control des p-values
obj$genes$DetectionPValue <- detectionPValues(obj)

#décryptage des fichier idat avec illuminao
#mydata = list()
#for (i in 1:24){
#  mydata[[nameCol[i]]] <- readIDAT(idatFiles[i])
#}

#targetData <- readTargets(file = "../Data/Annot.txt")
###############################################################################
## Control Data                                                              ##
###############################################################################
# source(file = "qa_metrics.R")
controlData <- obj[obj$genes$Status != "regular",]
bruteData <- obj[obj$genes$Status == "regular",]
source(file = "controlData.R")
# #estimation de proportion des sondes exprimées
# plot(propexpr(obj),type = "b", main = "Estimation de proportion des sondes exprimées",ylab = "",xlab = "")
# axis(1,1:24,nameCol , las = 2)
# abline(v=12.5, col="red")
###############################################################################
## Preprocessing Transformation Normalization Filtrage                                                              ##
###############################################################################
#correction de bruit de fond et normalisaion (par quantille) avec neqc (ajustement des paramètres au fur et à mesure)
# à tester la fonction backgroundCorrect()
dCorect <- neqc(obj)
# plotMD(obj)
# plotMD(dCorect)
# plotDensities(obj$E[,1:8],legend = "topright",main = "Densité de la population 1 avant la normalisation")
# plotDensities(dCorect$E[,1:8],legend = "topright",main = "Densité de la population 1 après la normalisation")
# boxplot(dCorect$E[,1:8],main="Box plot des signaux de la population 1 après normalisation", las=2)
################################
# Calcul des moyennes des P-Values
################################
# tab = c()
# for(i in c(1:nrow(dCorect$E))){
#   tab <- c(tab,mean(dCorect$genes$DetectionPValue[i,1:8]))
# }
# dCorect$genes$mPValuePop1 <- matrix(tab,nrow = nrow(dCorect$E), ncol = 1)
# tab = c()
# for(i in c(1:nrow(dCorect$E))){
#   tab <- c(tab,mean(dCorect$genes$DetectionPValue[i,9:16]))
# }
# dCorect$genes$mPValuePop2 <- matrix(tab,nrow = nrow(dCorect$E), ncol = 1)
# tab = c()
# for(i in c(1:nrow(dCorect$E))){
#   tab <- c(tab,mean(dCorect$genes$DetectionPValue[i,17:24]))
# }
# dCorect$genes$mPValuePop3 <- matrix(tab,nrow = nrow(dCorect$E), ncol = 1)
# boxplot(dCorect$genes$mPValuePop1)
# abline(h=0.01, col="red")
# abline(h=0.05, col="red")
################################
# Analyse différentielle regression linéaire
################################

#tester la relation entre var
summary(lm(formula = P1_GAL ~ P1_G,data=as.data.frame(dCorect$E[,c("P1_GAL","P1_G")])))
summary(lm(formula = P1_M ~ P1_G,data=as.data.frame(dCorect$E[,c("P1_M","P1_G")])))

#P1_G vs P1_GAL
design <- model.matrix(~1+nameCol)
colnames(design) <- nameCol
contMatrix <- makeContrasts(GvsGAL=P1_G-P1_GAL,MvsG=P1_M-P1_G,levels = design)
fit <- lmFit(dCorect, design)
fitC <- contrasts.fit(fit, contMatrix)
fitC <- eBayes(fitC)

designGvsGAL <- cbind(G=as.numeric(nameCol == "P1_G"),GAL=as.numeric(nameCol == "P1_GAL"))
p1GvsGAL <- lmFit(dCorect,designGvsGAL)
p1GvsGAL <- eBayes(p1GvsGAL)
resultat1 <- toptable(p1GvsGAL,number = nrow(dCorect))

design2GvsGAL <- cbind(G=as.numeric(nameCol == "P2_G"),GAL=as.numeric(nameCol == "P2_GAL"))
p2GvsGAL <- lmFit(dCorect,design2GvsGAL)
p2GvsGAL <- eBayes(p2GvsGAL)
resultat2 <- toptable(p2GvsGAL,number = nrow(dCorect))

design3GvsGAL <- cbind(G=as.numeric(nameCol == "P3_G"),GAL=as.numeric(nameCol == "P3_GAL"))
p3GvsGAL <- lmFit(dCorect,design3GvsGAL)
p3GvsGAL <- eBayes(p3GvsGAL)
resultat3 <- toptable(p3GvsGAL,number = nrow(dCorect))
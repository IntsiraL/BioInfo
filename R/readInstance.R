library(limma)
library(illuminaio)

idatfilesPath<-getwd()
#######################################################################################################
#######################################################################################################
############################## LECTURE DES FICHIERS IDATS
#######################################################################################################
#######################################################################################################
idatFiles44 <- list.files(paste(idatfilesPath,"/../Data/200729890044/",sep = ""),all.files=FALSE,pattern=".idat")
idatFiles64 <- list.files(paste(idatfilesPath,"/../Data/200796240064/",sep = ""),all.files=FALSE,pattern=".idat")
bgxfile=paste(idatfilesPath,"/../Data/HumanHT-12_V4_0_R2_15002873_B.bgx",sep = "")

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
nameCol <- c("P1_G_Co","P1_M_AL","P1_G_LPS","P1_G_AL_LPS","P1_M_Co","P1_G_AL","P1_M_LPS","P1_M_AL_LPS")

  for (i in c("Co","AL","LPS","AL_LPS")){
    nameCol <- c(nameCol, paste("P2","G",i,sep = "_"))
  }

  for (i in c("Co","AL","LPS","AL_LPS")){
    nameCol <- c(nameCol, paste("P2","M",i,sep = "_"))
  }

  for (i in c("Co","AL","LPS","AL_LPS")){
    nameCol <- c(nameCol, paste("P3","G",i,sep = "_"))
  }

  for (i in c("Co","AL","LPS","AL_LPS")){
    nameCol <- c(nameCol, paste("P3","M",i,sep = "_"))
  }
colnames(obj$E) <- nameCol
colnames(obj$other$NumBeads) <- nameCol
colnames(obj$other$STDEV) <- nameCol

#juste un petit control des p-values
obj$genes$DetectionPValue <- detectionPValues(obj)


#correction de bruit de fond et normalisaion (par quantille) avec neqc (ajustement des paramètres au fur et à mesure)
# à tester la fonction backgroundCorrect()
#dCorect <- neqc(obj)


#décryptage des fichier idat avec illuminao
mydata = list()
for (i in 1:24){
  mydata[[nameCol[i]]] <- readIDAT(idatFiles[i])
}
#print(mydata$P1_G_Co)

targetData <- readTargets(file = "../Data/Annot.txt")
#######################################################################################################
#######################################################################################################
############################## CONTROL DATA
#######################################################################################################
#######################################################################################################
controlData <- obj[obj$genes$Status != "regular",]
bruteData <- obj[obj$genes$Status == "regular",]
###############################################################################
## Hybridation controls                                                      ##
###############################################################################
hybC <- obj[obj$genes$Status == "cy3_hyb",]
for(i in c(1:24)){
  hybC$E[,i] <- sort(hybC$E[,i])
}
low = c()
for(i in c(1:24)){
  low <- c(low,mean(c(hybC$E[1,i],hybC$E[2,i])))
}
medium = c()
for(i in c(1:24)){
  medium <- c(medium,mean(c(hybC$E[3,i],hybC$E[4,i])))
}
high = c()
for(i in c(1:24)){
  high <- c(high,mean(c(hybC$E[5,i],hybC$E[6,i])))
}
vv <- c(mean(low),mean(medium),mean(high))
names(vv) <- c("Low","Medium", "High")
barplot(vv)
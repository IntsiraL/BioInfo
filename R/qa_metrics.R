par(mfcol=c(2,1))
ht12metrics  <- read.table(paste(idatfilesPath,"/../Data/200729890044/Metrics.txt",sep = ""), sep = "\t", header = TRUE ,as.is = TRUE)
ht12snr  <- ht12metrics$P95Grn/ht12metrics$P05Grn
labs  <- paste(ht12metrics[, 2],  ht12metrics[, 3], sep = "_")
par(mai = c(1.5, 0.8, 0.3, 0.1))
plot (1:12 , type = "b",  ht12snr , pch = 19, ylab = "P95 / P05", xlab = "",main = "Signal -to-noise  ratio  for  200729890044  data", axes = FALSE ,
      frame.plot = TRUE)
axis (2)
axis(1, 1:12, nameCol[1:12] , las = 2)
abline(h=10, col="red")
#########################################################
ht12metricsB  <- read.table(paste(idatfilesPath,"/../Data/200796240064/Metrics.txt",sep = ""), sep = "\t", header = TRUE ,as.is = TRUE)
ht12snrB  <- ht12metricsB$P95Grn/ht12metricsB$P05Grn
par(mai = c(1.5, 0.8, 0.3, 0.1))
plot (1:12 ,  ht12snrB ,type = "b", pch = 19, ylab = "P95 / P05", xlab = "",main = "Signal -to-noise  ratio  for  200796240064  data", axes = FALSE ,
      frame.plot = TRUE)
axis (2)
axis(1, 1:12, nameCol[13:24] , las = 2)
abline(h=10, col="red")
par(mfcol=c(1,1))
#########################################
summary(ht12snr)#Signal -to-noise  ratio  for  200729890044  data
summary(ht12snrB)#Signal -to-noise  ratio  for  200796240064  data

# Region<- c("PA","PA","PA","RA","RA","RA","PA","PA","PA","RA","RA","RA")
# Labor <- c("TNL","TNL","TNL","TNL","TNL","TNL","TIL","TIL","TIL","TIL","TIL","TIL")
# REG <- factor(Region)
# dddd<- model.matrix(~0+REG)
# colnames(dddd) <- levels(REG)
#sampleGene  <- read.table(paste(idatfilesPath,"/../Data/dataTest_A.csv",sep = ""), sep = "\t", header = TRUE ,as.is = TRUE)
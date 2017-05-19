path <- getwd()

nameCol <- c("P1_G_Co","P1_M_AL","P1_G_LPS","P1_G_AL_LPS","P1_M_Co","P1_G_AL","P1_M_LPS","P1_M_AL_LPS")

for (i in c("Co","AL","LPS","AL_LPS")){
  nameCol <- c(nameCol, paste("P2","G",i,sep = "_"))
}

ht12metrics  <- read.table(paste(path,"/../Data/200729890044/Metrics.txt",sep = ""), sep = "\t", header = TRUE ,as.is = TRUE)
ht12snr  <- ht12metrics$P95Grn/ht12metrics$P05Grn
labs  <- paste(ht12metrics[, 2],  ht12metrics[, 3], sep = "_")
par(mai = c(1.5, 0.8, 0.3, 0.1))
plot (1:12 , type = "b",  ht12snr , pch = 19, ylab = "P95 / P05", xlab = "",main = "Signal -to-noise  ratio  for  200729890044  data", axes = FALSE ,
      frame.plot = TRUE)
axis (2)
axis(1, 1:12, nameCol , las = 2)
abline(h=10, col="red")
#########################################################
nameCol <- c()
for (i in c("Co","AL","LPS","AL_LPS")){
  nameCol <- c(nameCol, paste("P2","M",i,sep = "_"))
}

for (i in c("Co","AL","LPS","AL_LPS")){
  nameCol <- c(nameCol, paste("P3","G",i,sep = "_"))
}

for (i in c("Co","AL","LPS","AL_LPS")){
  nameCol <- c(nameCol, paste("P3","M",i,sep = "_"))
}
ht12metrics  <- read.table(paste(path,"/../Data/200796240064/Metrics.txt",sep = ""), sep = "\t", header = TRUE ,as.is = TRUE)
ht12snr  <- ht12metrics$P95Grn/ht12metrics$P05Grn
par(mai = c(1.5, 0.8, 0.3, 0.1))
plot (1:12 ,  ht12snr ,type = "b", pch = 19, ylab = "P95 / P05", xlab = "",main = "Signal -to-noise  ratio  for  200796240064  data", axes = FALSE ,
      frame.plot = TRUE)
axis (2)
axis(1, 1:12, nameCol , las = 2)
abline(h=10, col="red")
#######################################################################################################
#######################################################################################################
############################## CONTROL DATA
#######################################################################################################
#######################################################################################################
controlData <- obj[obj$genes$Status != "regular",]
bruteData <- obj[obj$genes$Status == "regular",]
#function to add the error bar representing the confidence interval
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
###############################################################################
## Hybridation controls                                                      ##
###############################################################################
low = c()
for(i in c(1:24)){
  low <- c(low,mean(controlData$E[controlProfil$Reporter_Group_id=="phage_lambda_genome:low",i]))
}
medium = c()
for(i in c(1:24)){
  medium <- c(medium,mean(controlData$E[controlProfil$Reporter_Group_id=="phage_lambda_genome:med",i]))
}
high = c()
for(i in c(1:24)){
  high <- c(high,mean(controlData$E[controlProfil$Reporter_Group_id=="phage_lambda_genome:high",i]))
}
vv <- c(mean(low),mean(medium),mean(high))
vvSD <- c(sd(low),sd(medium),sd(high))
names(vv) <- c("Low","Medium", "High")
par(mfcol=c(2,2))
hybC_Bar <- barplot(vv,ylim = c(0,1.2*max(vv)), ylab = "Signal",main = "Hybridation Controls")
error.bar(hybC_Bar,vv, vvSD)
plot(low, type = "b", ylab = "Signal",xaxt="n", xlab = "",main = "QC Report Hybridation Controls: low")
axis(1, at=c(1:24), labels = nameCol, las=2)
plot(medium, type = "b", ylab = "Signal",xaxt="n", xlab = "",main = "QC Report Hybridation Controls: medium")
axis(1, at=c(1:24), labels = nameCol, las=2)
plot(high, type = "b", ylab = "Signal",xaxt="n", xlab = "",main = "QC Report Hybridation Controls: high")
axis(1, at=c(1:24), labels = nameCol, las=2)
###############################################################################
## Negative controls                                                         ##
###############################################################################
negC <- obj[obj$genes$Status == "negative",]
background = c()
for(i in c(1:24)){
  background <- c(background,mean(c(negC$E[,i])))
}
noise = c()
for(i in c(1:24)){
  noise <- c(noise,sd(c(negC$E[,i])))
}
nn <- c(mean(background),mean(noise))
nnSD <- c(sd(background),sd(noise))
names(nn) <- c("Background", "Noise")
par(mfcol=c(2,2))
negC_Bar <- barplot(nn,ylim = c(0,1.2*max(nn)), ylab = "Signal",main = "Negative Controls")
error.bar(negC_Bar,nn, nnSD)
plot(background, type = "b", ylab = "Signal",xaxt="n", xlab = "",main = "QC Report background")
axis(1, at=c(1:24), labels = nameCol, las=2)
plot(noise, type = "b", ylab = "Signal",xaxt="n", xlab = "",main = "QC Report noise")
axis(1, at=c(1:24), labels = nameCol, las=2)
###############################################################################
## Biotin & High Stringency   | Low Stringency                               ##
###############################################################################
biotC <- obj[obj$genes$Status == "biotin",]
#stringC <- obj[obj$genes$Status == "low_stringency_hyb",]
#for(i in c(1:24)){
#  stringC$E[,i] <- sort(stringC$E[,i])
#}
biot = c()
for(i in c(1:24)){
  biot <- c(biot,mean(c(biotC$E[,i])))
}
#A voir
highStingc = rep(55,24)
#for(i in c(1:24)){
#  highStingc <- c(highStingc,mean(c(stringC$E[8,i])))
#}
mm2 = c()
for(i in c(1:24)){
  mm2 <- c(mm2,mean(controlData$E[controlProfil$Reporter_Group_id=="phage_lambda_genome:mm2",i]))
}
pm = c()
for(i in c(1:24)){
  pm <- c(pm,mean(controlData$E[controlProfil$Reporter_Group_id=="phage_lambda_genome:pm",i]))
}
ss <- c(mean(biot),mean(highStingc))
ssSD <- c(sd(biot),sd(highStingc))
names(ss) <- c("Biotin", "High stringency")

lss <- c(mean(mm2),mean(pm))
lssSD <- c(sd(mm2),sd(pm))
names(lss) <- c("mm2", "pm")

par(mfcol=c(2,2))
bhC_Bar <- barplot(ss,ylim = c(0,1.2*max(ss)), ylab = "Signal",main = "Biotin & High Stringency")
error.bar(bhC_Bar,ss, ssSD)
plot(biot, type = "b", ylab = "Signal",xaxt="n", xlab = "",main = "QC Report biotin")
axis(1, at=c(1:24), labels = nameCol, las=2)
plot(highStingc, type = "b", ylab = "Signal",xaxt="n", xlab = "",main = "QC Report high stringency")
axis(1, at=c(1:24), labels = nameCol, las=2)

par(mfcol=c(2,2))
lsC_Bar <- barplot(lss,ylim = c(0,1.2*max(lss)), ylab = "Signal",main = "Low Stringency")
error.bar(lsC_Bar,lss, lssSD)
plot(mm2, type = "b", ylab = "Signal",xaxt="n", xlab = "",main = "QC Report mm2")
axis(1, at=c(1:24), labels = nameCol, las=2)
plot(pm, type = "b", ylab = "Signal",xaxt="n", xlab = "",main = "QC Report pm")
axis(1, at=c(1:24), labels = nameCol, las=2)
###############################################################################
## Gene Intensity                                                            ##
###############################################################################
houseKC <- obj[obj$genes$Status == "housekeeping",]
housk = c()
for(i in c(1:24)){
  housk <- c(housk,mean(c(houseKC$E[,i])))
}
allge = c()
for(i in c(1:24)){
  allge <- c(allge,mean(c(bruteData$E[,i])))
}
gi <- c(mean(housk),mean(allge))
giSD <- c(sd(housk),sd(allge))
names(gi) <- c("Housekeeping genes", "all genes")
par(mfcol=c(2,2))
giC_Bar <- barplot(gi, ylim = c(0,1.2*max(gi)), ylab = "Signal",main = "Gene Intensity")
error.bar(giC_Bar,gi, giSD)
plot(housk, type = "b", ylab = "Signal",xaxt="n", xlab = "",main = "QC Report Housekeeping")
axis(1, at=c(1:24), labels = nameCol, las=2)
plot(allge, type = "b", ylab = "Signal",xaxt="n", xlab = "",main = "QC Report gene")
axis(1, at=c(1:24), labels = nameCol, las=2)
###############################################################################
## Labeling & Background  | Control Summary                                  ##
###############################################################################
labelC <- obj[obj$genes$Status == "labeling",]
labl = c()
for(i in c(1:24)){
  labl <- c(labl,mean(c(labelC$E[,i])))
}
labk <- c(mean(labl),mean(background))
labkSD <- c(sd(labl),sd(background))
names(labk) <- c("Labeling","Background")
par(mfcol=c(3,2))
hybC_Bar <- barplot(vv,ylim = c(0,1.2*max(vv)), ylab = "Signal",main = "Hybridation Controls")
error.bar(hybC_Bar,vv, vvSD)
bhC_Bar <- barplot(ss,ylim = c(0,1.2*max(ss)), ylab = "Signal",main = "Biotin & High Stringency")
error.bar(bhC_Bar,ss, ssSD)
lsC_Bar <- barplot(lss,ylim = c(0,1.2*max(lss)), ylab = "Signal",main = "Low Stringency")
error.bar(lsC_Bar,lss, lssSD)
negC_Bar <- barplot(nn,ylim = c(0,1.2*max(nn)), ylab = "Signal",main = "Negative Controls")
error.bar(negC_Bar,nn, nnSD)
giC_Bar <- barplot(gi, ylim = c(0,1.2*max(gi)), ylab = "Signal",main = "Gene Intensity")
error.bar(giC_Bar,gi, giSD)
lbkC_Bar <- barplot(labk, ylim = c(0,1.2*max(labk)), ylab = "Signal",main = "Labeling & Background")
error.bar(lbkC_Bar,labk, labkSD)
par(mfcol=c(1,1))
boxplot(log2(bruteData$E), las= 2)
path <- getwd()
controlProb  <- read.table(paste(path,"/../Data/ControlProbeprof.txt",sep = ""), sep = "\t",dec = ",", header = TRUE ,as.is = TRUE)
controlGene  <- read.table(paste(path,"/../Data/ControlGeneprof.txt",sep = ""), sep = "\t",dec = ",", header = TRUE ,as.is = TRUE)
controlSample  <- read.table(paste(path,"/../Data/TableSampleControl.txt",sep = ""), sep = "\t",dec = ",", header = TRUE ,as.is = TRUE)
hyb <- controlProb[controlProb$TargetID=="CY3_HYB",]
biot <- controlProb[controlProb$TargetID=="BIOTIN",]
background <- controlProb[controlProb$TargetID == "NEGATIVE",]
pvalue <- function(x){
  ((length(background$X200729890044_A.AVG_Signal[background$X200729890044_A.AVG_Signal == x])/2)+length(background$X200729890044_A.AVG_Signal[background$X200729890044_A.AVG_Signal > x]))/ length(background$X200729890044_A.AVG_Signal)
}
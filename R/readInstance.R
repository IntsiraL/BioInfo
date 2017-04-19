library(limma)
idatFiles= file.choose()
bgxfile=file.choose()
obj<-read.idat(idatFiles, bgxfile, dateinfo=FALSE, tolerance=0)


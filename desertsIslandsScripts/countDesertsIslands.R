options(stringsAsFactors = F)

popList <- read.table("popList_all.txt", header=F)
populations <- popList$V1
regTypes <- c("Deserts", "Islands")

#popName <- "CHPL2021"
#regType <- "Deserts"
for (popName in populations) {
  for (regType in regTypes) {
    regData <- read.table(paste(popName, "_minorParent", regType, "_AIMs_cMwins_merged.txt", sep=""), header=T)
    #print(paste("popName", "regType", "cutoff", "filter", "nRegs", "medcMlen", "meancMlen", "meanPropcM", "medAIMlen", "meanAIMen", "meanPropAIM", "meanDist", "medDist", sep="    "))
    print(paste(popName, regType, nrow(regData), median(regData$cMlen, na.rm = T), mean(regData$cMlen, na.rm = T), mean(regData$cMpropAIMs, na.rm = T), median(regData$aimLen, na.rm = T), mean(regData$aimLen, na.rm = T), mean(regData$propAIMs, na.rm = T), median(regData$distToPrev, na.rm=T), mean(regData$distToPrev, na.rm=T), sep="    "))
  }
}


popList <- read.table("popList_thinned_xcorxbir.txt", header=F)
populations <- popList$V1
regTypes <- c("Deserts", "Islands")

#popName <- "CHPL2021"
#regType <- "Deserts"
for (popName in populations) {
  for (regType in regTypes) {
    regData <- read.table(paste(popName, "_minorParent", regType, "_AIMs_cMwins_merged.txt", sep=""), header=T)
    #print(paste("popName", "regType", "cutoff", "filter", "nRegs", "medcMlen", "meancMlen", "meanPropcM", "medAIMlen", "meanAIMen", "meanPropAIM", "meanDist", "medDist", sep="    "))
    print(paste(popName, regType, nrow(regData), median(regData$cMlen, na.rm = T), mean(regData$cMlen, na.rm = T), mean(regData$cMpropAIMs, na.rm = T), median(regData$aimLen, na.rm = T), mean(regData$aimLen, na.rm = T), mean(regData$propAIMs, na.rm = T), median(regData$distToPrev, na.rm=T), mean(regData$distToPrev, na.rm=T), sep="    "))
  }
}



popList <- read.table("popList_thinned_xmalxbir.txt", header=F)
populations <- popList$V1
regTypes <- c("Deserts", "Islands")

#popName <- "CHPL2021"
#regType <- "Deserts"
for (popName in populations) {
  for (regType in regTypes) {
    regData <- read.table(paste(popName, "_minorParent", regType, "_AIMs_cMwins_merged.txt", sep=""), header=T)
    #print(paste("popName", "regType", "cutoff", "filter", "nRegs", "medcMlen", "meancMlen", "meanPropcM", "medAIMlen", "meanAIMen", "meanPropAIM", "meanDist", "medDist", sep="    "))
    print(paste(popName, regType, nrow(regData), median(regData$cMlen, na.rm = T), mean(regData$cMlen, na.rm = T), mean(regData$cMpropAIMs, na.rm = T), median(regData$aimLen, na.rm = T), mean(regData$aimLen, na.rm = T), mean(regData$propAIMs, na.rm = T), median(regData$distToPrev, na.rm=T), mean(regData$distToPrev, na.rm=T), sep="    "))
  }
}
options(stringAsFactors=FALSE)
options(scipen=999)
library(data.table)

inputParams <- fread("simParams_ABC_2M_initprop_0.5-1_mig.csv", header=F, col.names = c("seed","popsize","gen","init_prop","mig1","mig2"), sep=",")

#popParams <- read.table("allPops_chr2_popParams.txt", header=T)
#popParams <- read.table("pacbio_chr2_popParams.txt", header=T)
popParams <- read.table("pacbio2023_chr2_filteredSegs_popParams.txt", header=T)
uniPops <- unique(popParams$population)
for (popName in uniPops) {
  #popName = "CHPL2006"
  if (file.exists(paste("summary_params_", popName, ".txt", sep=""))) {
    print(popName)
    print(subset(popParams, population==popName))
    popTractLen <- subset(popParams, population==popName)$medianTractLen
    popAnc <- subset(popParams, population==popName)$meanAnc
    popSD <- subset(popParams, population==popName)$ancSD
    if (popAnc<0.5) {
      print(popName)
      popAnc <- 1-popAnc
      print(paste("popAnc=", popAnc))
      popCov <- popSD/popAnc
      print(paste("popCov=", popCov))
    } else {
      popCov <- subset(popParams, population==popName)$ancCvar
    }
    
    #readInName <- paste("summary_params_", popName, ".txt", sep="")
    #popSum <- fread(paste("summary_params_", popName, ".txt", sep=""), header=F, col.names = c("simID", "numberOfTracts", "medianLengthOfMinorTracts", "chrWideAvgAnc", "chrWideVarAnc", "chrWideCoV"), sep=" ")
    popSum <- fread(paste("grep -v '^Fatal' summary_params_", popName, ".txt", sep=""), header=F, col.names = c("simID", "numberOfTracts", "medianLengthOfMinorTracts", "chrWideAvgAnc", "chrWideVarAnc", "chrWideCoV"), sep=" ")
        
    popRev <- popSum[nrow(popSum):1,]
    #Then remove duplicates
    dim(popRev[duplicated(popRev$simID),])
    popRev <- popRev[!duplicated(popRev$simID),]
    #Then put it back accending order
    popOut <-popRev[order(popRev$simID),]
    #fwrite(chplOut, "summary_params_chplLikeUnique_nov2.txt", row.names = F, col.names = F, quote = F, sep = "\t")
    
    popSort <- popOut[order(popOut$simID),]
    print("N sims finished")
    print(dim(popSort))
    #fwrite(popSort, paste("clean_summary_params_", popName, "_", Sys.Date(), ".txt", sep=""), row.names=F, quote = F, sep="\t")
    inputIn <- inputParams[which(inputParams$seed %in% popSort$simID),]
    popFull <- cbind(inputIn, popSort)
    print(paste("First sim = ", min(popFull$seed)))
    print(paste("Last sim = ", max(popFull$seed)))
    #fwrite(popFull, paste("input_summary_params_", popName, "_", Sys.Date(), ".txt", sep=""), row.names=F, quote = F, sep="\t")
    toOutput <- data.frame(popsize=popFull$popsize, gen=popFull$gen, init_prop=popFull$init_prop, mig1=popFull$mig1, mig2=popFull$mig2, medTractLen=popFull$medianLengthOfMinorTracts, meanAnc=popFull$chrWideAvgAnc, ancCvar=popFull$chrWideCoV)
    fwrite(toOutput, paste(popName, "_priors.txt", sep=""), row.names = F, col.names = F, quote = F, sep="\t")
    
    pass.tract <- subset(popFull, (medianLengthOfMinorTracts<(popTractLen*1.05) & medianLengthOfMinorTracts>(popTractLen*0.95)))
    print("N tract pass")
    print(nrow(pass.tract))
    pass.popAnc <- subset(popFull, (chrWideAvgAnc<(popAnc*1.05) & chrWideAvgAnc>(popAnc*0.95)))
    print("N anc pass")
    print(nrow(pass.popAnc))
    #pass.popSD <- subset(popFull, (chrWideVarAnc<(popSD*1.05) & chrWideVarAnc>(popSD*0.95)))
    #print("N SD pass")
    #print(nrow(pass.popSD))
    pass.popCov <- subset(popFull, (chrWideCoV<(popCov*1.05) & chrWideCoV>(popCov*0.95)))
    print("N CoV pass")
    print(nrow(pass.popCov))
    #print("N tract and anc pass")
    #print(nrow(subset(popFull, (medianLengthOfMinorTracts<(popTractLen*1.05) & medianLengthOfMinorTracts>(popTractLen*0.95)) & (chrWideAvgAnc<(popAnc*1.05) & chrWideAvgAnc>(popAnc*0.95)))))
    #print("N tract and cov anc pass")
    #print(nrow(subset(popFull, (medianLengthOfMinorTracts<(popTractLen*1.05) & medianLengthOfMinorTracts>(popTractLen*0.95)) & (chrWideCoV<(popCov*1.05) & chrWideCoV>(popCov*0.95)))))
    #print("N anc and cov anc pass")
    #print(nrow(subset(popFull, (chrWideAvgAnc<(popAnc*1.05) & chrWideAvgAnc>(popAnc*0.95)) & (chrWideCoV<(popCov*1.05) & chrWideCoV>(popCov*0.95)))))
    
    #5%
    pass.Chr2 <- subset(popFull, (medianLengthOfMinorTracts<(popTractLen*1.05) & medianLengthOfMinorTracts>(popTractLen*0.95)) & (chrWideAvgAnc<(popAnc*1.05) & chrWideAvgAnc>(popAnc*0.95)) & (chrWideCoV<(popCov*1.05) & chrWideCoV>(popCov*0.95)))
    print("N 5% pass")
    print(nrow(pass.Chr2))
    if (nrow(pass.Chr2)>0) {
      print(summary(pass.Chr2))
    }
    if (nrow(pass.Chr2)>20) {
      fwrite(pass.Chr2, paste("abcOut_", popName, "_pass3_", Sys.Date(), ".txt", sep=""), row.names = F, quote=F, sep="\t")
    }
  } 
}


passMollyFormat <- pass.Chr2[,which(colnames(pass.Chr2) != "chrWideVarAnc")]

fwrite(passMollyFormat, paste("abcOut_", popName, "_pass3_", Sys.Date(), ".txt", sep=""), row.names = F, col.names = F, quote=F, sep="\t")




simNums <- inputParams$seed

for (i in 101:300) {
  start <- i*1000
  end <- (i+1)*1000
  seeds <- simNums[which(simNums>start & simNums<=end)]
  write.table(seeds, paste("simSeeds/", i, "-", i+1, "K_simSeeds.txt", sep=""), col.names = F, row.names = F, quote = F)
}


##Getting missing seeds
inputParams35 <- subset(inputParams, seed<=35000)
missingSeedsData <- inputParams35[which(!(inputParams35$seed %in% popFull$seed)),]
missingSeedsPopSize <- subset(missingSeedsData, popsize>25)
missingSeeds <- missingSeedsPopSize$seed
#write.table(missingSeeds, "seedsMissing_0-150K.txt", row.names = F, col.names = F, quote = F)



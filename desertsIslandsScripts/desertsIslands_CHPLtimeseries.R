options(stringsAsFactors = F)
library("ggplot2")
library(RColorBrewer)

setwd("~/Documents/SchumerLab/Xcortezi_Xbirchmanni_hybridization/desertsIslandsPacBio2023/")
tranInv <- read.table("../bir_cor_translocations_inversions.txt", header=T)

CHPLyears <- c("CHPL2021", "CHPL2017", "CHPL2006", "CHPL2003")
quants <- read.table("minorParentQuantiles_thinnedAIMs_allPops_0.05cMwins.txt", header=T)
quants$popName <- quants$pop
focalQuants <- subset(quants, pop %in% CHPLyears)

reData <- read.table("xbir_2023pacbio_genome_0.05cm_windows_recRate_codingBP_xcorxbir_minPar_deserts_island.txt", header=T)
aimReg <- read.table("xcorxbir_minPar_deserts_island_2023pacbioXbir.txt", header = T)
regSum <- read.table("xcorxbir_minPar_deserts_island_2023pacbioXbir_summary.txt", header=T)

chplWins <- subset(reData, popName %in% CHPLyears)
write.table(chplWins, "CHPLyears_0.05cm_windows_recRate_codingBP_xcorxbir_minPar_deserts_island.txt", row.names = F, quote=F, sep="\t")
chplRegs <- subset(aimReg, popName %in% CHPLyears)
chplSum <- subset(regSum, pop %in% CHPLyears)

chplSharedLong <- data.frame()
chplSharedShort <- data.frame()
for (chrID in unique(chplRegs$chrName)) {
  print(chrID)
  chrDes <- subset(chplRegs, chrName==chrID & regType=="Desert")
  chrIsle <- subset(chplRegs, chrName==chrID & regType=="Island")
  if (nrow(chrDes)>0) {
    newChrDesertsShort <- data.frame()
    newChrDesertsLong <- data.frame()
    chrDesSort <- chrDes[order(chrDes$cMstart),]
    chrDesSort$rows <- 1:nrow(chrDesSort)
    alreadyMerged <- c()
    for (i in 1:nrow(chrDesSort)) {
      if (!(i %in% alreadyMerged)) {
        toMerge <- chrDesSort[i,]
        alreadyMerged <- c(alreadyMerged, i)
        for (j in (i+1):nrow(chrDesSort)) {
          if (j<=nrow(chrDesSort)) {
            cMmin <- min(toMerge$cMstart)
            cMmax <- max(toMerge$cMend)
            curReg <- chrDesSort[j,]
            startMin <- curReg$cMstart-cMmin
            startMax <- curReg$cMstart-cMmax
            endMin <- curReg$cMend-cMmin
            endMax <- curReg$cMend-cMmin
            if (startMin>=0 & startMax<=0 & endMin>=0 & endMax<=0) {
              toMerge <- rbind(toMerge, curReg)
              alreadyMerged <- c(alreadyMerged, j)
            } else if (startMin>=0 & startMax<=0 & endMin>=0 & endMax>=0) {
              toMerge <- rbind(toMerge, curReg)
              alreadyMerged <- c(alreadyMerged, j)
            }
          }
        }
        cMmin <- min(toMerge$cMstart)
        cMmax <- max(toMerge$cMend)
        regWins <- subset(chplWins, chr==chrID & start>=cMmin & end<=cMmax)
        tmpDF <- data.frame(regionType="Desert", chrName=chrID, cMstart=cMmin, cMend=cMmax, nSNPs=sum(regWins$nSNPs, na.rm=T), nAIMs=sum(regWins$numAIMs, na.rm=T), nWins=length(unique(regWins$start)))
        tmpShort <- data.frame(regionType="Desert", chrName=chrID, cMstart=cMmin, cMend=cMmax, nSNPs=sum(regWins$nSNPs, na.rm=T), nAIMs=sum(regWins$numAIMs, na.rm=T), nWins=length(unique(regWins$start)))
        shortHeader <- colnames(tmpShort)
        tmpLong <- data.frame()
        for (popID in CHPLyears) {
          popWins <- subset(regWins, popName==popID)
          nDes <- nrow(subset(popWins, regType=="Desert"))
          nIsle <- nrow(subset(popWins, regType=="Island"))
          meanAnc <- mean(popWins$minParAnc, na.rm=T)
          tmpShort$nRegions <- nDes
          #tmpShort$nDesert <- nDes
          #tmpShort$nIsle <- nIsle
          tmpShort$minParAnc <- meanAnc
          shortHeader <- c(shortHeader, paste(popID, "nRegions", sep=""), paste(popID, "anc", sep=""))
          colnames(tmpShort) <- shortHeader
          tmpDF$popName <- popID
          tmpDF$nRegions <- nDes
          #tmpDF$nDeserts <- nDes
          #tmpDF$nIslands <- nIsle
          tmpDF$minParAnc <- meanAnc
          tmpLong <- rbind(tmpLong, tmpDF)
        }
        newChrDesertsShort <- rbind(newChrDesertsShort, tmpShort)
        newChrDesertsLong <- rbind(newChrDesertsLong, tmpLong)
      }
    }
    chplSharedLong <- rbind(chplSharedLong, newChrDesertsLong)
    chplSharedShort <- rbind(chplSharedShort, newChrDesertsShort)
  }
  if (nrow(chrIsle)>0) {
    newChrIslandsShort <- data.frame()
    newChrIslandsLong <- data.frame()
    chrIsleSort <- chrIsle[order(chrIsle$cMstart),]
    chrIsleSort$rows <- 1:nrow(chrIsleSort)
    alreadyMerged <- c()
    for (i in 1:nrow(chrIsleSort)) {
      if (!(i %in% alreadyMerged)) {
        toMerge <- chrIsleSort[i,]
        alreadyMerged <- c(alreadyMerged, i)
        for (j in (i+1):nrow(chrIsleSort)) {
          if (j<=nrow(chrIsleSort)) {
            cMmin <- min(toMerge$cMstart)
            cMmax <- max(toMerge$cMend)
            curReg <- chrIsleSort[j,]
            startMin <- curReg$cMstart-cMmin
            startMax <- curReg$cMstart-cMmax
            endMin <- curReg$cMend-cMmin
            endMax <- curReg$cMend-cMmin
            if (startMin>=0 & startMax<=0 & endMin>=0 & endMax<=0) {
              toMerge <- rbind(toMerge, curReg)
              alreadyMerged <- c(alreadyMerged, j)
            } else if (startMin>=0 & startMax<=0 & endMin>=0 & endMax>=0) {
              toMerge <- rbind(toMerge, curReg)
              alreadyMerged <- c(alreadyMerged, j)
            }
          }
        }
        cMmin <- min(toMerge$cMstart)
        cMmax <- max(toMerge$cMend)
        regWins <- subset(chplWins, chr==chrID & start>=cMmin & end<=cMmax)
        tmpDF <- data.frame(regionType="Island", chrName=chrID, cMstart=cMmin, cMend=cMmax,  nSNPs=sum(regWins$nSNPs, na.rm=T), nAIMs=sum(regWins$numAIMs, na.rm=T), nWins=length(unique(regWins$start)))
        tmpShort <- data.frame(regionType="Island", chrName=chrID, cMstart=cMmin, cMend=cMmax,  nSNPs=sum(regWins$nSNPs, na.rm=T), nAIMs=sum(regWins$numAIMs, na.rm=T), nWins=length(unique(regWins$start)))
        shortHeader <- colnames(tmpShort)
        tmpLong <- data.frame()
        for (popID in CHPLyears) {
          popWins <- subset(regWins, popName==popID)
          nDes <- nrow(subset(popWins, regType=="Desert"))
          nIsle <- nrow(subset(popWins, regType=="Island"))
          meanAnc <- mean(popWins$minParAnc, na.rm=T)
          tmpShort$nRegions <- nIsle
          #tmpShort$nDesert <- nDes
          #tmpShort$nIsle <- nIsle
          tmpShort$minParAnc <- meanAnc
          shortHeader <- c(shortHeader, paste(popID, "nRegions", sep=""), paste(popID, "anc", sep=""))
          colnames(tmpShort) <- shortHeader
          tmpDF$popName <- popID
          tmpDF$nRegions <- nIsle
          #tmpDF$nDeserts <- nDes
          #tmpDF$nIslands <- nIsle
          tmpDF$minParAnc <- meanAnc
          tmpLong <- rbind(tmpLong, tmpDF)
        }
        newChrIslandsShort <- rbind(newChrIslandsShort, tmpShort)
        newChrIslandsLong <- rbind(newChrIslandsLong, tmpLong)
      }
    }
    chplSharedLong <- rbind(chplSharedLong, newChrIslandsLong)
    chplSharedShort <- rbind(chplSharedShort, newChrIslandsShort)
  }
}

chplSharedLong$propWins <- chplSharedLong$nRegions/chplSharedLong$nWins

write.table(chplSharedShort, "shared_CHPLregions_timeseries_short.txt", row.names = F, quote = F, sep="\t")
write.table(chplSharedLong, "shared_CHPLregions_timeseries_long.txt", row.names = F, quote = F, sep="\t")

dim(chplRegs)
dim(chplSharedShort)
summary(chplSharedShort)
nrow(subset(chplSharedShort, chrName=="chr-20"))
summary(subset(chplSharedShort, regionType=="Island"))
regFillScale <- scale_fill_manual(values=c("Desert"="tan", "Island"="lightpink"))

chrID <- "chr-24"
pdf("CHPLsharedRegions.pdf")
for (chrID in sort(unique(chplRegs$chrName))) {
  chrWinData <- subset(chplWins, chr==chrID)
  chrReg <- subset(chplSharedLong, chrName==chrID)
  
  print(ggplot() + 
    geom_rect(data=chrReg, aes(xmin=cMstart, xmax=cMend, ymin=0, ymax=0.99, fill=regionType, alpha=propWins)) + 
    geom_rect(data=chrWinData, aes(xmin=start, xmax=end, ymin=minParAnc-0.01, ymax=minParAnc+0.01, color=popName), fill=NA, linewidth=1) + 
    geom_hline(data=focalQuants, aes(yintercept=meanMinPar), color="grey40", linetype="dashed") + 
    geom_hline(data=focalQuants, aes(yintercept=lower10), color="grey20", linetype="dotted") + 
    geom_hline(data=focalQuants, aes(yintercept=upper10), color="grey20", linetype="dotted") + 
    theme_bw() + 
    ylab("Minor Parent Ancestry") + 
    geom_hline(yintercept=0, color="grey20") + regFillScale + xlab(chrID) + 
    facet_wrap(~popName, ncol=1, strip.position="right"))
}
dev.off()


chplSharedShort$regLen <- chplSharedShort$cMend-chplSharedShort$cMstart
chplSharedShort$distToNext <- NA
for (i in 1:nrow(chplSharedShort)) {
  curReg <- chplSharedShort[i,]
  nextReg <- chplSharedShort[i+1,]
  if (curReg$chrName==nextReg$chrName & curReg$regionType==nextReg$regionType) {
    chplSharedShort$distToNext[i] <- nextReg$cMstart-curReg$cMend
  }
}

chrID <- "chr-13"
chrWinData <- subset(chplWins, chr==chrID)
chrReg <- subset(chplSharedLong, chrName==chrID)

print(ggplot() + 
        geom_rect(data=chrReg, aes(xmin=cMstart, xmax=cMend, ymin=0, ymax=0.99, fill=regionType, alpha=propWins)) + 
        geom_rect(data=chrWinData, aes(xmin=start, xmax=end, ymin=minParAnc-0.01, ymax=minParAnc+0.01, color=popName), fill=NA, linewidth=1) + 
        geom_hline(data=focalQuants, aes(yintercept=meanMinPar), color="grey40", linetype="dashed") + 
        geom_hline(data=focalQuants, aes(yintercept=lower10), color="grey20", linetype="dotted") + 
        geom_hline(data=focalQuants, aes(yintercept=upper10), color="grey20", linetype="dotted") + 
        theme_bw() + 
        ylab("Minor Parent Ancestry") + 
        geom_hline(yintercept=0, color="grey20") + regFillScale + xlab(chrID) + 
        facet_wrap(~popName, ncol=1, strip.position="right") )


##Merge any within 50000
##Filter less than 10000 
##Filter less than 10 SNPs and 10 AIM

head(chplSharedShort)
head(chplSharedLong)

chplSharedMerged <- data.frame()
for (chrID in sort(unique(chplRegs$chrName))) {
  print(chrID)
  chrDes <- subset(chplSharedShort, chrName==chrID & regionType=="Desert")
  chrIsles <- subset(chplSharedShort, chrName==chrID & regionType=="Island")
  print(summary(chrDes$distToNext))
  if (nrow(chrDes)>0 & min(chrDes$distToNext, na.rm=T)<50000) {
    chrDes$toMerge <- chrDes$distToNext<50000
    chrDes$row <- 1:nrow(chrDes)
    desNew <- data.frame()
    i<-1
    while (i <= nrow(chrDes)) {
      curReg <- chrDes[i,]
      if (curReg$toMerge==F | is.na(curReg$toMerge)) {
        desNew <- rbind(desNew, curReg[1:4])
        i <- i+1
      } else if (curReg$toMerge==T) {
        mergeDF <- curReg
        j <- i+1
        while (j <= nrow(chrDes)) {
          nextReg <- chrDes[j,]
          mergeDF <- rbind(mergeDF, nextReg)
          if (nextReg$toMerge==F | is.na(nextReg$toMerge)) {
            newMerge <- data.frame(regionType="Desert", chrName=chrID, cMstart=min(mergeDF$cMstart), cMend=max(mergeDF$cMend))
            desNew <- rbind(desNew, newMerge)
            i <- j +1
            j <- nrow(chrDes) + 1
          } else if (nextReg$toMerge==T) {
            j <- j + 1
          }
        }
      }
    }
  } else if (nrow(chrDes)>0) {
    desNew <- chrDes[,1:4]
  } else {
    print(chrID)
  }
  if (min(chrIsles$distToNext, na.rm=T)<50000) {
    chrIsles$toMerge <- chrIsles$distToNext<50000
    chrIsles$row <- 1:nrow(chrIsles)
    isleNew <- data.frame()
    i<-1
    while (i <= nrow(chrIsles)) {
      curReg <- chrIsles[i,]
      if (curReg$toMerge==F | is.na(curReg$toMerge)) {
        isleNew <- rbind(isleNew, curReg[1:4])
        i <- i+1
      } else if (curReg$toMerge==T) {
        mergeDF <- curReg
        j <- i+1
        while (j <= nrow(chrIsles)) {
          nextReg <- chrIsles[j,]
          mergeDF <- rbind(mergeDF, nextReg)
          if (nextReg$toMerge==F | is.na(nextReg$toMerge)) {
            newMerge <- data.frame(regionType="Island", chrName=chrID, cMstart=min(mergeDF$cMstart), cMend=max(mergeDF$cMend))
            isleNew <- rbind(isleNew, newMerge)
            i <- j +1
            j <- nrow(chrIsles) + 1
          } else if (nextReg$toMerge==T) {
            j <- j + 1
          }
        }
      }
    }
  } else {
    isleNew <- chrIsles[,1:4]
  }
  chplSharedMerged <- rbind(chplSharedMerged, desNew, isleNew)
}

chplSharedMerged$regLen <- chplSharedMerged$cMend-chplSharedMerged$cMstart
chplSharedMerged$distToNext <- NA
for (i in 1:nrow(chplSharedMerged)) {
  curReg <- chplSharedMerged[i,]
  nextReg <- chplSharedMerged[i+1,]
  if (curReg$chrName==nextReg$chrName & curReg$regionType==nextReg$regionType) {
    chplSharedMerged$distToNext[i] <- nextReg$cMstart-curReg$cMend
  }
}

nrow(subset(chplSharedMerged, distToNext<50000))
subset(chplSharedMerged, regLen<10000)

chplSharedMerged <- subset(chplSharedMerged, regLen>10000)

chplSharedMergedShort <- data.frame()
chplSharedMergedLong <- data.frame()
for (i in 1:nrow(chplSharedMerged)) {
  curReg <- chplSharedMerged[i,]
  curType <- curReg$regionType
  chrID <- curReg$chrName
  cMmin <- curReg$cMstart
  cMmax <- curReg$cMend
  cMmidP <- cMmin + ((cMmax-cMmin)/2)
  regWins <- subset(chplWins, chr==chrID & start>=cMmin & end<=cMmax)
  tmpDF <- data.frame(regionType=curType, chrName=chrID, cMstart=cMmin, cMend=cMmax, regLen=curReg$regLen, distToNext=curReg$distToNext, nSNPs=sum(regWins$nSNPs, na.rm=T), nAIMs=sum(regWins$numAIMs, na.rm=T), nWins=length(unique(regWins$start)))
  tmpShort <- data.frame(regionType=curType, chrName=chrID, cMstart=cMmin, cMend=cMmax, regLen=curReg$regLen, distToNext=curReg$distToNext, nSNPs=sum(regWins$nSNPs, na.rm=T), nAIMs=sum(regWins$numAIMs, na.rm=T), nWins=length(unique(regWins$start)))
  shortHeader <- colnames(tmpShort)
  tmpLong <- data.frame()
  for (popID in CHPLyears) {
    popWins <- subset(regWins, popName==popID)
    nRegs <- nrow(subset(popWins, regType==curType))
    meanAnc <- mean(popWins$minParAnc, na.rm=T)
    midWinAnc <- subset(popWins, start<cMmidP & end>cMmidP)$minParAnc
    midPass <- NA
    if (nRegs>0) {
      popDes <- subset(focalQuants, pop==popID)$lower10
      popIsle <- subset(focalQuants, pop==popID)$upper10
      if (curType=="Desert") {
        if (is.na(midWinAnc)) {midPass <- NA}
        else if (midWinAnc<popDes) {midPass <- T}
        else {midPass <- F}
      }
      if (curType=="Island") {
        if (is.na(midWinAnc)) {midPass <- NA}
        else if (midWinAnc>popIsle) {midPass <- T}
        else {midPass <- F}
      }
    }
    tmpShort$nRegions <- nRegs
    tmpShort$minParAnc <- meanAnc
    tmpShort$midPass <- midPass
    shortHeader <- c(shortHeader, paste(popID, "nRegions", sep=""), paste(popID, "anc", sep=""), paste(popID, "cMmidPass", sep=""))
    colnames(tmpShort) <- shortHeader
    tmpDF$popName <- popID
    tmpDF$nRegions <- nRegs
    tmpDF$minParAnc <- meanAnc
    tmpDF$midPass <- midPass
    tmpLong <- rbind(tmpLong, tmpDF)
  }
  chplSharedMergedShort <- rbind(chplSharedMergedShort, tmpShort)
  chplSharedMergedLong <- rbind(chplSharedMergedLong, tmpLong)
}

##Filter less than 10 SNPs and 10 AIM
#something's wrong with the nAIMs because there are windows with no AIMs but sill an average ancestry...
#nrow(subset(chplSharedMergedShort, nAIMs<10))
#subset(chplSharedMergedShort, nAIMs<10)

chplSharedMergedLong$propWins <- chplSharedMergedLong$nRegions/chplSharedMergedLong$nWins

write.table(chplSharedMergedShort, "shared_CHPLregionsMerged_timeseries_short.txt", row.names = F, quote = F, sep="\t")
write.table(chplSharedMergedLong, "shared_CHPLregionsMerged_timeseries_long.txt", row.names = F, quote = F, sep="\t")


chrID <- "chr-06"
pdf("CHPLsharedRegionsMerged.pdf")
for (chrID in sort(unique(chplRegs$chrName))) {
  chrWinData <- subset(chplWins, chr==chrID)
  chrReg <- subset(chplSharedMergedLong, chrName==chrID)
  
  print(ggplot() + 
          geom_rect(data=chrReg, aes(xmin=cMstart, xmax=cMend, ymin=0, ymax=0.99, fill=regionType, alpha=propWins)) + 
          geom_rect(data=chrWinData, aes(xmin=start, xmax=end, ymin=minParAnc-0.01, ymax=minParAnc+0.01, color=popName), fill=NA, linewidth=1) + 
          geom_hline(data=focalQuants, aes(yintercept=meanMinPar), color="grey40", linetype="dashed") + 
          geom_hline(data=focalQuants, aes(yintercept=lower10), color="grey20", linetype="dotted") + 
          geom_hline(data=focalQuants, aes(yintercept=upper10), color="grey20", linetype="dotted") + 
          theme_bw() + 
          ylab("Minor Parent Ancestry") + 
          geom_hline(yintercept=0, color="grey20") + regFillScale + xlab(chrID) + 
          facet_wrap(~popName, ncol=1, strip.position="right"))
}
dev.off()

chplSharedMergedShort$cMpassAll <- NA
chplSharedMergedShort$nShared <- NA
for (i in 1:nrow(chplSharedMergedShort)) {
  curReg <- chplSharedMergedShort[i,]
  cMpass <- sum(curReg$CHPL2021cMmidPass, curReg$CHPL2017cMmidPass, curReg$CHPL2006cMmidPass, curReg$CHPL2003cMmidPass, na.rm=T)
  chplSharedMergedShort$cMpassAll[i] <- cMpass
  nShared <- 0
  if (curReg$CHPL2021nRegions>0) {nShared <- nShared + 1}
  if (curReg$CHPL2017nRegions>0) {nShared <- nShared + 1}
  if (curReg$CHPL2006nRegions>0) {nShared <- nShared + 1}
  if (curReg$CHPL2003nRegions>0) {nShared <- nShared + 1}
  chplSharedMergedShort$nShared[i] <- nShared
}

write.table(chplSharedMergedShort, "shared_CHPLregionsMerged_timeseries_short.txt", row.names = F, quote = F, sep="\t")

nrow(subset(chplSharedMergedShort, nShared==4))
nrow(subset(chplSharedMergedShort, nShared==3))
nrow(subset(chplSharedMergedShort, nShared==2))
nrow(subset(chplSharedMergedShort, nShared==1))

nrow(subset(chplSharedMergedShort, cMpassAll==0))
nrow(subset(chplSharedMergedShort, cMpassAll==4))

nrow(subset(chplSharedMergedShort, nShared==4 & regionType=="Desert"))
nrow(subset(chplSharedMergedShort, nShared==4 & regionType=="Island"))

ggplot(chplSharedMergedShort) + geom_histogram(aes(x=regLen, fill=as.character(nShared)), binwidth=10000) + facet_wrap(~regionType, nrow=2) + theme_bw()
ggplot(chplSharedMergedShort) + geom_violin(aes(y=regLen, x=as.character(nShared))) + facet_wrap(~regionType, nrow=2) + theme_bw()


head(chplSharedMergedLong)
chplSharedMergedLong$nShared <- NA
for (i in 1:nrow(chplSharedMergedShort)) {
  chplSharedMergedLong$nShared[which(chplSharedMergedLong$chrName==chplSharedMergedShort$chrName[i] & chplSharedMergedLong$cMstart==chplSharedMergedShort$cMstart[i])] <- chplSharedMergedShort$nShared[i]
}

chplSharedMergedShort$change2003.2006 <- chplSharedMergedShort$CHPL2006anc - chplSharedMergedShort$CHPL2003anc
chplSharedMergedShort$change2006.2017 <- chplSharedMergedShort$CHPL2017anc - chplSharedMergedShort$CHPL2006anc
chplSharedMergedShort$change2017.2021 <- chplSharedMergedShort$CHPL2021anc - chplSharedMergedShort$CHPL2017anc

summary(subset(chplSharedMergedShort, regionType=="Desert"))
summary(subset(chplSharedMergedShort, regionType=="Island"))

write.table(chplSharedMergedShort, "shared_CHPLregionsMerged_timeseries_short.txt", row.names = F, quote = F, sep="\t")

chplSharedMergedLong$quant <- "mid"
chplSharedMergedLong$quant[which(chplSharedMergedLong$nRegions>0)] <- chplSharedMergedLong$regionType[which(chplSharedMergedLong$nRegions>0)]

chplSharedMergedLong$year <- NA
chplSharedMergedLong$year[which(chplSharedMergedLong$popName=="CHPL2003")] <- 2003
chplSharedMergedLong$year[which(chplSharedMergedLong$popName=="CHPL2006")] <- 2006
chplSharedMergedLong$year[which(chplSharedMergedLong$popName=="CHPL2017")] <- 2017
chplSharedMergedLong$year[which(chplSharedMergedLong$popName=="CHPL2021")] <- 2021

write.table(chplSharedMergedLong, "shared_CHPLregionsMerged_timeseries_long.txt", row.names = F, quote = F, sep="\t")

max(chplSharedMergedLong$nWins)

chrData <- subset(chplSharedMergedLong, chrName=="chr-01")

ggplot(chrData, aes(x=year, y=minParAnc)) + geom_line() + geom_point(aes(color=quant)) + ylim(c(0,1)) + theme_bw() + facet_wrap(regionType~cMstart)

pdf("CHPLsharedRegionsTrajectories.pdf")
for (chrID in sort(unique(chplRegs$chrName))) {
  chrData <- subset(chplSharedMergedLong, chrName==chrID)
  print(ggplot(chrData, aes(x=year, y=minParAnc)) + geom_line() + geom_point(aes(color=quant)) + ylim(c(0,1)) + theme_bw() + facet_wrap(regionType~regLen) + ggtitle(chrID))
}  
dev.off()



chplSharedMergedLong$regID <- paste(chplSharedMergedLong$chrName, chplSharedMergedLong$cMstart, sep="_")
chplSharedMergedLong$diffFrom2003 <- NA
chplSharedMergedLong$diffFrom2021 <- NA
for (region in unique(chplSharedMergedLong$regID)) {
  curReg <- subset(chplSharedMergedLong, regID==region)
  reg2021 <- subset(curReg, popName=="CHPL2021")$minParAnc
  reg2003 <- subset(curReg, popName=="CHPL2003")$minParAnc
  for (yearID in unique(curReg$year)) {
    chplSharedMergedLong$diffFrom2003[which(chplSharedMergedLong$regID==region & chplSharedMergedLong$year==yearID)] <- reg2003 - curReg$minParAnc[which(curReg$year==yearID)]
    chplSharedMergedLong$diffFrom2021[which(chplSharedMergedLong$regID==region & chplSharedMergedLong$year==yearID)] <- reg2021 - curReg$minParAnc[which(curReg$year==yearID)]
  }
}

chplSharedMergedLong <- read.table("shared_CHPLregionsMerged_timeseries_wPercentile.txt", header=T) 

ggplot(chplSharedMergedLong, aes(x=year, y=minParAnc)) + geom_path(alpha=0.6, color="grey20", aes(group=regID)) + geom_point(aes(color=quant, size=regLen)) + facet_grid(chrName~regionType, scales = "free_y") + theme_bw() 
ggsave("desertIslandsTrajectoryByYearByChr.pdf", height=24, width=8)

ggplot(chplSharedMergedLong, aes(x=year, y=minParAnc)) + geom_path(alpha=0.6, color="grey20", aes(alpha=nShared, group=regID)) + geom_point(aes(color=quant, alpha=nShared, size=nShared)) + facet_grid(~regionType, scales = "free_y") + theme_bw() 
ggsave("desertIslandsTrajectoryByYear.pdf", height=12, width=12)

ggplot(chplSharedMergedLong, aes(x=year, y=minParAnc)) + geom_path(color="grey20", aes(alpha=regLen, group=regID)) + geom_point(aes(color=quant, alpha=regLen, size=nWins)) + facet_grid(nShared~regionType) + theme_bw() 
ggsave("desertIslandsTrajectoryByYearByNumberShared.pdf", height=18, width=9)


ggplot(chplSharedMergedLong, aes(x=year, y=diffFrom2003)) + geom_path(alpha=0.6, color="grey20", aes(group=regID)) + geom_point(aes(color=quant, alpha=0.8)) + facet_grid(nShared~regionType) + theme_bw() 

ggplot(chplSharedMergedLong, aes(x=year, y=diffFrom2021)) + geom_path(alpha=0.6, color="grey20", aes(group=regID)) + geom_point(aes(color=quant, alpha=0.8)) + facet_grid(nShared~regionType) + theme_bw() 


ecdf(subset(chplWins, popName=="CHPL2021")$minparAnc)(0.136608658)

chpl2021anc <- na.omit(subset(chplWins, popName=="CHPL2021")$minParAnc)
chpl2021regs <- subset(chplSharedMergedLong, popName=="CHPL2021")
chpl2021regs$percentile <- ecdf(chpl2021anc)(chpl2021regs$minParAnc)

chpl2017anc <- na.omit(subset(chplWins, popName=="CHPL2017")$minParAnc)
chpl2017regs <- subset(chplSharedMergedLong, popName=="CHPL2017")
chpl2017regs$percentile <- ecdf(chpl2017anc)(chpl2017regs$minParAnc)

chpl2006anc <- na.omit(subset(chplWins, popName=="CHPL2006")$minParAnc)
chpl2006regs <- subset(chplSharedMergedLong, popName=="CHPL2006")
chpl2006regs$percentile <- ecdf(chpl2006anc)(chpl2006regs$minParAnc)

chpl2003anc <- na.omit(subset(chplWins, popName=="CHPL2003")$minParAnc)
chpl2003regs <- subset(chplSharedMergedLong, popName=="CHPL2003")
chpl2003regs$percentile <- ecdf(chpl2003anc)(chpl2003regs$minParAnc)

chplSharedPercByYear <- rbind(chpl2021regs, chpl2017regs, chpl2006regs, chpl2003regs)

chplSharedPerc <- data.frame()
chplSharedPercFoundIn <- data.frame()
for(regName in unique(chplSharedPercByYear$regID)) {
  regData <- subset(chplSharedPercByYear, regID==regName) 
  chplSharedPerc <- rbind(chplSharedPerc, regData)
  if (regData$nRegions[which(regData$popName=="CHPL2003")]>0) {
    regTmp <- regData
    regTmp$foundIn <- "CHPL2003"
    chplSharedPercFoundIn <- rbind(chplSharedPercFoundIn, regTmp)
  }
  if (regData$nRegions[which(regData$popName=="CHPL2006")]>0) {
    regTmp <- regData
    regTmp$foundIn <- "CHPL2006"
    chplSharedPercFoundIn <- rbind(chplSharedPercFoundIn, regTmp)
  }
  if (regData$nRegions[which(regData$popName=="CHPL2017")]>0) {
    regTmp <- regData
    regTmp$foundIn <- "CHPL2017"
    chplSharedPercFoundIn <- rbind(chplSharedPercFoundIn, regTmp)
  }
  if (regData$nRegions[which(regData$popName=="CHPL2021")]>0) {
    regTmp <- regData
    regTmp$foundIn <- "CHPL2021"
    chplSharedPercFoundIn <- rbind(chplSharedPercFoundIn, regTmp)
  }
}

chplSharedPercFoundIn <- read.table("shared_CHPLregionsMerged_timeseries_byYearFoundIn.txt", header=T)


ggplot(chplSharedPerc, aes(x=year, y=percentile)) + geom_path(color="grey20", aes(alpha=regLen, group=regID)) + geom_point(aes(color=quant, alpha=regLen, size=nWins)) + facet_grid(nShared~regionType) + theme_bw() 
ggsave("desertIslandsPercentilTrajectoryByYearByNumberShared.pdf", height=12, width=8)

ggplot(chplSharedPercFoundIn, aes(x=year, y=percentile)) + geom_path(color="grey20", aes(alpha=nWins, group=regID)) + geom_point(aes(color=quant, alpha=propWins, size=nShared)) + facet_grid(foundIn~regionType) + theme_bw() 
ggsave("desertIslandsPercentilTrajectoryByYearByYearFoundIn.pdf", height=12, width=8)

ggplot(chplSharedPercFoundIn, aes(x=year, y=minParAnc)) + geom_path(color="grey20", aes(alpha=nWins, group=regID)) + geom_point(aes(color=quant, alpha=propWins, size=nShared)) + facet_grid(foundIn~regionType) + theme_bw() 
ggsave("desertIslandsMinParAncTrajectoryByYearByYearFoundIn.pdf", height=12, width=8)

write.table(chplSharedPerc, "shared_CHPLregionsMerged_timeseries_wPercentile.txt", row.names = F, quote = F, sep="\t")
write.table(chplSharedPercFoundIn, "shared_CHPLregionsMerged_timeseries_byYearFoundIn.txt", row.names = F, quote = F, sep="\t")

chplSharedPerc <- read.table("shared_CHPLregionsMerged_timeseries_wPercentile.txt", header=T)
head(chplSharedPerc)
summary(subset(chplSharedPerc, regionType=="Island"))

justIslands <- subset(chplSharedPerc, regionType=="Island")
islandsLM <- lm(minParAnc~year, justIslands)
summary(islandsLM)
lmdfIslands <- data.frame()
for (i in 1:length(unique(justIslands$regID))) {
  curRegID <- unique(justIslands$regID)[i]
  curRegData <- subset(justIslands, regID==curRegID)
  islandLM <- lm(minParAnc~year, curRegData)
  #print(curRegID)
  #print(summary(islandLM))
  tmpDF <- data.frame(regID = curRegID, lmEstYear = summary(islandLM)$coefficients[,1][2], lmPvalueYear = summary(islandLM)$coefficients[,4][2])
  lmdfIslands <- rbind(lmdfIslands, tmpDF)
}

sigIslandRun <- subset(lmdfIslands, lmPvalueYear<0.05)
sigIslandIDs <- sigIslandRun$regID
sigIslands <- subset(justIslands, regID%in%sigIslandIDs)

ggplot(sigIslands, aes(x=year, y=minParAnc)) + geom_path(color="grey20", aes(alpha=nWins, group=regID)) + geom_point(aes(color=quant, alpha=propWins, size=nShared)) + theme_bw() 
ggplot(sigIslands, aes(x=year, y=percentile)) + geom_path(color="grey20", aes(alpha=nWins, group=regID)) + geom_point(aes(color=quant, alpha=propWins, size=nShared)) + theme_bw() 

negIslandID <- sigIslands$regID[sigIslands$minParAnc<0.6]
negIsland <- subset(sigIslands, regID==negIslandID)
ggplot(negIsland, aes(x=year, y=minParAnc)) + geom_path(color="grey20", aes(alpha=nWins, group=regID)) + geom_point(aes(color=quant, alpha=propWins, size=nShared)) + theme_bw() 

posIslands <- subset(sigIslands, !(regID %in% negIslandID))
ggplot(posIslands, aes(x=year, y=minParAnc)) + geom_path(color="grey20", aes(alpha=nWins, group=regID)) + geom_point(aes(color=quant, alpha=propWins, size=nShared)) + theme_bw() 


justDeserts <- subset(chplSharedPerc, regionType=="Desert")
desertsLM <- lm(minParAnc~year, justDeserts)
summary(desertsLM)
lmdfDeserts <- data.frame()
for (i in 1:length(unique(justDeserts$regID))) {
  curRegID <- unique(justDeserts$regID)[i]
  curRegData <- subset(justDeserts, regID==curRegID)
  desertLM <- lm(minParAnc~year, curRegData)
  #print(curRegID)
  #print(summary(desertLM))
  tmpDF <- data.frame(regID = curRegID, lmEstYear = summary(desertLM)$coefficients[,1][2], lmPvalueYear = summary(desertLM)$coefficients[,4][2])
  lmdfDeserts <- rbind(lmdfDeserts, tmpDF)
}

sigDesertLM <- subset(lmdfDeserts, lmPvalueYear<0.05)
sigDesertIDs <- sigDesertLM$regID
sigDeserts <- subset(justDeserts, regID%in%sigDesertIDs)

ggplot(sigDeserts, aes(x=year, y=minParAnc)) + geom_path(color="grey20", aes(alpha=nWins, group=regID)) + geom_point(aes(color=quant, alpha=propWins, size=nShared)) + theme_bw() 
ggplot(sigDeserts, aes(x=year, y=percentile)) + geom_path(color="grey20", aes(alpha=nWins, group=regID)) + geom_point(aes(color=quant, alpha=propWins, size=nShared)) + theme_bw() 

negDesertIDs <- sigDesertLM$regID[which(sigDesertLM$lmEstYear<0)]
negDesert <- subset(justDeserts, regID%in%negDesertIDs)
ggplot(negDesert, aes(x=year, y=minParAnc)) + geom_path(color="grey20", aes(alpha=nWins, group=regID)) + geom_point(aes(color=quant, alpha=propWins, size=nShared)) + theme_bw() 

posDeserts <- subset(sigDeserts, !(regID %in% negDesertIDs))
ggplot(posDeserts, aes(x=year, y=minParAnc)) + geom_path(color="grey20", aes(alpha=nWins, group=regID)) + geom_point(aes(color=quant, alpha=propWins, size=nShared)) + theme_bw() 

dim(posIslands)
unique(posIslands$chrName)


chplSharedMergedShort <- read.table("shared_CHPLregionsMerged_timeseries_short.txt", header=T)
chplSharedMergedShort$regID <- paste(chplSharedMergedShort$chrName, chplSharedMergedShort$cMstart, sep="_")
sigIslandShort <- subset(chplSharedMergedShort, regID %in% unique(posIslands$regID))

write.table(sigIslandShort, "CHPL_significantIslands_timeseries_short.txt", row.names = F, quote = F, sep="\t")

pdf("CHPL_significantIslandsByChr.pdf")
for (chrID in sort(unique(posIslands$chrName))) {
  chrWinData <- subset(chplWins, chr==chrID)
  chrReg <- subset(posIslands, chrName==chrID)
  
  print(ggplot() + 
          geom_rect(data=chrReg, aes(xmin=cMstart, xmax=cMend, ymin=0, ymax=0.99, fill=regionType)) + 
          geom_rect(data=chrWinData, aes(xmin=start, xmax=end, ymin=minParAnc-0.01, ymax=minParAnc+0.01, color=popName), fill=NA, linewidth=1) + 
          geom_hline(data=focalQuants, aes(yintercept=meanMinPar), color="grey40", linetype="dashed") + 
          geom_hline(data=focalQuants, aes(yintercept=lower10), color="grey20", linetype="dotted") + 
          geom_hline(data=focalQuants, aes(yintercept=upper10), color="grey20", linetype="dotted") + 
          theme_bw() + 
          ylab("Minor Parent Ancestry") + 
          geom_hline(yintercept=0, color="grey20") + regFillScale + xlab(chrID) + 
          facet_wrap(~popName, ncol=1, strip.position="right"))
}
dev.off()



head(chplSharedMergedShort)
nrow(subset(chplSharedMergedShort, regionType=="Desert"))
nrow(subset(chplSharedMergedShort, regionType=="Desert" & CHPL2021nRegions>0))
nrow(subset(chplSharedMergedShort, regionType=="Desert" & CHPL2017nRegions>0))
nrow(subset(chplSharedMergedShort, regionType=="Desert" & CHPL2006nRegions>0))
nrow(subset(chplSharedMergedShort, regionType=="Desert" & CHPL2003nRegions>0))

nrow(subset(chplSharedMergedShort, regionType=="Island"))
nrow(subset(chplSharedMergedShort, regionType=="Island" & CHPL2021nRegions>0))
nrow(subset(chplSharedMergedShort, regionType=="Island" & CHPL2017nRegions>0))
nrow(subset(chplSharedMergedShort, regionType=="Island" & CHPL2006nRegions>0))
nrow(subset(chplSharedMergedShort, regionType=="Island" & CHPL2003nRegions>0))

ggplot(chplSharedMergedShort) + geom_violin(aes(x=as.character(nShared), y=regLen, color=as.character(nShared)), draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(~regionType) + theme_bw()

ggplot(chplSharedMergedShort) + geom_violin(aes(x=as.character(nShared), y=regLen, color=as.character(nShared)), draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(~regionType) + theme_bw()

chplSharedMergedShort$CHPL2021propPass <- chplSharedMergedShort$CHPL2021nRegions/chplSharedMergedShort$nWins
chplSharedMergedShort$CHPL2017propPass <- chplSharedMergedShort$CHPL2017nRegions/chplSharedMergedShort$nWins
chplSharedMergedShort$CHPL2006propPass <- chplSharedMergedShort$CHPL2006nRegions/chplSharedMergedShort$nWins
chplSharedMergedShort$CHPL2003propPass <- chplSharedMergedShort$CHPL2003nRegions/chplSharedMergedShort$nWins

ggplot(chplSharedMergedShort) + geom_violin(aes(x=as.character(nShared), y=CHPL2021propPass, color=as.character(nShared)), draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(~regionType) + theme_bw()
ggplot(chplSharedMergedShort) + geom_violin(aes(x=as.character(nShared), y=CHPL2017propPass, color=as.character(nShared)), draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(~regionType) + theme_bw()
ggplot(chplSharedMergedShort) + geom_violin(aes(x=as.character(nShared), y=CHPL2006propPass, color=as.character(nShared)), draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(~regionType) + theme_bw()
ggplot(chplSharedMergedShort) + geom_violin(aes(x=as.character(nShared), y=CHPL2003propPass, color=as.character(nShared)), draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(~regionType) + theme_bw()


chplSharedPercFoundInUnique <- chplSharedPercFoundIn[which(chplSharedPercFoundIn$popName==chplSharedPercFoundIn$foundIn),]

ggplot(chplSharedPercFoundInUnique) + geom_density(aes(x=regLen, color=foundIn)) + facet_wrap(~regionType) + theme_bw()

ggplot(chplSharedPercFoundInUnique) + geom_violin(aes(x=as.character(nShared), y=percentile, color=as.character(nShared)), draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(regionType~foundIn) + theme_bw()

ggplot(chplSharedPercFoundInUnique) + geom_violin(aes(x=as.character(nShared), y=regLen, color=as.character(nShared)), draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(regionType~foundIn, scales="free_y") + theme_bw()


head(chplSharedMergedShort)
nrow(subset(chplSharedMergedShort, regionType=="Desert"))
nrow(subset(chplSharedMergedShort, regionType=="Desert" & nShared==4))
nrow(subset(chplSharedMergedShort, regionType=="Desert" & nShared==3))
nrow(subset(chplSharedMergedShort, regionType=="Desert" & nShared==2))
nrow(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1))
nrow(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1 & CHPL2021nRegions>0))
nrow(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1 & CHPL2017nRegions>0))
nrow(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1 & CHPL2006nRegions>0))
nrow(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1 & CHPL2003nRegions>0))


nrow(subset(chplSharedMergedShort, regionType=="Island"))
nrow(subset(chplSharedMergedShort, regionType=="Island" & nShared==4))
nrow(subset(chplSharedMergedShort, regionType=="Island" & nShared==3))
nrow(subset(chplSharedMergedShort, regionType=="Island" & nShared==2))
nrow(subset(chplSharedMergedShort, regionType=="Island" & nShared==1))
nrow(subset(chplSharedMergedShort, regionType=="Island" & nShared==1 & CHPL2021nRegions>0))
nrow(subset(chplSharedMergedShort, regionType=="Island" & nShared==1 & CHPL2017nRegions>0))
nrow(subset(chplSharedMergedShort, regionType=="Island" & nShared==1 & CHPL2006nRegions>0))
nrow(subset(chplSharedMergedShort, regionType=="Island" & nShared==1 & CHPL2003nRegions>0))


head(chplSharedMergedShort)
summary(subset(chplSharedMergedShort, regionType=="Desert")$regLen)
summary(subset(chplSharedMergedShort, regionType=="Desert" & CHPL2021nRegions>0)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Desert" & CHPL2017nRegions>0)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Desert" & CHPL2006nRegions>0)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Desert" & CHPL2003nRegions>0)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==4)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==3)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==2)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1 & CHPL2021nRegions>0)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1 & CHPL2017nRegions>0)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1 & CHPL2006nRegions>0)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1 & CHPL2003nRegions>0)$regLen)

summary(subset(chplSharedMergedShort, regionType=="Island")$regLen)
summary(subset(chplSharedMergedShort, regionType=="Island" & CHPL2021nRegions>0)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Island" & CHPL2017nRegions>0)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Island" & CHPL2006nRegions>0)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Island" & CHPL2003nRegions>0)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==4)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==3)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==2)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==1)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==1 & CHPL2021nRegions>0)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==1 & CHPL2017nRegions>0)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==1 & CHPL2006nRegions>0)$regLen)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==1 & CHPL2003nRegions>0)$regLen)

summary(subset(chplSharedMergedShort, regionType=="Desert")$nWins)
summary(subset(chplSharedMergedShort, regionType=="Desert" & CHPL2021nRegions>0)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Desert" & CHPL2017nRegions>0)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Desert" & CHPL2006nRegions>0)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Desert" & CHPL2003nRegions>0)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==4)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==3)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==2)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1 & CHPL2021nRegions>0)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1 & CHPL2017nRegions>0)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1 & CHPL2006nRegions>0)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1 & CHPL2003nRegions>0)$nWins)

summary(subset(chplSharedMergedShort, regionType=="Island")$nWins)
summary(subset(chplSharedMergedShort, regionType=="Island" & CHPL2021nRegions>0)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Island" & CHPL2017nRegions>0)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Island" & CHPL2006nRegions>0)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Island" & CHPL2003nRegions>0)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==4)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==3)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==2)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==1)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==1 & CHPL2021nRegions>0)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==1 & CHPL2017nRegions>0)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==1 & CHPL2006nRegions>0)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==1 & CHPL2003nRegions>0)$nWins)

chplSharedMergedShort$meanNregWin <- NA
for (i in 1:nrow(chplSharedMergedShort)) {
  curReg <- chplSharedMergedShort[i,]
  nReg <- 0
  if (curReg$CHPL2021nRegions>0) {nReg <- nReg + curReg$CHPL2021nRegions}
  if (curReg$CHPL2017nRegions>0) {nReg <- nReg + curReg$CHPL2017nRegions}
  if (curReg$CHPL2006nRegions>0) {nReg <- nReg + curReg$CHPL2006nRegions}
  if (curReg$CHPL2003nRegions>0) {nReg <- nReg + curReg$CHPL2003nRegions}
  chplSharedMergedShort$meanNregWin[i] <- nReg/curReg$nShared
}

summary(subset(chplSharedMergedShort, regionType=="Desert")$meanNregWin)
summary(subset(chplSharedMergedShort, regionType=="Desert" & CHPL2021nRegions>0)$CHPL2021nRegions)
summary(subset(chplSharedMergedShort, regionType=="Desert" & CHPL2017nRegions>0)$CHPL2017nRegions)
summary(subset(chplSharedMergedShort, regionType=="Desert" & CHPL2006nRegions>0)$CHPL2006nRegions)
summary(subset(chplSharedMergedShort, regionType=="Desert" & CHPL2003nRegions>0)$nWins)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==4)$meanNregWin)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==3)$meanNregWin)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==2)$meanNregWin)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1)$meanNregWin)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1 & CHPL2021nRegions>0)$CHPL2021nRegions)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1 & CHPL2017nRegions>0)$CHPL2017nRegions)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1 & CHPL2006nRegions>0)$CHPL2006nRegions)
summary(subset(chplSharedMergedShort, regionType=="Desert" & nShared==1 & CHPL2003nRegions>0)$CHPL2003nRegions)

summary(subset(chplSharedMergedShort, regionType=="Island")$meanNregWin)
summary(subset(chplSharedMergedShort, regionType=="Island" & CHPL2021nRegions>0)$CHPL2021nRegions)
summary(subset(chplSharedMergedShort, regionType=="Island" & CHPL2017nRegions>0)$CHPL2017nRegions)
summary(subset(chplSharedMergedShort, regionType=="Island" & CHPL2006nRegions>0)$CHPL2006nRegions)
summary(subset(chplSharedMergedShort, regionType=="Island" & CHPL2003nRegions>0)$CHPL2003nRegions)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==4)$meanNregWin)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==3)$meanNregWin)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==2)$meanNregWin)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==1)$meanNregWin)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==1 & CHPL2021nRegions>0)$CHPL2021nRegions)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==1 & CHPL2017nRegions>0)$CHPL2017nRegions)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==1 & CHPL2006nRegions>0)$CHPL2006nRegions)
summary(subset(chplSharedMergedShort, regionType=="Island" & nShared==1 & CHPL2003nRegions>0)$CHPL2003nRegions)

write.table(chplSharedMergedShort, "shared_CHPLregionsMerged_timeseries_short.txt", row.names = F, quote = F, sep="\t")

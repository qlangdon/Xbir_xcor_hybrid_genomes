options(stringAsFactors=FALSE)
library(ggplot2)

setwd("~/Documents/SchumerLab/Xcortezi_Xbirchmanni_hybridization/desertsIslandsPacBio2023/")

desOverlap <- read.table("deserts_shared_CHPL_STAC_allReps_allSelection.txt", header=T)

selectValues <- c("s0.1", "s0.05", "s0.025", "s0.01")
desSum <- data.frame()
selSum <- data.frame()
for (selID in selectValues) {
  selData <- subset(desOverlap, selStrength == selID)
  desNoSel <- subset(selData, selected == F)
  for (repID in 1:105) {
    repData <- subset(desNoSel, repNum==repID)
    if (nrow(repData)>0) {
      desCount <- nrow(repData)
      desUniChr <- length(unique(repData$chrName))
      meanWinNum <- mean(repData$nWin)
      medianWinNum <- median(repData$nWin)
      sdWinNum <- sd(repData$nWin)
      meanLen <- mean(repData$cMlen)
      medianLen <- median(repData$cMlen)
      sdLen <- sd(repData$cMlen)
      meanChpl <- mean(repData$nChpl)
      meanStac <- mean(repData$nStac)
      tmpDF <- data.frame(selStrength=selID, repNum=repID, nDes=desCount, nUniDes=desUniChr, meanWinNum=meanWinNum, medianWinNum=medianWinNum, sdWinNum=sdWinNum, meanLen=meanLen, medianLen=medianLen, nChpl=meanChpl, nStac=meanStac)
    } else {
      tmpDF <- data.frame(selStrength=selID, repNum=repID, nDes=0, nUniDes=0, meanWinNum=0, medianWinNum=0, sdWinNum=0, meanLen=0, medianLen=0, nChpl=0, nStac=0)
    }
    desSum <- rbind(desSum, tmpDF)
  }
  selCur <- subset(selData, selected==T)
  selSumCur <-data.frame(selStrength=selID, nSelShared=nrow(selCur), meanWinNum=mean(selCur$nWin), medianWinNum=median(selCur$nWin), sdWinNum=sd(selCur$nWin), meanLen=mean(selCur$cMlen), medianLen=median(selCur$cMlen), nChpl=mean(selCur$nChpl), nStac=mean(selCur$nStac))
  selSum <- rbind(selSum, selSumCur)
}
write.table(selSum, "summary_shared_simulated_selected_deserts.txt", row.names = F, quote = F , sep="\t")
write.table(desSum, "summary_shared_simulated_non-selected_deserts.txt", row.names = F, quote = F , sep="\t")

ggplot(desSum) + geom_bar(aes(x=nUniDes, fill=selStrength), color="grey50", position = "dodge") + scale_x_continuous(breaks = seq(0, 24, 1), lim = c(NA, 24)) + theme_bw() + theme(panel.grid.minor = element_line(colour = "white")) + xlab("Number of non-selected chromosomes with shared deserts")
ggsave("simulated_non-selected_shared_deserts.pdf")

ggplot(desSum) + geom_violin(aes(x=selStrength, y=nDes), draw_quantiles = c(0.25,0.5,0.75)) + theme_bw()
ggplot(desSum) + geom_density(aes(color=selStrength, x=nDes)) + theme_bw()

ggplot(desSum) + geom_violin(aes(x=selStrength, y=nUniDes), draw_quantiles = c(0.25,0.5,0.75)) + theme_bw()
ggplot(desSum) + geom_density(aes(color=selStrength, x=nUniDes), bw=0.5) + xlim(c(0,24)) + theme_bw()
ggplot(desSum) + geom_histogram(aes(x=nUniDes, fill=selStrength)) + facet_wrap(~selStrength) + theme_bw()
ggplot(desSum) + geom_histogram(aes(x=nUniDes, fill=selStrength), color="grey50", position = "dodge") + scale_x_continuous(breaks = seq(0, 24, 1), lim = c(0, 24)) + theme_bw() + theme(panel.grid.minor = element_line(colour = "white"))
ggplot(desSum) + geom_histogram(aes(x=nUniDes, fill=selStrength), color="grey50", position = "stack") + scale_x_continuous(breaks = seq(0, 24, 1), lim = c(0, 24)) + theme_bw() + theme(panel.grid.minor = element_line(colour = "white"))

ggplot() + geom_violin(data=desSum, aes(x=selStrength, y=medianLen), color="grey20", draw_quantiles=c(0.25,0.5,0.75)) + geom_point(data=selSum, aes(x=selStrength, y=medianLen), color="coral", size=3) + geom_errorbar(data=selSum, aes(x=selStrength, ymin=medianLen-sdLen, ymax=medianLen+sdLen), color="coral") + theme_bw()

ggplot() + geom_violin(data=desSum, aes(x=selStrength, y=meanLen), color="grey20", draw_quantiles=c(0.25,0.5,0.75)) + geom_point(data=selSum, aes(x=selStrength, y=meanLen), color="coral", size=3) + geom_errorbar(data=selSum, aes(x=selStrength, ymin=meanLen-sdLen, ymax=meanLen+sdLen), color="coral") + theme_bw()

ggplot() + geom_violin(data=desSum, aes(x=selStrength, y=meanWinNum), color="grey20", draw_quantiles=c(0.25,0.5,0.75)) + geom_point(data=selSum, aes(x=selStrength, y=meanWinNum), color="coral", size=3) + geom_errorbar(data=selSum, aes(x=selStrength, ymin=meanWinNum-sdWinNum, ymax=meanWinNum+sdWinNum), color="coral") + theme_bw()


doptions(stringAsFactors=FALSE)
### On Sherlock getting the overlapping deserts and tagging as encompassing or not the selected site
#repNum <- "86"
repsCrossShared <- data.frame()
for (repNum in 1:105) {
  chplDes <- read.table(paste("CHPL_minorParentDeserts_AIMs_cMwins_merged_", repNum, ".txt", sep=""), header=T)
  stacDes <- read.table(paste("STAC_minorParentDeserts_AIMs_cMwins_merged_", repNum, ".txt", sep=""), header=T)
  selSite <- read.table(paste("selected_sites_", repNum, sep=""), header = F)
  chplWin <- read.table("average_average_ancestry_by_site_STAC_sel_generation86_ancestry_xbir_genome_0.05cM_windows_admixem.bed_WG", header=F)
  
  chplWin <- read.table(paste("average_average_ancestry_by_site_CHPL_sel_generation", repNum, "_ancestry_xbir_genome_0.05cM_windows_admixem.bed_WG", sep=""), header=F, col.names = c("chrID", "start", "end", "hybrid_index", "nInds", "nAims"))
  chplWin$chrPos <- paste(chplWin$chrID, chplWin$start, sep="_")
  chplWin$chplAnc <- chplWin$hybrid_index
  chplWin$desert <- NA
  
  stacWin <- read.table(paste("average_average_ancestry_by_site_STAC_sel_generation", repNum, "_ancestry_xbir_genome_0.05cM_windows_admixem.bed_WG", sep=""), header=F, col.names = c("chrID", "start", "end", "hybrid_index", "nInds", "nAims"))
  stacWin$chrPos <- paste(stacWin$chrID, stacWin$start, sep="_")
  stacWin$stacAnc <- stacWin$hybrid_index
  stacWin$desert <- NA
  
  if (nrow(chplDes)>0) {
    chplDes$nWins <- NA
    for (i in 1:nrow(chplDes)) {
      curDes <- chplDes[i,]
      desWins <- subset(chplWin, chrID==as.character(curDes$chrName) & start>=curDes$cMstart & end<=curDes$cMend)
      chplDes$nWins[i] <- nrow(desWins)
      chplWin$desert[which(chplWin$chrPos %in% desWins$chrPos)] <- TRUE
    }
  }
  
  if (nrow(stacDes)>0) {
    stacDes$nWins <- NA
    for (i in 1:nrow(stacDes)) {
      curDes <- stacDes[i,]
      desWins <- subset(stacWin, chrID==as.character(curDes$chrName) & start>=curDes$cMstart & end<=curDes$cMend)
      stacDes$nWins[i] <- nrow(desWins)
      stacWin$desert[which(stacWin$chrPos %in% desWins$chrPos)] <- TRUE
    }
  }
  
  if (nrow(chplDes)>0) {
    chplDes$pop <- "CHPL"
    chplDes$crossN <- NA
    for (i in 1:nrow(chplDes)) {
      curDes <- chplDes[i,]
      crossWins <- subset(stacWin, chrID==as.character(curDes$chrName) & start>=curDes$cMstart & end<=curDes$cMend)
      chplDes$crossN[i] <- nrow(subset(crossWins, desert==T))
    }
  }
  
  if (nrow(stacDes)>0) {
    stacDes$pop <- "STAC"
    stacDes$crossN <- NA
    for (i in 1:nrow(stacDes)) {
      curDes <- stacDes[i,]
      crossWins <- subset(chplWin, chrID==as.character(curDes$chrName) & start>=curDes$cMstart & end<=curDes$cMend)
      stacDes$crossN[i] <- nrow(subset(crossWins, desert==T))
    }
  }
 
  comboShare <- data.frame()
  if (nrow(chplDes)>0 & nrow(stacDes)>0) {
    chplShared <- subset(chplDes, crossN>0)
    stacShared <- subset(stacDes, crossN>0)
    desertShare <- rbind(stacShared, chplShared)
    
    for (chr in unique(desertShare$chrName)) {
      chrDesert <- subset(desertShare, chrName==chr)
      chrShare <- data.frame()
      if (nrow(chrDesert)>0) {
        desSort <- chrDesert[order(chrDesert$cMstart),]
        desSort$chrRow <- seq(1,nrow(desSort))
        curRow <- 1
        desChrComboCombo <- data.frame()
        while (curRow <= nrow(desSort)) {
          curReg <- desSort[curRow,]
          minStart <- curReg$cMstart
          maxStart <- curReg$cMend
          curOver <- subset(desSort, cMstart>=minStart & cMstart<=maxStart)
          toComb <- nrow(curOver)
          minStart <- min(curOver$cMstart)
          maxStart <- max(curOver$cMend)
          newOver <- subset(desSort, cMstart>=minStart & cMstart<=maxStart)
          while (nrow(newOver) > toComb) {
            toComb <- nrow(newOver)
            minStart <- min(newOver$cMstart)
            maxStart <- max(newOver$cMend)
            newOver <- subset(desSort, cMstart>=minStart & cMstart<=maxStart)
          }
          newStart <- min(newOver$cMstart)
          newEnd <- max(newOver$cMend)
          newLen <- newEnd - newStart
          newReg <- data.frame(chrName=chr, cMstart=newStart, cMend=newEnd, cMmid=newStart+newLen/2, cMlen=newLen, aimStart=min(newOver$aimStart), aimEnd=max(newOver$aimEnd), aimLen=max(newOver$aimEnd)-min(newOver$aimStart), regType="Desert")
          chplNewWins <- subset(chplWin, chrID==as.character(newReg$chrName) & start>=newStart & end<=newEnd)
          stacNewWins <- subset(stacWin, chrID==as.character(newReg$chrName) & start>=newStart & end<=newEnd)
          #winNew <- subset(stacChpl, chr==newReg$chrName & start>=newStart & end<=newEnd)
          nWin <- nrow(chplNewWins)/2
          nChpl <- nrow(subset(chplNewWins, desert==T))
          nStac <- nrow(subset(stacNewWins, desert==T))
          #propTotal <- (sum(chplNewWins$desert))/nrow(winNew)
          newReg$nWin <- nWin
          newReg$nChpl <- nChpl
          newReg$nStac <- nStac
          newReg$selected <- FALSE
          if (newReg$chrName==selSite$V1) { 
            if (selSite$V2>=newReg$cMstart & selSite$V2<=newReg$cMend) {
              newReg$selected <- TRUE
            } else {
              newReg$selected <- "OnChr"
            }
          }
          #newReg$propPass <- propTotal
          desChrComboCombo <- rbind(desChrComboCombo, newReg)
          curRow <- max(newOver$chrRow) + 1
        }
        desChrComboCombo$distToPrev <- desChrComboCombo$cMstart-c(0,desChrComboCombo$cMend[-nrow(desChrComboCombo)])
        chrShare <- rbind(chrShare, desChrComboCombo)
      }
      chrShareSort <- chrShare[order(chrShare$cMstart),]
      comboShare <- rbind(comboShare, chrShareSort)
    } 
  }
  
  print(paste("On rep:", repNum))
  if (nrow(comboShare)>0) {
    comboShare$repNum <- repNum
    repsCrossShared <- rbind(repsCrossShared, comboShare)
    print(nrow(comboShare))
  }
}

#repsCrossShared$selStrength <- "s0.1"
#repsCrossShared$selStrength <- "s0.05"
#repsCrossShared$selStrength <- "s0.025"
repsCrossShared$selStrength <- "s0.01"

write.table(repsCrossShared, "deserts_shared_CHPL_STAC_allReps.txt", row.names = F, quote = F, sep = "\t")

length(unique(repsCrossShared$repNum))
nrow(subset(repsCrossShared, selected==T))
summary(subset(repsCrossShared, selected==T)$nWin)
summary(subset(repsCrossShared, selected==T)$cMlen/1000000)

nrow(subset(repsCrossShared, selected=="OnChr"))
summary(subset(repsCrossShared, selected=="OnChr")$nWin)
summary(subset(repsCrossShared, selected=="OnChr")$cMlen/1000000)

nrow(subset(repsCrossShared, selected==F))
summary(subset(repsCrossShared, selected==F)$nWin)
summary(subset(repsCrossShared, selected==F)$cMlen/1000000)


sharedDeserts <- data.frame()
selectValues <- c("s0.1", "s0.05", "s0.025", "s0.01")
for (selVal in selectValues) {
  selShared <- read.table(paste("selection_simulations_", selVal, "/deserts_shared_CHPL_STAC_allReps.txt", sep=""), header=T)
  sharedDeserts <- rbind(sharedDeserts, selShared)
}

write.table(sharedDeserts, "deserts_shared_CHPL_STAC_allReps_allSelection.txt", row.names = F, quote = F, sep = "\t")


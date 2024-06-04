options(stringAsFactors=FALSE)
library(ggplot2)
setwd("~/Documents/SchumerLab/Xcortezi_Xbirchmanni_hybridization/desertsIslandsPacBio2023/")

allSites <- read.table("selected_sites_quantInfo_allReps_allS_allDeserts.txt", header=T)

#ggplot(allSites) + geom_bar(aes(x=selection, y=pass5, fill=pop), color="black", stat="identity", position="dodge") +theme_bw() + facet_grid(~desType)  

pops <- c("CHPL", "STAC")
desTypes <- c("merged", "raw")
selectValues <- c("s0.1", "s0.05", "s0.025", "s0.01")
passingCounts <- data.frame()
mergeCounts <- data.frame()
for (popID in pops) {
  popData <- subset(allSites, pop==popID)
  for (desID in desTypes) {
    popDes <- subset(popData, desType==desID)
    for (selID in selectValues) {
      popDesSel <- subset(popDes, selection==selID)
      tmpDF <- data.frame(popName=popID, desType=desID, selection=selID, totQuant5=sum(popDesSel$pass2.5,na.rm=T), propQuant5=sum(popDesSel$pass2.5,na.rm=T)/nrow(subset(popDesSel, !is.na(pass2.5))), totcMquant10=sum(popDesSel$pass5,na.rm=T), propcMquant10=sum(popDesSel$pass5,na.rm=T)/nrow(subset(popDesSel, !is.na(pass5))), totWinDes=sum(popDesSel$withinDesert,na.rm=T), propWinDesert=sum(popDesSel$withinDesert,na.rm=T)/nrow(subset(popDesSel, !is.na(withinDesert))))
      passingCounts <- rbind(passingCounts, tmpDF)
      if (desID == "merged") {
        newTmp <- data.frame(popName=popID, desType=desID, selection=selID, passCatagory=c("AIMquant5", "cMwinQuant10", "withinDesert"), propPass=c(sum(popDesSel$pass2.5,na.rm=T)/nrow(subset(popDesSel, !is.na(pass2.5))), sum(popDesSel$pass5,na.rm=T)/nrow(subset(popDesSel, !is.na(pass5))), sum(popDesSel$withinDesert,na.rm=T)/nrow(subset(popDesSel, !is.na(withinDesert)))), totPass=c(sum(popDesSel$pass2.5,na.rm=T), sum(popDesSel$pass5,na.rm=T), sum(popDesSel$withinDesert,na.rm=T)))
        mergeCounts <- rbind(mergeCounts, newTmp)
      }
    }
  }
}

ggplot(passingCounts) + geom_bar(aes(x=selection, y=propQuant5, fill=popName), color="black", stat="identity", position="dodge") +theme_bw() + facet_grid(~desType) + ylim(c(0,1)) 

ggplot(mergeCounts) + geom_bar(aes(x=selection, y=totPass, fill=passCatagory), color="black", stat="identity", position="dodge") +theme_bw() + facet_grid(~popName) + ylim(c(0,105)) 
ggplot(mergeCounts) + geom_bar(aes(x=selection, y=propPass, fill=passCatagory), color="black", stat="identity", position="dodge") +theme_bw() + facet_grid(~popName) + ylim(c(0,1)) 

desertsFound <- subset(mergeCounts, passCatagory=="withinDesert")
newColors <- c("CHPL"="#45b9bc", "STAC"="#308385")

ggplot(desertsFound) + geom_bar(aes(x=selection, y=propPass, fill=popName), color="black", stat="identity", position="dodge") +theme_bw() + ylim(c(0,1)) + scale_fill_manual(values=newColors, name="Population") + ylab("Proportion of Deserts Found")
ggsave("deserts_power.pdf")
write.table(desertsFound, "deserts_power_counts.txt", row.names = F, quote = F, sep = "\t")

##Combining on Sherlock
options(stringAsFactors=FALSE)

pops <- c("CHPL", "STAC")
desTypes <- c("merged", "raw")
selectValues <- c("s0.1", "s0.05", "s0.025", "s0.01")
allOut <- data.frame()
for (popID in pops) {
  for (desID in desTypes) {
    for (selID in selectValues) {
      powerDF <- read.table(paste("selection_simulations_", selID, "/selected_sites_allReps_", desID, "Deserts_", popID, ".txt", sep=""), header = F)
      if (desID=="merged") {
        colnames(powerDF) <- c("chrName", "start", "end", "distToAim", "ancestry", "pass2.5", "pass5", "withinDesert", "closestDesertDist", "closestDesertLength", "closestDesertName", "closestDesertDistMB", "repNum")
      } else if (desID=="raw") {
        tmpPower <- powerDF
        colnames(tmpPower) <- c("chrName", "start", "end", "distToAim", "ancestry", "pass2.5", "pass5", "withinDesert", "closestDesertDist", "closestDesertLength", "closestDesertName", "repNum")
        powerDF <- cbind(tmpPower[,1:11], data.frame(closestDesertDistMB=NA), repNum=tmpPower[,12])
        powerDF$closestDesertDistMB <- powerDF$closestDesertDist/1000000
      }
      powerDF$desType <- desID
      powerDF$pop <- popID
      powerDF$selection <- selID
      print(desID)
      print(popID)
      print(selID)
      print(summary(powerDF$withinDesert))
      allOut <- rbind(allOut, powerDF)
    }
  }
}

write.table(allOut, "selected_sites_quantInfo_allReps_allS_allDeserts.txt", row.names = F, quote = F, sep = "\t")



##Checking if the selected site was captured
options(stringAsFactors=FALSE)

#popName <- "CHPL"
popName <- "STAC"
print(popName)
#print(repNum)

popRawSummary <- data.frame()
popMergedSummary <- data.frame()
for (j in 1:100) {
  repNum <- j
  print(repNum)
  aimData <- read.table(paste("average_ancestry_by_site_", popName, "_sel_generation", repNum, sep=""), header=T)
  meanHI <- mean(aimData$hybrid_index, na.rm=T)
  if (meanHI > 0.5) {
    aimData$majPar <- aimData$hybrid_index
    aimData$minPar <- 1-aimData$hybrid_index
  } else {
    aimData$majPar <- 1-aimData$hybrid_index
    aimData$minPar <- aimData$hybrid_index
  }
  meanMajAnc <- mean(aimData$majPar, na.rm = T)
  meanMinAnc <- mean(aimData$minPar, na.rm = T)
  quant99 <- quantile(aimData$minPar, 0.99, na.rm = T)
  quant97.5 <- quantile(aimData$minPar, 0.975, na.rm = T)
  quant95 <- quantile(aimData$minPar, 0.95, na.rm = T)
  quant90 <- quantile(aimData$minPar, 0.9, na.rm = T)
  quant1 <- quantile(aimData$minPar, 0.01, na.rm = T)
  quant2.5 <- quantile(aimData$minPar, 0.025, na.rm = T)
  quant5 <- quantile(aimData$minPar, 0.05, na.rm = T)
  quant10 <- quantile(aimData$minPar, 0.1, na.rm = T)
  
  cMdata <- read.table(paste("average_average_ancestry_by_site_", popName, "_sel_generation", repNum, "_ancestry_xbir_genome_0.05cM_windows_admixem.bed_WG", sep=""), header=F, col.names = c("chrID", "start", "end", "hybrid_index", "nInds", "nAims"))
  meanHI <- mean(cMdata$hybrid_index, na.rm=T)
  if (meanHI > 0.5) {
    cMdata$majPar <- cMdata$hybrid_index
    cMdata$minPar <- 1-cMdata$hybrid_index
  } else {
    cMdata$majPar <- 1-cMdata$hybrid_index
    cMdata$minPar <- cMdata$hybrid_index
  }
  cMdata$SNPs <- 15
  
  cMmean <- mean(cMdata$minPar, na.rm=T)
  cMupper10 <- quantile(cMdata$minPar, 0.9, na.rm=T)
  cMlower10 <- quantile(cMdata$minPar, 0.1, na.rm=T)
  cMupper5 <- quantile(cMdata$minPar, 0.95, na.rm=T)
  cMlower5 <- quantile(cMdata$minPar, 0.05, na.rm=T)
  
  cMdata$win10quant <- "mid"
  cMdata$win10quant[which(cMdata$minPar<=cMlower10)] <- "lower"
  cMdata$win10quant[which(cMdata$minPar>=cMupper10)] <- "upper"
  cMdata$win5quant <- "mid"
  cMdata$win5quant[which(cMdata$minPar<=cMlower5)] <- "lower"
  cMdata$win5quant[which(cMdata$minPar>=cMupper5)] <- "upper"
  
  selectedRegs <- read.table(paste("selected_sites_", repNum, sep=""), header = F)
  colnames(selectedRegs) <- c("chrName", "start", "end")
  
  selectedRegs$distToAim <- NA
  selectedRegs$ancestry <- NA
  selectedRegs$quant5 <- NA
  selectedRegs$quant10 <- NA
  selectedRegs$cMquant5 <- NA
  selectedRegs$cMquant10 <- NA
  for (i in 1:nrow(selectedRegs)) {
    regChr <- selectedRegs$chrName[i]
    regPos <- selectedRegs$start[i]
    chrData <- subset(aimData, as.character(group)==regChr)
    chrData$distToSel <- chrData$position-regPos
    closestAIM <- chrData[which(abs(chrData$distToSel)==min(abs(chrData$distToSel))),]
    selectedRegs$ancestry[i] <- closestAIM$minPar
    selectedRegs$distToAim[i] <- closestAIM$distToSel
    selectedRegs$quant10[i] <- closestAIM$minPar<=quant10
    selectedRegs$quant5[i] <- closestAIM$minPar<=quant5
    site.cMwin <- subset(cMdata, chrID==as.character(regChr) & start<=regPos & end>=regPos)
    if (nrow(site.cMwin)==1) {
      if (site.cMwin$win10quant=="lower") {
        selectedRegs$cMquant10[i] <- TRUE
      } else {
        selectedRegs$cMquant10[i] <- FALSE
      }
      if (site.cMwin$win5quant=="lower") {
        selectedRegs$cMquant5[i] <- TRUE
      } else {
        selectedRegs$cMquant5[i] <- FALSE
      }
    }
  }
  selectedRegsMerge <- selectedRegs
  
  selectedRegs$withinDesert <- NA
  selectedRegs$closestDesertDist <- NA
  selectedRegs$closestDesertLength <- NA
  selectedRegs$closestDesertName <- NA
  selectedRegs$closestDesertDistMB <- NA
  if (file.exists(paste(popName, "_minorParentDeserts_AIMs_cMwins_raw_", repNum, ".txt", sep=""))) {
    desertDF <- read.table(paste(popName, "_minorParentDeserts_AIMs_cMwins_raw_", repNum, ".txt", sep=""), header=T)
    for (i in 1:nrow(selectedRegs)) {
      regChr <- selectedRegs$chrName[i]
      regPos <- selectedRegs$start[i]
      desertRegs <- subset(desertDF, as.character(chrName)==regChr & cMstart<=regPos & cMend>=regPos)
      if (nrow(desertRegs)>0) {
        selectedRegs$withinDesert[i] <- TRUE
        selectedRegs$closestDesertDist[i] <- 0
        selectedRegs$closestDesertLength[i] <- desertRegs$cMlen
        selectedRegs$closestDesertName[i] <- paste(desertRegs$chrName, desertRegs$cMstart, sep="_")
      } else {
        selectedRegs$withinDesert[i] <- FALSE
        chrDeserts <- subset(desertDF, as.character(chrName)==regChr)
        if (nrow(chrDeserts) > 0) {
          chrDeserts$midDistToSel <- chrDeserts$cMmid-regPos
          chrDeserts$midDistToSel[which(is.na(chrDeserts$midDistToSel))] <- max(chrDeserts$midDistToSel, na.rm = T)
          closestDesert <- chrDeserts[which(abs(chrDeserts$midDistToSel)==min(abs(chrDeserts$midDistToSel))),]
          selectedRegs$closestDesertDist[i] <- closestDesert$midDistToSel
          selectedRegs$closestDesertLength[i] <-  closestDesert$cMlen
          selectedRegs$closestDesertName[i] <- paste(closestDesert$chrName, closestDesert$cMstart, sep="_")
        } else {
          selectedRegs$closestDesertDist[i] <- NA
          selectedRegs$closestDesertLength[i] <- NA
          selectedRegs$closestDesertName[i] <- NA
        }
      }
    }
    selectedRegs$closestDesertDistMB <- selectedRegs$closestDesertDist/1000000
  }
  selectedRegs$repNum <- repNum
  write.table(selectedRegs, paste("selected_sites_", repNum, "_rawDeserts_", popName, ".txt", sep=""), quote = F, row.names = F, sep="\t")
  
  selectedRegsMerge$withinDesert <- NA
  selectedRegsMerge$closestDesertDist <- NA
  selectedRegsMerge$closestDesertLength <- NA
  selectedRegsMerge$closestDesertName <- NA
  selectedRegsMerge$closestDesertDistMB <- NA
  #selectedRegsMerge$desert <- NA
  if (file.exists(paste(popName, "_minorParentDeserts_AIMs_cMwins_merged_", repNum, ".txt", sep=""))) {
    mergeDeserts <- read.table(paste(popName, "_minorParentDeserts_AIMs_cMwins_merged_", repNum, ".txt", sep=""), header = T)
    for (i in 1:nrow(selectedRegsMerge)) {
      regChr <- selectedRegsMerge$chrName[i]
      regPos <- selectedRegsMerge$start[i]
      desertRegs <- subset(mergeDeserts, as.character(chrName)==regChr & cMstart<=regPos & cMend>=regPos)
      if (nrow(desertRegs)>0) {
        selectedRegsMerge$withinDesert[i] <- TRUE
        selectedRegsMerge$closestDesertDist[i] <- 0
        selectedRegsMerge$closestDesertLength[i] <- desertRegs$cMlen
        selectedRegsMerge$closestDesertName[i] <- paste(desertRegs$chrName, desertRegs$cMstart, sep="_")
      } else {
        selectedRegsMerge$withinDesert[i] <- FALSE
        chrDeserts <- subset(mergeDeserts, as.character(chrName)==regChr)
        if (nrow(chrDeserts)>0) {
          chrDeserts$midDistToSel <- chrDeserts$cMmid-regPos
          chrDeserts$midDistToSel[which(is.na(chrDeserts$midDistToSel))] <- max(chrDeserts$midDistToSel, na.rm = T)
          closestDesert <- chrDeserts[which(abs(chrDeserts$midDistToSel)==min(abs(chrDeserts$midDistToSel))),]
          selectedRegsMerge$closestDesertDist[i] <- closestDesert$midDistToSel
          selectedRegsMerge$closestDesertLength[i] <- closestDesert$cMlen
          selectedRegsMerge$closestDesertName[i] <- paste(closestDesert$chrName, closestDesert$cMstart, sep="_")
        }
      }
    }
    selectedRegsMerge$closestDesertDistMB <- selectedRegsMerge$closestDesertDist/1000000
  }
  selectedRegsMerge$repNum <- repNum
  write.table(selectedRegsMerge, paste("selected_sites_", repNum, "_desertsMerged_", popName, ".txt", sep=""), quote = F, row.names = F, sep="\t")
  selectedRegsMerge
  
  selectedRegs$popName <- popName
  popRawSummary <- rbind(popRawSummary, selectedRegs)
  selectedRegsMerge$popName <- popName
  popMergedSummary <- rbind(popMergedSummary, selectedRegsMerge)
  
}

write.table(popRawSummary, paste("selected_sites_quantInfo_allReps_desertsRaw_", popName, ".txt", sep=""), quote = F, row.names = F, sep="\t")
write.table(popMergedSummary, paste("selected_sites_quantInfo_allReps_desertsMerged_", popName, ".txt", sep=""), quote = F, row.names = F, sep="\t")

#popMergedSummaryCHPL <- popMergedSummary
#popRawSummaryCHPL <- popRawSummary

popRawSummary$desType <- "raw"
popRawSummaryCHPL$desType <- "raw"
popRawSummary$selection <- "s0.01"
popRawSummaryCHPL$selection <- "s0.01"
summary(popRawSummary)
summary(popRawSummaryCHPL)

rawSum <- rbind(popRawSummaryCHPL, popRawSummary)
write.table(rawSum, paste("selected_sites_quantInfo_allReps_s0.01_desertsRaw.txt", sep=""), quote = F, row.names = F, sep="\t")

popMergedSummary$desType <- "merged"
popMergedSummaryCHPL$desType <- "merged"
popMergedSummary$selection <- "s0.01"
popMergedSummaryCHPL$selection <- "s0.01"
summary(popMergedSummary)
summary(popMergedSummaryCHPL)

mergedSum <- rbind(popMergedSummaryCHPL, popMergedSummary)
write.table(mergedSum, paste("selected_sites_quantInfo_allReps_s0.01_desertsMerged.txt", sep=""), quote = F, row.names = F, sep="\t")

allDeserts <- rbind(mergedSum, rawSum)
write.table(allDeserts, paste("selected_sites_quantInfo_allReps_s0.01_allDeserts.txt", sep=""), quote = F, row.names = F, sep="\t")


des.s0.1 <- read.table("selection_simulations_s0.1/selected_sites_quantInfo_allReps_s0.1_allDeserts.txt", header = T)
des.s0.05 <- read.table("selection_simulations_s0.05/selected_sites_quantInfo_allReps_s0.05_allDeserts.txt", header = T)
des.s0.025 <- read.table("selection_simulations_s0.025/selected_sites_quantInfo_allReps_s0.025_allDeserts.txt", header = T)
des.s0.01 <- read.table("selection_simulations_s0.01/selected_sites_quantInfo_allReps_s0.01_allDeserts.txt", header = T)

comboDeserts <- rbind(des.s0.1, des.s0.05, des.s0.025, des.s0.01)
write.table(comboDeserts, "selected_sites_quantInfo_allReps_allS_allDeserts.txt", row.names = F, quote = F, sep = "\t")

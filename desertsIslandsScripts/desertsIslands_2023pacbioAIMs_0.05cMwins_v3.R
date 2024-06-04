options(stringAsFactors=FALSE)
args <- commandArgs(TRUE)
popName <- args[1]
comboPopPrefix <- args[2]
print(popName)

aimData <- read.table(paste("average_ancestry_by_site_thinnedAIMs_", popName, ".txt", sep=""), header=T)
meanHI <- mean(aimData$hybrid_index, na.rm=T)
if (meanHI > 0.5) {
  aimData$majPar <- aimData$hybrid_index
  aimData$minPar <- 1-aimData$hybrid_index
} else {
  aimData$majPar <- 1-aimData$hybrid_index
  aimData$minPar <- aimData$hybrid_index
}

allPops.cM <- read.table(paste("xbir_pacbio2023_allChrs_0.05cM_windows_recRate_codingBP_conservedBP_thinnedAIMs_", comboPopPrefix, "_minPar.txt", sep=""), header=T)
popToGrep <- paste("^", popName, "$", sep="")
popCol <- grep(popName, colnames(allPops.cM))
minParAnc <- allPops.cM[,popCol]
cMdata <- data.frame(chrID=allPops.cM$chrID, start=allPops.cM$start, end=allPops.cM$end, cM=allPops.cM$recRate, SNPs=allPops.cM$nSNPs, minPar=minParAnc)

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

cMmean <- mean(cMdata$minPar, na.rm=T)
cMupper10 <- quantile(cMdata$minPar, 0.9, na.rm=T)
cMlower10 <- quantile(cMdata$minPar, 0.1, na.rm=T)
cMupper5 <- quantile(cMdata$minPar, 0.95, na.rm=T)
cMlower5 <- quantile(cMdata$minPar, 0.05, na.rm=T)

aimData$aim5quant <- "mid"
aimData$aim5quant[which(aimData$minPar<=quant5)] <- "lower"
aimData$aim5quant[which(aimData$minPar>=quant95)] <- "upper"

aimData$aim10quant <- "mid"
aimData$aim10quant[which(aimData$minPar<=quant10)] <- "lower"
aimData$aim10quant[which(aimData$minPar>=quant90)] <- "upper"

cMdata$win10quant <- "mid"
cMdata$win10quant[which(cMdata$minPar<=cMlower10)] <- "lower"
cMdata$win10quant[which(cMdata$minPar>=cMupper10)] <- "upper"
cMdata$win5quant <- "mid"
cMdata$win5quant[which(cMdata$minPar<=cMlower5)] <- "lower"
cMdata$win5quant[which(cMdata$minPar>=cMupper5)] <- "upper"

#chr="ScPPXeE-107-HRSCAF-158"

desertDF <- data.frame()
uniChrs <- unique(cMdata$chrID)
for (chr in uniChrs) {
  chrAIMdata <- subset(aimData, group==chr)
  chrAIMdata$chrRow <- seq(1,nrow(chrAIMdata))
  chr.cMdata <- subset(cMdata, chrID==chr)
  chr.cMdata$chrRow <- seq(1,nrow(chr.cMdata))
  aimRow <- min(chrAIMdata$chrRow[which(chrAIMdata$aim5quant=="lower")])
  while (aimRow < nrow(chrAIMdata)) {
    #chrAIMdata[aimRow,]
    aimInfo <- chrAIMdata[aimRow,]
    aimPos <- aimInfo$position
    aim.cMwin <- subset(chr.cMdata, start<=aimInfo$position & end>=aimInfo$position)
    if (nrow(aim.cMwin)==1) {
      if ((aim.cMwin$win10quant=="lower" | is.na(aim.cMwin$minPar))) {
        #print(aimRow)
        cmRow <- aim.cMwin$chrRow
        ##Something to check if previous window is actually in the lower10 (unlikely)
        preWin <- chr.cMdata[cmRow-1,]
        if (nrow(preWin)>0) {
          if (preWin$win10quant=="lower") {
            checkPreWin <- T
            while (checkPreWin == T) {
              if (preWin$win10quant=="lower") {
                preWinRow <- preWin$chrRow
                #print(prePos)
                if (preWinRow > 1) {
                  preWin <- chr.cMdata[preWinRow-1,]
                } else {
                  checkPreWin <- F
                }
              } else {
                checkPreWin <- F
              }
              #print(paste("Check upstream:", preAIM$chrRow))
            }
            cmReg <- chr.cMdata[preWinRow+1,]
            cmRow <- chr.cMdata$chrRow[preWinRow+1]
          } else {
            cmReg <- aim.cMwin
          }
        } else {
          cmReg <- aim.cMwin
        }
        #print(paste("Previous window passes:", aim.cMwin$chrRow))
        ##Increase by cM wins until they stop being in the lower10
        nextWinQuant <- NA
        checkRow <- cmRow+1
        if (chr.cMdata$win10quant[checkRow]=="lower" | is.na(chr.cMdata$minPar[checkRow])) {
          nextWinQuant <- T
        } else {
          nextWinQuant <- F
        }
        while (nextWinQuant == T) {
          nextWin <- chr.cMdata[checkRow,]
          cmReg <- rbind(cmReg, nextWin)
          checkRow <- checkRow + 1
          if (checkRow <= max(chr.cMdata$chrRow))  {
            if (chr.cMdata$win10quant[checkRow]=="lower" | is.na(chr.cMdata$minPar[checkRow])) {
              nextWinQuant <- T
            } else {
              nextWinQuant <- F
            }
          } else {
            nextWinQuant <- F
          }
        }
        ##Check the AIMs just outside of bounds of cM regions
        firstAIMpos <- NA
        lastAIMpos <- NA
        if (aimRow > 1) {
          preAIM <- tail(subset(chrAIMdata, position <= min(cmReg$start)), n=1)
          if (nrow(preAIM)>0) {
            if (preAIM$aim10quant=="lower") {
              checkPre <- T
              while (checkPre == T) {
                if (preAIM$aim10quant=="lower") {
                  prePos <- preAIM$position
                  #print(prePos)
                  preAIM <- tail(subset(chrAIMdata, position < prePos), n=1)
                  if (nrow(preAIM) == 0) {
                    checkPre <- F
                  }
                } else {
                  checkPre <- F
                }
                #print(paste("Check upstream:", preAIM$chrRow))
              }
              firstAIMpos <- prePos
            }
          }
        }
        if (nrow(subset(chrAIMdata, position >= max(cmReg$end)))>0) {
          postAIM <- subset(chrAIMdata, position >= max(cmReg$end))[1,]
          if (postAIM$aim10quant=="lower") {
            checkPost <- T
            while (checkPost == T) {
              if (postAIM$aim10quant=="lower") {
                postPos <- postAIM$position
                #print(postPos)
                if (nrow(subset(chrAIMdata, position > postPos))>0) {
                  postAIM <- subset(chrAIMdata, position > postPos)[1,]
                } else {
                  checkPost <- F
                }
              } else {
                checkPost <- F
              }
              #print(paste("Check upstream:", preAIM$chrRow))
            }
            lastAIMpos <- postPos
            #print(paste("Check downstream:", postAIM$chrRow))
          }
        }
        ##trim cM windows where minPar is NA if they are the first of last
        if (is.na(cmReg$minPar[1])) {
          cmReg <- cmReg[2:nrow(cmReg),] 
          if (is.na(cmReg$minPar[1])) {
            print(paste("Need to check upstream NAs:", cmReg$chrRow[1]))
          }
        } else if (is.na(cmReg$minPar[nrow(cmReg)])) {
          checking <- T 
          cmReg <- cmReg[1:nrow(cmReg)-1,] 
          while (checking == T) {
            if (is.na(cmReg$minPar[nrow(cmReg)])) {
              cmReg <- cmReg[1:nrow(cmReg)-1,]
            } else {
              checking <- F
            }
          }
        }
        ##Get all AIMs within the final region
        aimsInReg <- subset(chrAIMdata, position >= min(cmReg$start)  & position <= max(cmReg$end))
        if (is.na(firstAIMpos)) {
          firstAIMpos <- min(aimsInReg$position[which(aimsInReg$aim10quant=="lower")])
        }
        if (is.na(lastAIMpos)) {
          lastAIMpos <- max(aimsInReg$position[which(aimsInReg$aim10quant=="lower")])
        }
        aimSubset <- subset(aimsInReg, position >= firstAIMpos  & position <= lastAIMpos)
        cMstart <- cmReg$start[1]
        cMend <- cmReg$end[nrow(cmReg)]
        cMlen <- cMend - cMstart
        cMmid <- cMstart + (cMlen/2) 
        cMsnps <- sum(cmReg$SNPs)
        cMpropAIMs <- nrow(subset(aimsInReg, aim10quant=="lower"))/nrow(aimsInReg)
        minAnc <- mean(aimsInReg$minPar)
        aimLen <- lastAIMpos- firstAIMpos
        aimMid <- firstAIMpos + (aimLen/2)
        propAIMs <- nrow(subset(aimSubset, aim10quant=="lower"))/nrow(aimSubset) 
        desertInfo <- data.frame(chrName=chr, cMstart=cMstart, cMend=cMend, cMmid=cMmid, cMlen=cMlen, cMsnps=cMsnps, cMpropAIMs=cMpropAIMs, minAnc=minAnc, aimStart=firstAIMpos, aimEnd=lastAIMpos, aimLen=aimLen, aimMid=aimMid, nAIMs=nrow(aimSubset), propAIMs=propAIMs)
        #print(desertInfo)
        desertDF <- rbind(desertDF, desertInfo)
        if (nrow(chrAIMdata[which(chrAIMdata$aim5quant=="lower" & chrAIMdata$position>max(cMend, lastAIMpos)),])>0) {
          aimRowTest <- min(chrAIMdata$chrRow[which(chrAIMdata$aim5quant=="lower" & chrAIMdata$position>max(cMend, lastAIMpos))])
          if (aimRowTest > aimRow) {
            aimRow <-aimRowTest
          } else  {
            aimRow <- min(subset(chrAIMdata, aim5quant=="lower" &  chrRow>aimRowTest)$chrRow)
          }
        } else {
          aimRow <- nrow(chrAIMdata) + 1
        }
      }
      else {
        if (nrow(chrAIMdata[which(chrAIMdata$aim5quant=="lower" & chrAIMdata$chrRow>aimRow),])>0) {
          aimRow <- min(chrAIMdata$chrRow[which(chrAIMdata$aim5quant=="lower" & chrAIMdata$chrRow>aimRow)])
        } else {
          aimRow <- nrow(chrAIMdata) + 1
        }
      }
    } else {
      if (nrow(chrAIMdata[which(chrAIMdata$aim5quant=="lower" & chrAIMdata$chrRow>aimRow),])>0) {
        aimRow <- min(chrAIMdata$chrRow[which(chrAIMdata$aim5quant=="lower" & chrAIMdata$chrRow>aimRow)])
      } else {
        aimRow <- nrow(chrAIMdata) + 1
      }
    }
  }
  print(chr)
  #print(dim(desertDF))
}
print("Deserts done")
print(dim(desertDF))
write.table(desertDF, paste(popName, "_minorParentDeserts_AIMs_cMwins_raw.txt", sep=""), quote=F, row.names = F, sep="\t")

islandDF <- data.frame()
uniChrs <- unique(cMdata$chrID)
for (chr in uniChrs) {
  chrAIMdata <- subset(aimData, group==chr)
  chrAIMdata$chrRow <- seq(1,nrow(chrAIMdata))
  chr.cMdata <- subset(cMdata, chrID==chr)
  chr.cMdata$chrRow <- seq(1,nrow(chr.cMdata))
  aimRow <- min(chrAIMdata$chrRow[which(chrAIMdata$aim5quant=="upper")])
  while (aimRow < nrow(chrAIMdata)) {
    #print(aimRow)
    #chrAIMdata[aimRow,]
    aimInfo <- chrAIMdata[aimRow,]
    aimPos <- aimInfo$position
    aim.cMwin <- subset(chr.cMdata, start<=aimInfo$position & end>=aimInfo$position)
    if (nrow(aim.cMwin)==1) {
      if ((aim.cMwin$win10quant=="upper" | is.na(aim.cMwin$minPar))) {
        #print(aimRow)
        cmRow <- aim.cMwin$chrRow
        ##Something to check if previous window is actually in the upper10 (unlikely)
        preWin <- chr.cMdata[cmRow-1,]
        if (nrow(preWin)>0) {
          if (preWin$win10quant=="upper") {
            checkPreWin <- T
            while (checkPreWin == T) {
              if (nrow(preWin)==0) {
                checkPreWin <- F
              } else if (preWin$win10quant=="upper") {
                preWinRow <- preWin$chrRow
                #print(prePos)
                preWin <- chr.cMdata[preWinRow-1,]
              } else {
                checkPreWin <- F
              }
              #print(paste("Check upstream:", preAIM$chrRow))
            }
            cmReg <- chr.cMdata[preWinRow+1,]
            cmRow <- chr.cMdata$chrRow[preWinRow+1]
          } else {
            cmReg <- aim.cMwin
          }
        } else {
          cmReg <- aim.cMwin
        }
        #print(paste("Previous window passes:", aim.cMwin$chrRow))
        #check if it's alreay covered by an island
        lastIsland <- islandDF[nrow(islandDF),]
        if (nrow(lastIsland)==0) {
          lastIsland <- data.frame(cMend=0, aimEnd=0) 
        } else if (lastIsland!=chr) {
          lastIsland <- data.frame(cMend=0, aimEnd=0) 
        }
        if (cmReg$start <= lastIsland$cMend | cmReg$start <= lastIsland$aimEnd) {
          if (nrow(chrAIMdata[which(chrAIMdata$aim5quant=="upper" & chrAIMdata$chrRow>aimRow),])>0) {
            aimRow <- min(chrAIMdata$chrRow[which(chrAIMdata$aim5quant=="upper" & chrAIMdata$chrRow>aimRow)])
          } else {
            aimRow <- nrow(chrAIMdata) + 1
          }
        } else {
          nextWinQuant <- NA
          checkRow <- cmRow+1
          if (chr.cMdata$win10quant[checkRow]=="upper" | is.na(chr.cMdata$minPar[checkRow])) {
            nextWinQuant <- T
          } else {
            nextWinQuant <- F
          }
          while (nextWinQuant == T) {
            nextWin <- chr.cMdata[checkRow,]
            cmReg <- rbind(cmReg, nextWin)
            checkRow <- checkRow + 1
            if (checkRow <= max(chr.cMdata$chrRow))  {
              if (chr.cMdata$win10quant[checkRow]=="upper" | is.na(chr.cMdata$minPar[checkRow])) {
                nextWinQuant <- T
              } else {
                nextWinQuant <- F
              }
            } else {
              nextWinQuant <- F
            }
          }
          ##Check the AIMs just outside of bounds of cM regions
          firstAIMpos <- NA
          lastAIMpos <- NA
          if (aimRow > 1) {
            preAIM <- tail(subset(chrAIMdata, position <= min(cmReg$start)), n=1)
            if (nrow(preAIM)>0) {
              if (preAIM$aim10quant=="upper") {
                checkPre <- T
                while (checkPre == T) {
                  if (preAIM$aim10quant=="upper") {
                    prePos <- preAIM$position
                    #print(prePos)
                    preAIM <- tail(subset(chrAIMdata, position < prePos), n=1)
                    if (nrow(preAIM) == 0) {
                      checkPre <- F
                    }
                  } else {
                    checkPre <- F
                  }
                  #print(paste("Check upstream:", preAIM$chrRow))
                }
                firstAIMpos <- prePos
              }
            }
          }
          if (nrow(subset(chrAIMdata, position >= max(cmReg$end)))>0) {
            postAIM <- subset(chrAIMdata, position >= max(cmReg$end))[1,]
            if (postAIM$aim10quant=="upper") {
              checkPost <- T
              while (checkPost == T) {
                if (postAIM$aim10quant=="upper") {
                  postPos <- postAIM$position
                  #print(postPos)
                  if (nrow(subset(chrAIMdata, position > postPos))>0) {
                    postAIM <- subset(chrAIMdata, position > postPos)[1,]
                  } else {
                    checkPost <- F
                  }
                } else {
                  checkPost <- F
                }
                #print(paste("Check upstream:", preAIM$chrRow))
              }
              lastAIMpos <- postPos
              #print(paste("Check downstream:", postAIM$chrRow))
            }
          }
          ##trim cM windows where minPar is NA if they are the first of last
          if (is.na(cmReg$minPar[1])) {
            cmReg <- cmReg[2:nrow(cmReg),] 
            if (is.na(cmReg$minPar[1])) {
              print(paste("Need to check upstream NAs:", cmReg$chrRow[1]))
            }
          } else if (is.na(cmReg$minPar[nrow(cmReg)])) {
            checking <- T 
            cmReg <- cmReg[1:nrow(cmReg)-1,] 
            while (checking == T) {
              if (is.na(cmReg$minPar[nrow(cmReg)])) {
                cmReg <- cmReg[1:nrow(cmReg)-1,]
              } else {
                checking <- F
              }
            }
          }
          ##Get all AIMs within the final region
          aimsInReg <- subset(chrAIMdata, position >= min(cmReg$start)  & position <= max(cmReg$end))
          if (is.na(firstAIMpos)) {
            firstAIMpos <- min(aimsInReg$position[which(aimsInReg$aim10quant=="upper")])
          }
          if (is.na(lastAIMpos)) {
            lastAIMpos <- max(aimsInReg$position[which(aimsInReg$aim10quant=="upper")])
          }
          aimSubset <- subset(aimsInReg, position >= firstAIMpos  & position <= lastAIMpos)
          cMstart <- cmReg$start[1]
          cMend <- cmReg$end[nrow(cmReg)]
          cMlen <- cMend - cMstart
          cMmid <- cMstart + (cMlen/2) 
          cMsnps <- sum(cmReg$SNPs)
          cMpropAIMs <- nrow(subset(aimsInReg, aim10quant=="upper"))/nrow(aimsInReg)
          minAnc <- mean(aimsInReg$minPar)
          aimLen <- lastAIMpos- firstAIMpos
          aimMid <- firstAIMpos + (aimLen/2)
          propAIMs <- nrow(subset(aimSubset, aim10quant=="upper"))/nrow(aimSubset) 
          islandInfo <- data.frame(chrName=chr, cMstart=cMstart, cMend=cMend, cMmid=cMmid, cMlen=cMlen, cMsnps=cMsnps, cMpropAIMs=cMpropAIMs, minAnc=minAnc, aimStart=firstAIMpos, aimEnd=lastAIMpos, aimLen=aimLen, aimMid=aimMid, nAIMs=nrow(aimSubset), propAIMs=propAIMs)
          #print(islandInfo)
          #print(islandInfo)
          islandDF <- rbind(islandDF, islandInfo)
          if (nrow(chrAIMdata[which(chrAIMdata$aim5quant=="upper" & chrAIMdata$position>max(cMend, lastAIMpos)),])>0) {
            aimRowTest <- min(chrAIMdata$chrRow[which(chrAIMdata$aim5quant=="upper" & chrAIMdata$position>max(cMend, lastAIMpos))])
            if (aimRowTest > aimRow) {
              aimRow <-aimRowTest
            } else  {
              aimRow <- min(subset(chrAIMdata, aim5quant=="upper" &  chrRow>aimRowTest)$chrRow)
            }
          } else {
            aimRow <- nrow(chrAIMdata) + 1
          }
        }
        ##Increase by cM wins until they stop being in the upper10
      }
      else {
        if (nrow(chrAIMdata[which(chrAIMdata$aim5quant=="upper" & chrAIMdata$chrRow>aimRow),])>0) {
          aimRow <- min(chrAIMdata$chrRow[which(chrAIMdata$aim5quant=="upper" & chrAIMdata$chrRow>aimRow)])
        } else {
          aimRow <- nrow(chrAIMdata) + 1
        }
      }
    } else {
      if (nrow(chrAIMdata[which(chrAIMdata$aim5quant=="upper" & chrAIMdata$chrRow>aimRow),])>0) {
        aimRow <- min(chrAIMdata$chrRow[which(chrAIMdata$aim5quant=="upper" & chrAIMdata$chrRow>aimRow)])
      } else {
        aimRow <- nrow(chrAIMdata) + 1
      }
    }
  }
  print(chr)
  #print(dim(islandDF))
}
print("Islands done")
print(dim(islandDF))
write.table(islandDF, paste(popName, "_minorParentIslands_AIMs_cMwins_raw.txt", sep=""), quote=F, row.names = F, sep="\t")


##Merge any within 50000
##Filter less than 10000 
##Filter less than 10 SNPs and 10 AIM

regions <- desertDF

regions$distToPrev <- NA
regions$distToNext <- NA
uniChrs <- unique(regions$chrName)
for (chr in uniChrs) {
  chrData <- subset(regions, chrName==chr)
  preEnd <- NA
  nextStart <- NA
  for (i in 1:nrow(chrData)) {
    nextStart <- chrData$aimStart[i+1]
    regInfo <- chrData[i,]
    regStart <- regInfo$aimStart
    regEnd <- regInfo$aimEnd
    preDist <- regStart-preEnd
    nextDist <- nextStart-regEnd
    regions$distToPrev[which(regions$chrName==chr & regions$aimStart==regStart & regions$aimEnd==regEnd)] <- preDist
    regions$distToNext[which(regions$chrName==chr & regions$aimStart==regStart & regions$aimEnd==regEnd)] <- nextDist
    preEnd <- regEnd
  }
}

mergeDF <- data.frame()
i = 1
while(i < nrow(regions)) {
  regDist <- regions$distToNext[i]
  j <- i
  if (!is.na(regDist) & regDist < 50000) {
    #print(i)
    mergeRegs <- regions[i,]
    nextRegDist <- regions$distToNext[i+(j-i)]
    while (!is.na(nextRegDist) & nextRegDist < 50000) {
      #print(nextRegDist)
      j <- j + 1
      nextRegDist <- regions$distToNext[i+(j-i)]
      mergeRegs <- rbind(mergeRegs, regions[i+(j-i),])
    }
    #print(paste(as.character(regions$regName[i]), "to merge", (j-i), sep=" "))
    mergeChr <- as.character(mergeRegs$chrName[1])
    new.cMstart <- min(mergeRegs$cMstart)
    new.cMend <- max(mergeRegs$cMend)
    new.cMlen <- new.cMend-new.cMstart
    new.cMmid <- as.integer(new.cMstart+(new.cMlen/2))
    new.cMsnps <- sum(mergeRegs$cMsnps)
    new.cMaims <- subset(aimData, group==as.character(mergeChr) & position>=new.cMstart & position<=new.cMend)
    new.cMpropAims <- nrow(subset(new.cMaims, aim10quant=="lower"))/nrow(new.cMaims)
    new.minAnc <- mean(new.cMaims$minPar, na.rm=T)
    #newRegName <- paste(mergeChr, ":", newStart, "_", newEnd, sep="")
    new.aimStart <- min(mergeRegs$aimStart)
    new.aimEnd <- max(mergeRegs$aimEnd)
    new.aimLen <- new.aimEnd-new.aimStart
    new.aimMid <- as.integer(new.aimStart+(new.aimLen/2))
    newRegAIMs <- subset(aimData, group==as.character(mergeChr) & position>=new.aimStart & position<=new.aimEnd)
    #newHI <- mean(newRegAIMs$hybrid_index, na.rm=T)
    newNumAIMs <- nrow(newRegAIMs)
    newPropAIMs <- nrow(subset(newRegAIMs, aim10quant=="lower"))/nrow(newRegAIMs)
    newPreDist <- NA
    if (i!=1) {
      if (as.character(mergeChr) == regions$chrName[i-1]) {
        newPreDist <- new.aimStart - regions$aimEnd[i-1]
        if (newPreDist==0) {
          #print(regions[i-1,])
        }
      }
    }
    newNextDist <- NA
    if (j+1 <= nrow(regions)) {
      if (as.character(mergeChr) == regions$chrName[j+1]) {
        newNextDist <- regions$aimStart[j+1] - new.aimEnd
      }
    }
    newRegDF <- data.frame(chrName=as.character(mergeChr), cMstart=new.cMstart, cMend=new.cMend, cMmid=new.cMmid, cMlen=new.cMlen, cMsnps=new.cMsnps, cMpropAIMs=new.cMpropAims, minAnc=new.minAnc, aimStart=new.aimStart, aimEnd=new.aimEnd, aimLen=new.aimLen, aimMid=new.aimMid, nAIMs=newNumAIMs, propAIMs=newPropAIMs, distToPrev=newPreDist, distToNext=newNextDist)
    mergeDF <- rbind(mergeDF, newRegDF)
    i = j + 1
  } else {
    mergeDF <- rbind(mergeDF, regions[i,])
    unMergeDist <- regions$distToPrev[i]
    if (!is.na(unMergeDist) & unMergeDist==0) {
      print(i)
    }
    i = i+1
  }
}

toKeep <- subset(mergeDF, aimLen>10000)
toKeepAIMs <- subset(toKeep, nAIMs>=10)
toKeepSNPs <- subset(toKeepAIMs, cMsnps>=10)
write.table(toKeepSNPs, paste(popName, "_minorParentDeserts_AIMs_cMwins_merged.txt", sep=""), row.names=F, quote=F, sep="\t")
print(popName)
print("Deserts")
print(dim(regions))
print(dim(mergeDF))
print(dim(toKeepSNPs))
print(summary(toKeepSNPs))


###And again for islands
regions <- islandDF

regions$distToPrev <- NA
regions$distToNext <- NA
uniChrs <- unique(regions$chrName)
for (chr in uniChrs) {
  chrData <- subset(regions, chrName==chr)
  preEnd <- NA
  nextStart <- NA
  for (i in 1:nrow(chrData)) {
    nextStart <- chrData$aimStart[i+1]
    regInfo <- chrData[i,]
    regStart <- regInfo$aimStart
    regEnd <- regInfo$aimEnd
    preDist <- regStart-preEnd
    nextDist <- nextStart-regEnd
    regions$distToPrev[which(regions$chrName==chr & regions$aimStart==regStart & regions$aimEnd==regEnd)] <- preDist
    regions$distToNext[which(regions$chrName==chr & regions$aimStart==regStart & regions$aimEnd==regEnd)] <- nextDist
    preEnd <- regEnd
  }
}

mergeDF <- data.frame()
i = 1
while(i < nrow(regions)) {
  regDist <- regions$distToNext[i]
  j <- i
  if (!is.na(regDist) & regDist < 50000) {
    #print(i)
    mergeRegs <- regions[i,]
    nextRegDist <- regions$distToNext[i+(j-i)]
    while (!is.na(nextRegDist) & nextRegDist < 50000) {
      #print(nextRegDist)
      j <- j + 1
      nextRegDist <- regions$distToNext[i+(j-i)]
      mergeRegs <- rbind(mergeRegs, regions[i+(j-i),])
    }
    #print(paste(as.character(regions$regName[i]), "to merge", (j-i), sep=" "))
    mergeChr <- as.character(mergeRegs$chrName[1])
    new.cMstart <- min(mergeRegs$cMstart)
    new.cMend <- max(mergeRegs$cMend)
    new.cMlen <- new.cMend-new.cMstart
    new.cMmid <- as.integer(new.cMstart+(new.cMlen/2))
    new.cMsnps <- sum(mergeRegs$cMsnps)
    new.cMaims <- subset(aimData, group==as.character(mergeChr) & position>=new.cMstart & position<=new.cMend)
    new.cMpropAims <- nrow(subset(new.cMaims, aim10quant=="upper"))/nrow(new.cMaims)
    new.minAnc <- mean(new.cMaims$minPar, na.rm=T)
    #newRegName <- paste(mergeChr, ":", newStart, "_", newEnd, sep="")
    new.aimStart <- min(mergeRegs$aimStart)
    new.aimEnd <- max(mergeRegs$aimEnd)
    new.aimLen <- new.aimEnd-new.aimStart
    new.aimMid <- as.integer(new.aimStart+(new.aimLen/2))
    newRegAIMs <- subset(aimData, group==as.character(mergeChr) & position>=new.aimStart & position<=new.aimEnd)
    #newHI <- mean(newRegAIMs$hybrid_index, na.rm=T)
    newNumAIMs <- nrow(newRegAIMs)
    newPropAIMs <- nrow(subset(newRegAIMs, aim10quant=="upper"))/nrow(newRegAIMs)
    newPreDist <- NA
    if (i!=1) {
      if (as.character(mergeChr) == regions$chrName[i-1]) {
        newPreDist <- new.aimStart - regions$aimEnd[i-1]
        if (newPreDist==0) {
          #print(regions[i-1,])
        }
      }
    }
    newNextDist <- NA
    if (j+1 <= nrow(regions)) {
      if (as.character(mergeChr) == regions$chrName[j+1]) {
        newNextDist <- regions$aimStart[j+1] - new.aimEnd
      }
    }
    newRegDF <- data.frame(chrName=as.character(mergeChr), cMstart=new.cMstart, cMend=new.cMend, cMmid=new.cMmid, cMlen=new.cMlen, cMsnps=new.cMsnps, cMpropAIMs=new.cMpropAims, minAnc=new.minAnc, aimStart=new.aimStart, aimEnd=new.aimEnd, aimLen=new.aimLen, aimMid=new.aimMid, nAIMs=newNumAIMs, propAIMs=newPropAIMs, distToPrev=newPreDist, distToNext=newNextDist)
    mergeDF <- rbind(mergeDF, newRegDF)
    i = j + 1
  } else {
    mergeDF <- rbind(mergeDF, regions[i,])
    unMergeDist <- regions$distToPrev[i]
    if (!is.na(unMergeDist) & unMergeDist==0) {
      print(i)
    }
    i = i+1
  }
}

toKeep <- subset(mergeDF, aimLen>10000)
toKeepAIMs <- subset(toKeep, nAIMs>=10)
toKeepSNPs <- subset(toKeepAIMs, cMsnps>=10)
write.table(toKeepSNPs, paste(popName, "_minorParentIslands_AIMs_cMwins_merged.txt", sep=""), row.names=F, quote=F, sep="\t")
print(popName)
print("Islands")
print(dim(regions))
print(dim(mergeDF))
print(dim(toKeepSNPs))
print(summary(toKeepSNPs))



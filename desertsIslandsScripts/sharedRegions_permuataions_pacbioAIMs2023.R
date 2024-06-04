options(stringsAsFactors = F)

setwd("~/Documents/SchumerLab/Xcortezi_Xbirchmanni_hybridization/desertsIslandsPacBio2023/")
popsXcorXbir <- c("CHPL2021", "STAC2020")
popsXmalXbir <- c("AGZC", "ACUA2018", "TLMC2017")
#chrKey <- read.table("~/Documents/SchumerLab/match_chromnames_xbir-10x_xbir-pacbio_xmal-10x_chroms.txt", header=T)

##Thinned
quants <- read.table("minorParentQuantiles_thinnedAIMs_allPops_0.05cMwins.txt", header=T)
#quants$pop <- data.frame(do.call('rbind', strsplit(quants$pop, "_")))$X2

winCM <- read.table("xbir_pacbio2023_allChrs_0.05cM_windows_recRate_codingBP_conservedBP_thinnedAIMs_everything_minPar.txt", header=T)
#winCM$CHPL2021 <- winCM$thinned_CHPL2021
#winCM$STAC2020 <- winCM$thinned_STAC2020
#winCM$AGZC <- winCM$thinned_AGZC
#winCM$ACUA2018 <- winCM$thinned_ACUA2018
#winCM$TLMC2017 <- winCM$thinned_TLMC2017

chplWin <- data.frame(chrID=winCM$chrID, start=winCM$start, end=winCM$end, chrPos=paste(winCM$chrID, winCM$start, sep="_"), minAnc=winCM$CHPL2021)
chplWin$regType <- NA

chplDes <- read.table("CHPL2021_minorParentDeserts_AIMs_cMwins_merged.txt", header=T)
chplDes$regType <- "Desert"
chplDes$cMmidChrPos <- NA
chplDes$cMmidPass <- NA
chplDes$nWin <- NA
chplDes$propWin <- NA
for (i in 1:nrow(chplDes)) {
  regWins <- subset(chplWin, chrID==chplDes$chrName[i] & start>=chplDes$cMstart[i] & end<=chplDes$cMend[i])
  chplWin$regType[which(chplWin$chrID==chplDes$chrName[i] & chplWin$start>=chplDes$cMstart[i] & chplWin$end<=chplDes$cMend[i])] <- "Desert"
  midWin <- subset(regWins, start<=chplDes$cMmid[i] & end>=chplDes$cMmid[i])
  cutoff <- subset(quants, pop=="CHPL2021")$lower10
  chplDes$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  } else if (midWin$minAnc<=cutoff) {
    pass <- TRUE
  } else if (midWin$minAnc>cutoff) {
    pass <- FALSE
  }
  chplDes$cMmidPass[i] <- pass
  chplDes$nWin[i] <- nrow(regWins)
  chplDes$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc<=cutoff))
}
dim(chplDes)
nrow(subset(chplDes, is.na(cMmidPass)))
sum(chplDes$cMmidPass, na.rm=T)

chplIsle <- read.table("CHPL2021_minorParentIslands_AIMs_cMwins_merged.txt", header=T)
chplIsle$regType <- "Island"
chplIsle$cMmidChrPos <- NA
chplIsle$cMmidPass <- NA
chplIsle$nWin <- NA
chplIsle$propWin <- NA
for (i in 1:nrow(chplIsle)) {
  regWins <- subset(chplWin, chrID==chplIsle$chrName[i] & start>=chplIsle$cMstart[i] & end<=chplIsle$cMend[i])
  chplWin$regType[which(chplWin$chrID==chplIsle$chrName[i] & chplWin$start>=chplIsle$cMstart[i] & chplWin$end<=chplIsle$cMend[i])] <- "Island"
  midWin <- subset(regWins,  start<=chplIsle$cMmid[i] & end>=chplIsle$cMmid[i])
  cutoff <- subset(quants, pop=="CHPL2021")$upper10
  chplIsle$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  } else if (midWin$minAnc>=cutoff) {
    pass <- TRUE
  } else if (midWin$minAnc<cutoff) {
    pass <- FALSE
  }
  chplIsle$cMmidPass[i] <- pass
  chplIsle$nWin[i] <- nrow(regWins)
  chplIsle$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc>=cutoff))
}
dim(chplIsle)
nrow(subset(chplIsle, is.na(cMmidPass)))
sum(chplIsle$cMmidPass, na.rm=T)


stacWin <- data.frame(chrID=winCM$chrID, start=winCM$start, end=winCM$end, chrPos=paste(winCM$chrID, winCM$start, sep="_"), minAnc=winCM$STAC2020, regType=NA)

stacDes <- read.table("STAC2020_minorParentDeserts_AIMs_cMwins_merged.txt", header=T)
stacDes$regType <- "Desert"
stacDes$cMmidChrPos <- NA
stacDes$cMmidPass <- NA
stacDes$nWin <- NA
stacDes$propWin <- NA
for (i in 1:nrow(stacDes)) {
  regWins <- subset(stacWin, chrID==stacDes$chrName[i] & start>=stacDes$cMstart[i] & end<=stacDes$cMend[i])
  stacWin$regType[which(stacWin$chrID==stacDes$chrName[i] & stacWin$start>=stacDes$cMstart[i] & stacWin$end<=stacDes$cMend[i])] <- "Desert"
  midWin <- subset(regWins,  start<=stacDes$cMmid[i] & end>=stacDes$cMmid[i])
  cutoff <- subset(quants, pop=="STAC2020")$lower10
  stacDes$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  } else if (midWin$minAnc<=cutoff) {
    pass <- TRUE
  } else if (midWin$minAnc>cutoff) {
    pass <- FALSE
  }
  stacDes$cMmidPass[i] <- pass
  stacDes$nWin[i] <- nrow(regWins)
  stacDes$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc<=cutoff))
}
dim(stacDes)
nrow(subset(stacDes, is.na(cMmidPass)))
sum(stacDes$cMmidPass, na.rm=T)

stacIsle <- read.table("STAC2020_minorParentIslands_AIMs_cMwins_merged.txt", header=T)
stacIsle$regType <- "Island"
stacIsle$cMmidChrPos <- NA
stacIsle$cMmidPass <- NA
stacIsle$nWin <- NA
stacIsle$propWin <- NA
for (i in 1:nrow(stacIsle)) {
  regWins <- subset(stacWin, chrID==stacIsle$chrName[i] & start>=stacIsle$cMstart[i] & end<=stacIsle$cMend[i])
  stacWin$regType[which(stacWin$chrID==stacIsle$chrName[i] & stacWin$start>=stacIsle$cMstart[i] & stacWin$end<=stacIsle$cMend[i])] <- "Island"
  midWin <- subset(regWins,  start<=stacIsle$cMmid[i] & end>=stacIsle$cMmid[i])
  cutoff <- subset(quants, pop=="STAC2020")$upper10
  stacIsle$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  } else if (midWin$minAnc>=cutoff) {
    pass <- TRUE
  } else if (midWin$minAnc<cutoff) {
    pass <- FALSE
  }
  stacIsle$cMmidPass[i] <- pass
  stacIsle$nWin[i] <- nrow(regWins)
  stacIsle$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc>=cutoff))
}
dim(stacIsle)
nrow(subset(stacIsle, is.na(cMmidPass)))
sum(stacIsle$cMmidPass, na.rm=T)


agzcWin <- data.frame(chrID=winCM$chrID, start=winCM$start, end=winCM$end, chrPos=paste(winCM$chrID, winCM$start, sep="_"), minAnc=winCM$AGZC, regType=NA)

agzcDes <- read.table("AGZC_minorParentDeserts_AIMs_cMwins_merged.txt", header=T)
agzcDes$regType <- "Desert"
agzcDes$cMmidChrPos <- NA
agzcDes$cMmidPass <- NA
agzcDes$nWin <- NA
agzcDes$propWin <- NA
for (i in 1:nrow(agzcDes)) {
  regWins <- subset(agzcWin, chrID==agzcDes$chrName[i] & start>=agzcDes$cMstart[i] & end<=agzcDes$cMend[i])
  agzcWin$regType[which(agzcWin$chrID==agzcDes$chrName[i] & agzcWin$start>=agzcDes$cMstart[i] & agzcWin$end<=agzcDes$cMend[i])] <- "Desert"
  midWin <- subset(regWins,  start<=agzcDes$cMmid[i] & end>=agzcDes$cMmid[i])
  cutoff <- subset(quants, pop=="AGZC")$lower10
  agzcDes$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  } else if (midWin$minAnc<=cutoff) {
    pass <- TRUE
  } else if (midWin$minAnc>cutoff) {
    pass <- FALSE
  }
  agzcDes$cMmidPass[i] <- pass
  agzcDes$nWin[i] <- nrow(regWins)
  agzcDes$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc<=cutoff))
}
dim(agzcDes)
nrow(subset(agzcDes, is.na(cMmidPass)))
sum(agzcDes$cMmidPass, na.rm=T)

agzcIsle <- read.table("AGZC_minorParentIslands_AIMs_cMwins_merged.txt", header=T)
agzcIsle$regType <- "Island"
agzcIsle$cMmidChrPos <- NA
agzcIsle$cMmidPass <- NA
agzcIsle$nWin <- NA
agzcIsle$propWin <- NA
for (i in 1:nrow(agzcIsle)) {
  regWins <- subset(agzcWin, chrID==agzcIsle$chrName[i] & start>=agzcIsle$cMstart[i] & end<=agzcIsle$cMend[i])
  agzcWin$regType[which(agzcWin$chrID==agzcIsle$chrName[i] & agzcWin$start>=agzcIsle$cMstart[i] & agzcWin$end<=agzcIsle$cMend[i])] <- "Island"
  midWin <- subset(regWins,  start<=agzcIsle$cMmid[i] & end>=agzcIsle$cMmid[i])
  cutoff <- subset(quants, pop=="AGZC")$upper10
  agzcIsle$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  } else if (midWin$minAnc>=cutoff) {
    pass <- TRUE
  } else if (midWin$minAnc<cutoff) {
    pass <- FALSE
  }
  agzcIsle$cMmidPass[i] <- pass
  agzcIsle$nWin[i] <- nrow(regWins)
  agzcIsle$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc>=cutoff))
}
dim(agzcIsle)
nrow(subset(agzcIsle, is.na(cMmidPass)))
sum(agzcIsle$cMmidPass, na.rm=T)


acuaWin <- data.frame(chrID=winCM$chrID, start=winCM$start, end=winCM$end, chrPos=paste(winCM$chrID, winCM$start, sep="_"), minAnc=winCM$ACUA2018, regType=NA)

acuaDes <- read.table("ACUA2018_minorParentDeserts_AIMs_cMwins_merged.txt", header=T)
acuaDes$regType <- "Desert"
acuaDes$cMmidChrPos <- NA
acuaDes$cMmidPass <- NA
acuaDes$nWin <- NA
acuaDes$propWin <- NA
for (i in 1:nrow(acuaDes)) {
  regWins <- subset(acuaWin, chrID==acuaDes$chrName[i] & start>=acuaDes$cMstart[i] & end<=acuaDes$cMend[i])
  acuaWin$regType[which(acuaWin$chrID==acuaDes$chrName[i] & acuaWin$start>=acuaDes$cMstart[i] & acuaWin$end<=acuaDes$cMend[i])] <- "Desert"
  midWin <- subset(regWins,  start<=acuaDes$cMmid[i] & end>=acuaDes$cMmid[i])
  cutoff <- subset(quants, pop=="ACUA2018")$lower10
  acuaDes$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  } else if (midWin$minAnc<=cutoff) {
    pass <- TRUE
  } else if (midWin$minAnc>cutoff) {
    pass <- FALSE
  }
  acuaDes$cMmidPass[i] <- pass
  acuaDes$nWin[i] <- nrow(regWins)
  acuaDes$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc<=cutoff))
}
dim(acuaDes)
nrow(subset(acuaDes, is.na(cMmidPass)))
sum(acuaDes$cMmidPass, na.rm=T)

acuaIsle <- read.table("ACUA2018_minorParentIslands_AIMs_cMwins_merged.txt", header=T)
acuaIsle$regType <- "Island"
acuaIsle$cMmidChrPos <- NA
acuaIsle$cMmidPass <- NA
acuaIsle$nWin <- NA
acuaIsle$propWin <- NA
for (i in 1:nrow(acuaIsle)) {
  regWins <- subset(acuaWin, chrID==acuaIsle$chrName[i] & start>=acuaIsle$cMstart[i] & end<=acuaIsle$cMend[i])
  acuaWin$regType[which(acuaWin$chrID==acuaIsle$chrName[i] & acuaWin$start>=acuaIsle$cMstart[i] & acuaWin$end<=acuaIsle$cMend[i])] <- "Island"
  midWin <- subset(regWins,  start<=acuaIsle$cMmid[i] & end>=acuaIsle$cMmid[i])
  cutoff <- subset(quants, pop=="ACUA2018")$upper10
  acuaIsle$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  } else if (midWin$minAnc>=cutoff) {
    pass <- TRUE
  } else if (midWin$minAnc<cutoff) {
    pass <- FALSE
  }
  acuaIsle$cMmidPass[i] <- pass
  acuaIsle$nWin[i] <- nrow(regWins)
  acuaIsle$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc>=cutoff))
}
dim(acuaIsle)
nrow(subset(acuaIsle, is.na(cMmidPass)))
sum(acuaIsle$cMmidPass, na.rm=T)


tlmcWin <- data.frame(chrID=winCM$chrID, start=winCM$start, end=winCM$end, chrPos=paste(winCM$chrID, winCM$start, sep="_"), minAnc=winCM$TLMC2017, regType=NA)

tlmcDes <- read.table("TLMC2017_minorParentDeserts_AIMs_cMwins_merged.txt", header=T)
tlmcDes$regType <- "Desert"
tlmcDes$cMmidChrPos <- NA
tlmcDes$cMmidPass <- NA
tlmcDes$nWin <- NA
tlmcDes$propWin <- NA
for (i in 1:nrow(tlmcDes)) {
  regWins <- subset(tlmcWin, chrID==tlmcDes$chrName[i] & start>=tlmcDes$cMstart[i] & end<=tlmcDes$cMend[i])
  tlmcWin$regType[which(tlmcWin$chrID==tlmcDes$chrName[i] & tlmcWin$start>=tlmcDes$cMstart[i] & tlmcWin$end<=tlmcDes$cMend[i])] <- "Desert"
  midWin <- subset(regWins,  start<=tlmcDes$cMmid[i] & end>=tlmcDes$cMmid[i])
  cutoff <- subset(quants, pop=="TLMC2017")$lower10
  tlmcDes$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  } else if (midWin$minAnc<=cutoff) {
    pass <- TRUE
  } else if (midWin$minAnc>cutoff) {
    pass <- FALSE
  }
  tlmcDes$cMmidPass[i] <- pass
  tlmcDes$nWin[i] <- nrow(regWins)
  tlmcDes$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc<=cutoff))
}
dim(tlmcDes)
nrow(subset(tlmcDes, is.na(cMmidPass)))
sum(tlmcDes$cMmidPass, na.rm=T)

tlmcIsle <- read.table("TLMC2017_minorParentIslands_AIMs_cMwins_merged.txt", header=T)
tlmcIsle$regType <- "Island"
tlmcIsle$cMmidChrPos <- NA
tlmcIsle$cMmidPass <- NA
tlmcIsle$nWin <- NA
tlmcIsle$propWin <- NA
for (i in 1:nrow(tlmcIsle)) {
  regWins <- subset(tlmcWin, chrID==tlmcIsle$chrName[i] & start>=tlmcIsle$cMstart[i] & end<=tlmcIsle$cMend[i])
  tlmcWin$regType[which(tlmcWin$chrID==tlmcIsle$chrName[i] & tlmcWin$start>=tlmcIsle$cMstart[i] & tlmcWin$end<=tlmcIsle$cMend[i])] <- "Island"
  midWin <- subset(regWins,  start<=tlmcIsle$cMmid[i] & end>=tlmcIsle$cMmid[i])
  cutoff <- subset(quants, pop=="TLMC2017")$upper10
  tlmcIsle$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  } else if (midWin$minAnc>=cutoff) {
    pass <- TRUE
  } else if (midWin$minAnc<cutoff) {
    pass <- FALSE
  }
  tlmcIsle$cMmidPass[i] <- pass
  tlmcIsle$nWin[i] <- nrow(regWins)
  tlmcIsle$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc>=cutoff))
}
dim(tlmcIsle)
nrow(subset(tlmcIsle, is.na(cMmidPass)))
sum(tlmcIsle$cMmidPass, na.rm=T)

agzcDes$propPass <- agzcDes$propWin/agzcDes$nWin
agzcIsle$propPass <- agzcIsle$propWin/agzcIsle$nWin
stacDes$propPass <- stacDes$propWin/stacDes$nWin
stacIsle$propPass <- stacIsle$propWin/stacIsle$nWin
agzcDes$propPass <- agzcDes$propWin/agzcDes$nWin
agzcIsle$propPass <- agzcIsle$propWin/agzcIsle$nWin
acuaDes$propPass <- acuaDes$propWin/acuaDes$nWin
acuaIsle$propPass <- acuaIsle$propWin/acuaIsle$nWin
tlmcDes$propPass <- tlmcDes$propWin/tlmcDes$nWin
tlmcIsle$propPass <- tlmcIsle$propWin/tlmcIsle$nWin

summary(chplDes)
summary(chplIsle)
summary(stacDes)
summary(stacIsle)
summary(agzcDes)
summary(agzcIsle)
summary(acuaDes)
summary(acuaIsle)
summary(tlmcDes)
summary(tlmcIsle)

agzcDesPass <- subset(agzcDes, cMmidPass==T)
agzcIslePass <- subset(agzcIsle, cMmidPass==T)
stacDesPass <- subset(stacDes, cMmidPass==T)
stacIslePass <- subset(stacIsle, cMmidPass==T)

sum(agzcDes$cMmidChrPos %in% stacDes$cMmidChrPos)
sum(agzcIsle$cMmidChrPos %in% stacIsle$cMmidChrPos)

sum(stacDes$cMmidChrPos %in% agzcDes$cMmidChrPos)
sum(stacDes$cMmidChrPos %in% agzcDes$cMmidChrPos)
sum(stacDes$cMmidChrPos %in% acuaDes$cMmidChrPos)
sum(stacDes$cMmidChrPos %in% tlmcDes$cMmidChrPos)
sum(stacIsle$cMmidChrPos %in% agzcIsle$cMmidChrPos)
sum(stacIsle$cMmidChrPos %in% agzcIsle$cMmidChrPos)
sum(stacIsle$cMmidChrPos %in% acuaIsle$cMmidChrPos)
sum(stacIsle$cMmidChrPos %in% tlmcIsle$cMmidChrPos)

sum(agzcDes$cMmidChrPos %in% agzcDes$cMmidChrPos)
sum(agzcDes$cMmidChrPos %in% acuaDes$cMmidChrPos)
sum(agzcDes$cMmidChrPos %in% tlmcDes$cMmidChrPos)
sum(agzcIsle$cMmidChrPos %in% agzcIsle$cMmidChrPos)
sum(agzcIsle$cMmidChrPos %in% acuaIsle$cMmidChrPos)
sum(agzcIsle$cMmidChrPos %in% tlmcIsle$cMmidChrPos)

chplDes$stacN <- NA
for (i in 1:nrow(chplDes)) {
  stacRegWins <- subset(stacWin, chrID==chplDes$chrName[i] & start>=chplDes$cMstart[i] & end<=chplDes$cMend[i])
  chplDes$stacN[i] <- nrow(subset(stacRegWins, regType=="Desert"))
}
nrow(chplDes)
nrow(subset(chplDes, stacN>0))
sum(subset(chplDes, stacN>0)$cMmidPass, na.rm=T)

chplIsle$stacN <- NA
for (i in 1:nrow(chplIsle)) {
  stacRegWins <- subset(stacWin, chrID==chplIsle$chrName[i] & start>=chplIsle$cMstart[i] & end<=chplIsle$cMend[i])
  chplIsle$stacN[i] <- nrow(subset(stacRegWins, regType=="Island"))
}
nrow(chplIsle)
nrow(subset(chplIsle, stacN>0))
sum(subset(chplIsle, stacN>0)$cMmidPass, na.rm=T)

stacDes$chplN <- NA
for (i in 1:nrow(stacDes)) {
  chplRegWins <- subset(chplWin, chrID==stacDes$chrName[i] & start>=stacDes$cMstart[i] & end<=stacDes$cMend[i])
  stacDes$chplN[i] <- nrow(subset(chplRegWins, regType=="Desert"))
}
nrow(stacDes)
nrow(subset(stacDes, chplN>0))
sum(subset(stacDes, chplN>0)$cMmidPass, na.rm=T)

stacIsle$chplN <- NA
for (i in 1:nrow(stacIsle)) {
  chplRegWins <- subset(chplWin, chrID==stacIsle$chrName[i] & start>=stacIsle$cMstart[i] & end<=stacIsle$cMend[i])
  stacIsle$chplN[i] <- nrow(subset(chplRegWins, regType=="Island"))
}
nrow(stacIsle)
nrow(subset(stacIsle, chplN>0))
head(subset(stacIsle, chplN>0))
sum(subset(stacIsle, chplN>0)$cMmidPass, na.rm=T)

chplWin$quant <- NA
chplWin$quant[which(chplWin$minAnc<=quants$lower10[which(quants$pop=="CHPL2021")])] <- "lower10"
chplWin$quant[which(chplWin$minAnc>=quants$upper10[which(quants$pop=="CHPL2021")])] <- "upper10"

stacWin$quant <- NA
stacWin$quant[which(stacWin$minAnc<=quants$lower10[which(quants$pop=="STAC2020")])] <- "lower10"
stacWin$quant[which(stacWin$minAnc>=quants$upper10[which(quants$pop=="STAC2020")])] <- "upper10"

agzcWin$quant <- NA
agzcWin$quant[which(agzcWin$minAnc<=quants$lower10[which(quants$pop=="AGZC")])] <- "lower10"
agzcWin$quant[which(agzcWin$minAnc>=quants$upper10[which(quants$pop=="AGZC")])] <- "upper10"

acuaWin$quant <- NA
acuaWin$quant[which(acuaWin$minAnc<=quants$lower10[which(quants$pop=="ACUA2018")])] <- "lower10"
acuaWin$quant[which(acuaWin$minAnc>=quants$upper10[which(quants$pop=="ACUA2018")])] <- "upper10"

tlmcWin$quant <- NA
tlmcWin$quant[which(tlmcWin$minAnc<=quants$lower10[which(quants$pop=="TLMC2017")])] <- "lower10"
tlmcWin$quant[which(tlmcWin$minAnc>=quants$upper10[which(quants$pop=="TLMC2017")])] <- "upper10"

stacDesRows <- stacWin[which(stacWin$chrPos %in% stacDes$cMmidChrPos[which(stacDes$cMmidPass=="TRUE" & stacDes$chplN>0)]),]
focal<-stacDesRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-chplWin
  pop_null$quant<-sample(chplWin$quant)
  pass<-subset(pop_null,pop_null$quant=="lower10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
stacDesPermut <- data.frame(popComp="STAC_Desert", focalPop="STAC", regType="Desert", compPop='CHPL', permutations=null_positives)
print(summary(stacDesPermut))

stacIsleRows <- stacWin[which(stacWin$chrPos %in% stacIsle$cMmidChrPos[which(stacIsle$cMmidPass=="TRUE" & stacIsle$chplN>0)]),]
focal<-stacIsleRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-chplWin
  pop_null$quant<-sample(chplWin$quant)
  pass<-subset(pop_null,pop_null$quant=="upper10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
stacIslePermut <- data.frame(popComp="STAC_Island", focalPop="STAC", regType="Island", compPop='CHPL', permutations=null_positives)
print(summary(stacIslePermut))


chplDesRows <- chplWin[which(chplWin$chrPos %in% chplDes$cMmidChrPos[which(chplDes$cMmidPass=="TRUE" & chplDes$stacN>0)]),]
focal<-chplDesRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-stacWin
  pop_null$quant<-sample(stacWin$quant)
  pass<-subset(pop_null,pop_null$quant=="lower10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
chplDesPermut <- data.frame(popComp="CHPL_Desert", focalPop="CHPL", regType="Desert", compPop='STAC', permutations=null_positives)
print(summary(chplDesPermut))

chplIsleRows <- chplWin[which(chplWin$chrPos %in% chplIsle$cMmidChrPos[which(chplIsle$cMmidPass=="TRUE" & chplIsle$stacN>0)]),]
focal<-chplIsleRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-stacWin
  pop_null$quant<-sample(stacWin$quant)
  pass<-subset(pop_null,pop_null$quant=="upper10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
chplIslePermut <- data.frame(popComp="CHPL_Island", focalPop="CHPL", regType="Island", compPop='STAC', permutations=null_positives)
print(summary(chplIslePermut))

comboShare <- read.table("combinedShared_CHPL_STAC_crossType_desertsIslands.txt", header=T)
comboShare$cMmidChrPos <- NA
comboShare$stacMidPass <- NA
comboShare$chplMidPass <- NA
for (i in 1:nrow(comboShare)) {
  stacMidWin <- subset(stacWin, chrID==comboShare$chrName[i] & start<=comboShare$cMmid[i] & end>=comboShare$cMmid[i])
  chplMidWin <- subset(chplWin, chrID==comboShare$chrName[i] & start<=comboShare$cMmid[i] & end>=comboShare$cMmid[i])
  comboShare$cMmidChrPos[i] <- paste(stacMidWin$chrID, stacMidWin$start, sep="_")
  stacPass <- NA
  if (is.na(stacMidWin$minAnc)) {
    stacPass <- NA
  } else if (comboShare$regType[i]=="Desert" & stacMidWin$quant=="lower10" & !is.na(stacMidWin$quant)) {
    stacPass <- TRUE
  } else if (comboShare$regType[i]=="Island" & stacMidWin$quant=="upper10" & !is.na(stacMidWin$quant)) {
    stacPass <- TRUE
  } else if (comboShare$regType[i]=="Desert" & (stacMidWin$quant!="lower10" | is.na(stacMidWin$quant))) {
    stacPass <- FALSE
  } else if (comboShare$regType[i]=="Island" & (stacMidWin$quant!="upper10" | is.na(stacMidWin$quant))) {
    stacPass <- FALSE
  }
  chplPass <- NA
  if (is.na(chplMidWin$minAnc)) {
    chplPass <- NA
  } else if (comboShare$regType[i]=="Desert" & chplMidWin$quant=="lower10" & !is.na(chplMidWin$quant)) {
    chplPass <- TRUE
  } else if (comboShare$regType[i]=="Island" & chplMidWin$quant=="upper10" & !is.na(chplMidWin$quant)) {
    chplPass <- TRUE
  } else if (comboShare$regType[i]=="Desert" & (chplMidWin$quant!="lower10" | is.na(chplMidWin$quant))) {
    chplPass <- FALSE
  } else if (comboShare$regType[i]=="Island" & (chplMidWin$quant!="upper10" | is.na(chplMidWin$quant))) {
    chplPass <- FALSE
  }
  comboShare$stacMidPass[i] <- stacPass
  comboShare$chplMidPass[i] <- chplPass
}

comboDes <- subset(comboShare, regType=="Desert")
comboIsle <- subset(comboShare, regType=="Island")
write.table(comboShare, "thinnedCombinedShared_CHPL_STAC_desertsIslands_midCMwin.txt", row.name=F, quote=F, sep="\t")

comboDesPass <- subset(comboDes, stacMidPass==T & chplMidPass==T)
comboIslePass <- subset(comboIsle, stacMidPass==T & chplMidPass==T)

nrow(subset(comboDesPass, nAGZC>0))
nrow(subset(comboDesPass, nACUA>0))
nrow(subset(comboDesPass, nTLMC>0))
nrow(subset(comboDesPass, nAGZC>0 & nACUA>0))
nrow(subset(comboDesPass, nAGZC>0 & nACUA>0 & nTLMC>0))
nrow(subset(comboDesPass, nAGZC>0 | nACUA>0 | nTLMC>0))

nrow(subset(comboIslePass, nAGZC>0))
nrow(subset(comboIslePass, nACUA>0))
nrow(subset(comboIslePass, nTLMC>0))
nrow(subset(comboIslePass, nAGZC>0 & nACUA>0))
nrow(subset(comboIslePass, nAGZC>0 & nACUA>0 & nTLMC>0))
nrow(subset(comboIslePass, nAGZC>0 | nACUA>0 | nTLMC>0))


comboAgzcDesRows <- stacWin[which(stacWin$chrPos %in% comboDesPass$cMmidChrPos[which(comboDesPass$nAGZC>0)]),]
focal<-comboAgzcDesRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-agzcWin
  pop_null$quant<-sample(agzcWin$quant)
  pass<-subset(pop_null,pop_null$quant=="lower10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAgzcDesPermut <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop='AGZC', permutations=null_positives)
print(summary(comboAgzcDesPermut))

comboAgzcIsleRows <- stacWin[which(stacWin$chrPos %in% comboIslePass$cMmidChrPos[which(comboIslePass$nAGZC>0)]),]
focal<-comboAgzcIsleRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-agzcWin
  pop_null$quant<-sample(agzcWin$quant)
  pass<-subset(pop_null,pop_null$quant=="upper10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAgzcIslePermut <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop='AGZC', permutations=null_positives)
print(summary(comboAgzcIslePermut))

comboAcuaDesRows <- stacWin[which(stacWin$chrPos %in% comboDesPass$cMmidChrPos[which(comboDesPass$nACUA>0)]),]
focal<-comboAcuaDesRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-acuaWin
  pop_null$quant<-sample(acuaWin$quant)
  pass<-subset(pop_null,pop_null$quant=="lower10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAcuaDesPermut <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop='ACUA', permutations=null_positives)
print(summary(comboAcuaDesPermut))

comboAcuaIsleRows <- stacWin[which(stacWin$chrPos %in% comboIslePass$cMmidChrPos[which(comboIslePass$nACUA>0)]),]
focal<-comboAcuaIsleRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-acuaWin
  pop_null$quant<-sample(acuaWin$quant)
  pass<-subset(pop_null,pop_null$quant=="upper10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAcuaIslePermut <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop='ACUA', permutations=null_positives)
print(summary(comboAcuaIslePermut))


comboTlmcDesRows <- stacWin[which(stacWin$chrPos %in% comboDesPass$cMmidChrPos[which(comboDesPass$nTLMC>0)]),]
focal<-comboTlmcDesRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-tlmcWin
  pop_null$quant<-sample(tlmcWin$quant)
  pass<-subset(pop_null,pop_null$quant=="lower10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboTlmcDesPermut <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop='TLMC', permutations=null_positives)
print(summary(comboTlmcDesPermut))

comboTlmcIsleRows <- stacWin[which(stacWin$chrPos %in% comboIslePass$cMmidChrPos[which(comboIslePass$nTLMC>0)]),]
focal<-comboTlmcIsleRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-tlmcWin
  pop_null$quant<-sample(tlmcWin$quant)
  pass<-subset(pop_null,pop_null$quant=="upper10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboTlmcIslePermut <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop='TLMC', permutations=null_positives)
print(summary(comboTlmcIslePermut))


comboAgzcAcuaDesRows <- stacWin[which(stacWin$chrPos %in% comboDesPass$cMmidChrPos[which(comboDesPass$nACUA>0 & comboDesPass$nAGZC>0)]),]
focal<-comboAgzcAcuaDesRows
null_positives<-{}
for(x in 1:1000){
  acua_null<-acuaWin
  agzc_null<-agzcWin
  acua_null$quant<-sample(acuaWin$quant)
  agzc_null$quant<-sample(agzcWin$quant)
  pass<-subset(acua_null,acua_null$quant=="lower10"  & agzc_null$quant=="lower10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAgzcAcuaDesPermut <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop='AGZC+ACUA', permutations=null_positives)
print(summary(comboAgzcAcuaDesPermut))

comboAgzcAcuaIsleRows <- stacWin[which(stacWin$chrPos %in% comboIslePass$cMmidChrPos[which(comboIslePass$nACUA>0 & comboIslePass$nAGZC>0)]),]
focal<-comboAgzcAcuaIsleRows
null_positives<-{}
for(x in 1:1000){
  acua_null<-acuaWin
  agzc_null<-agzcWin
  acua_null$quant<-sample(acuaWin$quant)
  agzc_null$quant<-sample(agzcWin$quant)
  pass<-subset(acua_null,acua_null$quant=="upper10"  & agzc_null$quant=="upper10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAgzcAcuaIslePermut <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop='AGZC+ACUA', permutations=null_positives)
print(summary(comboAgzcAcuaIslePermut))

comboAgzcAcuaTlmcDesRows <- stacWin[which(stacWin$chrPos %in% comboDesPass$cMmidChrPos[which(comboDesPass$nACUA>0 & comboDesPass$nAGZC>0 & comboDesPass$nTLMC>0)]),]
focal<-comboAgzcAcuaTlmcDesRows
null_positives<-{}
for(x in 1:1000){
  acua_null<-acuaWin
  agzc_null<-agzcWin
  tlmc_null<-tlmcWin
  acua_null$quant<-sample(acuaWin$quant)
  agzc_null$quant<-sample(agzcWin$quant)
  tlmc_null$quant<-sample(tlmcWin$quant)
  pass<-subset(acua_null,acua_null$quant=="lower10" & agzc_null$quant=="lower10" & tlmc_null$quant=="lower10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAgzcAcuaTlmcDesPermut <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop='AGZC+ACUA+TLMC', permutations=null_positives)
print(summary(comboAgzcAcuaTlmcDesPermut))

comboAgzcAcuaTlmcIsleRows <- stacWin[which(stacWin$chrPos %in% comboIslePass$cMmidChrPos[which(comboIslePass$nACUA>0 & comboIslePass$nAGZC>0& comboIslePass$nTLMC>0)]),]
focal<-comboAgzcAcuaTlmcIsleRows
null_positives<-{}
for(x in 1:1000){
  acua_null<-acuaWin
  agzc_null<-agzcWin
  tlmc_null<-tlmcWin
  acua_null$quant<-sample(acuaWin$quant)
  agzc_null$quant<-sample(agzcWin$quant)
  tlmc_null$quant<-sample(tlmcWin$quant)
  pass<-subset(acua_null,acua_null$quant=="upper10"  & agzc_null$quant=="upper10"  & tlmc_null$quant=="upper10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAgzcAcuaTlmcIslePermut <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop='AGZC+ACUA+TLMC', permutations=null_positives)
print(summary(comboAgzcAcuaTlmcIslePermut))

comboAnyAgzcAcuaTlmcDesRows <- stacWin[which(stacWin$chrPos %in% comboDesPass$cMmidChrPos[which(comboDesPass$nACUA>0 | comboDesPass$nAGZC>0 | comboDesPass$nTLMC>0)]),]
focal<-comboAnyAgzcAcuaTlmcDesRows
null_positives<-{}
for(x in 1:1000){
  acua_null<-acuaWin
  agzc_null<-agzcWin
  tlmc_null<-tlmcWin
  acua_null$quant<-sample(acuaWin$quant)
  agzc_null$quant<-sample(agzcWin$quant)
  tlmc_null$quant<-sample(tlmcWin$quant)
  pass<-subset(acua_null,acua_null$quant=="lower10" | agzc_null$quant=="lower10" | tlmc_null$quant=="lower10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAnyAgzcAcuaTlmcDesPermut <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop='AGZC_ACUA_TLMC', permutations=null_positives)
print(summary(comboAnyAgzcAcuaTlmcDesPermut))

comboAnyAgzcAcuaTlmcIsleRows <- stacWin[which(stacWin$chrPos %in% comboIslePass$cMmidChrPos[which(comboIslePass$nACUA>0 | comboIslePass$nAGZC>0 | comboIslePass$nTLMC>0)]),]
focal<-comboAnyAgzcAcuaTlmcIsleRows
null_positives<-{}
for(x in 1:1000){
  acua_null<-acuaWin
  agzc_null<-agzcWin
  tlmc_null<-tlmcWin
  acua_null$quant<-sample(acuaWin$quant)
  agzc_null$quant<-sample(agzcWin$quant)
  tlmc_null$quant<-sample(tlmcWin$quant)
  pass<-subset(acua_null,acua_null$quant=="upper10"  | agzc_null$quant=="upper10"  | tlmc_null$quant=="upper10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAnyAgzcAcuaTlmcIslePermut <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop='AGZC_ACUA_TLMC', permutations=null_positives)
print(summary(comboAnyAgzcAcuaTlmcIslePermut))

permutations <- rbind(stacDesPermut, stacIslePermut, chplDesPermut, chplIslePermut, comboAgzcDesPermut, comboAgzcIslePermut, comboAcuaDesPermut, comboAcuaIslePermut, comboTlmcDesPermut, comboTlmcIslePermut, comboAgzcAcuaDesPermut, comboAgzcAcuaIslePermut, comboAgzcAcuaTlmcDesPermut, comboAgzcAcuaTlmcIslePermut, comboAnyAgzcAcuaTlmcDesPermut, comboAnyAgzcAcuaTlmcIslePermut)
write.table(permutations, "thinned_desert_islands_shared_permutations.txt", row.name=F, sep="\t", quote=F)

###Jackknife
chunkSize <- 200

winCM$chrPos <- paste(winCM$chrID, winCM$start, sep="_")
uniChrs <- unique(winCM$chrID)

stacDesCross <- subset(stacDes, chplN>0)
stacDesCrossCount <- nrow(stacDesCross)
stacDesJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- stacDesCross[which(!(stacDesCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="STAC_Desert", focalPop="STAC", regType="Desert", compPop="CHPL", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    popCount <- nrow(subset(sharedOutside, chplN>0))
    tempDF <- cbind(tempDF, popCount=popCount)
    stacDesJackknifeDF <- rbind(stacDesJackknifeDF, tempDF)
    start = start+chunkSize
  }
}
summary(stacDesJackknifeDF)

stacIsleCross <- subset(stacIsle, chplN>0 & cMmidPass==T)
stacIsleCrossCount <- nrow(stacIsleCross)
stacIsleJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- stacIsleCross[which(!(stacIsleCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="STAC_Island", focalPop="STAC", regType="Island", compPop="CHPL", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    popCount <- nrow(subset(sharedOutside, chplN>0))
    tempDF <- cbind(tempDF, popCount=popCount)
    stacIsleJackknifeDF <- rbind(stacIsleJackknifeDF, tempDF)
    start = start+chunkSize
  }
}

chplDesCross <- subset(chplDes, stacN>0)
chplDesCrossCount <- nrow(chplDesCross)
chplDesJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- chplDesCross[which(!(chplDesCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="CHPL_Desert", focalPop="CHPL", regType="Desert", compPop="STAC", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    popCount <- nrow(subset(sharedOutside, stacN>0))
    tempDF <- cbind(tempDF, popCount=popCount)
    chplDesJackknifeDF <- rbind(chplDesJackknifeDF, tempDF)
    start = start+chunkSize
  }
}

chplIsleCross <- subset(chplIsle, stacN>0 & cMmidPass==T)
chplIsleCrossCount <- nrow(chplIsleCross)
chplIsleJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- chplIsleCross[which(!(chplIsleCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="CHPL_Island", focalPop="CHPL", regType="Island", compPop="STAC", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    popCount <- nrow(subset(sharedOutside, stacN>0))
    tempDF <- cbind(tempDF, popCount=popCount)
    chplIsleJackknifeDF <- rbind(chplIsleJackknifeDF, tempDF)
    start = start+chunkSize
  }
}

comboAgzcDesCross <- subset(comboDesPass, nAGZC>0)
comboAgzcDesCrossCount <- nrow(comboAgzcDesCross)
comboAgzcDesJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- comboAgzcDesCross[which(!(comboAgzcDesCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop="AGZC", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    agzcCount <- nrow(subset(sharedOutside, nAGZC>0))
    tempDF <- cbind(tempDF, popCount=agzcCount)
    comboAgzcDesJackknifeDF <- rbind(comboAgzcDesJackknifeDF, tempDF)
    start = start+chunkSize
  }
}

comboAgzcIsleCross <- subset(comboIslePass, nAGZC>0)
comboAgzcIsleCrossCount <- nrow(comboAgzcIsleCross)
comboAgzcIsleJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- comboAgzcIsleCross[which(!(comboAgzcIsleCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop="AGZC", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    agzcCount <- nrow(subset(sharedOutside, nAGZC>0))
    tempDF <- cbind(tempDF, popCount=agzcCount)
    comboAgzcIsleJackknifeDF <- rbind(comboAgzcIsleJackknifeDF, tempDF)
    start = start+chunkSize
  }
}

comboAcuaDesCross <- subset(comboDesPass, nACUA>0)
comboAcuaDesCrossCount <- nrow(comboAcuaDesCross)
comboAcuaDesJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- comboAcuaDesCross[which(!(comboAcuaDesCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop="ACUA", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    acuaCount <- nrow(subset(sharedOutside, nACUA>0))
    tempDF <- cbind(tempDF, popCount=acuaCount)
    comboAcuaDesJackknifeDF <- rbind(comboAcuaDesJackknifeDF, tempDF)
    start = start+chunkSize
  }
}
summary(comboAcuaDesJackknifeDF)

comboAcuaIsleCross <- subset(comboIslePass, nACUA>0)
comboAcuaIsleCrossCount <- nrow(comboAcuaIsleCross)
comboAcuaIsleJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- comboAcuaIsleCross[which(!(comboAcuaIsleCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop="ACUA", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    acuaCount <- nrow(subset(sharedOutside, nACUA>0))
    tempDF <- cbind(tempDF, popCount=acuaCount)
    comboAcuaIsleJackknifeDF <- rbind(comboAcuaIsleJackknifeDF, tempDF)
    start = start+chunkSize
  }
}
summary(comboAcuaIsleJackknifeDF)

comboTlmcDesCross <- subset(comboDesPass, nTLMC>0)
comboTlmcDesCrossCount <- nrow(comboTlmcDesCross)
comboTlmcDesJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- comboTlmcDesCross[which(!(comboTlmcDesCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop="TLMC", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    tlmcCount <- nrow(subset(sharedOutside, nTLMC>0))
    tempDF <- cbind(tempDF, popCount=tlmcCount)
    comboTlmcDesJackknifeDF <- rbind(comboTlmcDesJackknifeDF, tempDF)
    start = start+chunkSize
  }
}
summary(comboTlmcDesJackknifeDF)

comboTlmcIsleCross <- subset(comboIslePass, nTLMC>0)
comboTlmcIsleCrossCount <- nrow(comboTlmcIsleCross)
comboTlmcIsleJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- comboTlmcIsleCross[which(!(comboTlmcIsleCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop="TLMC", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    tlmcCount <- nrow(subset(sharedOutside, nTLMC>0))
    tempDF <- cbind(tempDF, popCount=tlmcCount)
    comboTlmcIsleJackknifeDF <- rbind(comboTlmcIsleJackknifeDF, tempDF)
    start = start+chunkSize
  }
}
summary(comboTlmcIsleJackknifeDF)


comboAnyDesCross <- subset(comboDesPass, nAGZC>0 | nACUA>0 | nTLMC>0)
comboAnyDesCrossCount <- nrow(comboAnyDesCross)
comboAnyDesJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- comboAnyDesCross[which(!(comboAnyDesCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop="Any", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    anyCount <- nrow(subset(sharedOutside, nAGZC>0 | nACUA>0 | nTLMC>0))
    tempDF <- cbind(tempDF, popCount=anyCount)
    comboAnyDesJackknifeDF <- rbind(comboAnyDesJackknifeDF, tempDF)
    start = start+chunkSize
  }
}
summary(comboAnyDesJackknifeDF)

comboAnyIsleCross <- subset(comboIslePass, nAGZC>0 | nACUA>0 | nTLMC>0)
comboAnyIsleCrossCount <- nrow(comboAnyIsleCross)
comboAnyIsleJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- comboAnyIsleCross[which(!(comboAnyIsleCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop="Any", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    anyCount <- nrow(subset(sharedOutside, nAGZC>0 | nACUA>0 | nTLMC>0))
    tempDF <- cbind(tempDF, popCount=anyCount)
    comboAnyIsleJackknifeDF <- rbind(comboAnyIsleJackknifeDF, tempDF)
    start = start+chunkSize
  }
}
summary(comboAnyIsleJackknifeDF)

allJackknifeDF <- rbind(stacDesJackknifeDF, stacIsleJackknifeDF, chplDesJackknifeDF, chplIsleJackknifeDF, comboAgzcDesJackknifeDF, comboAgzcIsleJackknifeDF, comboAcuaDesJackknifeDF, comboAcuaIsleJackknifeDF, comboTlmcDesJackknifeDF, comboTlmcIsleJackknifeDF, comboAnyDesJackknifeDF, comboAnyIsleJackknifeDF)
write.table(allJackknifeDF, "thinned_desert_island_shared_jackknifed.txt", row.names = F, sep="\t", quote=F)

###Unthinned
winCM <- read.table("../windowedData/xbir_pacbio2018_allChrs_0.05cM_windows_recRate_codingBP_xcorxbir_xmalxbir_minPar.txt", header=T)

chplWin <- data.frame(chrID=winCM$chrID, start=winCM$start, end=winCM$end, chrPos=paste(winCM$chrID, winCM$start, sep="_"), minAnc=winCM$CHPL2021)
chplWin$regType <- NA

chplDes <- read.table("CHPL2021_minorParentDeserts_AIMs_cMwins_merged.txt", header=T)
chplDes$regType <- "Desert"
chplDes$cMmidChrPos <- NA
chplDes$cMmidPass <- NA
chplDes$nWin <- NA
chplDes$propWin <- NA
for (i in 1:nrow(chplDes)) {
  regWins <- subset(chplWin, chrID==chplDes$chrName[i] & start>=chplDes$cMstart[i] & end<=chplDes$cMend[i])
  chplWin$regType[which(chplWin$chrID==chplDes$chrName[i] & chplWin$start>=chplDes$cMstart[i] & chplWin$end<=chplDes$cMend[i])] <- "Desert"
  midWin <- subset(regWins, start<=chplDes$cMmid[i] & end>=chplDes$cMmid[i])
  cutoff <- subset(quants, popName=="CHPL2021")$lower10
  chplDes$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  }
  else if (midWin$minAnc<=cutoff) {
    pass <- TRUE
  }
  else if (midWin$minAnc>cutoff) {
    pass <- FALSE
  }
  chplDes$cMmidPass[i] <- pass
  chplDes$nWin[i] <- nrow(regWins)
  chplDes$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc<=cutoff))
}
dim(chplDes)
nrow(subset(chplDes, is.na(cMmidPass)))
sum(chplDes$cMmidPass, na.rm=T)

chplIsle <- read.table("CHPL2021_minorParentIslands_AIMs_cMwins_merged.txt", header=T)
chplIsle$regType <- "Island"
chplIsle$cMmidChrPos <- NA
chplIsle$cMmidPass <- NA
chplIsle$nWin <- NA
chplIsle$propWin <- NA
for (i in 1:nrow(chplIsle)) {
  regWins <- subset(chplWin, chrID==chplIsle$chrName[i] & start>=chplIsle$cMstart[i] & end<=chplIsle$cMend[i])
  chplWin$regType[which(chplWin$chrID==chplIsle$chrName[i] & chplWin$start>=chplIsle$cMstart[i] & chplWin$end<=chplIsle$cMend[i])] <- "Island"
  midWin <- subset(regWins,  start<=chplIsle$cMmid[i] & end>=chplIsle$cMmid[i])
  cutoff <- subset(quants, popName=="CHPL2021")$upper10
  chplIsle$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  }
  else if (midWin$minAnc>=cutoff) {
    pass <- TRUE
  }
  else if (midWin$minAnc<cutoff) {
    pass <- FALSE
  }
  chplIsle$cMmidPass[i] <- pass
  chplIsle$nWin[i] <- nrow(regWins)
  chplIsle$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc>=cutoff))
}
dim(chplIsle)
nrow(subset(chplIsle, is.na(cMmidPass)))
sum(chplIsle$cMmidPass, na.rm=T)


stacWin <- data.frame(chrID=winCM$chrID, start=winCM$start, end=winCM$end, chrPos=paste(winCM$chrID, winCM$start, sep="_"), minAnc=winCM$STAC2020, regType=NA)

stacDes <- read.table("STAC2020_minorParentDeserts_AIMs_cMwins_merged.txt", header=T)
stacDes$regType <- "Desert"
stacDes$cMmidChrPos <- NA
stacDes$cMmidPass <- NA
stacDes$nWin <- NA
stacDes$propWin <- NA
for (i in 1:nrow(stacDes)) {
  regWins <- subset(stacWin, chrID==stacDes$chrName[i] & start>=stacDes$cMstart[i] & end<=stacDes$cMend[i])
  stacWin$regType[which(stacWin$chrID==stacDes$chrName[i] & stacWin$start>=stacDes$cMstart[i] & stacWin$end<=stacDes$cMend[i])] <- "Desert"
  midWin <- subset(regWins,  start<=stacDes$cMmid[i] & end>=stacDes$cMmid[i])
  cutoff <- subset(quants, popName=="STAC2020")$lower10
  stacDes$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  }
  else if (midWin$minAnc<=cutoff) {
    pass <- TRUE
  }
  else if (midWin$minAnc>cutoff) {
    pass <- FALSE
  }
  stacDes$cMmidPass[i] <- pass
  stacDes$nWin[i] <- nrow(regWins)
  stacDes$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc<=cutoff))
}
dim(stacDes)
nrow(subset(stacDes, is.na(cMmidPass)))
sum(stacDes$cMmidPass, na.rm=T)

stacIsle <- read.table("STAC2020_minorParentIslands_AIMs_cMwins_merged.txt", header=T)
stacIsle$regType <- "Island"
stacIsle$cMmidChrPos <- NA
stacIsle$cMmidPass <- NA
stacIsle$nWin <- NA
stacIsle$propWin <- NA
for (i in 1:nrow(stacIsle)) {
  regWins <- subset(stacWin, chrID==stacIsle$chrName[i] & start>=stacIsle$cMstart[i] & end<=stacIsle$cMend[i])
  stacWin$regType[which(stacWin$chrID==stacIsle$chrName[i] & stacWin$start>=stacIsle$cMstart[i] & stacWin$end<=stacIsle$cMend[i])] <- "Island"
  midWin <- subset(regWins,  start<=stacIsle$cMmid[i] & end>=stacIsle$cMmid[i])
  cutoff <- subset(quants, popName=="STAC2020")$upper10
  stacIsle$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  }
  else if (midWin$minAnc>=cutoff) {
    pass <- TRUE
  }
  else if (midWin$minAnc<cutoff) {
    pass <- FALSE
  }
  stacIsle$cMmidPass[i] <- pass
  stacIsle$nWin[i] <- nrow(regWins)
  stacIsle$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc>=cutoff))
}
dim(stacIsle)
nrow(subset(stacIsle, is.na(cMmidPass)))
sum(stacIsle$cMmidPass, na.rm=T)


agzcWin <- data.frame(chrID=winCM$chrID, start=winCM$start, end=winCM$end, chrPos=paste(winCM$chrID, winCM$start, sep="_"), minAnc=winCM$AGZC, regType=NA)

agzcDes <- read.table("AGZC_minorParentDeserts_AIMs_cMwins_merged.txt", header=T)
agzcDes$regType <- "Desert"
agzcDes$cMmidChrPos <- NA
agzcDes$cMmidPass <- NA
agzcDes$nWin <- NA
agzcDes$propWin <- NA
for (i in 1:nrow(agzcDes)) {
  regWins <- subset(agzcWin, chrID==agzcDes$chrName[i] & start>=agzcDes$cMstart[i] & end<=agzcDes$cMend[i])
  agzcWin$regType[which(agzcWin$chrID==agzcDes$chrName[i] & agzcWin$start>=agzcDes$cMstart[i] & agzcWin$end<=agzcDes$cMend[i])] <- "Desert"
  midWin <- subset(regWins,  start<=agzcDes$cMmid[i] & end>=agzcDes$cMmid[i])
  cutoff <- subset(quants, popName=="AGZC")$lower10
  agzcDes$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  }
  else if (midWin$minAnc<=cutoff) {
    pass <- TRUE
  }
  else if (midWin$minAnc>cutoff) {
    pass <- FALSE
  }
  agzcDes$cMmidPass[i] <- pass
  agzcDes$nWin[i] <- nrow(regWins)
  agzcDes$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc<=cutoff))
}
dim(agzcDes)
nrow(subset(agzcDes, is.na(cMmidPass)))
sum(agzcDes$cMmidPass, na.rm=T)

agzcIsle <- read.table("AGZC_minorParentIslands_AIMs_cMwins_merged.txt", header=T)
agzcIsle$regType <- "Island"
agzcIsle$cMmidChrPos <- NA
agzcIsle$cMmidPass <- NA
agzcIsle$nWin <- NA
agzcIsle$propWin <- NA
for (i in 1:nrow(agzcIsle)) {
  regWins <- subset(agzcWin, chrID==agzcIsle$chrName[i] & start>=agzcIsle$cMstart[i] & end<=agzcIsle$cMend[i])
  agzcWin$regType[which(agzcWin$chrID==agzcIsle$chrName[i] & agzcWin$start>=agzcIsle$cMstart[i] & agzcWin$end<=agzcIsle$cMend[i])] <- "Island"
  midWin <- subset(regWins,  start<=agzcIsle$cMmid[i] & end>=agzcIsle$cMmid[i])
  cutoff <- subset(quants, popName=="AGZC")$upper10
  agzcIsle$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  }
  else if (midWin$minAnc>=cutoff) {
    pass <- TRUE
  }
  else if (midWin$minAnc<cutoff) {
    pass <- FALSE
  }
  agzcIsle$cMmidPass[i] <- pass
  agzcIsle$nWin[i] <- nrow(regWins)
  agzcIsle$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc>=cutoff))
}
dim(agzcIsle)
nrow(subset(agzcIsle, is.na(cMmidPass)))
sum(agzcIsle$cMmidPass, na.rm=T)


acuaWin <- data.frame(chrID=winCM$chrID, start=winCM$start, end=winCM$end, chrPos=paste(winCM$chrID, winCM$start, sep="_"), minAnc=winCM$ACUA2018, regType=NA)

acuaDes <- read.table("ACUA2018_minorParentDeserts_AIMs_cMwins_merged.txt", header=T)
acuaDes$regType <- "Desert"
acuaDes$cMmidChrPos <- NA
acuaDes$cMmidPass <- NA
acuaDes$nWin <- NA
acuaDes$propWin <- NA
for (i in 1:nrow(acuaDes)) {
  regWins <- subset(acuaWin, chrID==acuaDes$chrName[i] & start>=acuaDes$cMstart[i] & end<=acuaDes$cMend[i])
  acuaWin$regType[which(acuaWin$chrID==acuaDes$chrName[i] & acuaWin$start>=acuaDes$cMstart[i] & acuaWin$end<=acuaDes$cMend[i])] <- "Desert"
  midWin <- subset(regWins,  start<=acuaDes$cMmid[i] & end>=acuaDes$cMmid[i])
  cutoff <- subset(quants, popName=="ACUA2018")$lower10
  acuaDes$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  }
  else if (midWin$minAnc<=cutoff) {
    pass <- TRUE
  }
  else if (midWin$minAnc>cutoff) {
    pass <- FALSE
  }
  acuaDes$cMmidPass[i] <- pass
  acuaDes$nWin[i] <- nrow(regWins)
  acuaDes$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc<=cutoff))
}
dim(acuaDes)
nrow(subset(acuaDes, is.na(cMmidPass)))
sum(acuaDes$cMmidPass, na.rm=T)

acuaIsle <- read.table("ACUA2018_minorParentIslands_AIMs_cMwins_merged.txt", header=T)
acuaIsle$regType <- "Island"
acuaIsle$cMmidChrPos <- NA
acuaIsle$cMmidPass <- NA
acuaIsle$nWin <- NA
acuaIsle$propWin <- NA
for (i in 1:nrow(acuaIsle)) {
  regWins <- subset(acuaWin, chrID==acuaIsle$chrName[i] & start>=acuaIsle$cMstart[i] & end<=acuaIsle$cMend[i])
  acuaWin$regType[which(acuaWin$chrID==acuaIsle$chrName[i] & acuaWin$start>=acuaIsle$cMstart[i] & acuaWin$end<=acuaIsle$cMend[i])] <- "Island"
  midWin <- subset(regWins,  start<=acuaIsle$cMmid[i] & end>=acuaIsle$cMmid[i])
  cutoff <- subset(quants, popName=="ACUA2018")$upper10
  acuaIsle$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  }
  else if (midWin$minAnc>=cutoff) {
    pass <- TRUE
  }
  else if (midWin$minAnc<cutoff) {
    pass <- FALSE
  }
  acuaIsle$cMmidPass[i] <- pass
  acuaIsle$nWin[i] <- nrow(regWins)
  acuaIsle$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc>=cutoff))
}
dim(acuaIsle)
nrow(subset(acuaIsle, is.na(cMmidPass)))
sum(acuaIsle$cMmidPass, na.rm=T)


tlmcWin <- data.frame(chrID=winCM$chrID, start=winCM$start, end=winCM$end, chrPos=paste(winCM$chrID, winCM$start, sep="_"), minAnc=winCM$TLMC2017, regType=NA)

tlmcDes <- read.table("TLMC2017_minorParentDeserts_AIMs_cMwins_merged.txt", header=T)
tlmcDes$regType <- "Desert"
tlmcDes$cMmidChrPos <- NA
tlmcDes$cMmidPass <- NA
tlmcDes$nWin <- NA
tlmcDes$propWin <- NA
for (i in 1:nrow(tlmcDes)) {
  regWins <- subset(tlmcWin, chrID==tlmcDes$chrName[i] & start>=tlmcDes$cMstart[i] & end<=tlmcDes$cMend[i])
  tlmcWin$regType[which(tlmcWin$chrID==tlmcDes$chrName[i] & tlmcWin$start>=tlmcDes$cMstart[i] & tlmcWin$end<=tlmcDes$cMend[i])] <- "Desert"
  midWin <- subset(regWins,  start<=tlmcDes$cMmid[i] & end>=tlmcDes$cMmid[i])
  cutoff <- subset(quants, popName=="TLMC2017")$lower10
  tlmcDes$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  }
  else if (midWin$minAnc<=cutoff) {
    pass <- TRUE
  }
  else if (midWin$minAnc>cutoff) {
    pass <- FALSE
  }
  tlmcDes$cMmidPass[i] <- pass
  tlmcDes$nWin[i] <- nrow(regWins)
  tlmcDes$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc<=cutoff))
}
dim(tlmcDes)
nrow(subset(tlmcDes, is.na(cMmidPass)))
sum(tlmcDes$cMmidPass, na.rm=T)

tlmcIsle <- read.table("TLMC2017_minorParentIslands_AIMs_cMwins_merged.txt", header=T)
tlmcIsle$regType <- "Island"
tlmcIsle$cMmidChrPos <- NA
tlmcIsle$cMmidPass <- NA
tlmcIsle$nWin <- NA
tlmcIsle$propWin <- NA
for (i in 1:nrow(tlmcIsle)) {
  regWins <- subset(tlmcWin, chrID==tlmcIsle$chrName[i] & start>=tlmcIsle$cMstart[i] & end<=tlmcIsle$cMend[i])
  tlmcWin$regType[which(tlmcWin$chrID==tlmcIsle$chrName[i] & tlmcWin$start>=tlmcIsle$cMstart[i] & tlmcWin$end<=tlmcIsle$cMend[i])] <- "Island"
  midWin <- subset(regWins,  start<=tlmcIsle$cMmid[i] & end>=tlmcIsle$cMmid[i])
  cutoff <- subset(quants, popName=="TLMC2017")$upper10
  tlmcIsle$cMmidChrPos[i] <- paste(midWin$chrID, midWin$start, sep="_")
  pass <- NA
  if (is.na(midWin$minAnc)) {
    pass <- NA
  }
  else if (midWin$minAnc>=cutoff) {
    pass <- TRUE
  }
  else if (midWin$minAnc<cutoff) {
    pass <- FALSE
  }
  tlmcIsle$cMmidPass[i] <- pass
  tlmcIsle$nWin[i] <- nrow(regWins)
  tlmcIsle$propWin[i] <- nrow(subset(regWins, !is.na(minAnc) & minAnc>=cutoff))
}
dim(tlmcIsle)
nrow(subset(tlmcIsle, is.na(cMmidPass)))
sum(tlmcIsle$cMmidPass, na.rm=T)

agzcDes$propPass <- agzcDes$propWin/agzcDes$nWin
agzcIsle$propPass <- agzcIsle$propWin/agzcIsle$nWin
stacDes$propPass <- stacDes$propWin/stacDes$nWin
stacIsle$propPass <- stacIsle$propWin/stacIsle$nWin
agzcDes$propPass <- agzcDes$propWin/agzcDes$nWin
agzcIsle$propPass <- agzcIsle$propWin/agzcIsle$nWin
acuaDes$propPass <- acuaDes$propWin/acuaDes$nWin
acuaIsle$propPass <- acuaIsle$propWin/acuaIsle$nWin
tlmcDes$propPass <- tlmcDes$propWin/tlmcDes$nWin
tlmcIsle$propPass <- tlmcIsle$propWin/tlmcIsle$nWin

summary(chplDes)
summary(chplIsle)
summary(stacDes)
summary(stacIsle)
summary(agzcDes)
summary(agzcIsle)
summary(acuaDes)
summary(acuaIsle)
summary(tlmcDes)
summary(tlmcIsle)

agzcDesPass <- subset(agzcDes, cMmidPass==T)
agzcIslePass <- subset(agzcIsle, cMmidPass==T)
stacDesPass <- subset(stacDes, cMmidPass==T)
stacIslePass <- subset(stacIsle, cMmidPass==T)

sum(agzcDes$cMmidChrPos %in% stacDes$cMmidChrPos)
sum(agzcIsle$cMmidChrPos %in% stacIsle$cMmidChrPos)

sum(stacDes$cMmidChrPos %in% agzcDes$cMmidChrPos)
sum(stacDes$cMmidChrPos %in% agzcDes$cMmidChrPos)
sum(stacDes$cMmidChrPos %in% acuaDes$cMmidChrPos)
sum(stacDes$cMmidChrPos %in% tlmcDes$cMmidChrPos)
sum(stacIsle$cMmidChrPos %in% agzcIsle$cMmidChrPos)
sum(stacIsle$cMmidChrPos %in% agzcIsle$cMmidChrPos)
sum(stacIsle$cMmidChrPos %in% acuaIsle$cMmidChrPos)
sum(stacIsle$cMmidChrPos %in% tlmcIsle$cMmidChrPos)

sum(agzcDes$cMmidChrPos %in% agzcDes$cMmidChrPos)
sum(agzcDes$cMmidChrPos %in% acuaDes$cMmidChrPos)
sum(agzcDes$cMmidChrPos %in% tlmcDes$cMmidChrPos)
sum(agzcIsle$cMmidChrPos %in% agzcIsle$cMmidChrPos)
sum(agzcIsle$cMmidChrPos %in% acuaIsle$cMmidChrPos)
sum(agzcIsle$cMmidChrPos %in% tlmcIsle$cMmidChrPos)

chplDes$stacN <- NA
for (i in 1:nrow(chplDes)) {
  stacRegWins <- subset(stacWin, chrID==chplDes$chrName[i] & start>=chplDes$cMstart[i] & end<=chplDes$cMend[i])
  chplDes$stacN[i] <- nrow(subset(stacRegWins, regType=="Desert"))
}
nrow(chplDes)
nrow(subset(chplDes, stacN>0))
sum(subset(chplDes, stacN>0)$cMmidPass, na.rm=T)

chplIsle$stacN <- NA
for (i in 1:nrow(chplIsle)) {
  stacRegWins <- subset(stacWin, chrID==chplIsle$chrName[i] & start>=chplIsle$cMstart[i] & end<=chplIsle$cMend[i])
  chplIsle$stacN[i] <- nrow(subset(stacRegWins, regType=="Island"))
}
nrow(chplIsle)
nrow(subset(chplIsle, stacN>0))
sum(subset(chplIsle, stacN>0)$cMmidPass, na.rm=T)

stacDes$chplN <- NA
for (i in 1:nrow(stacDes)) {
  chplRegWins <- subset(chplWin, chrID==stacDes$chrName[i] & start>=stacDes$cMstart[i] & end<=stacDes$cMend[i])
  stacDes$chplN[i] <- nrow(subset(chplRegWins, regType=="Desert"))
}
nrow(stacDes)
nrow(subset(stacDes, chplN>0))
sum(subset(stacDes, chplN>0)$cMmidPass, na.rm=T)

stacIsle$chplN <- NA
for (i in 1:nrow(stacIsle)) {
  chplRegWins <- subset(chplWin, chrID==stacIsle$chrName[i] & start>=stacIsle$cMstart[i] & end<=stacIsle$cMend[i])
  stacIsle$chplN[i] <- nrow(subset(chplRegWins, regType=="Island"))
}
nrow(stacIsle)
nrow(subset(stacIsle, chplN>0))
head(subset(stacIsle, chplN>0))
sum(subset(stacIsle, chplN>0)$cMmidPass, na.rm=T)

chplWin$quant <- NA
chplWin$quant[which(chplWin$minAnc<=quants$lower10[which(quants$pop=="CHPL2021")])] <- "lower10"
chplWin$quant[which(chplWin$minAnc>=quants$upper10[which(quants$pop=="CHPL2021")])] <- "upper10"

stacWin$quant <- NA
stacWin$quant[which(stacWin$minAnc<=quants$lower10[which(quants$pop=="STAC2020")])] <- "lower10"
stacWin$quant[which(stacWin$minAnc>=quants$upper10[which(quants$pop=="STAC2020")])] <- "upper10"

agzcWin$quant <- NA
agzcWin$quant[which(agzcWin$minAnc<=quants$lower10[which(quants$pop=="AGZC")])] <- "lower10"
agzcWin$quant[which(agzcWin$minAnc>=quants$upper10[which(quants$pop=="AGZC")])] <- "upper10"

acuaWin$quant <- NA
acuaWin$quant[which(acuaWin$minAnc<=quants$lower10[which(quants$pop=="ACUA2018")])] <- "lower10"
acuaWin$quant[which(acuaWin$minAnc>=quants$upper10[which(quants$pop=="ACUA2018")])] <- "upper10"

tlmcWin$quant <- NA
tlmcWin$quant[which(tlmcWin$minAnc<=quants$lower10[which(quants$pop=="TLMC2017")])] <- "lower10"
tlmcWin$quant[which(tlmcWin$minAnc>=quants$upper10[which(quants$pop=="TLMC2017")])] <- "upper10"

stacDesRows <- stacWin[which(stacWin$chrPos %in% stacDes$cMmidChrPos[which(stacDes$cMmidPass=="TRUE" & stacDes$chplN>0)]),]
focal<-stacDesRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-chplWin
  pop_null$quant<-sample(chplWin$quant)
  pass<-subset(pop_null,pop_null$quant=="lower10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
stacDesPermut <- data.frame(popComp="STAC_Desert", focalPop="STAC", regType="Desert", compPop='CHPL', permutations=null_positives)
print(summary(stacDesPermut))

stacIsleRows <- stacWin[which(stacWin$chrPos %in% stacIsle$cMmidChrPos[which(stacIsle$cMmidPass=="TRUE" & stacIsle$chplN>0)]),]
focal<-stacIsleRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-chplWin
  pop_null$quant<-sample(chplWin$quant)
  pass<-subset(pop_null,pop_null$quant=="upper10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
stacIslePermut <- data.frame(popComp="STAC_Island", focalPop="STAC", regType="Island", compPop='CHPL', permutations=null_positives)
print(summary(stacIslePermut))


chplDesRows <- chplWin[which(chplWin$chrPos %in% chplDes$cMmidChrPos[which(chplDes$cMmidPass=="TRUE" & chplDes$stacN>0)]),]
focal<-chplDesRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-stacWin
  pop_null$quant<-sample(stacWin$quant)
  pass<-subset(pop_null,pop_null$quant=="lower10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
chplDesPermut <- data.frame(popComp="CHPL_Desert", focalPop="CHPL", regType="Desert", compPop='STAC', permutations=null_positives)
print(summary(chplDesPermut))

chplIsleRows <- chplWin[which(chplWin$chrPos %in% chplIsle$cMmidChrPos[which(chplIsle$cMmidPass=="TRUE" & chplIsle$stacN>0)]),]
focal<-chplIsleRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-stacWin
  pop_null$quant<-sample(stacWin$quant)
  pass<-subset(pop_null,pop_null$quant=="upper10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
chplIslePermut <- data.frame(popComp="CHPL_Island", focalPop="CHPL", regType="Island", compPop='STAC', permutations=null_positives)
print(summary(chplIslePermut))

comboShare <- read.table("combinedShared_CHPL_STAC_desertsIslands.txt", header=T)
comboShare$cMmidChrPos <- NA
comboShare$stacMidPass <- NA
comboShare$chplMidPass <- NA
for (i in 1:nrow(comboShare)) {
  stacMidWin <- subset(stacWin, chrID==comboShare$chrName[i] & start<=comboShare$cMmid[i] & end>=comboShare$cMmid[i])
  chplMidWin <- subset(chplWin, chrID==comboShare$chrName[i] & start<=comboShare$cMmid[i] & end>=comboShare$cMmid[i])
  comboShare$cMmidChrPos[i] <- paste(stacMidWin$chrID, stacMidWin$start, sep="_")
  stacPass <- NA
  if (is.na(stacMidWin$minAnc)) {
    stacPass <- NA
  } else if (comboShare$regType[i]=="Desert" & stacMidWin$quant=="lower10" & !is.na(stacMidWin$quant)) {
    stacPass <- TRUE
  } else if (comboShare$regType[i]=="Island" & stacMidWin$quant=="upper10" & !is.na(stacMidWin$quant)) {
    stacPass <- TRUE
  } else if (comboShare$regType[i]=="Desert" & (stacMidWin$quant!="lower10" | is.na(stacMidWin$quant))) {
    stacPass <- FALSE
  } else if (comboShare$regType[i]=="Island" & (stacMidWin$quant!="upper10" | is.na(stacMidWin$quant))) {
    stacPass <- FALSE
  }
  chplPass <- NA
  if (is.na(chplMidWin$minAnc)) {
    chplPass <- NA
  } else if (comboShare$regType[i]=="Desert" & chplMidWin$quant=="lower10" & !is.na(chplMidWin$quant)) {
    chplPass <- TRUE
  } else if (comboShare$regType[i]=="Island" & chplMidWin$quant=="upper10" & !is.na(chplMidWin$quant)) {
    chplPass <- TRUE
  } else if (comboShare$regType[i]=="Desert" & (chplMidWin$quant!="lower10" | is.na(chplMidWin$quant))) {
    chplPass <- FALSE
  } else if (comboShare$regType[i]=="Island" & (chplMidWin$quant!="upper10" | is.na(chplMidWin$quant))) {
    chplPass <- FALSE
  }
  comboShare$stacMidPass[i] <- stacPass
  comboShare$chplMidPass[i] <- chplPass
}

comboDes <- subset(comboShare, regType=="Desert")
comboIsle <- subset(comboShare, regType=="Island")
write.table(comboShare, "combinedShared_CHPL_STAC_desertsIslands_midCMwin.txt", row.name=F, quote=F, sep="\t")

comboDesPass <- subset(comboDes, stacMidPass==T & chplMidPass==T)
comboIslePass <- subset(comboIsle, stacMidPass==T & chplMidPass==T)

comboShare$nAGZC <- NA
for (i in 1:nrow(comboShare)) {
  agzcRegWins <- subset(agzcWin, chrID==comboShare$chrName[i] & start>=comboShare$cMstart[i] & end<=comboShare$cMend[i])
  comboShare$nAGZC[i] <- nrow(subset(agzcRegWins, regType==comboShare$regType[i]))
}

comboShare$nTLMC <- NA
for (i in 1:nrow(comboShare)) {
  tlmcRegWins <- subset(tlmcWin, chrID==comboShare$chrName[i] & start>=comboShare$cMstart[i] & end<=comboShare$cMend[i])
  comboShare$nTLMC[i] <- nrow(subset(tlmcRegWins, regType==comboShare$regType[i]))
}

comboShare$nACUA <- NA
for (i in 1:nrow(comboShare)) {
  acuaRegWins <- subset(acuaWin, chrID==comboShare$chrName[i] & start>=comboShare$cMstart[i] & end<=comboShare$cMend[i])
  comboShare$nACUA[i] <- nrow(subset(acuaRegWins, regType==comboShare$regType[i]))
}

comboDes <- subset(comboShare, regType=="Desert")
comboIsle <- subset(comboShare, regType=="Island")
write.table(comboShare, "combinedShared_CHPL_STAC_desertsIslands_midCMwin_crossInfo.txt", row.name=F, quote=F, sep="\t")

comboDesPass <- subset(comboDes, stacMidPass==T & chplMidPass==T)
comboIslePass <- subset(comboIsle, stacMidPass==T & chplMidPass==T)

nrow(subset(comboDesPass, nAGZC>0))
nrow(subset(comboDesPass, nACUA>0))
nrow(subset(comboDesPass, nTLMC>0))
nrow(subset(comboDesPass, nAGZC>0 & nACUA>0))
nrow(subset(comboDesPass, nAGZC>0 & nACUA>0 & nTLMC>0))

nrow(subset(comboIslePass, nAGZC>0))
nrow(subset(comboIslePass, nACUA>0))
nrow(subset(comboIslePass, nTLMC>0))
nrow(subset(comboIslePass, nAGZC>0 & nACUA>0))
nrow(subset(comboIslePass, nAGZC>0 & nACUA>0 & nTLMC>0))

comboAgzcDesRows <- stacWin[which(stacWin$chrPos %in% comboDesPass$cMmidChrPos[which(comboDesPass$nAGZC>0)]),]
focal<-comboAgzcDesRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-agzcWin
  pop_null$quant<-sample(agzcWin$quant)
  pass<-subset(pop_null,pop_null$quant=="lower10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAgzcDesPermut <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop='AGZC', permutations=null_positives)
print(summary(comboAgzcDesPermut))

comboAgzcIsleRows <- stacWin[which(stacWin$chrPos %in% comboIslePass$cMmidChrPos[which(comboIslePass$nAGZC>0)]),]
focal<-comboAgzcIsleRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-agzcWin
  pop_null$quant<-sample(agzcWin$quant)
  pass<-subset(pop_null,pop_null$quant=="upper10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAgzcIslePermut <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop='AGZC', permutations=null_positives)
print(summary(comboAgzcIslePermut))

comboAcuaDesRows <- stacWin[which(stacWin$chrPos %in% comboDesPass$cMmidChrPos[which(comboDesPass$nACUA>0)]),]
focal<-comboAcuaDesRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-acuaWin
  pop_null$quant<-sample(acuaWin$quant)
  pass<-subset(pop_null,pop_null$quant=="lower10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAcuaDesPermut <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop='ACUA', permutations=null_positives)
print(summary(comboAcuaDesPermut))

comboAcuaIsleRows <- stacWin[which(stacWin$chrPos %in% comboIslePass$cMmidChrPos[which(comboIslePass$nACUA>0)]),]
focal<-comboAcuaIsleRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-acuaWin
  pop_null$quant<-sample(acuaWin$quant)
  pass<-subset(pop_null,pop_null$quant=="upper10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAcuaIslePermut <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop='ACUA', permutations=null_positives)
print(summary(comboAcuaIslePermut))


comboTlmcDesRows <- stacWin[which(stacWin$chrPos %in% comboDesPass$cMmidChrPos[which(comboDesPass$nTLMC>0)]),]
focal<-comboTlmcDesRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-tlmcWin
  pop_null$quant<-sample(tlmcWin$quant)
  pass<-subset(pop_null,pop_null$quant=="lower10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboTlmcDesPermut <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop='TLMC', permutations=null_positives)
print(summary(comboTlmcDesPermut))

comboTlmcIsleRows <- stacWin[which(stacWin$chrPos %in% comboIslePass$cMmidChrPos[which(comboIslePass$nTLMC>0)]),]
focal<-comboTlmcIsleRows
null_positives<-{}
for(x in 1:1000){
  pop_null<-tlmcWin
  pop_null$quant<-sample(tlmcWin$quant)
  pass<-subset(pop_null,pop_null$quant=="upper10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboTlmcIslePermut <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop='TLMC', permutations=null_positives)
print(summary(comboTlmcIslePermut))


comboAgzcAcuaDesRows <- stacWin[which(stacWin$chrPos %in% comboDesPass$cMmidChrPos[which(comboDesPass$nACUA>0 & comboDesPass$nAGZC>0)]),]
focal<-comboAgzcAcuaDesRows
null_positives<-{}
for(x in 1:1000){
  acua_null<-acuaWin
  agzc_null<-agzcWin
  acua_null$quant<-sample(acuaWin$quant)
  agzc_null$quant<-sample(agzcWin$quant)
  pass<-subset(acua_null,acua_null$quant=="lower10"  & agzc_null$quant=="lower10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAgzcAcuaDesPermut <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop='AGZC+ACUA', permutations=null_positives)
print(summary(comboAgzcAcuaDesPermut))

comboAgzcAcuaIsleRows <- stacWin[which(stacWin$chrPos %in% comboIslePass$cMmidChrPos[which(comboIslePass$nACUA>0 & comboIslePass$nAGZC>0)]),]
focal<-comboAgzcAcuaIsleRows
null_positives<-{}
for(x in 1:1000){
  acua_null<-acuaWin
  agzc_null<-agzcWin
  acua_null$quant<-sample(acuaWin$quant)
  agzc_null$quant<-sample(agzcWin$quant)
  pass<-subset(acua_null,acua_null$quant=="upper10"  & agzc_null$quant=="upper10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAgzcAcuaIslePermut <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop='AGZC+ACUA', permutations=null_positives)
print(summary(comboAgzcAcuaIslePermut))

comboAgzcAcuaTlmcDesRows <- stacWin[which(stacWin$chrPos %in% comboDesPass$cMmidChrPos[which(comboDesPass$nACUA>0 & comboDesPass$nAGZC>0 & comboDesPass$nTLMC>0)]),]
focal<-comboAgzcAcuaTlmcDesRows
null_positives<-{}
for(x in 1:1000){
  acua_null<-acuaWin
  agzc_null<-agzcWin
  tlmc_null<-tlmcWin
  acua_null$quant<-sample(acuaWin$quant)
  agzc_null$quant<-sample(agzcWin$quant)
  tlmc_null$quant<-sample(tlmcWin$quant)
  pass<-subset(acua_null,acua_null$quant=="lower10" & agzc_null$quant=="lower10" & tlmc_null$quant=="lower10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAgzcAcuaTlmcDesPermut <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop='AGZC+ACUA+TLMC', permutations=null_positives)
print(summary(comboAgzcAcuaTlmcDesPermut))

comboAgzcAcuaTlmcIsleRows <- stacWin[which(stacWin$chrPos %in% comboIslePass$cMmidChrPos[which(comboIslePass$nACUA>0 & comboIslePass$nAGZC>0& comboIslePass$nTLMC>0)]),]
focal<-comboAgzcAcuaTlmcIsleRows
null_positives<-{}
for(x in 1:1000){
  acua_null<-acuaWin
  agzc_null<-agzcWin
  tlmc_null<-tlmcWin
  acua_null$quant<-sample(acuaWin$quant)
  agzc_null$quant<-sample(agzcWin$quant)
  tlmc_null$quant<-sample(tlmcWin$quant)
  pass<-subset(acua_null,acua_null$quant=="upper10"  & agzc_null$quant=="upper10"  & tlmc_null$quant=="upper10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAgzcAcuaTlmcIslePermut <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop='AGZC+ACUA+TLMC', permutations=null_positives)
print(summary(comboAgzcAcuaTlmcIslePermut))

comboAnyAgzcAcuaTlmcDesRows <- stacWin[which(stacWin$chrPos %in% comboDesPass$cMmidChrPos[which(comboDesPass$nACUA>0 | comboDesPass$nAGZC>0 | comboDesPass$nTLMC>0)]),]
focal<-comboAnyAgzcAcuaTlmcDesRows
null_positives<-{}
for(x in 1:1000){
  acua_null<-acuaWin
  agzc_null<-agzcWin
  tlmc_null<-tlmcWin
  acua_null$quant<-sample(acuaWin$quant)
  agzc_null$quant<-sample(agzcWin$quant)
  tlmc_null$quant<-sample(tlmcWin$quant)
  pass<-subset(acua_null,acua_null$quant=="lower10" | agzc_null$quant=="lower10" | tlmc_null$quant=="lower10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAnyAgzcAcuaTlmcDesPermut <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop='AGZC_ACUA_TLMC', permutations=null_positives)
print(summary(comboAnyAgzcAcuaTlmcDesPermut))

comboAnyAgzcAcuaTlmcIsleRows <- stacWin[which(stacWin$chrPos %in% comboIslePass$cMmidChrPos[which(comboIslePass$nACUA>0 | comboIslePass$nAGZC>0 | comboIslePass$nTLMC>0)]),]
focal<-comboAnyAgzcAcuaTlmcIsleRows
null_positives<-{}
for(x in 1:1000){
  acua_null<-acuaWin
  agzc_null<-agzcWin
  tlmc_null<-tlmcWin
  acua_null$quant<-sample(acuaWin$quant)
  agzc_null$quant<-sample(agzcWin$quant)
  tlmc_null$quant<-sample(tlmcWin$quant)
  pass<-subset(acua_null,acua_null$quant=="upper10"  | agzc_null$quant=="upper10"  | tlmc_null$quant=="upper10")
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
comboAnyAgzcAcuaTlmcIslePermut <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop='AGZC_ACUA_TLMC', permutations=null_positives)
print(summary(comboAnyAgzcAcuaTlmcIslePermut))

permutations <- rbind(stacDesPermut, stacIslePermut, chplDesPermut, chplIslePermut, comboAgzcDesPermut, comboAgzcIslePermut, comboAcuaDesPermut, comboAcuaIslePermut, comboTlmcDesPermut, comboTlmcIslePermut, comboAgzcAcuaDesPermut, comboAgzcAcuaIslePermut, comboAgzcAcuaTlmcDesPermut, comboAgzcAcuaTlmcIslePermut, comboAnyAgzcAcuaTlmcDesPermut, comboAnyAgzcAcuaTlmcIslePermut)
write.table(permutations, "desert_islands_shared_permutations_Nov14.txt", row.name=F, sep="\t", quote=F)

###Jackknife
chunkSize <- 200

winCM$chrPos <- paste(winCM$chrID, winCM$start, sep="_")
uniChrs <- unique(winCM$chrID)

#stacDesCross <- subset(stacDes, chplN>0 & cMmidPass==T)
stacDesCross <- comboDesPass
stacDesCrossCount <- nrow(stacDesCross)
stacDesJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- stacDesCross[which(!(stacDesCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="STAC_Desert", focalPop="STAC", regType="Desert", compPop="CHPL", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    popCount <- nrow(subset(sharedOutside, nChpl>0))
    tempDF <- cbind(tempDF, popCount=popCount)
    stacDesJackknifeDF <- rbind(stacDesJackknifeDF, tempDF)
    start = start+chunkSize
  }
}
summary(stacDesJackknifeDF)

#stacIsleCross <- subset(stacIsle, chplN>0 & cMmidPass==T)
stacIsleCross <- comboIslePass
stacIsleCrossCount <- nrow(stacIsleCross)
stacIsleJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- stacIsleCross[which(!(stacIsleCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="STAC_Island", focalPop="STAC", regType="Island", compPop="CHPL", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    popCount <- nrow(subset(sharedOutside, nChpl>0))
    tempDF <- cbind(tempDF, popCount=popCount)
    stacIsleJackknifeDF <- rbind(stacIsleJackknifeDF, tempDF)
    start = start+chunkSize
  }
}

#chplDesCross <- subset(chplDes, stacN>0 & cMmidPass==T)
chplDesCross <- comboDesPass
chplDesCrossCount <- nrow(chplDesCross)
chplDesJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- chplDesCross[which(!(chplDesCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="CHPL_Desert", focalPop="CHPL", regType="Desert", compPop="STAC", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    popCount <- nrow(subset(sharedOutside, nStac>0))
    tempDF <- cbind(tempDF, popCount=popCount)
    chplDesJackknifeDF <- rbind(chplDesJackknifeDF, tempDF)
    start = start+chunkSize
  }
}

#chplIsleCross <- subset(chplIsle, stacN>0 & cMmidPass==T)
chplIsleCross <- comboIslePass
chplIsleCrossCount <- nrow(chplIsleCross)
chplIsleJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- chplIsleCross[which(!(chplIsleCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="CHPL_Island", focalPop="CHPL", regType="Island", compPop="STAC", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    popCount <- nrow(subset(sharedOutside, nStac>0))
    tempDF <- cbind(tempDF, popCount=popCount)
    chplIsleJackknifeDF <- rbind(chplIsleJackknifeDF, tempDF)
    start = start+chunkSize
  }
}

comboAgzcDesCross <- subset(comboDesPass, nAGZC>0)
comboAgzcDesCrossCount <- nrow(comboAgzcDesCross)
comboAgzcDesJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- comboAgzcDesCross[which(!(comboAgzcDesCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop="AGZC", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    agzcCount <- nrow(subset(sharedOutside, nAGZC>0))
    tempDF <- cbind(tempDF, popCount=agzcCount)
    comboAgzcDesJackknifeDF <- rbind(comboAgzcDesJackknifeDF, tempDF)
    start = start+chunkSize
  }
}

comboAgzcIsleCross <- subset(comboIslePass, nAGZC>0)
comboAgzcIsleCrossCount <- nrow(comboAgzcIsleCross)
comboAgzcIsleJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- comboAgzcIsleCross[which(!(comboAgzcIsleCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop="AGZC", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    agzcCount <- nrow(subset(sharedOutside, nAGZC>0))
    tempDF <- cbind(tempDF, popCount=agzcCount)
    comboAgzcIsleJackknifeDF <- rbind(comboAgzcIsleJackknifeDF, tempDF)
    start = start+chunkSize
  }
}

comboAcuaDesCross <- subset(comboDesPass, nACUA>0)
comboAcuaDesCrossCount <- nrow(comboAcuaDesCross)
comboAcuaDesJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- comboAcuaDesCross[which(!(comboAcuaDesCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop="ACUA", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    acuaCount <- nrow(subset(sharedOutside, nACUA>0))
    tempDF <- cbind(tempDF, popCount=acuaCount)
    comboAcuaDesJackknifeDF <- rbind(comboAcuaDesJackknifeDF, tempDF)
    start = start+chunkSize
  }
}
summary(comboAcuaDesJackknifeDF)

comboAcuaIsleCross <- subset(comboIslePass, nACUA>0)
comboAcuaIsleCrossCount <- nrow(comboAcuaIsleCross)
comboAcuaIsleJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- comboAcuaIsleCross[which(!(comboAcuaIsleCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop="ACUA", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    acuaCount <- nrow(subset(sharedOutside, nACUA>0))
    tempDF <- cbind(tempDF, popCount=acuaCount)
    comboAcuaIsleJackknifeDF <- rbind(comboAcuaIsleJackknifeDF, tempDF)
    start = start+chunkSize
  }
}
summary(comboAcuaIsleJackknifeDF)

comboTlmcDesCross <- subset(comboDesPass, nTLMC>0)
comboTlmcDesCrossCount <- nrow(comboTlmcDesCross)
comboTlmcDesJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- comboTlmcDesCross[which(!(comboTlmcDesCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop="TLMC", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    tlmcCount <- nrow(subset(sharedOutside, nTLMC>0))
    tempDF <- cbind(tempDF, popCount=tlmcCount)
    comboTlmcDesJackknifeDF <- rbind(comboTlmcDesJackknifeDF, tempDF)
    start = start+chunkSize
  }
}
summary(comboTlmcDesJackknifeDF)

comboTlmcIsleCross <- subset(comboIslePass, nTLMC>0)
comboTlmcIsleCrossCount <- nrow(comboTlmcIsleCross)
comboTlmcIsleJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- comboTlmcIsleCross[which(!(comboTlmcIsleCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop="TLMC", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    tlmcCount <- nrow(subset(sharedOutside, nTLMC>0))
    tempDF <- cbind(tempDF, popCount=tlmcCount)
    comboTlmcIsleJackknifeDF <- rbind(comboTlmcIsleJackknifeDF, tempDF)
    start = start+chunkSize
  }
}
summary(comboTlmcIsleJackknifeDF)


comboAnyDesCross <- subset(comboDesPass, nAGZC>0 | nACUA>0 | nTLMC>0)
comboAnyDesCrossCount <- nrow(comboAnyDesCross)
comboAnyDesJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- comboAnyDesCross[which(!(comboAnyDesCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="combo_Desert", focalPop="combo", regType="Desert", compPop="Any", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    anyCount <- nrow(subset(sharedOutside, nAGZC>0 | nACUA>0 | nTLMC>0))
    tempDF <- cbind(tempDF, popCount=anyCount)
    comboAnyDesJackknifeDF <- rbind(comboAnyDesJackknifeDF, tempDF)
    start = start+chunkSize
  }
}
summary(comboAnyDesJackknifeDF)

comboAnyIsleCross <- subset(comboIslePass, nAGZC>0 | nACUA>0 | nTLMC>0)
comboAnyIsleCrossCount <- nrow(comboAnyIsleCross)
comboAnyIsleJackknifeDF <- data.frame()
for (chrName in uniChrs) {
  chrWins <- subset(winCM, chrID==chrName)
  chrStartRow <- as.numeric(row.names(chrWins)[1])
  chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
  start <- 1
  while(start<nrow(chrWins)) {
    chunkData <- chrWins[start:(start+chunkSize),]
    sharedOutside <- comboAnyIsleCross[which(!(comboAnyIsleCross$cMmidChrPos %in% chunkData$chrPos)),]
    outsiteTotal <- nrow(sharedOutside)
    #outsidePair <- sum(sharedOutside$pairPops)
    #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
    tempDF <- data.frame(popComp="combo_Island", focalPop="combo", regType="Island", compPop="Any", winName=chunkData$chrPos[1], crossPopTotal=outsiteTotal)
    anyCount <- nrow(subset(sharedOutside, nAGZC>0 | nACUA>0 | nTLMC>0))
    tempDF <- cbind(tempDF, popCount=anyCount)
    comboAnyIsleJackknifeDF <- rbind(comboAnyIsleJackknifeDF, tempDF)
    start = start+chunkSize
  }
}
summary(comboAnyIsleJackknifeDF)

allJackknifeDF <- rbind(stacDesJackknifeDF, stacIsleJackknifeDF, chplDesJackknifeDF, chplIsleJackknifeDF, comboAgzcDesJackknifeDF, comboAgzcIsleJackknifeDF, comboAcuaDesJackknifeDF, comboAcuaIsleJackknifeDF, comboTlmcDesJackknifeDF, comboTlmcIsleJackknifeDF, comboAnyDesJackknifeDF, comboAnyIsleJackknifeDF)
write.table(allJackknifeDF, "desert_island_shared_jackknifed_Nov14.txt", row.names = F, sep="\t", quote=F)


options(stringAsFactors=FALSE)
library(ggplot2)
library(gridExtra)

setwd("~/Documents/SchumerLab/Xcortezi_Xbirchmanni_hybridization/desertsIslandsPacBio2023/")

empRegionsAll <- read.table("combinedShared_CHPL_STAC_desertsIslands_midCMwin_crossInfo.txt", header = T)
#empRegions$cMwin <- paste(empRegions$chrName, empRegions$cMmid, sep="_")
#

allPopData <- read.table("xbir_pacbio2023_allChrs_0.05cM_windows_recRate_codingBP_conservedBP_thinnedAIMs_allPops_minPar.txt", header=T)
allPopData$cMwin <- paste(allPopData$chr, allPopData$start, sep="_")

#empRegions <- subset(empRegions, popName=="STAC2020")
empRegions <- empRegionsAll
empRegions$cMwin <- NA
for (i in 1:nrow(empRegions)) {
  midWin <- subset(allPopData, allPopData$chrID==empRegions$chrName[i] & start<=empRegions$cMmid[i] & end>=empRegions$cMmid[i])
  empRegions$cMwin[i] <- midWin$cMwin
}

empDesertsAll <- subset(empRegions, regType=="Desert")
empDeserts <- subset(empDesertsAll, stacMidPass==T & chplMidPass==T)
desertRows <- allPopData[which(allPopData$cMwin %in% empDeserts$cMwin),]

stac <- data.frame(chr=allPopData$chr, start=allPopData$start, end=allPopData$end, cMwin=allPopData$cMwin, ancestry=allPopData$STAC2020)
chpl <- data.frame(chr=allPopData$chr, start=allPopData$start, end=allPopData$end, cMwin=allPopData$cMwin, ancestry=allPopData$CHPL2021)

quant=0.1
crossPopPass<-subset(stac, stac$ancestry<=quantile(stac$ancestry,quant,na.rm=TRUE) & (chpl$ancestry<=quantile(chpl$ancestry,quant,na.rm=TRUE)))
crossPopDes <- subset(crossPopPass, crossPopPass$cMwin%in%empDeserts$cMwin)

#shift two out of three by different step sizes (~5e6)
###By pop
step1=1
step2=500
step_size=250

pop1<-stac
pop2<-chpl

quant=0.1
counts<-{}
overlap<-{}

#I think it should actually be 125 to do the full shuffle
for (k in 1:132){
  
  total=length(pop1$chr)
  
  pop1_null<-pop1
  
  s1<-total-step1
  s2<-1
  s3<-s1-1
  pop1_null$chr<-c(pop1$chr[s1:total],pop1$chr[s2:s3])
  pop1_null<-pop1_null[order(pop1_null[,1],pop1_null[,2]),]
  rownames(pop1_null)<-1:length(pop1_null[,1])
  
  step1=step1+step_size
  
  pop2_null<-pop2
  
  s1<-total-step2
  s2<-1
  s3<-s1-1
  pop2_null$chr<-c(pop2$chr[s1:total],pop2$chr[s2:s3])
  pop2_null<-pop2_null[order(pop2_null[,1],pop1_null[,2]),]
  rownames(pop2_null)<-1:length(pop2_null[,1])
  
  step2=step2+step_size
  
  low_both<-subset(cbind(pop1_null,pop2_null),pop1_null$ancestry<=quantile(pop1_null$ancestry,quant,na.rm=TRUE) & pop2_null$ancestry<=quantile(pop2_null$ancestry,quant,na.rm=TRUE) & pop1_null$chr==pop2_null$chr & as.numeric(rownames(pop1_null))==as.numeric(rownames(pop2_null)))
  print(length(low_both[,1]))
  counts<-c(counts,length(low_both[,1]))
  
  empOverlap <- subset(low_both, low_both$cMwin%in%empDeserts$cMwin)
  print(length(empOverlap[,1]))
  overlap<-c(overlap,length(empOverlap[,1]))
  
}

trueDF<-subset(cbind(pop1,pop2),(pop1$ancestry<=quantile(pop1$ancestry,quant,na.rm=TRUE) & pop2$ancestry<=quantile(pop2$ancestry,quant,na.rm=TRUE) & pop1$chr==pop2$chr & pop1$start==pop2$start))
#trueDF<-subset(cbind(pop1,pop2),((pop1$ancestry<=quantile(pop1$ancestry,quant,na.rm=TRUE) | pop2$ancestry<=quantile(pop2$ancestry,quant,na.rm=TRUE)) & pop1$chr==pop2$chr & pop1$start==pop2$start))
true <- length(trueDF$chr)
true
mean(counts)
true/mean(counts)

overlapDF <- subset(trueDF, trueDF$cMwin%in%empDeserts$cMwin)
trueOverlap <- nrow(overlapDF)
trueOverlap
summary(overlap)
trueOverlap/mean(overlap)

stac_chpl_comp<-counts

stac_chpl_true<-true

stac_chpl_overlap<-overlap

stac_chpl_trueOverlap<-trueOverlap

desertShuffle <- rbind(data.frame(count=stac_chpl_trueOverlap,popComp="STAC/CHPL",comparison="observed",regType="Deserts"),data.frame(count=stac_chpl_overlap,popComp="STAC/CHPL",comparison="shuffled",regType="Deserts"))
desertShuffle25cM <- desertShuffle
#desertShuffle50cM <- desertShuffle

ggplot() + geom_boxplot(data=subset(desertShuffle, comparison=="shuffled"), aes(x=popComp, y=count), color="grey50", alpha=0.5) + geom_jitter(data=subset(desertShuffle, comparison=="shuffled"), aes(x=popComp, y=count), color="grey50", alpha=0.5) + geom_point(data=subset(desertShuffle, comparison=="observed"), aes(x=popComp, y=count), color="black", alpha=1, size=3) + theme_bw()


##Islands
empIslandsAll <- subset(empRegions, regType=="Island")
empIslands <- subset(empIslandsAll, stacMidPass==T & chplMidPass==T)

quant=0.9
crossPopPass<-subset(stac, stac$ancestry>=quantile(stac$ancestry,quant,na.rm=TRUE) & (chpl$ancestry>=quantile(chpl$ancestry,quant,na.rm=TRUE)))
crossPopIsle <- subset(crossPopPass, crossPopPass$cMwin%in%empIslands$cMwin)


###By pop
step1=1
step2=500
step_size=250

pop1<-stac
pop2<-chpl

quant=0.9
counts<-{}
overlap<-{}

#I think it should actually be 125 to do the full shuffle
for (k in 1:132){
  
  total=length(pop1$chr)
  
  pop1_null<-pop1
  
  s1<-total-step1
  s2<-1
  s3<-s1-1
  pop1_null$chr<-c(pop1$chr[s1:total],pop1$chr[s2:s3])
  pop1_null<-pop1_null[order(pop1_null[,1],pop1_null[,2]),]
  rownames(pop1_null)<-1:length(pop1_null[,1])
  
  step1=step1+step_size
  
  pop2_null<-pop2
  
  s1<-total-step2
  s2<-1
  s3<-s1-1
  pop2_null$chr<-c(pop2$chr[s1:total],pop2$chr[s2:s3])
  pop2_null<-pop2_null[order(pop2_null[,1],pop1_null[,2]),]
  rownames(pop2_null)<-1:length(pop2_null[,1])
  
  step2=step2+step_size
  
  low_both<-subset(cbind(pop1_null,pop2_null),pop1_null$ancestry>=quantile(pop1_null$ancestry,quant,na.rm=TRUE) & pop2_null$ancestry>=quantile(pop2_null$ancestry,quant,na.rm=TRUE) & pop1_null$chr==pop2_null$chr & as.numeric(rownames(pop1_null))==as.numeric(rownames(pop2_null)))
  print(length(low_both[,1]))
  counts<-c(counts,length(low_both[,1]))
  
  empOverlap <- subset(low_both, low_both$cMwin%in%empIslands$cMwin)
  print(length(empOverlap[,1]))
  overlap<-c(overlap,length(empOverlap[,1]))
  
}

trueDF<-subset(cbind(pop1,pop2),(pop1$ancestry>=quantile(pop1$ancestry,quant,na.rm=TRUE) & pop2$ancestry>=quantile(pop2$ancestry,quant,na.rm=TRUE) & pop1$chr==pop2$chr & pop1$start==pop2$start))
true <- length(trueDF$chr)
true
mean(counts)
true/mean(counts)

overlapDF <- subset(trueDF, trueDF$cMwin%in%empIslands$cMwin)
trueOverlap <- nrow(overlapDF)
trueOverlap
summary(overlap)
trueOverlap/mean(overlap)

stac_chpl_comp<-counts

stac_chpl_true<-true

stac_chpl_overlap<-overlap

stac_chpl_trueOverlap<-trueOverlap

islandShuffle <- rbind(data.frame(count=stac_chpl_trueOverlap,popComp="STAC/CHPL",comparison="observed",regType="Islands"),data.frame(count=stac_chpl_overlap,popComp="STAC/CHPL",comparison="shuffled",regType="Islands"))
islandShuffle25cM <- islandShuffle

ggplot() + geom_boxplot(data=subset(islandShuffle, comparison=="shuffled"), aes(x=popComp, y=count), color="grey50", alpha=0.5, width=0.5) + geom_jitter(data=subset(islandShuffle, comparison=="shuffled"), aes(x=popComp, y=count), color="grey50", alpha=0.5) + geom_point(data=subset(islandShuffle, comparison=="observed"), aes(x=popComp, y=count), color="black", alpha=1, size=3) + theme_bw()

shuffledData <- rbind(desertShuffle, islandShuffle)

simpleColors <- c("shuffled"="grey50", "observed"="black")
simpleColor <- scale_color_manual(values=simpleColors, name="Comparison")

ggplot() + geom_boxplot(data=subset(shuffledData, comparison=="shuffled"), aes(x=popComp, y=count, color=comparison), alpha=0.5, width=0.5) + geom_jitter(data=subset(shuffledData, comparison=="shuffled"), aes(x=popComp, y=count, color=comparison), alpha=0.5) + geom_point(data=subset(shuffledData, comparison=="observed"), aes(x=popComp, y=count, color=comparison), alpha=1, size=4) + simpleColor + theme_bw() + facet_wrap(~regType) + xlab("Comparison") + ylab("Shared Regions Detected")
ggsave("desertIslands_2023_CHPL_STAC.pdf")

write.table(shuffledData, file="STAC_CHPL_desertsIslandCounts_shuffled.txt",row.names=FALSE,quote=FALSE,sep="\t")













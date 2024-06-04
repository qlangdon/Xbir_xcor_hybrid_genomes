options(stringAsFactors=FALSE)
library(ggplot2)
setwd("~/Documents/SchumerLab/Xcortezi_Xbirchmanni_hybridization/ABCdemography/ABCregPower/")

chplTrue <- read.table("holdouts_CHPL2021_fullInfo.txt", header=T)
chplTrue$popsize95 <- NA
chplTrue$popsize50 <- NA
chplTrue$gen95 <- NA
chplTrue$gen50 <- NA
chplTrue$init_prop95 <- NA
chplTrue$init_prop50 <- NA
chplTrue$mig195 <- NA
chplTrue$mig150 <- NA
chplTrue$mig295 <- NA
chplTrue$mig250 <- NA

chplTrueLong <- data.frame()
for (i in 0:99) {
  curTrue <- chplTrue[i+1,]
  curABCreg <- read.table(paste("output.CHPL.holds.05." ,i, ".tangent.post", sep=""), header = F, col.names = c("popsize", "gen", "init_prop", "mig1", "mig2"))
  popsize95 <- curTrue$popsize>=quantile(curABCreg$popsize, c(0.025)) & curTrue$popsize<=quantile(curABCreg$popsize, c(0.975))
  popsize50 <- curTrue$popsize>=quantile(curABCreg$popsize, c(0.25)) & curTrue$popsize<=quantile(curABCreg$popsize, c(0.75))
  gen95 <- curTrue$gen>=quantile(curABCreg$gen, c(0.025)) & curTrue$gen<=quantile(curABCreg$gen, c(0.975))
  gen50 <- curTrue$gen>=quantile(curABCreg$gen, c(0.25)) & curTrue$gen<=quantile(curABCreg$gen, c(0.75)) 
  init_prop95 <- curTrue$init_prop>=quantile(curABCreg$init_prop, c(0.025)) & curTrue$init_prop<=quantile(curABCreg$init_prop, c(0.975))
  init_prop50 <- curTrue$init_prop>=quantile(curABCreg$init_prop, c(0.25)) & curTrue$init_prop<=quantile(curABCreg$init_prop, c(0.75))
  mig195 <- curTrue$mig1>=quantile(curABCreg$mig1, c(0.025)) & curTrue$mig1<=quantile(curABCreg$mig1, c(0.975))
  mig150 <- curTrue$mig1>=quantile(curABCreg$mig1, c(0.25)) & curTrue$mig1<=quantile(curABCreg$mig1, c(0.75))
  mig295 <- curTrue$mig2>=quantile(curABCreg$mig2, c(0.025)) & curTrue$mig2<=quantile(curABCreg$mig2, c(0.975))
  mig250 <- curTrue$mig2>=quantile(curABCreg$mig2, c(0.25)) & curTrue$mig2<=quantile(curABCreg$mig2, c(0.75))
  
  curTrueSize <- cbind(curTrue, prior="popsize", pass95=popsize95, pass50=popsize50)
  curTrueGen <- cbind(curTrue, prior="gen", pass95=gen95, pass50=gen50)
  curTrueInit <- cbind(curTrue, prior="init_prop", pass95=init_prop95, pass50=init_prop50)
  curTrueMig1 <- cbind(curTrue, prior="mig1", pass95=mig195, pass50=mig150)
  curTrueMig2 <- cbind(curTrue, prior="mig2", pass95=mig295, pass50=mig250)
  chplTrueLong <- rbind(chplTrueLong, curTrueSize, curTrueGen, curTrueInit, curTrueMig1, curTrueMig2)
  
  chplTrue$popsize95[i+1] <- popsize95
  chplTrue$popsize50[i+1] <- popsize50
  chplTrue$gen95[i+1] <- gen95
  chplTrue$gen50[i+1] <- gen50
  chplTrue$init_prop95[i+1] <- init_prop95
  chplTrue$init_prop50[i+1] <- init_prop50
  chplTrue$mig195[i+1] <- mig195
  chplTrue$mig150[i+1] <- mig150
  chplTrue$mig295[i+1] <- mig295
  chplTrue$mig250[i+1] <- mig250
  
}

summary(chplTrue)
write.table(chplTrue, "holdouts_CHPL2021_fullInfo_wPass.txt", row.names = F, quote = F, sep="\t")

chplSum <- data.frame(pop="CHPL", prior=c("popsize", "popsize", "gen", "gen", "init_prop", "init_prop", "mig1", "mig1", "mig2", "mig2"), criteria=c("inner95", "inner50", "inner95", "inner50", "inner95", "inner50", "inner95", "inner50", "inner95", "inner50"), propPass=c(sum(chplTrue$popsize95)/100, sum(chplTrue$popsize50)/100, sum(chplTrue$gen95)/100, sum(chplTrue$gen50)/100, sum(chplTrue$init_prop95)/100, sum(chplTrue$init_prop50)/100, sum(chplTrue$mig195)/100, sum(chplTrue$mig150)/100, sum(chplTrue$mig295)/100, sum(chplTrue$mig250)/100))

ggplot(chplSum) + geom_bar(aes(x=prior, y=propPass, fill=prior), color="black", stat="identity", position="dodge") + theme_bw() + facet_wrap(~criteria) + guides(fill="none")

#Same with STAC
stacTrue <- read.table("holdouts_STAC2020_fullInfo.txt", header=T)
stacTrue$popsize95 <- NA
stacTrue$popsize50 <- NA
stacTrue$gen95 <- NA
stacTrue$gen50 <- NA
stacTrue$init_prop95 <- NA
stacTrue$init_prop50 <- NA
stacTrue$mig195 <- NA
stacTrue$mig150 <- NA
stacTrue$mig295 <- NA
stacTrue$mig250 <- NA

stacTrueLong <- data.frame()
for (i in 0:92) {
  curTrue <- stacTrue[i+1,]
  curABCreg <- read.table(paste("output.CHPL.holds.05." ,i, ".tangent.post", sep=""), header = F, col.names = c("popsize", "gen", "init_prop", "mig1", "mig2"))
  popsize95 <- curTrue$popsize>=quantile(curABCreg$popsize, c(0.025)) & curTrue$popsize<=quantile(curABCreg$popsize, c(0.975))
  popsize50 <- curTrue$popsize>=quantile(curABCreg$popsize, c(0.25)) & curTrue$popsize<=quantile(curABCreg$popsize, c(0.75))
  gen95 <- curTrue$gen>=quantile(curABCreg$gen, c(0.025)) & curTrue$gen<=quantile(curABCreg$gen, c(0.975))
  gen50 <- curTrue$gen>=quantile(curABCreg$gen, c(0.25)) & curTrue$gen<=quantile(curABCreg$gen, c(0.75)) 
  init_prop95 <- curTrue$init_prop>=quantile(curABCreg$init_prop, c(0.025)) & curTrue$init_prop<=quantile(curABCreg$init_prop, c(0.975))
  init_prop50 <- curTrue$init_prop>=quantile(curABCreg$init_prop, c(0.25)) & curTrue$init_prop<=quantile(curABCreg$init_prop, c(0.75))
  mig195 <- curTrue$mig1>=quantile(curABCreg$mig1, c(0.025)) & curTrue$mig1<=quantile(curABCreg$mig1, c(0.975))
  mig150 <- curTrue$mig1>=quantile(curABCreg$mig1, c(0.25)) & curTrue$mig1<=quantile(curABCreg$mig1, c(0.75))
  mig295 <- curTrue$mig2>=quantile(curABCreg$mig2, c(0.025)) & curTrue$mig2<=quantile(curABCreg$mig2, c(0.975))
  mig250 <- curTrue$mig2>=quantile(curABCreg$mig2, c(0.25)) & curTrue$mig2<=quantile(curABCreg$mig2, c(0.75))
  
  curTrueSize <- cbind(curTrue, prior="popsize", pass95=popsize95, pass50=popsize50)
  curTrueGen <- cbind(curTrue, prior="gen", pass95=gen95, pass50=gen50)
  curTrueInit <- cbind(curTrue, prior="init_prop", pass95=init_prop95, pass50=init_prop50)
  curTrueMig1 <- cbind(curTrue, prior="mig1", pass95=mig195, pass50=mig150)
  curTrueMig2 <- cbind(curTrue, prior="mig2", pass95=mig295, pass50=mig250)
  stacTrueLong <- rbind(stacTrueLong, curTrueSize, curTrueGen, curTrueInit, curTrueMig1, curTrueMig2)
  
  stacTrue$popsize95[i+1] <- popsize95
  stacTrue$popsize50[i+1] <- popsize50
  stacTrue$gen95[i+1] <- gen95
  stacTrue$gen50[i+1] <- gen50
  stacTrue$init_prop95[i+1] <- init_prop95
  stacTrue$init_prop50[i+1] <- init_prop50
  stacTrue$mig195[i+1] <- mig195
  stacTrue$mig150[i+1] <- mig150
  stacTrue$mig295[i+1] <- mig295
  stacTrue$mig250[i+1] <- mig250
  
}

summary(stacTrue)
write.table(stacTrue, "holdouts_STAC2020_fullInfo_wPass.txt", row.names = F, quote = F, sep="\t")

stacSum <- data.frame(pop="STAC", prior=c("popsize", "popsize", "gen", "gen", "init_prop", "init_prop", "mig1", "mig1", "mig2", "mig2"), criteria=c("inner95", "inner50", "inner95", "inner50", "inner95", "inner50", "inner95", "inner50", "inner95", "inner50"), propPass=c(sum(stacTrue$popsize95)/94, sum(stacTrue$popsize50)/94, sum(stacTrue$gen95)/94, sum(stacTrue$gen50)/94, sum(stacTrue$init_prop95)/94, sum(stacTrue$init_prop50)/94, sum(stacTrue$mig195)/94, sum(stacTrue$mig150)/94, sum(stacTrue$mig295)/94, sum(stacTrue$mig250)/94))

ggplot(stacSum) + geom_bar(aes(x=prior, y=propPass, fill=prior), color="black", stat="identity", position="dodge") + theme_bw() + facet_wrap(~criteria) + guides(fill="none")

popsSum <- rbind(chplSum, stacSum)
ggplot(popsSum) + geom_bar(aes(x=prior, y=propPass, fill=prior), color="black", stat="identity", position="dodge") + theme_bw() + facet_grid(criteria~pop) + guides(fill="none")
ggsave("ABCreg_power.pdf")

write.table(popsSum, "ABCreg_power_passingSummary.txt", row.names=F, quote=F, sep="\t")


options(stringAsFactors=FALSE)

seedsToCheck <- read.table("seedsToTest.txt", header=F)
stacPriors <- read.table("STAC2019_2020_priors.txt", header=F, col.names = c("popsize", "gen", "init_prop", "mig1", "mig2", "tractlen", "meanAnc", "covAnc"))
stacSum <- read.csv("summary_params_STAC2019_2020.txt", header=F, sep="", col.names = c("seed", "numtracts", "tracklen", "meanAnc", "varAnc", "covAnc"))
stacPriors$seed <- NA
stacPriors$seed <- stacSum$seed[match(stacPriors$covAnc, stacSum$covAnc)]

chplPriors <- read.table("CHPL2021_priors.txt", header=F, col.names = c("popsize", "gen", "init_prop", "mig1", "mig2", "tractlen", "meanAnc", "covAnc"))
chplSum <- read.csv("summary_params_CHPL2021.txt", header=F, sep="", col.names = c("seed", "numtracts", "tracklen", "meanAnc", "varAnc", "covAnc"))
chplPriors$seed <- NA
chplPriors$seed <- chplSum$seed[match(chplPriors$covAnc, chplSum$covAnc)]

sharedSeeds <- unique(chplPriors$seed, stacPriors$seed)

newSeeds <- sample(sharedSeeds, 100)
write.table(newSeeds, "seedsToTest.txt", row.names=F, quote=F, sep="\t")

chplHoldouts <- subset(chplPriors, seed %in% newSeeds)
stacHoldouts <- subset(stacPriors, seed %in% newSeeds)

write.table(chplHoldouts, "holdouts_CHPL2021_fullInfo.txt", row.names=F, quote=F, sep="\t")
write.table(stacHoldouts, "holdouts_STAC2020_fullInfo.txt", row.names=F, quote=F, sep="\t")

chplHoldPost <- cbind(chplHoldouts$tractlen, chplHoldouts$meanAnc, chplHoldouts$covAnc)
stacHoldPost <- cbind(stacHoldouts$tractlen, stacHoldouts$meanAnc, stacHoldouts$covAnc)

write.table(chplHoldPost, "holdouts_CHPL2021_data.txt", row.names=F, quote=F, sep="\t", col.names = F)
write.table(stacHoldPost, "holdouts_STAC2020_data.txt", row.names=F, quote=F, sep="\t", col.names = F)

chplPriorsNoHold <- subset(chplPriors, !(seed %in% newSeeds))
stacPriorsNoHold <- subset(stacPriors, !(seed %in% newSeeds))

write.table(chplPriorsNoHold[,1:7], "CHPL2021_priors_noHoldout.txt", row.names=F, quote=F, sep="\t", col.names = F)
write.table(stacPriorsNoHold[,1:7], "STAC2020_priors_noHoldout.txt", row.names=F, quote=F, sep="\t", col.names = F)

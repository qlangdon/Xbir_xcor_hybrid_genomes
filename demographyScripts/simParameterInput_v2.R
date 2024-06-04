###########################
###new ABC - reviewer 1
##########################

#install.packages("KScorrect")
library("KScorrect")

sim<-1:2e6
popsize<-round(runif(2e6,2,10000))
gen<-round(runif(2e6,10,400))
admixprop<-runif(2e6,0.5,1)
mig1<-rlunif(2e6, 1e-6, 0.03, base = exp(10))
mig2<-rlunif(2e6, 1e-6, 0.03, base = exp(10))

params<-cbind(sim,popsize,gen,admixprop,mig1,mig2)

write.table(params,file="ABC_2M_initprop_0.5-1_sim_params_mig_v3.csv",sep=",",row.names=FALSE,col.names=FALSE,quote=FALSE)

library(data.table)
library(parallel)
library(tidyverse)


dams=fread("Data/illustration2_empirical.csv")

#--------------does diversification help?---------------------------

netposition_stationary=mean(dams$stationary)-dams$stationary
netposition_nonstationary=mean(dams$non_stationary)-dams$non_stationary

#insurer has similar coverage in multiple markets facing same distribution of risk but independent from each other
nmarkets=1:10

nsamps=1000000

divers_function=function(markets,risk=netposition_nonstationary,samples=nsamps){
  samps=split(sample(risk,samples*markets,replace=TRUE),1:samples)
  return(unlist(lapply(samps,FUN="mean")))
}

divers_function_stationary=function(markets,risk=netposition_stationary,samples=nsamps){
  samps=split(sample(risk,samples*markets,replace=TRUE),1:samples)
  return(unlist(lapply(samps,FUN="mean")))
}

cl=makeCluster(10)
clusterExport(cl,list("divers_function","netposition_nonstationary","nsamps","nmarkets","netposition_stationary"))

diversification=parSapply(cl,nmarkets,FUN="divers_function")
diversification_stationary=parSapply(cl,nmarkets,FUN="divers_function_stationary")

stopCluster(cl)

#variance in net position - diversification does decrease variance
#normalize by variance by variance in one market under stationarity
norm=var(diversification_stationary[,1])
x11()
par(mar=c(5,5,3,2))
plot(x=nmarkets,y=apply(diversification,MARGIN=2,FUN=function(x) var(x)/norm),ylim=c(0,3.5),type="o",col="#FFE74C",pch=20,xlab="Number of Markets",ylab="Variance Insurer Position\n(Relative to One Market Under Stationarity)",las=1,cex=1.4)
points(x=nmarkets,y=apply(diversification_stationary,MARGIN=2,FUN=function(x) var(x)/norm),col="#44BBA4",pch=20)
lines(x=nmarkets,y=apply(diversification_stationary,MARGIN=2,FUN=function(x) var(x)/norm),col="#44BBA4")
legend("topright",col=c("#44BBA4","#FFE74C"),legend=c("Assumed Stationarity","Potential Non-Stationarity"),bty="n",cex=1.4,lwd=1,pch=20)

#probability of bankruptcy - diversification increases probability of net losses
plot(x=nmarkets,y=apply(diversification,MARGIN=2,function(x) {sum(x<0)/length(x)}))

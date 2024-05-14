library(data.table)
library(parallel)
library(tidyverse)

dams=fread("Data/illustration2_empirical.csv")

#make table comparing elements of damage distributions
anydamage=round(sum(dams$non_stationary!=0)/sum(dams$stationary!=0)*100)
mean=round(mean(dams$non_stationary)/mean(dams$stationary)*100)
flooddams_stationary=dams$stationary[which(dams$stationary>0)];flooddams_nonstationary=dams$non_stationary[which(dams$non_stationary>0)]
percentiles_flooddams=round(quantile(flooddams_nonstationary,probs=c(0.25,0.5,0.75,0.9,0.95,0.975,0.99,0.995))/quantile(flooddams_stationary,probs=c(0.25,0.5,0.75,0.9,0.95,0.975,0.99,0.995))*100)

sumtable=append(append(anydamage,mean),percentiles_flooddams)

#insurers get expected losses in premiums but face long-tailed damages

#distribution of claims minus premiums

netposition_stationary=mean(dams$stationary)-dams$stationary
netposition_nonstationary=mean(dams$non_stationary)-dams$non_stationary

#-----------probability of bankruptcy - inability to pay claims--------------
pb_stationary=sum(netposition_stationary<0)/length(netposition_stationary)*100
pb_nonstationary=sum(netposition_nonstationary<0)/length(netposition_nonstationary)*100

#---------reinsurance costs -------------------------------


#an insurer / regulator is willing to accept a certain risk of bankruptcy. Up to that point, an insurer must pay a reinsurer to take on risks
#model as an insurance-linked security - insurer pays a return on the capital to the investor, if a covered loss event occurs, then investor loses capital to insurer
#reinsurance levels are determined by a naive regulator - set at a standard amount based on max of historic damage distribution

#read in daily NY Central Park rainfall (from https://www.ncdc.noaa.gov/cdo-web/)
dat=read.csv("data/CentralParkNY_DailyRainfall.csv")

#find annual maxima
dat$year=year(dat$DATE)
yearmax=dat%>%
  group_by(year)%>%
  dplyr::summarise(maxpprcp=max(PRCP))

#empirical damage function
load("Data/NYCdamagefunc.Rdat")
damagefunc=function(rain,predictfunction=predictfunc,scaling=500000000){return(predictfunction(rain)/scaling)}

threshold=-1*max(damagefunc(yearmax$maxpprcp))


#insurer must pay a return to an investor - typical return on a safe asset plus compensation for risk of loss

securityloss_stationary=ifelse(netposition_stationary>0,0,ifelse(netposition_stationary<threshold,threshold,netposition_stationary))
securityloss_nonstationary=ifelse(netposition_nonstationary>0,0,ifelse(netposition_nonstationary<threshold,threshold,netposition_nonstationary))

#change in expected loss as fraction of total
el_stationary=(sum(securityloss_stationary)/length(securityloss_stationary))/abs(threshold);el_nonstationary=sum(securityloss_nonstationary)/length(securityloss_nonstationary)/abs(threshold)

#expected loss increases from 2.5% to 3.5%

#tail value at risk99
var99_stationary=quantile(securityloss_stationary,0.01);var99_nonstationary=quantile(securityloss_nonstationary,0.01)
tvar99_stationary=sum(securityloss_stationary[which(securityloss_stationary<var99_stationary)])/length(securityloss_stationary[which(securityloss_stationary<var99_stationary)])/abs(threshold)
tvar99_nonstationary=sum(securityloss_nonstationary[which(securityloss_nonstationary<var99_nonstationary)])/length(securityloss_nonstationary[which(securityloss_nonstationary<var99_nonstationary)])/abs(threshold)

#using ILS pricing model from Land and Mahul - premium spread = expected_loss + 0.054 * TVaR99
spread_stationary=abs(el_stationary)+0.054*abs(tvar99_stationary);spread_nonstationary=abs(el_nonstationary)+0.054*abs(tvar99_nonstationary)

#spread increases from 5.6% to 8.0% under non-stationarity

#figure showing claim variance
hist(securityloss_stationary/abs(threshold)*100,col=rgb(t(col2rgb("#44BBA4"))/255,alpha=0.5),xlab="Collateral Loss (%)", yaxt="n",ylab="",main="")
hist(securityloss_nonstationary/abs(threshold)*100, col=rgb(t(col2rgb("#FFE74C"))/255,alpha=0.5),add=TRUE)
legend("topleft",fill=c("#44BBA4","#FFE74C"),legend=c("Assumed Stationarity","Potential Non-Stationarity"),bty="n",cex=1.4)

#### old version
# #consider an indemnity security - insurer pays a return on an amount of capital provided by investors. This capital is released to insurer in the event of losses
# 
# bankrupt_tolerance=0.005 #accept 0.5% chance of insurer bankruptcy
# bankrupt_threshold_stationary=quantile(netposition_stationary,probs = bankrupt_tolerance)
# bankrupt_threshold_nonstationary=quantile(netposition_nonstationary,probs = bankrupt_tolerance)
# 
# #threshold almost doubles from stationary to non-stationary
# 
# #insurer must compensates investors based on expected loss and a risk premium based on TvaR99
# 
# securityloss_stationary=ifelse(netposition_stationary>0,0,ifelse(netposition_stationary<bankrupt_threshold_stationary,bankrupt_threshold_stationary,netposition_stationary))
# securityloss_nonstationary=ifelse(netposition_nonstationary>0,0,ifelse(netposition_nonstationary<bankrupt_threshold_nonstationary,bankrupt_threshold_nonstationary,netposition_nonstationary))
# 
# #expected loss as fraction of total capital - very similar in % across stationary and nonstationary, but very different magnitudes 
# el_stationary=prob_anyloss_stationary*sum(lossfrac_stationary)/length(lossfrac_stationary)
# el_nonstationary=prob_anyloss_nonstationary*sum(lossfrac_nonstationary)/length(lossfrac_nonstationary)
# 
# #tail value at risk 99 - average loss in worst 1% of outcomes
# var99_stationary=quantile(securityloss_stationary,0.01);var99_nonstationary=quantile(securityloss_nonstationary,0.01)
# tvar99_stationary=-1*sum(securityloss_stationary[which(securityloss_stationary<var99_stationary)])/length(securityloss_stationary[which(securityloss_stationary<var99_stationary)])
# tvar99_nonstationary=-1*sum(securityloss_nonstationary[which(securityloss_nonstationary<var99_nonstationary)])/length(securityloss_nonstationary[which(securityloss_nonstationary<var99_nonstationary)])
# 
# tvar99frac_stationary=tvar99_stationary/abs(bankrupt_threshold_stationary)*100
# tvar99frac_nonstationary=tvar99_nonstationary/abs(bankrupt_threshold_nonstationary)*100

library(tidyverse)
library(fitdistrplus)
library(invgamma)
library(patchwork)
library(tidyverse)
library(parallel)
library(data.table)
library(ggpattern)

#---------Set Up: Dataset and Damage Function ------------
#assume rainfall events drawn from a Weibull distribution 

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


#---------Illustration 2: Known Non-stationarity -----------

datlength=30
nyear=nrow(yearmax)

climatology=yearmax[(nyear-datlength+1):nyear,]
climdist=fitdist(climatology$maxpprcp,"weibull",method="mle",lower=c(0,0))
climshape=climdist$estimate[1] #known shape parameter

#---establish initial period prior over shape parameter ------
priorclim=yearmax[(nyear-2*datlength+1):(nyear-datlength),]
priorpriorclim=yearmax[(nyear-3*datlength+1):(nyear-2*datlength),]

#find initial period posterior over scale parameter based on prior 30 year climatology
prior_shape=1.5
priorclimdist=fitdist(priorpriorclim$maxpprcp,"weibull",method="mle",lower=c(0,0))
prior_rate=(priorclimdist$estimate[2]^priorclimdist$estimate[1])*(prior_shape-1)

#posterior given period 30 year period
n=datlength;dat=priorclim$maxpprcp
#posterior distribution is an inverse gamma distribution with the following parameters (based on https://www.johndcook.com/CompendiumOfConjugatePriors.pdf)
post_shape=prior_shape+n+1; post_rate=prior_rate+sum(dat^climshape)
#create prior for learning based on posterior from prior 30 year period
theta_postdist=rinvgamma(100000,post_shape,post_rate)
scale_postdist=theta_postdist^(1/climshape)

#----------Modify Weather Data to Show Effect of Extreme Events on Expected Loss ---------
#set most recent year to be an unprecedented rain event in the climatology
climatology$maxpprcp[nrow(climatology)]=max(climatology$maxpprcp)+0.3

#------Time Trend Learning Model ------------

#priors over time trends are normally-distributed, centered at zero
#standard deviation is such that 95% CI includes increase or decrease of 1 in scale parameter over 30 years, compared to current 3.7
priortrend_mean=0;priortrend_sd=0.5/datlength

nsamps=4000 #number of samples of trend and initial scale parameter

samps=data.frame(trend_samp=rep(rnorm(nsamps,priortrend_mean,priortrend_sd),each=nsamps),scale_samps=sample(scale_postdist,nsamps^2,replace=TRUE))

bayeslikelihood=function(sampdat,dat=climatology$maxpprcp,shape=climshape){
  trend=sampdat[1];scale=sampdat[2]
  years=1:length(dat);scaletrend=scale+(years-1)*trend
  scaletrend=ifelse(scaletrend<=0,0.001,scaletrend)
  #probability of data given scale parameter
  pdata=prod(dweibull(dat,shape,scaletrend))
  return(c(trend,scale,pdata))
}

cl=makeCluster(12)

clusterExport(cl,list("bayeslikelihood","samps","climatology","climshape"))

probs=t(parApply(cl,samps,MARGIN=1,bayeslikelihood))

stopCluster(cl)

#probs=t(apply(samps,MARGIN=1,bayeslikelihood))
probs=as.data.frame(probs);colnames(probs)[3]="probs"
probs$probs=probs$probs/sum(probs$probs) #normalize to sum to 1

#--------Posterior Distribution Over 2024 Scale Parameter------------
postnsamps=10000

#resample joint distribution over trend and initial scale using posterior probabilities as weights
postsamp=sample(1:nrow(probs),postnsamps,prob=probs$probs,replace=TRUE)

#calculate 2024 scale distribution based on joint trend and initial scale 
scale_2024=probs$scale_samps[postsamp]+probs$trend_samp[postsamp]*datlength

#compare to scale distribution under known stationarity
#use posterior from prior 30 years as prior for current
n=datlength;dat=climatology$maxpprcp
post_shape_stationary=post_shape+n+1; post_rate_stationary=post_rate+sum(dat^climshape)
theta_postdist_stationary=rinvgamma(postnsamps,post_shape_stationary,post_rate_stationary)
scale_postdist_stationary=theta_postdist_stationary^(1/climshape)
probs_stationary=dinvgamma(theta_postdist_stationary,post_shape_stationary,post_rate_stationary)

#calculate expected damages in 2024 with and without stationarity
rainsamp=10000
rain_nonstationary=as.vector(mapply(rweibull,n=rainsamp,shape=climshape,scale=scale_2024))
rain_stationary=as.vector(mapply(rweibull,n=rainsamp,shape=climshape,scale=scale_postdist_stationary))

dams_nonstationary=damagefunc(rain_nonstationary);dams_stationary=damagefunc(rain_stationary)

premiums=data.frame(type=c("Assumed Stationarity","Possible Non-Stationarity"),meandams=c(mean(dams_stationary),mean(dams_nonstationary)),simulation=rep("Extreme",2))

dams_original=fread("illustration2_empirical.csv")

premiums=rbind(premiums,data.frame(type=c("Assumed Stationarity","Possible Non-Stationarity"),meandams=c(mean(dams_original$stationary),mean(dams_original$non_stationary)),simulation=rep("Original",2)))
#normalize by stationary original premiums
premiums$meandams=premiums$meandams/premiums$meandams[which(premiums$type=="Assumed Stationarity"&premiums$simulation=="Original")]

a=ggplot(premiums,aes(x=type,y=meandams,group=simulation,fill=type,pattern=simulation))+geom_bar_pattern(stat="identity",position="dodge",color = "black", pattern_fill = "black", pattern_angle = 45,pattern_density = 0.1,  pattern_spacing = 0.025, pattern_key_scale_factor = 0.6)+theme_classic()+
  labs(y="Change Relative to Original Data Under Stationarity",x="")+scale_fill_manual(values=c("#44BBA4","#FFE74C"),guide=NULL)+theme(text=element_text(size=14))+
  scale_pattern_manual(values=c(Extreme="stripe",Original="none"),labels=c("Additional Extreme\nObservation","Original Data"),name="")+guides(pattern = guide_legend(override.aes = list(fill = "white")))


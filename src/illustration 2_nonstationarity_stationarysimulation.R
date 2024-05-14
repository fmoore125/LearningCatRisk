library(tidyverse)
library(fitdistrplus)
library(invgamma)
library(patchwork)
library(tidyverse)
library(parallel)
library(data.table)

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

#---------Illustration 2: Possible Non-stationarity -----------

datlength=30
nyear=nrow(yearmax)

climatology=yearmax[(nyear-datlength+1):nyear,]
climdist=fitdist(climatology$maxpprcp,"weibull",method="mle",lower=c(0,0))
climshape=climdist$estimate[1] #known shape parameter

#create artificial climatology from stationary distribution
climatology$maxpprcp=rweibull(n=datlength,shape=climshape,scale=climdist$estimate[2])

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

quants=as.data.frame(rbind(quantile(scale_2024,probs=c(0.025,0.975)),quantile(scale_postdist_stationary,probs=c(0.025,0.975))))
quants$type=c("Non-Stationary","Stationary")
quants=pivot_longer(quants,cols=1:2,names_to="Quantile")

#figure of posterior scale parameters
a=ggplot(data=data.frame(scale=c(scale_2024,scale_postdist_stationary),type=c(rep("Non-Stationary",postnsamps),rep("Stationary",postnsamps))))+geom_density(aes(x=scale,group=type,col=type),lwd=0.8)+
  theme_classic()+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=unit(c(0.25,0,0.25,0.25),"inches"))+
  labs(y="Probability Density",x="Posterior Scale Parameter")+scale_color_manual(values=c("#FFE74C","#44BBA4"),labels=c("Potential Non-Stationarity","Assumed Stationarity"),name="")+
  coord_cartesian(xlim = c(2.9, 5.2), clip = "off")+
  annotate("text",x=2.9,y=1.9,label="a)",size=8)+
  geom_vline(data=quants,aes(xintercept=value,col=type),lty=2)

fwrite(data.frame(scale=c(scale_2024,scale_postdist_stationary),type=c(rep("Non-Stationary",postnsamps),rep("Stationary",postnsamps)),climatology="Constructed_Stationary"),file="Data/constructedstationaryposterior.csv")

   



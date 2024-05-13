library(tidyverse)
library(fitdistrplus)
library(invgamma)
library(patchwork)

#assume rainfall events drawn from a weibull distribution 

#read in daily NY Central Park rainfall (from https://www.ncdc.noaa.gov/cdo-web/)
dat=read.csv("data/CentralParkNY_DailyRainfall.csv")

#find annual maxima
dat$year=year(dat$DATE)
yearmax=dat%>%
  group_by(year)%>%
  dplyr::summarise(maxpprcp=max(PRCP))

#hypothetical damage function
damagefunc=function(rain,thresh=threshold,jump=0.03){
  dams=ifelse(rain<threshold,0,jump+exp(rain-1.85*threshold))
  return(dams)
}

##---------Illustration 1: Effect of Record Length Under Stationarity ------

#find top 5% of whole record for threshold
threshold=quantile(yearmax$maxpprcp,0.95)

#get posterior for multiple record length
lengths=c(10,30,100)

#chose prior shape parameter with longest tail with finite variance
prior_shape=1.5

#prior climatology (use first 30 years to establish )

posterior=matrix(nrow=length(lengths),ncol=2)
damages=numeric(length=length(lengths))

nsamp=10000 #number of samples for each posterior distribution to sample damage function for integration

for(i in 1:length(lengths)){
  print(i)
  climatology=yearmax[(155-lengths[i]):155,]
  climdist=fitdist(climatology$maxpprcp,"weibull",method="mle",lower=c(0,0))
  
  #shape parameter is known and scale parameter must be inferred
  climshape=climdist$estimate[1] #known shape parameter
  
  #establish priors based on previous years data
  if(lengths[i]<77) priorclim=yearmax[(155-2*lengths[i]):(155-lengths[i]-1),]
  if(lengths[i]>=77) priorclim=yearmax[1:(155-lengths[i]-1),]
  priorclimdist=fitdist(priorclim$maxpprcp,"weibull",method="mle",lower=c(0,0))
  #chose prior scale to match scale parameter of prior climate distribution
  #note that conjugate prior distribution is for theta = scale^shape - mean of prior should match mean of prior climate for this variable
  prior_rate=(priorclimdist$estimate[2]^priorclimdist$estimate[1])*(prior_shape-1)
  
  #data
  n=length(climatology$year);dat=climatology$maxpprcp
  
  #posterior distribution is an inverse gamma distribution with the following parameters (based on https://www.johndcook.com/CompendiumOfConjugatePriors.pdf)
  post_shape=prior_shape+n+1; post_rate=prior_rate+sum(dat^climshape)
  posterior[i,]=c(post_shape,post_rate)
  
  #convert posterior density over theta back to posterior density over Weibull scale parameter
  theta_postdist=rinvgamma(nsamp,post_shape,post_rate)
  scale_postdist=theta_postdist^(1/climshape)
  
  #computationally integrate posterior distribution over damage function to get expected damages
  means=numeric(length=nsamp)#means of damage distributions for each draw from the posterior
  for(j in 1:nsamp){
    #get damages for sampled posterior distribution
    dams=damagefunc(rweibull(nsamp,climshape,scale_postdist[j]))
    #mean of the means from each sample is equal to the mean over the whole (since each is equal sized samples)
    means[j]=mean(dams)
  }
  damages[i]=mean(means)
}

#plot uppper 10 percentile of posterior distributions under different record lengths

#show expected damages relative to 30 year, stationary

##---------Illustration 2: Hierarchical Bayesian Learning About a Trend ------







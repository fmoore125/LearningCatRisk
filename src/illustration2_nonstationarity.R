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

posteriors_real=data.frame(scale=c(scale_2024,scale_postdist_stationary),type=c(rep("Non-Stationary",postnsamps),rep("Stationary",postnsamps)),climatology="Real")

#add posteriors from constructed stationary distribution for comparison
posteriors=fread(file="Data/constructedstationaryposterior.csv")
posteriors=rbind(posteriors,posteriors_real)

#figure of posterior scale parameters
a=ggplot(data=posteriors)+geom_density(aes(x=scale,group=interaction(type,climatology),col=type,lty=climatology),lwd=0.8)+
  theme_classic()+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=unit(c(0.25,0,0.25,0.25),"inches"))+
  scale_linetype_manual(values=c(3,1),labels=c("Constructed Stationary","Real"),name="Climatology")+
  labs(y="Probability Density",x="Posterior Scale Parameter")+scale_color_manual(values=c("#FFE74C","#44BBA4"),labels=c("Potential Non-Stationarity","Assumed Stationarity"),name="Learning Model")+
  coord_cartesian(xlim = c(2.9, 5.2), clip = "off")+
  #annotate("text",x=2.9,y=1.9,label="a)",size=8)+
  geom_vline(data=quants,aes(xintercept=value,col=type),lty=1)

# #show upper 95th percentile of the posterior rainfall distribution over the damage function
# rain=seq(0,9,length.out=1000)
# upper95_nonstationary=dweibull(rain,shape=climshape,scale=quantile(scale_2024,0.975))
# upper95_stationary=dweibull(rain,shape=climshape,scale=quantile(scale_postdist_stationary,0.975))
# 
# upper95=data.frame(type=c(rep("Assumed Stationarity",length(rain)),rep("Potential Non-Stationarity",length(rain))),dist=c(upper95_stationary,upper95_nonstationary),rain=rep(rain,2))
# b=ggplot(upper95,aes(x=rain,y=dist,group=type,col=type))+geom_line()+theme_classic()+
#   theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=unit(c(0.25,0,0.25,0.25),"inches"))+
#   scale_color_manual(values=c("#FFE74C","#44BBA4"),name="")+labs(y="",x="Daily Rainfall (Annual Max, inches)")+
#   geom_line(data=data.frame(rain=rain,dam=damagefunc(rain)),aes(x=rain,y=dam),inherit.aes = FALSE,col="tomato3")+
#   annotate("text",x=c(0.2,7.5),y=c(0.2,0.3),label=c("b)","Damage\nFunction"),col=c("black","tomato3"),size=c(8,4))
#   
  
#calculate expected damages in 2024 with and without stationarity
rainsamp=10000
rain_nonstationary=as.vector(mapply(rweibull,n=rainsamp,shape=climshape,scale=scale_2024))
rain_stationary=as.vector(mapply(rweibull,n=rainsamp,shape=climshape,scale=scale_postdist_stationary))

dams_nonstationary=damagefunc(rain_nonstationary);dams_stationary=damagefunc(rain_stationary)

fwrite(data.frame(non_stationary=dams_nonstationary,stationary=dams_stationary),file="Data/illustration2_empirical.csv")

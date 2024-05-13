library(tidyverse)
library(fitdistrplus)
library(invgamma)
library(patchwork)

#assume rainfall events drawn from a gamma distribution 

#read in daily NY Central Park rainfall (from https://www.ncdc.noaa.gov/cdo-web/)
dat=read.csv("data/CentralParkNY_DailyRainfall.csv")

#find annual maxima
dat$year=year(dat$DATE)
yearmax=dat%>%
  group_by(year)%>%
  dplyr::summarise(maxpprcp=max(PRCP))


#split into first and last half

earlydat=yearmax%>%filter(year<=1945)
latedat=yearmax%>%filter(year>1945)

#fit a gamma distribution to early and late
startparams=list(shape=0.1,rate=0.1)
earlydist=fitdist(earlydat$maxpprcp,"weibull",method="mle",lower=c(0,0))
latedist=fitdist(latedat$maxpprcp,"weibull",method="mle",lower=c(0,0))

#visualize early vs late distributions
rainfall=seq(0.1,10,length.out=10000)
xearly=dweibull(rainfall,shape=earlydist$estimate[1],scale=earlydist$estimate[2])
xlate=dweibull(rainfall,shape=latedist$estimate[1],scale=latedist$estimate[2])

plot(x=rainfall,y=xearly,type="l",lwd=0.5,col="#FFB833",las=0,xlab="Annual Maximum Daily Rainfall (inches)",ylab="",yaxt="n",main="Rainfall Probability and Damage Function")
lines(x=rainfall,y=xlate,lwd=0.5,col="#2C6C8C")
legend(x=6.8,y=0.07,legend=c("1869-1945","1945-2023"),lwd=0.5,col=c("#FFB833","#2C6C8C"),bty="n")

#add threshold for 1 in 20 in early distribution
threshold=qweibull(0.95,shape=earlydist$estimate[1],scale=earlydist$estimate[2])
lines(x=rep(threshold,2),y=c(-0.05,0.2),lty=3)
text(x=threshold+0.5,y=0.22,labels="5% Rain Event \n Early Period")

#overlay damage function - damages are 0 below threshold and increase exponentially above
damagefunc=function(rain,thresh=threshold,jump=0.03){
  if(rain<thresh) return(0)
  if(rain>=thresh) return(jump+exp(rain-1.85*threshold))
}
damages=numeric(length=length(rainfall));for(i in 1:length(rainfall)) damages[i]=damagefunc(rainfall[i])
lines(x=rainfall,y=damages,lwd=1,col="#EB708D")
text(x=9.1,y=0.27,labels="Damage\nFunction",col="#EB708D")

#now make plot of damage distributions
nsamp=10000
earlysamp=rweibull(nsamp,shape=earlydist$estimate[1],scale=earlydist$estimate[2]);latesamp=rweibull(nsamp,shape=latedist$estimate[1],scale=latedist$estimate[2])
dearly=numeric(length=nsamp);for(i in 1:nsamp) dearly[i]=damagefunc(earlysamp[i])
dlate=numeric(length=nsamp);for(i in 1:nsamp) dlate[i]=damagefunc(latesamp[i])


########------------------Inference Problem -------------------------
#plot max rainfall over climatology - last 30 years
climatology=yearmax%>%filter(year>1993)
plot(climatology$year,climatology$maxpprcp,type="p", col="#2B3A67",pch=18,las=1,xlab="",ylab="Annual Maximum Daily Rainfall (inches)")

#distribution over climatology
climdist=fitdist(climatology$maxpprcp,"weibull",method="mle",lower=c(0,0))

#inference problem- assume shape parameter is known and scale parameter must be inferred
climshape=climdist$estimate[1] #known shape parameter

#conjugate prior of the shape parameter is the inverse gamma distribution
#constrain prior using previous climatology:
priorclim=yearmax%>%filter(1963<year&year<=1993)
priorclimdist=fitdist(priorclim$maxpprcp,"weibull",method="mle",lower=c(0,0))

#chose prior shape parameter with longest tail with finite variance
prior_shape=1.5

#chose prior scale to match scale parameter of prior climate distribution
#note that conjugate prior distribution is for theta = scale^shape - mean of prior should match mean of prior climate for this variable
prior_rate=(priorclimdist$estimate[2]^priorclimdist$estimate[1])*(prior_shape-1)

#data
n=length(climatology$year);dat=climatology$maxpprcp

#posterior distribution is an inverse gamma distribution with the following parameters (based on https://www.johndcook.com/CompendiumOfConjugatePriors.pdf)
post_shape=prior_shape+n+1; post_rate=prior_rate+sum(dat^climshape)

#convert posterior density over theta back to posterior density over Weibull scale parameter
theta_postdist=rinvgamma(100000,post_shape,post_rate)
scale_postdist=theta_postdist^(1/climshape)

#(compare to prior distribution)
theta_priordist=rinvgamma(100000,prior_shape,prior_rate)
scale_priordist=theta_priordist^(1/climshape)

#plot Weibull distribution for 95% confidence interval of posterior distribution
post_lower=quantile(scale_postdist,0.025);post_upper=quantile(scale_postdist,0.975)

#re-do inference but with full climatology
rainfall=seq(0.1,10,length.out=10000)
xlower=dweibull(rainfall,climshape,post_lower)
xupper=dweibull(rainfall,climshape,post_upper)

plot(x=rainfall,y=xlower,type="l",lwd=1.25,lty=2,col="#2B3A67",las=0,xlab="Annual Maximum Daily Rainfall (inches)",ylab="",yaxt="n",main="Rainfall Probability Under Upper and Lower Posterior Parameters\n(with Damage Function Overlaid)")
lines(x=rainfall,y=xupper,lwd=1.25,lty=3,col="#2B3A67")
lines(x=rainfall,y=damages,lwd=1,col="#EB708D")
text(x=9.15,y=0.26,labels="Damage\nFunction",col="#EB708D")
text(x=c(4.8,8.5),y=c(0.26,0.05),labels=c("2.5th Posterior\nQuantile","97.5th Posterior\nQuantile"),col="#2B3A67")

#expected losses under damage function for upper and lower distributions
nsamp=10000
lowersamp=rweibull(nsamp,shape=climshape,scale=post_lower);uppersamp=rweibull(nsamp,shape=climshape,scale=post_upper)
dlower=numeric(length=nsamp);for(i in 1:nsamp) dlower[i]=damagefunc(lowersamp[i])
dupper=numeric(length=nsamp);for(i in 1:nsamp) dupper[i]=damagefunc(uppersamp[i])

########------------------Inference with different lengths of records -------------------------

#create an artificial record that is *actually* stationary in order to highlight effects of record length (based on 30 year climatology)
nrecord=155
record=rweibull(nrecord,climdist$estimate[1],climdist$estimate[2])

#loop through different lengths of record and find posterior distributions - use priors from previous distribution
lengths=c(10,30,155)

posts=matrix(nrow=length(lengths),ncol=3);posts[,1]=lengths
for(i in 1:length(lengths)){
  dat=record[1:lengths[i]]
  posts[i,2]=prior_shape+length(dat)+1
  posts[i,3]=prior_rate+sum(dat^climdist$estimate[1])
}

#plot histograms of theta posteriors
nsamp=10000
dists=data.frame(recordlength=c(rep(posts[,1],each=nsamp)),value=c(rinvgamma(nsamp,posts[1,2],posts[1,3]),rinvgamma(nsamp,posts[2,2],posts[2,3]),rinvgamma(nsamp,posts[3,2],posts[3,3])))

a=ggplot(dists,aes(x=value,fill=as.factor(recordlength)))+geom_histogram(alpha=0.6,position="identity",binwidth = 10)+theme_bw()+labs(x="Theta (Weibull Distribution)",y="Frequency",fill="Record Length")
a=a+scale_fill_manual(values=c("#FED766","#FE4A49","#009FB7"))

#estimate damage ranges under posterior distributions
posts=as.data.frame(posts);colnames(posts)=c("recordlength","shape","rate")
lower=numeric(length=dim(posts)[1]);upper=numeric(length=dim(posts)[1])
#get upper and lower 95th percentile of distribution and convert back to R's Weibull parameterization
for(i in 1:length(lengths)){
  lower[i]=qinvgamma(0.025,posts$shape[i],posts$rate[i])^(1/climdist$estimate[1])
  upper[i]=qinvgamma(0.975,posts$shape[i],posts$rate[i])^(1/climdist$estimate[1])
}
posts=cbind(posts,lower);posts=cbind(posts,upper)

damages=matrix(nrow=length(lengths),ncol=2)
for(i in 1:dim(posts)[1]){
  lowerdist=rweibull(nsamp,climdist$estimate[1],posts$lower[i]);upperdist=rweibull(nsamp,climdist$estimate[1],posts$upper[i])
  temp=numeric(length=nsamp)
  for(j in 1:nsamp) {temp[j]=damagefunc(lowerdist[j]);damages[i,1]=mean(temp)}
  for(j in 1:nsamp) {temp[j]=damagefunc(upperdist[j]);damages[i,2]=mean(temp)}
}
damages=as.data.frame(damages);colnames(damages)=c("lower","upper")
damages$recordlength=posts$recordlength

b=ggplot(damages,aes(y=as.factor(recordlength),col=as.factor(recordlength),xmin=lower,xmax=upper))+geom_linerange(lwd=4,alpha=0.6)+theme_bw()
b=b+labs(y="Length of Stationary Weather Record",x="Range of Expected Damages (2.5th - 97.5th Percentile Posterior Distribution)")
b=b+scale_color_manual(values=c("#FED766","#FE4A49","#009FB7"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

x11()
a/b

library(tidyverse)
library(lubridate)

dat=read.csv("data/CentralParkNY_DailyRainfall.csv")
dat$year=year(dat$DATE)

claims=read.csv("data/nfip_nydata.csv")

claims_yearly=claims%>%
  group_by(yearOfLoss)%>%
  dplyr::summarise(across(totalClaims:contents_coverage,~sum(.x,na.rm=T)))

#get some extreme rainfall statistics to merge in
probs=c(0.8,0.85,0.9,0.95,0.99)
thresholds=quantile(dat$PRCP,probs=probs)

potfunc=function(data,thresh=thresholds,probs_thresholds=probs){
  #calculate number of days over threshold level
  peaks=c(unlist(map(thresh,function(x) sum(data>=x))),max(data))
  names(peaks)=c(probs_thresholds,"max")
  return(peaks)
}

rain_extremes=dat%>%
  group_by(year)%>%
  group_modify(~
    potfunc(.x$PRCP)%>%
      enframe()
  )%>%
  pivot_wider(id_cols = year,names_from = name,values_from = value)

fulldat=merge(rain_extremes,claims_yearly,all=FALSE,by.x="year",by.y="yearOfLoss")

fulldat$paidclaims=fulldat$buildings_paid+fulldat$contents_paid
fulldat$coverage=fulldat$buildings_coverage+fulldat$contents_coverage

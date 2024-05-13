library(tidyverse)
library(lubridate)

dat=read.csv("Data/CentralParkNY_DailyRainfall.csv")
dat$year=year(dat$DATE)

claims=read.csv("Data/NFIP/NYC_NFIP_corrected.csv")

#drop unnecessary column and spurious zip code
claims=claims[,-which(colnames(claims)=="X")]
claims=claims[-which(claims$reportedZipCode==99999),]

#impute zeros for any missing zip-code / year combinations
claims=claims%>%
  complete(reportedZipCode,year,fill=list(claims_count=0,buildings_paid=0,contents_paid=0))

#adjust to real 2020 dollar values
cpi=read.csv("Data/US_CPI_2020.csv")
cpi$year=year(as.Date(cpi$DATE))

claims=merge(claims,cpi[,2:3])

claims=claims%>%
  mutate(across(c(buildings_paid,contents_paid,buildings_coverage,contents_coverage),.fns=~(.x*100/USACPIALLMINMEI_NBD20200101)))
  
claims_yearly=claims%>%
  group_by(year)%>%
  dplyr::summarise(across(claims_count:contents_coverage,~sum(.x,na.rm=T)))

claims_yearly$totalclaims=claims_yearly$buildings_paid+claims_yearly$contents_paid
claims_yearly$totalcoverage=claims_yearly$buildings_coverage+claims_yearly$contents_coverage

#get some extreme rainfall statistics to merge in
probs=c(0.8,0.85,0.9,0.95,0.99)
thresholds=quantile(dat$PRCP,probs=probs)

potfunc=function(data,thresh=thresholds,probs_thresholds=probs){
  #calculate number of days over threshold level
  peaks=c(unlist(map(thresh,function(x) sum(data>=x))),max(data,na.rm=T))
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

fulldat=merge(rain_extremes,claims_yearly,all=FALSE,by.x="year",by.y="year")

#merge in tide data for controls
tides=read.csv("./Data/TheBattery_TideData_19782023.csv")
tides$Date.Time=as.Date(tides$Date.Time)
tides$year=year(tides$Date.Time);tides$day=yday(tides$Date.Time)
#find daily maximum dite
tides_daily=tides%>%
  group_by(year,day)%>%
  summarize(dailymax=max(Water.Level,na.rm=T))

tidemax=tides_daily%>%
  group_by(year)%>%
  group_modify(~
                 potfunc(.x$dailymax)%>%
                 enframe()
  )%>%
  pivot_wider(id_cols = year,names_from = name,values_from = value)

colnames(tidemax)[2:7]=paste0(colnames(tidemax)[2:7],"_tide")

fulldat=merge(fulldat,tidemax)

#-----------Damage function Regression---------------------

#limit to 2009 data and later, when policy / coverage data is available
fulldat_2009=fulldat%>%filter(year>=2009)

mod=lm(I(log(totalclaims))~max+I(log(totalcoverage))+max_tide,data=fulldat_2009)

#create extreme rainfall damage function based on 2023 coverage levels and average tide value
#smearing term for non-parametric Duan smearing retransformation
smear=mean(exp(mod$residuals))

intercept=mod$coefficients[1]+mod$coefficients[3]*log(fulldat$totalcoverage[which(fulldat$year==2023)])+mod$coefficients[4]*mean(fulldat_2009$max_tide)

rain_coef=mod$coefficients[2]

predictfunc=function(rain,coef=rain_coef,interceptterm=intercept,smearterm=smear){
  return(exp(interceptterm+rain*rain_coef)*smearterm)
}

save(smear,intercept,rain_coef,predictfunc,mod,fulldat_2009,file="Data/NYCdamagefunc.Rdat")




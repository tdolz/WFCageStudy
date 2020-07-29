#Survival revisited 10/14/19

#Cloned to Github 7/29/2020

library("survival")
library("car")
library("ggplot2")
library("reshape2")
library("plyr")
library("moments")
library("reshape")
library("ggfortify")
library("survminer")
library("KMsurv")
library("dplyr")
library("tidyr")
library("survminer")
library("prediction")
library("ggthemes")
library("MuMIn")
library("corrplot")
library("MASS")

setwd("/Users//tdolan/Documents//R-Github/WFCageStudy")


svwk2016 <-read.csv("weeklysummarydata_2016_10-8-19.csv", header=TRUE) # the weekly summary data, covariates are unscaled
svwk2017 <-read.csv("weeklysummarydata_2017.csv", header=TRUE) # the weekly summary data, covariates are unscaled

#####univariate plots######
#######2016 univariate models #########
#site
summary(svwk2016$site)
cox.site <-coxph(Surv(start,end,event2)~ site + cluster(cage), na.action="na.fail", data=svwk2016)
summary(cox.site) # p= very sig
zp.site <-cox.zph(cox.site)
zp.site #only RS is time varying?
plt.zp.site <-ggcoxzph(zp.site) 
fit.cx.site <-survfit(cox.site,data.frame(list(site=c("CP","MD", "RS"))))
plt.cx.site <-ggsurvplot(fit.cx.site,svwk2016,conf.int=TRUE,censor=FALSE, surv.median.line = "hv", xlim=c(2,11), palette=c("green","red","blue"), legend.labs=c("CP","MD", "RS")) + ggtitle("Site")

#depth
summary(svwk2016$depth)
cox.depth <-coxph(Surv(start,end,event2)~ depth + cluster(cage), na.action="na.fail", data=svwk2016)
summary(cox.depth) # p= 0.597
zp.depth <-cox.zph(cox.depth)
zp.depth #time varying - which, interesting. 
plt.zp.depth <-ggcoxzph(zp.depth) #definitely time varying
fit.cx.depth <-survfit(cox.depth,data.frame(list(depth=c("deep","shallow"))))
plt.cx.depth <-ggsurvplot(fit.cx.depth,svwk2016,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",  xlim=c(2,11), palette=c("blue","goldenrod"), legend.labs=c("deep","shallow")) +ggtitle("Depth")

#station
summary(svwk2016$station)
cox.station <-coxph(Surv(start,end,event2)~ station + cluster(cage), na.action="na.fail", data=svwk2016)
summary(cox.station) # all are significantly different from CPD (except CPS), but are they sig different from each other? 
zp.station <-cox.zph(cox.station)
zp.station #time varying - which, interesting. 
plt.zp.station <-ggcoxzph(zp.station) 
fit.cx.station <-survfit(cox.station,data.frame(list(station=c("CPD", "CPS", "MDD", "MDS", "RSD", "RSS"))))
plt.cx.station <-ggsurvplot(fit.cx.station,svwk2016,conf.int=TRUE,censor=FALSE, surv.median.line = "hv", xlim=c(2,11), palette=c("darkgreen","lightgreen", "firebrick4", "red","blue","lightskyblue"), legend.labs=c("CPD", "CPS", "MDD", "MDS", "RSD", "RSS"))+ggtitle("Station")
#lots of overlap with the cI

#mindo.
summary(svwk2016$min.do)
cox.mindo <-coxph(Surv(start,end,event2)~ min.do + cluster(cage), na.action="na.fail", data=svwk2016)
summary(cox.mindo) # p= 0.09
zp.mindo <-cox.zph(cox.mindo)
zp.mindo # really time varying.
plt.zp.mindo <-ggcoxzph(zp.mindo)
svwk2016$min.do_cat =as.factor(cut(svwk2016$min.do, breaks = c(0,1.35,2.91,4.86,6.68)))
#discretize the svwk2016data
cx.mindo_cat <- coxph(Surv(start,end,event2)~min.do_cat, na.action='na.fail', data=svwk2016)
summary(cx.mindo_cat) #sig diff between categories.
fit.cx.mindo_cat <-survfit(cx.mindo_cat,data.frame(list(min.do_cat=c("(0,1.35]","(1.35,2.91]", "(2.91,4.86]", "(4.86,6.68]"))))
plt.cx.mindo_cat <-ggsurvplot(fit.cx.mindo_cat,svwk2016,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("Q1:(0,1.35)","Q2:(1.35,2.91)", "Q3:(2.91,4.86)", "Q4:(4.86,6.68)")) + ggtitle("Minimum Dissolved Oxygen (mg/L)")
#lower DO, lower survival. min.do stays in. 

#biom.
summary(svwk2016$biom)
cox.biom <-coxph(Surv(start,end,event2)~ biom + cluster(cage), na.action="na.fail", data=svwk2016)
summary(cox.biom) # p= 0.0195
zp.biom <-cox.zph(cox.biom)
zp.biom # really time varying. p=0.00148, but that makes sense.
plt.zp.biom <-ggcoxzph(zp.biom)
svwk2016$biom_cat =as.factor(cut(svwk2016$biom, breaks = c(120,1454,2167,3284,8202)))
levels(svwk2016$biom_cat)
#discretize the svwk2016data
cx.biom_cat <- coxph(Surv(start,end,event2)~biom_cat, na.action='na.fail', data=svwk2016)
summary(cx.biom_cat)
fit.cx.biom_cat <-survfit(cx.biom_cat,data.frame(list(biom_cat=c("(120,1.45e+03]","(1.45e+03,2.17e+03]", "(2.17e+03,3.28e+03]", "(3.28e+03,8.2e+03]"))))
plt.cx.biom_cat <-ggsurvplot(fit.cx.biom_cat,svwk2016,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("Q1:(120 - 1453)","Q2:(1454-2169)", "Q3:(2170-3283)", "Q4:(3284-8202)"))+ ggtitle("Biomass in cage (g)")
#high biomass, high survival, but this is kind of circular. 

#maxdo. - redone 10/8/19
summary(svwk2016$max.do)
cox.max.do <-coxph(Surv(start,end,event2)~ max.do + cluster(cage), na.action="na.fail", data=svwk2016)
summary(cox.max.do) # p= 0.85
zp.max.do <-cox.zph(cox.max.do)
zp.max.do # p= 0.06 ---- no longer sig. 
plt.zp.max.do <-ggcoxzph(zp.max.do) #looks time varying
#svwk2016$max.do_cat =as.factor(cut(svwk2016$max.do, breaks = c(117,158,190,228,316)))
svwk2016$max.do_cat =as.factor(cut(svwk2016$max.do, breaks = c(8.6,11.576,14.345,15.960,23.9)))
levels(svwk2016$max.do_cat)
#discretize the svwk2016data
cx.max.do_cat <- coxph(Surv(start,end,event2)~max.do_cat, na.action='na.fail', data=svwk2016)
summary(cx.max.do_cat) # not all are significantly different from the default
fit.cx.max.do_cat <-survfit(cx.max.do_cat,data.frame(list(max.do_cat=c("(8.6,11.6]", "(11.6,14.3]", "(14.3,16]", "(16,23.9]"))))
plt.cx.max.do_cat <-ggsurvplot(fit.cx.max.do_cat,svwk2016,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("Q1:8.6-11.6)", "Q2:11.6-14.3", "Q3:14.3-16.0", "Q4:16.0-23.9"))+ ggtitle("Maximum Dissolved Oxygen (mg/L)")
#mid to low max do has better survival. eliminated from the first pass model. 

#mean.do
summary(svwk2016$mean.do)
cox.mean.do <-coxph(Surv(start,end,event2)~ mean.do + cluster(cage), na.action="na.fail", data=svwk2016)
summary(cox.mean.do) # p= 0.0262
zp.mean.do <-cox.zph(cox.mean.do)
zp.mean.do # p= 0.7.62e-12 
plt.zp.mean.do <-ggcoxzph(zp.mean.do) #looks time varying
svwk2016$mean.do_cat =as.factor(cut(svwk2016$mean.do, breaks = c(3.7,7.4,8.0,8.5,10.25)))
levels(svwk2016$mean.do_cat)
#discretize the svwk2016data
cx.mean.do_cat <- coxph(Surv(start,end,event2)~mean.do_cat, na.action='na.fail', data=svwk2016)
summary(cx.mean.do_cat)#not all are sig different from the intercept
fit.cx.mean.do_cat <-survfit(cx.mean.do_cat,data.frame(list(mean.do_cat=c("(3.7,7.4]",  "(7.4,8]",    "(8,8.5]",    "(8.5,10.2]"))))
plt.cx.mean.do_cat <-ggsurvplot(fit.cx.mean.do_cat,svwk2016,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("Q1:3.7-7.4)", "Q2:7.4-8", "Q3:8-8.5", "Q4:8.5-10.2"))+ ggtitle("Mean Dissolved Oxygen (mg/L)")
#high mean do has better survival. makes sense. 

#min.temp
summary(svwk2016$min.temp)
cox.min.temp <-coxph(Surv(start,end,event2)~ min.temp + cluster(cage), na.action="na.fail", data=svwk2016)
summary(cox.min.temp) # p= 0.0149 *
zp.min.temp <-cox.zph(cox.min.temp)
zp.min.temp 
plt.zp.min.temp <-ggcoxzph(zp.min.temp) #looks time varying
svwk2016$min.temp_cat =as.factor(cut(svwk2016$min.temp, breaks = c(12.3,17.4,18.8,21.0,24)))
levels(svwk2016$min.temp_cat)
#discretize the svwk2016data
cx.min.temp_cat <- coxph(Surv(start,end,event2)~min.temp_cat, na.action='na.fail', data=svwk2016)
summary(cx.min.temp_cat) #not all categories are different from the intercept
fit.cx.min.temp_cat <-survfit(cx.min.temp_cat,data.frame(list(min.temp_cat=c("(12.3,17.4]", "(17.4,18.8]", "(18.8,21]", "(21,24]"))))
plt.cx.min.temp_cat <-ggsurvplot(fit.cx.min.temp_cat,svwk2016,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("Q1:12.3-17.4)", "Q2:17.4-18.8", "Q3:18.8-21.0", "Q4:21.0-24.0")) +ggtitle("Minimum Temperature (Celcius)")
#not an obvious intepretation. highest minimum temperature has the LOWEST survival. lowest min temperature has the HIGHEST temperature, but intermediate categories are not in order. 

#max.temp
summary(svwk2016$max.temp)
cox.max.temp <-coxph(Surv(start,end,event2)~ max.temp + cluster(cage), na.action="na.fail", data=svwk2016)
summary(cox.max.temp) # p= 5.04e-05
zp.max.temp <-cox.zph(cox.max.temp)
zp.max.temp 
plt.zp.max.temp <-ggcoxzph(zp.max.temp) #looks time varying
svwk2016$max.temp_cat =as.factor(cut(svwk2016$max.temp, breaks = c(22.1,24.7,26.1,27.6,30.2)))
levels(svwk2016$max.temp_cat)
#discretize the svwk2016data
cx.max.temp_cat <- coxph(Surv(start,end,event2)~max.temp_cat, na.action='na.fail', data=svwk2016)
summary(cx.max.temp_cat) #only one of these is significantly different from the intercept.
fit.cx.max.temp_cat <-survfit(cx.max.temp_cat,data.frame(list(max.temp_cat=c("(22.1,24.7]", "(24.7,26.1]", "(26.1,27.6]", "(27.6,30.2]"))))
plt.cx.max.temp_cat <-ggsurvplot(fit.cx.max.temp_cat,svwk2016,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("Q1:22.1-24.6)", "Q2:24.6-26.1", "Q3:26.1-27.6", "Q4:27.6-30.2"))+ggtitle("Max Temperature (Celcius)")
#high max temp has worse survival. but the other categories ton't follow so perfectly...

#mean.temp
summary(svwk2016$mean.temp)
cox.mean.temp <-coxph(Surv(start,end,event2)~ mean.temp + cluster(cage), na.action="na.fail", data=svwk2016)
summary(cox.mean.temp) # p= 2.56e-05
zp.mean.temp <-cox.zph(cox.mean.temp)
zp.mean.temp # p= 0.074 low evidence for time varying
plt.zp.mean.temp <-ggcoxzph(zp.mean.temp) #looks slightly time varying
svwk2016$mean.temp_cat =as.factor(cut(svwk2016$mean.temp, breaks = c(17.0,20.7,22.7,23.4,27.5)))
levels(svwk2016$mean.temp_cat)
#discretize the svwk2016data
cx.mean.temp_cat <- coxph(Surv(start,end,event2)~mean.temp_cat, na.action='na.fail', data=svwk2016)
summary(cx.mean.temp_cat) #some of the categories are significantly different from the intercept. 
fit.cx.mean.temp_cat <-survfit(cx.mean.temp_cat,data.frame(list(mean.temp_cat=c("(17,20.7]",   "(20.7,22.7]", "(22.7,23.4]", "(23.4,27.5]"))))
plt.cx.mean.temp_cat <-ggsurvplot(fit.cx.mean.temp_cat,svwk2016,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("Q1:17-20.7)", "Q2:20.7-22.7", "Q3:22.7-23.4", "Q4:23.4-27.5"))+ ggtitle("Mean Temperature (Celcius)")
#high mean temp has worse survival. makes sense but a lot of overlap which makes sense. 

#do.dur
summary(svwk2016$do.dur)
cox.do.dur <-coxph(Surv(start,end,event2)~ do.dur + cluster(cage), na.action="na.fail", data=svwk2016)
summary(cox.do.dur) # p= 0.518
zp.do.dur <-cox.zph(cox.do.dur)
zp.do.dur # p= 0.0962 low evidence for time varying
plt.zp.do.dur <-ggcoxzph(zp.do.dur) #this is just a weird variable. 
svwk2016$do.dur_cat =as.factor(cut(svwk2016$do.dur, breaks = c(0,45,2955)))
levels(svwk2016$do.dur_cat)
svwk2016 <-mutate(svwk2016, do.dur_cat2=ifelse(is.na(do.dur_cat),"(0,45]",do.dur_cat)) %>% mutate(do.dur_cat=ifelse(is.na(do.dur_cat),do.dur_cat2,do.dur_cat))
svwk2016 <-mutate(svwk2016, do.dur_cat=ifelse(do.dur_cat=="(0,45]",1,do.dur_cat2)) %>% select(-do.dur_cat2)
#discretize the svwk2016data
cx.do.dur_cat <- coxph(Surv(start,end,event2)~do.dur_cat, na.action='na.fail', data=svwk2016)
summary(cx.do.dur_cat) #sig difference
fit.cx.do.dur_cat <-survfit(cx.do.dur_cat,data.frame(list(do.dur_cat=c("1", "2"))))
plt.cx.do.dur_cat <-ggsurvplot(fit.cx.do.dur_cat,svwk2016,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("Q1,Q2,Q3:0-45", "Q4:45-2955"))+ ggtitle("Duration of dissolved oxygen < 2.5 mg/l (minutes)")
# data is highly zero dominated. 

#temp.dur
summary(svwk2016$temp.dur)
cox.temp.dur <-coxph(Surv(start,end,event2)~ temp.dur + cluster(cage), na.action="na.fail", data=svwk2016)
summary(cox.temp.dur) # p= 7.71e-07
zp.temp.dur <-cox.zph(cox.temp.dur)
zp.temp.dur # p= 0.844 low evidence for time varying
plt.zp.temp.dur <-ggcoxzph(zp.temp.dur) #does not look time varying.
svwk2016$temp.dur_cat =as.factor(cut(svwk2016$temp.dur, breaks = c(0,153,8910)))
levels(svwk2016$temp.dur_cat)
svwk2016 <-mutate(svwk2016, temp.dur_cat=ifelse(is.na(temp.dur_cat),"(0,153]",temp.dur_cat)) 
svwk2016 <-mutate(svwk2016, temp.dur_cat=ifelse(temp.dur_cat=="(0,153]",1,temp.dur_cat))
#discretize the svwk2016data
cx.temp.dur_cat <- coxph(Surv(start,end,event2)~temp.dur_cat, na.action='na.fail', data=svwk2016)
summary(cx.temp.dur_cat) #sig split
fit.cx.temp.dur_cat <-survfit(cx.temp.dur_cat,data.frame(list(temp.dur_cat=c("1", "2"))))
plt.cx.temp.dur_cat <-ggsurvplot(fit.cx.temp.dur_cat,svwk2016,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("Q1,Q2,Q3:0-45", "Q4:153-8910")) + ggtitle("Duration of temperature > 27 degrees C")
#highly zero dominated data but there appears to be a split. 

#mean.sal
summary(svwk2016$mean.sal)
cox.mean.sal <-coxph(Surv(start,end,event2)~ mean.sal + cluster(cage), na.action="na.fail", data=svwk2016)
summary(cox.mean.sal) # p= 0.00149 **
zp.mean.sal <-cox.zph(cox.mean.sal)
zp.mean.sal # p= 1 no evidence for time varying
plt.zp.mean.sal <-ggcoxzph(zp.mean.sal) #looks slightly time varying, but i think overwhelmed by variability?
svwk2016$mean.sal_cat =as.factor(cut(svwk2016$mean.sal, breaks = c(22.1,27.1,29.0,29.8,31.6)))
levels(svwk2016$mean.sal_cat)
#discretize the svwk2016data
cx.mean.sal_cat <- coxph(Surv(start,end,event2)~mean.sal_cat, na.action='na.fail', data=svwk2016)
summary(cx.mean.sal_cat)#differences
fit.cx.mean.sal_cat <-survfit(cx.mean.sal_cat,data.frame(list(mean.sal_cat=c("(22.1,27.1]", "(27.1,29]","(29,29.8]",  "(29.8,31.6]"))))
plt.cx.mean.sal_cat <-ggsurvplot(fit.cx.mean.sal_cat,svwk2016,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("Q1:22.1-27.1)", "Q2:27.1-29.0", "Q3:29.0-29.8", "Q4:29.8-31.6"))+ggtitle("Mean salinity (ppt)")
#not much of a pattern here. I imagine the sites closer to the inlet etc have higer salinity. 
#high salinity, higher survival, but weirdly second highest category has low survival. 
#I think we can eliminate this as a variable because of that. 

#ssat
summary(svwk2016$ssat)
cox.ssat <-coxph(Surv(start,end,event2)~ ssat + cluster(cage), na.action="na.fail", data=svwk2016)
summary(cox.ssat) # p= 2.71e-07 ***
zp.ssat <-cox.zph(cox.ssat)
zp.ssat # p=  0.172 no evidence for time varying
plt.zp.ssat <-ggcoxzph(zp.ssat) #looks slightly time varying, but i think overwhelmed by variability?
svwk2016$ssat_cat =as.factor(cut(svwk2016$ssat, breaks = c(134,2730,4080,5430,9555)))
levels(svwk2016$ssat_cat)
#discretize the svwk2016data
cx.ssat_cat <- coxph(Surv(start,end,event2)~ssat_cat, na.action='na.fail', data=svwk2016)
summary(cx.ssat_cat)
fit.cx.ssat_cat <-survfit(cx.ssat_cat,data.frame(list(ssat_cat=c("(134,2.73e+03]","(2.73e+03,4.08e+03]", "(4.08e+03,5.43e+03]", "(5.43e+03,9.56e+03]"))))
plt.cx.ssat_cat <-ggsurvplot(fit.cx.ssat_cat,svwk2016,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("Q1:135-2730)", "Q2:2730-4080", "Q3:4080-5430", "Q4:5430-9555"))+ggtitle("Duration of dissolved oxygen saturation > 115%")
#high duration of oversaturated DO is good for survival here, so im not sure that this is a good metric. 


###### 2017 univariate models #######

#site only - from now on, cluster on cage. 
cox.site <-coxph(Surv(start, end, event2)~ site + cluster(cage), na.action="na.fail", data=svwk2017)
summary(cox.site) #not significant even at 0.2
cx.surv.est <-survfit(cox.site, newdata=data.frame(list(site=c("M","N","S"))))
# do not add the se fit because when you use cluster the SE is computed differently. 
plot(cx.surv.est, col=c("purple","orange","blue"))
zp.site <-cox.zph(cox.site)
zp.site #only S is time varying?
plt.zp.site <-ggcoxzph(zp.site) 
fit.cx.site <-survfit(cox.site,data.frame(list(site=c("M","N", "S"))))
summary(fit.cx.site)
plt.cx.site <-ggsurvplot(fit.cx.site,svwk2017,conf.int=TRUE,censor=FALSE, surv.median.line = "hv", xlim=c(2,11), palette=c("purple","orange","blue"), legend.labs=c("M","N", "S")) + ggtitle("Site")
#so much overlap. 


#depth
summary(svwk2017$depth)
cox.depth <-coxph(Surv(start,end,event2)~ depth + cluster(cage), na.action="na.fail", data=svwk2017)
summary(cox.depth) #p =0.11 *passes first fit/
zp.depth <-cox.zph(cox.depth)
zp.depth #Very time varying 0.00599 
plt.zp.depth <-ggcoxzph(zp.depth) #definitely time varying
fit.cx.depth <-survfit(cox.depth,data.frame(list(depth=c("deep","shallow"))))
plt.cx.depth <-ggsurvplot(fit.cx.depth,svwk2017,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",  xlim=c(2,11), palette=c("blue","goldenrod"), legend.labs=c("deep","shallow")) +ggtitle("Depth")
#looks different but not significant. 

#mindo.
summary(svwk2017$min.do)
cox.mindo <-coxph(Surv(start,end,event2)~ min.do + cluster(cage), na.action="na.fail", data=svwk2017)
summary(cox.mindo) # p= 0.00944 very significant?
zp.mindo <-cox.zph(cox.mindo)
zp.mindo # not time varying 0.377
plt.zp.mindo <-ggcoxzph(zp.mindo)
svwk2017$min.do_cat =as.factor(cut(svwk2017$min.do, breaks = c(1.0,3.390,4.330,5.280,8.281)))
#discretize the svwk2017data
cx.mindo_cat <- coxph(Surv(start,end,event2)~min.do_cat, data=svwk2017)
summary(cx.mindo_cat) #not sig
fit.cx.mindo_cat <-survfit(cx.mindo_cat,data.frame(list(min.do_cat=c("(1,3.39]","(3.39,4.33]","(4.33,5.28]","(5.28,8.28]"))))
plt.cx.mindo_cat <-ggsurvplot(fit.cx.mindo_cat,svwk2017,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("Q1:1.0-3.93","Q2:3.93-4.33)", "Q3:4.33-5.28", "Q4:5.28-8.28")) + ggtitle("Minimum Dissolved Oxygen (mg/L)")
#trends make sense now!

#maxdo.
summary(svwk2017$max.do)
cox.max.do <-coxph(Surv(start,end,event2)~ max.do + cluster(cage), na.action="na.fail", data=svwk2017)
summary(cox.max.do) # p= 0.517 - do not include. 
zp.max.do <-cox.zph(cox.max.do)
zp.max.do # p= 0.105
plt.zp.max.do <-ggcoxzph(zp.max.do) #looks time varying but not significant. 
svwk2017$max.do_cat =as.factor(cut(svwk2017$max.do, breaks = c(8.7,10.47,11.87,13.65,20.0)))
levels(svwk2017$max.do_cat)
#discretize the svwk2017data
cx.max.do_cat <- coxph(Surv(start,end,event2)~max.do_cat, na.action='na.fail', data=svwk2017)
summary(cx.max.do_cat) # none are significantly different from the default
fit.cx.max.do_cat <-survfit(cx.max.do_cat,data.frame(list(max.do_cat=c("(8.7,10.5]", "(10.5,11.9]", "(11.9,13.7]", "(13.7,20]"))))
plt.cx.max.do_cat <-ggsurvplot(fit.cx.max.do_cat,svwk2017, conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("Q1:8.7-10.5)", "Q2:10.5-11.9", "Q3:11.9-13.7", "Q4:13.7-20"))+ ggtitle("Maximum Dissolved Oxygen (mg/L)")
#these values for Max.DO don't even make sense for DOmgl. they must be in % saturation. 
#total overlap. but sort of makes sense. 

#mean.do
summary(svwk2017$mean.do)
cox.mean.do <-coxph(Surv(start,end,event2)~ mean.do + cluster(cage), na.action="na.fail", data=svwk2017)
summary(cox.mean.do) # p= 0.84
zp.mean.do <-cox.zph(cox.mean.do)
zp.mean.do # p= 0.543
plt.zp.mean.do <-ggcoxzph(zp.mean.do) #looks time varying, but not significant.
svwk2017$mean.do_cat =as.factor(cut(svwk2017$mean.do, breaks = c(6.24,6.29,7.49,8.40,13.110)))
levels(svwk2017$mean.do_cat)
#discretize the svwk2017data
cx.mean.do_cat <- coxph(Surv(start,end,event2)~mean.do_cat, na.action='na.fail', data=svwk2017)
summary(cx.mean.do_cat)#none are sig different from the intercept
fit.cx.mean.do_cat <-survfit(cx.mean.do_cat,data.frame(list(mean.do_cat=c("(6.24,6.29]", "(6.29,7.49]", "(7.49,8.4]",  "(8.4,13.1]"))))
plt.cx.mean.do_cat <-ggsurvplot(fit.cx.mean.do_cat,svwk2017,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("Q1:6.24-6.29", "Q2:6.29-7.49", "Q3:7.49-8.40", "Q4:8.40-13.11"))+ ggtitle("Mean Dissolved Oxygen (mg/L)")
# i think something went wrong here. I looked at the data and most fish in these categories had high death weeks. 
# regardless, we're not including it. 

#min.temp
summary(svwk2017$min.temp)
cox.min.temp <-coxph(Surv(start,end,event2)~ min.temp + cluster(cage), na.action="na.fail", data=svwk2017)
summary(cox.min.temp) # p= 0.784 not sig. 
zp.min.temp <-cox.zph(cox.min.temp)
zp.min.temp # p= 0.886
plt.zp.min.temp <-ggcoxzph(zp.min.temp) #not time varying. 
svwk2017$min.temp_cat =as.factor(cut(svwk2017$min.temp, breaks = c(13.0,16.64,19.56,20.32,22.03)))
levels(svwk2017$min.temp_cat)
#discretize the svwk2017data
cx.min.temp_cat <- coxph(Surv(start,end,event2)~min.temp_cat, na.action='na.fail', data=svwk2017)
summary(cx.min.temp_cat) #no categories are different from the intercept
fit.cx.min.temp_cat <-survfit(cx.min.temp_cat,data.frame(list(min.temp_cat=c("(13,16.6]",   "(16.6,19.6]", "(19.6,20.3]", "(20.3,22]"))))
plt.cx.min.temp_cat <-ggsurvplot(fit.cx.min.temp_cat,svwk2017,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("Q1:13.0-16.6", "Q2:16.6-19.6", "Q3:19.6-20.3", "Q4:20.3-22.0")) +ggtitle("Minimum Temperature (Celcius)")
# definitely makes a lot of sense. 

#max.temp
summary(svwk2017$max.temp)
cox.max.temp <-coxph(Surv(start,end,event2)~ max.temp + cluster(cage), na.action="na.fail", data=svwk2017)
summary(cox.max.temp) # p= 0.569
zp.max.temp <-cox.zph(cox.max.temp)
zp.max.temp # p= 0.203
plt.zp.max.temp <-ggcoxzph(zp.max.temp) #looks time varying but not significant
svwk2017$max.temp_cat =as.factor(cut(svwk2017$max.temp, breaks = c(20.52,23.60,24.54,25.54,30.65)))
levels(svwk2017$max.temp_cat)
#discretize the svwk2017data
cx.max.temp_cat <- coxph(Surv(start,end,event2)~max.temp_cat, na.action='na.fail', data=svwk2017)
summary(cx.max.temp_cat) #none are significantly different from the intercept.
fit.cx.max.temp_cat <-survfit(cx.max.temp_cat,data.frame(list(max.temp_cat=c("(20.5,23.6]", "(23.6,24.5]", "(24.5,25.5]", "(25.5,30.6]"))))
plt.cx.max.temp_cat <-ggsurvplot(fit.cx.max.temp_cat,svwk2017,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("Q1:20.5-23.6", "Q2:23.6-24.5", "Q3:24.5-25.5", "Q4:25.5-30.6"))+ggtitle("Max Temperature (Celcius)")
#highest max temp has lowest survival, so that makes sense. but otherwise trends don't make a ton of sense. 

#mean.temp
summary(svwk2017$mean.temp)
cox.mean.temp <-coxph(Surv(start,end,event2)~ mean.temp + cluster(cage), na.action="na.fail", data=svwk2017)
summary(cox.mean.temp) # p= 0.485
zp.mean.temp <-cox.zph(cox.mean.temp)
zp.mean.temp # p= 0.234
plt.zp.mean.temp <-ggcoxzph(zp.mean.temp) #not time varying
svwk2017$mean.temp_cat =as.factor(cut(svwk2017$mean.temp, breaks = c(17.3,19.33,21.5,22.60,23.71)))
levels(svwk2017$mean.temp_cat)
#discretize the svwk2017data
cx.mean.temp_cat <- coxph(Surv(start,end,event2)~mean.temp_cat, na.action='na.fail', data=svwk2017)
summary(cx.mean.temp_cat) #middle category is significantly different 
fit.cx.mean.temp_cat <-survfit(cx.mean.temp_cat,data.frame(list(mean.temp_cat=c("(17.3,19.3]", "(19.3,21.5]", "(21.5,22.6]", "(22.6,23.7]"))))
plt.cx.mean.temp_cat <-ggsurvplot(fit.cx.mean.temp_cat,svwk2017,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("Q1:17.3-19.3", "Q2:19.3-21.5", "Q3:21.5-22.6", "Q4:22.6-23.7"))+ ggtitle("Mean Temperature (Celcius)")
# lowest mean temperature category looks very different from the others which makes sense. 

#do.dur
summary(svwk2017$do.dur) #because of the zero dominance, I will instead split the data 
cox.do.dur <-coxph(Surv(start,end,event2)~ do.dur + cluster(cage), na.action="na.fail", data=svwk2017) #wont properly converge because too zero dominated. maybe I need to scale it? 
summary(cox.do.dur) # p = 0.849
zp.do.dur <-cox.zph(cox.do.dur)
zp.do.dur # p= 0.941 
plt.zp.do.dur <-ggcoxzph(zp.do.dur) #this is just a weird variable. 
svwk2017$do.dur_cat =as.factor(cut(svwk2017$do.dur, breaks =c(0,1,195)))
levels(svwk2017$do.dur_cat)
svwk2017 <-mutate(svwk2017, do.dur_cat=ifelse(is.na(do.dur_cat),"(0,1]",do.dur_cat))# %>% mutate(do.dur_cat=ifelse(is.na(do.dur_cat),do.dur_cat2,do.dur_cat))
svwk2017 <-mutate(svwk2017, do.dur_cat=ifelse(do.dur_cat==2,"(1,195]",do.dur_cat)) 
#discretize the svwk2017data
cx.do.dur_cat <- coxph(Surv(start,end,event2)~do.dur_cat, na.action='na.fail', data=svwk2017) #won't converge
summary(cx.do.dur_cat) #0.79  
fit.cx.do.dur_cat <-survfit(cx.do.dur_cat,data.frame(list(do.dur_cat=c("(0,1]" ,  "(1,195]" ))))
plt.cx.do.dur_cat <-ggsurvplot(fit.cx.do.dur_cat,svwk2017,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("0", "> 0" ))+ ggtitle("Duration of dissolved oxygen < 2.0 mg/l (minutes)")
# data is highly zero dominated.
# this is reflected in the fact that there is 100% survival supposedly for those above 7.11, which brings back the spectre that the model doesn't understand the time structure. 
#do.dur, BUT SCALED has the same problem. it's too zero dominated. 
# maybe you should pick another threshold instead of two. 
filter(svwk2017, svwk2017$do.dur > 0) # only 75 fish-observations where do.dur > 0

#temp.dur
summary(svwk2017$temp.dur) # we're splitting by zero and not zero. 
cox.temp.dur <-coxph(Surv(start,end,event2)~ temp.dur + cluster(cage), na.action="na.fail", data=svwk2017)
summary(cox.temp.dur) # p= 0.0382 * significant
zp.temp.dur <-cox.zph(cox.temp.dur)
zp.temp.dur # p= 0.88 not time varying
plt.zp.temp.dur <-ggcoxzph(zp.temp.dur) #does not look time varying.
svwk2017$temp.dur_cat =as.factor(cut(svwk2017$temp.dur, breaks = c(0,1,195)))
levels(svwk2017$temp.dur_cat)
svwk2017 <-mutate(svwk2017, temp.dur_cat=ifelse(is.na(temp.dur_cat),"(0,1]",temp.dur_cat)) 
svwk2017 <-mutate(svwk2017, temp.dur_cat=ifelse(temp.dur_cat=="(0,1]",1,temp.dur_cat))
#discretize the svwk2017data
cx.temp.dur_cat <- coxph(Surv(start,end,event2)~temp.dur_cat, na.action='na.fail', data=svwk2017)
summary(cx.temp.dur_cat) #0.994
fit.cx.temp.dur_cat <-survfit(cx.temp.dur_cat,data.frame(list(temp.dur_cat=c("1", "2"))))
plt.cx.temp.dur_cat <-ggsurvplot(fit.cx.temp.dur_cat,svwk2017,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("0", ">0")) + ggtitle("Duration of temperature > 27 degrees C")
#highly zero dominated data and it appears like there's no difference. 

#mean.sal
summary(svwk2017$mean.sal)
cox.mean.sal <-coxph(Surv(start,end,event2)~ mean.sal + cluster(cage), na.action="na.fail", data=svwk2017)
summary(cox.mean.sal) # p= 0.193 #makes it in 
zp.mean.sal <-cox.zph(cox.mean.sal)
zp.mean.sal # p= 0.168 no evidence for time varying
plt.zp.mean.sal <-ggcoxzph(zp.mean.sal) #looks slightly time varying, but i think overwhelmed by variability?
svwk2017$mean.sal_cat =as.factor(cut(svwk2017$mean.sal, breaks = c(25.90,28.02,30.08,30.97,34.10)))
levels(svwk2017$mean.sal_cat)
#discretize the svwk2017data
cx.mean.sal_cat <- coxph(Surv(start,end,event2)~mean.sal_cat, na.action='na.fail', data=svwk2017)
summary(cx.mean.sal_cat)#differences!
fit.cx.mean.sal_cat <-survfit(cx.mean.sal_cat,data.frame(list(mean.sal_cat=c("(25.9,28]", "(28,30.1]", "(30.1,31]", "(31,34.1]"))))
plt.cx.mean.sal_cat <-ggsurvplot(fit.cx.mean.sal_cat,svwk2017,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("Q1:25.9-28.0", "Q2:28.0-30.1", "Q3:30.1-31.0", "Q4:31.0-34.1"))+ggtitle("Mean salinity (ppt)")
#lowest salnity has highest survival ok, highest salinity has lowest survival, but then the middle categories are swapped. I think this is like a depth effect. 

#ssat
summary(svwk2017$ssat)
cox.ssat <-coxph(Surv(start,end,event2)~ ssat + cluster(cage), na.action="na.fail", data=svwk2017)
summary(cox.ssat) # p= 0.86
zp.ssat <-cox.zph(cox.ssat)
zp.ssat # p=  0.0712  time varying
plt.zp.ssat <-ggcoxzph(zp.ssat) #looks slightly time varying, but i think overwhelmed by variability?
svwk2017$ssat_cat =as.factor(cut(svwk2017$ssat, breaks = c(0,525,1785,3405,11296)))
levels(svwk2017$ssat_cat)
svwk2017 <-mutate(svwk2017, ssat_cat=ifelse(is.na(ssat_cat),"(0,525]",ssat_cat)) 
svwk2017 <-mutate(svwk2017, ssat_cat=ifelse(ssat_cat=="(0,525]",1,ssat_cat)) 
#discretize the svwk2017data
cx.ssat_cat <- coxph(Surv(start,end,event2)~ssat_cat, na.action='na.fail', data=svwk2017)
summary(cx.ssat_cat) #one is signficantly different cat 2 
fit.cx.ssat_cat <-survfit(cx.ssat_cat,data.frame(list(ssat_cat=c("1","2","3", "4"))))
plt.cx.ssat_cat <-ggsurvplot(fit.cx.ssat_cat,svwk2017,conf.int=TRUE,censor=FALSE, surv.median.line = "hv",xlim=c(2,11), palette = "jco", legend.labs=c("Q1:0-525", "Q2:525-1785", "Q3:1785-3405", "Q4:3405-11296"))+ggtitle("Duration of dissolved oxygen saturation > 115%")
#low duration of oversaturation is good, but the second best is the highest duration. I think this is all confounded by site. 

########################

######### Model Selection ##########
#http://www.karlin.mff.cuni.cz/~pesta/NMFM404/ph.html#Model_selection
#I originally followed this model selection tutorial. However, this would result in different starting models for 2017 and 2016
# if we were to follow the tutorial and start by removing stuff based on a 0.25 pvalue cutoff from the univariate model
#for 2016, we would eliminate depth, max.do, and do.dur and start the full model from there. 
#for 2017, we would ONLY KEEP min.do, depth, temp.dur, and mean.salinity
# i think there's an argument that these years are fundamentally different because different location, but w/e

#VIF doesn't work for cox models because they lack an intercept. But a correlation analysis can be performed. 
#correlation plot
corr16 <-select(svwk2016,min.do,max.do,mean.do,min.temp,max.temp,mean.temp,temp.dur,do.dur,mean.sal,ssat)
M16 <-as.data.frame(cor(corr16))
colnames(M16) <-c("min.do16","max.do16","mean.do16","min.temp16","max.temp16","mean.temp16","temp.dur16","do.dur16","mean.sal16","ssat16")
#if we were to eliminate based on correlation in 2016 (cutoff = 0.7), we'd get rid of one of each of the following pairs:
## ssat & mean.do, mean.temp & min.temp, max.temp & mean.temp, 
cor16 <-coxph(Surv(start, end, event2)~ site + depth + min.do + max.do + mean.do + min.temp + max.temp  + do.dur + temp.dur + mean.sal + cluster(cage),na.action="na.fail", data=svwk2016)
dredge(cor16)
dcor16 <-coxph(Surv(start, end, event2)~ site + min.temp + cluster(cage),na.action="na.fail", data=svwk2016)
dcor16.2 <-coxph(Surv(start, end, event2)~ site + max.temp + cluster(cage),na.action="na.fail", data=svwk2016)

corr17 <-select(svwk2017,min.do,max.do,mean.do,min.temp,max.temp,mean.temp,temp.dur,do.dur,mean.sal,ssat)
M17 <-as.data.frame(cor(corr17))
colnames(M17) <- c("min.do17","max.do17","mean.do17","min.temp17","max.temp17","mean.temp17","temp.dur17","do.dur17","mean.sal17","ssat17")
#if we were to eliminate based on correlation in 2017 (cutoff = 0.7), we'd get rid of one of each of the following pairs:
## mean.do & max.do, ssat & max.do, ssat & mean.do, max.temp & min.temp, min.temp & mean.temp, max.temp & mean.temp, 
cor17 <-coxph(Surv(start, end, event2)~ site + depth + min.do  + mean.do + min.temp + max.temp  + do.dur + temp.dur + mean.sal +  cluster(cage),na.action="na.fail", data=svwk2017)
dredge(cor17)
dcor17 <- coxph(Surv(start, end, event2)~ depth + max.temp + mean.sal + temp.dur + cluster(cage),na.action="na.fail", data=svwk2017)
dcor17.2 <- coxph(Surv(start, end, event2)~ depth + min.temp + mean.sal + temp.dur + cluster(cage),na.action="na.fail", data=svwk2017)

#create the composite correlation matrix. 
#in this version we create a composite matrix that does correlation for each year separately then averages them. This seems better...
corav <-bind_cols(M16,M17) %>% transmute(min.do = (min.do16+min.do17)/2, max.do=(max.do16+max.do17)/2,mean.do=(mean.do16+mean.do17)/2,min.temp=(min.temp16+min.temp17)/2,max.temp=(max.temp16+max.temp17)/2,mean.temp=(mean.temp16+mean.temp17)/2,temp.dur=(temp.dur16+temp.dur17)/2,do.dur=(do.dur16+do.dur17)/2,mean.sal=(mean.sal16+mean.sal17)/2,ssat=(ssat16+ssat17)/2)
rownames(corav) <-c("min.do","max.do","mean.do","min.temp","max.temp","mean.temp","temp.dur","do.dur","mean.sal","ssat")
# the resulting pairs are: mean.do & ssat, min.temp & mean.temp, max.temp & mean.temp, 
# in choosing between eliminating max.temp vs. mean temp, max temp has more of a biological explanation and was more significant in the univariate models. 


#ssat,  min.temp, mean.temp  would do to remove based on both year correlations.
corboth16 <-coxph(Surv(start, end, event2)~site + depth + min.do + max.do + mean.do  + max.temp + do.dur + temp.dur + mean.sal  + cluster(cage),na.action="na.fail", data=svwk2016)
dredge(corboth16)
dcb16 <-coxph(Surv(start, end, event2)~  max.temp + cluster(cage),na.action="na.fail", data=svwk2016)
dcb16.2 <-coxph(Surv(start, end, event2)~  site + mean.sal + cluster(cage),na.action="na.fail", data=svwk2016)
corboth17 <-coxph(Surv(start, end, event2)~site + depth + min.do + max.do + mean.do  + max.temp + do.dur + temp.dur + mean.sal  + cluster(cage),na.action="na.fail", data=svwk2017)
dredge(corboth17)
dcb17 <-coxph(Surv(start, end, event2)~ depth + max.temp + mean.sal + temp.dur + cluster(cage),na.action="na.fail", data=svwk2017)
dcb17.2 <-coxph(Surv(start, end, event2)~ depth + max.temp + mean.sal + temp.dur + do.dur + cluster(cage),na.action="na.fail", data=svwk2017)


## previous best models for 2016 (from script surv_16_cox_10_9_19.R)
cox16.1 <-coxph(Surv(start, end, event2)~ site*mean.sal + mean.temp + cluster(cage),na.action="na.fail", data=svwk2016)
cox16.2 <-coxph(Surv(start, end, event2)~ site + mean.temp+ cluster(cage),na.action="na.fail", data=svwk2016)
cox16.3 <-coxph(Surv(start, end, event2)~ site + mean.temp + ssat + cluster(cage),na.action="na.fail", data=svwk2016)

cox16_all  <-coxph(Surv(start, end, event2)~ site + depth + min.do + max.do + mean.do + min.temp + max.temp + mean.temp + do.dur + temp.dur + mean.sal + ssat + cluster(cage),na.action="na.fail", data=svwk2016)
#d1 <-dredge(cox16_all)
dwin16cox <-coxph(Surv(start, end, event2)~ site + min.temp + cluster(cage),na.action="na.fail", data=svwk2016)
dwin16cox2 <-coxph(Surv(start, end, event2)~ site + min.temp + ssat + cluster(cage),na.action="na.fail", data=svwk2016)
#step <-stepAIC(cox16_all, direction="both", k=2)
# also steps to site + min temp, 2nd best is site + min.temp +ssat.

#all model 2017
cox17_all  <-coxph(Surv(start, end, event2)~ site + depth + min.do + max.do + mean.do + min.temp + max.temp + mean.temp + do.dur + temp.dur + mean.sal + ssat + cluster(cage),na.action="na.fail", data=svwk2017)
#d1 <-dredge(cox17_all)
dwincox17 <-coxph(Surv(start, end, event2)~depth + mean.sal + mean.temp + temp.dur + cluster(cage), data=svwk2017)
dwincox17.2 <-coxph(Surv(start, end, event2)~ depth + mean.do + mean.temp +  temp.dur + mean.sal + cluster(cage), data=svwk2017)
#step <-stepAIC(cox17_all, direction ="both",k=2)
step17cox <-coxph(Surv(start, end, event2)~depth + min.do + max.do + mean.do + min.temp + mean.temp + do.dur + temp.dur + mean.sal + cluster(cage), data=svwk2017)
# previous best model for 2017 (from script surv_17_cox_1_11_19.R)
cox_prev17 <-coxph(Surv(start, end, event2)~depth*mean.temp+temp.dur+mean.sal + cluster(cage), data=svwk2017)

  
faceoff16 <-as.data.frame(model.sel(dcor16,dcor16.2,dcb16,dcb16.2,cox16.1,cox16.2,cox16.3,dwin16cox,dwin16cox2))
faceoff17 <-as.data.frame(model.sel(dcor17,dcor17.2, dcb17,dcb17.2, dwincox17,dwincox17.2,step17cox,cox_prev17))

#write.csv(faceoff16,file="faceoff16.csv")
#write.csv(faceoff17,file="faceoff17.csv")

#Well that's fun.
# the model selection process is further detailed in the spreadsheet saved calculations_1_15_19
# this is only the model selection script, the graphics script will be called survgraphs10_17_19.R
### interpretation of coefficients http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Survival/BS704_Survival6.html
############

#2016 best models
dcor16 <-coxph(Surv(start, end, event2)~ site + min.temp + cluster(cage),na.action="na.fail", data=svwk2016)
cox16.1 <-coxph(Surv(start, end, event2)~ site*mean.sal + mean.temp + cluster(cage),na.action="na.fail", data=svwk2016)
dcor16.2 <-coxph(Surv(start, end, event2)~ site + max.temp + cluster(cage),na.action="na.fail", data=svwk2016)
win16 <- model.sel(dcor16,cox16.1,dcor16.2)

#2017 best models
cox_prev17 <-coxph(Surv(start, end, event2)~depth*mean.temp+temp.dur+mean.sal + cluster(cage),na.action ="na.fail", data=svwk2017)
dwincox17 <-coxph(Surv(start, end, event2)~depth + mean.sal + mean.temp + temp.dur + cluster(cage),na.action = "na.fail", data=svwk2017)
dwincox17.2 <-coxph(Surv(start, end, event2)~ depth + mean.do + mean.temp +  temp.dur + mean.sal + cluster(cage),na.action="na.fail", data=svwk2017)
win17 <-model.sel(cox_prev17,dwincox17,dwincox17.2)


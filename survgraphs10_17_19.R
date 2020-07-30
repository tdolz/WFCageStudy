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

setwd("/Users//tdolan/Documents//WIP research/Caging paper/data/heather meeting/caging_rproj")


svwk2016 <-read.csv("weeklysummarydata_2016_10-8-19.csv", header=TRUE) # the weekly summary data, covariates are unscaled
svwk2017 <-read.csv("weeklysummarydata_2017.csv", header=TRUE) # the weekly summary data, covariates are unscaled


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

###################graphics###################


###correlation plot ####
#change year as needed. 
csv <-dplyr::select(svwk2016,min.do,max.do,mean.do,min.temp,max.temp,mean.temp,temp.dur,do.dur,mean.sal,ssat)
library("corrplot")
M <-cor(csv)  # produces the correlation matrix. 
#These colnames are wrong depending on what year you use, so fix/pay attention
#colnames(M) <- c("min DO", "max DO", "min temp", "max temp", "mean temp", "minutes :> 27C","minutes :< 2.0 mg/L","mean salinity", "minutes :> 115:% DO" )
#rownames(M) <- c("min DO", "max DO", "min temp", "max temp", "mean temp", "minutes :> 27C","minutes :< 2.0 mg/L","mean salinity", "minutes :> 115:% DO" )
corrplot.mixed(M,method="circle", lower.col="black", number.cex=.7)
res1 <- cor.mtest(csv, conf.level = .95)
p.mat = res1$p

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M, method = "color", col = col(200),
         type = "upper", order = "hclust", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = FALSE)
######
######## Residuals Plots ###########
#ok I think we report the martingale residuals. 
#https://rpkgs.datanovia.com/survminer/survminer_cheatsheet.pdf

#remember to change the model and dataset you are calling
martres = as.vector(residuals(cox_prev17, type="martingale"))
martres2 <- cbind(svwk2017,martres)
ggplot(martres2,aes(x=site,y=martres))+ geom_violin()+ geom_boxplot(width=0.1)+ scale_fill_grey() + labs(x="site", y = "Score Residuals")+ theme_classic()

plot(martres2$min.temp, martres, xlab="Minimum Temperature", ylab="Martingale Residuals", col=2)
abline(h=0, lty=2)

#(1) Schoenfeld residual for checking the proportional hazards assumption for a covariate, 
#(2)Martingale residual used to examine overall test of the goodness-of-fit of a Cox model, 
#(3) Deviance residual for detection of poorly predicted observations, and
#(4) Score residual for determination of influential observations. 
fit0 <-coxph(Surv(start,end,event2)~1,data=svwk2017)

ggcoxdiagnostics(cox_prev17,type="martingale",ox.scale="linear.predictions") #observed # of deaths vs. predicted deaths for subject i between time 0 and Xi
ggcoxdiagnostics(dcor16,type="deviance",ox.scale="linear.predictions") # a derived martingale that is centered on 0. 
ggcoxdiagnostics(cox_prev17,type="score",ox.scale="linear.predictions") 
ggcoxdiagnostics(dcor16,type="schoenfeld",ox.scale="linear.predictions")
ggcoxdiagnostics(dcor16,type="scaledsch")

#if you calculate martingale or deviance residuals without any covariates and then plot them against covariates, you obtain a graphical relationship between the covariate and the hazard
mart0 <-residuals(fit0,type="deviance")
plot(svwk2016$min.temp, mart0, xlab="Minimum Temperature", ylab="Deviance Residuals")
lines(lowess(svwk2016$min.temp, mart0), col=2)

######
#### Forest plot ####

#have to make a no cluster version because it won't take a cluster. 
dcor16noclust <-coxph(Surv(start,end,event2)~site + min.temp, data = svwk2016, na.action = "na.fail")
ggforest(dcor16noclust, data=svwk2016)

cox_prev17noclust <-coxph(Surv(start,end,event2)~ depth * mean.temp + temp.dur + mean.sal, data = svwk2017, na.action = "na.fail")
ggforest(cox_prev17noclust, data=svwk2017)
############### 
##### Hazard, survivors, events plots ######

#create the ggsurv object. Change for each year. 
fitnew <-survfit(dcor16, data=svwk2016)

#hazard function
ggsurvplot(fitnew, fun="cumhaz", conf.int = TRUE, palette = "Dark2", 
           censor = FALSE, size=1, risk.table=FALSE, xlim=c(2,11), break.x.by=1)

#cumulative events
ggsurvplot(fitnew, fun="event", conf.int = TRUE, palette = "Dark2", 
           censor = FALSE, size=1, risk.table=FALSE, xlim=c(2,11), break.x.by=1)

#survivors plot
ggsurvplot(fitnew, conf.int = TRUE, palette = "Dark2", 
           censor = FALSE, fun="pct", size=1, surv.median.line = "hv", risk.table=TRUE, xlim=c(2,11), break.x.by=1)
############################
##### marginal effects ########

#2016#
#the effects for dcor16 are site, min.temp
av.mintemp <- mean(svwk2016$min.temp)

# looking at marginal effect of site, holding minimum temperature to it's average for the year. 
#average minimum temperature across whole year is 18.84001
fitnew2 <-survfit(dcor16, newdata=data.frame(site=c("CP","MD","RS"),min.temp=rep(av.mintemp),3), data=svwk2016)
ggsurvplot(fitnew2, conf.int = TRUE, palette = "Dark2", 
           censor = FALSE, fun="pct", size=1, surv.median.line = "hv", risk.table=FALSE,  legend.labs=c("CP","MD","RS"), xlim=c(2,11), break.x.by=1)
summary(fitnew2) # just printing the summary is the risk table
# The pval comes from the log rank test, so perform the log rank test to extrac the pval. The logrank test for the fitb model is p =3e-12.  


#Now these are the prediction graphs, we're looking at what would happen if minimum temperature was really high vs. really low. 
min.mintemp <-min(svwk2016$min.temp) #12.31
max.mintemp <-max(svwk2016$min.temp) #23.93 that is a high minimum temperature! 

#marginal effects with site.
fithigh <-survfit(dcor16, newdata=data.frame(site=c("CP","MD","RS"),min.temp=rep(max.mintemp),3), data=svwk2016)
fitlow <-survfit(dcor16, newdata=data.frame(site=c("CP","MD","RS"),min.temp=rep(min.mintemp),3), data=svwk2016)

#try without separating out site. 
#marginal effects with site.
fithigh <-survfit(dcor16, newdata=data.frame(site=c("CP","MD","RS"),min.temp=rep(max.mintemp),3), data=svwk2016)
fitlow <-survfit(dcor16, newdata=data.frame(site=c("CP","MD","RS"),min.temp=rep(min.mintemp),3), data=svwk2016)


#mycols <- c("#33a02c","#e31a1c","#1f78b4","#b2df8a","#fb9a99","#a6cee3")
mycols <- c("#8c510a", "#d8b365","#f6e8c3","#c7eae5","#5ab4ac","#01665e")
fit <- list(Hightemp = fithigh, Lowtemp = fitlow)
ggsurvplot(fit, data = svwk2016, combine = TRUE, # Combine curves
           risk.table = FALSE,                  # Add risk table
           conf.int = FALSE,                   # Add confidence interval
           conf.int.style = "ribbon",            # CI style, use "step" or "ribbon"
           censor = FALSE,                     # Remove censor points
           legend.labs=c("CP-high","MD-high","RS-high", "CP-low", "MD-low", "RS-low"),
           tables.theme = theme_cleantable(),  # Clean risk table
           xlim=c(2,11),
           break.x.by=1,
           palette = mycols)
#to me this plot suggests a significant site/min.temp interaction... 

#extract fitted values from the cox model, once you have decided.
cb <-prediction(dcor16, type="expected") %>% mutate(survprob=exp(-fitted)) %>% mutate(inst = 1-exp((1-survprob)/7))
cbgg <- ggplot(cb, aes(x=week, y=survprob, color=station))
cbgg + geom_line()+theme_bw()



##2017##
#the effects in 2017 are depth * mean.temp + temp.dur + mean.sal 
av.meantemp17 <-mean(svwk2017$mean.temp) #21.07824
av.tempdur17 <- mean(svwk2017$temp.dur)#8.764988
av.meansal17 <- mean(svwk2017$mean.sal)#29.8736

# looking at marginal effect of depth, holding salinity & temp & tempdur to mean. 
fitnew2 <-survfit(cox_prev17, newdata=data.frame(depth=c("deep","shallow"),mean.temp=rep(av.meantemp17),mean.sal=rep(av.meansal17),temp.dur=rep(av.tempdur17),2), data=svwk2017)
ggsurvplot(fitnew2, conf.int = TRUE, palette = "Dark2", 
           censor = FALSE, fun="pct", size=1, surv.median.line = "hv", risk.table=FALSE,  legend.labs=c("deep","shallow"), xlim=c(2,11), break.x.by=1)
summary(fitnew2) # just printing the summary is the risk table
# The pval comes from the log rank test, so perform the log rank test to extrac the pval. The logrank test for the fitb model is p =3e-12.  



#looking at marginal effect of mean.temp, holding everything else constant. 
#since temperature data changes day to day idk if this is appropriate.
max.meantemp17 <-max(svwk2017$mean.temp)#23.69685
min.meantemp17 <-min(svwk2017$mean.temp)#17.32684
fitmeanhigh <-survfit(cox_prev17, newdata=data.frame(depth=c("deep","shallow"),mean.temp=rep(max.meantemp17),mean.sal=rep(av.meansal17),temp.dur=rep(av.tempdur17),2), data=svwk2017)
fitmeanlow <-survfit(cox_prev17, newdata=data.frame(depth=c("deep","shallow"),mean.temp=rep(min.meantemp17),mean.sal=rep(av.meansal17),temp.dur=rep(av.tempdur17),2), data=svwk2017)

mycols <- c("#a6611a","#dfc27d","#018571","#80cdc1")

fit <- list(Hightemp = fitmeanhigh, Lowtemp = fitmeanlow)
ggsurvplot(fit, data = svwk2017, combine = TRUE, # Combine curves
           #risk.table = TRUE,                  # Add risk table
           #conf.int = TRUE,                    # Add confidence interval
           #conf.int.style = "step",            # CI style, use "step" or "ribbon"
           censor = FALSE,                     # Remove censor points
           legend.labs=c("deep-high","shallow-high","deep-low", "shallow-low"),
           tables.theme = theme_cleantable(),  # Clean risk table
           xlim=c(2,11),
           break.x.by=1,
           palette = mycols)

#looking at marginal effect of mean.sal, holding everything else constant. 
maxmeansal17 <-max(svwk2017$mean.sal)#34.09135
minmeansal17 <-min(svwk2017$mean.sal)#25.90873
fitmeanhigh <-survfit(cox_prev17, newdata=data.frame(depth=c("deep","shallow"),mean.temp=rep(av.meantemp17),mean.sal=rep(maxmeansal17),temp.dur=rep(av.tempdur17),2), data=svwk2017)
fitmeanlow <-survfit(cox_prev17, newdata=data.frame(depth=c("deep","shallow"),mean.temp=rep(av.meantemp17),mean.sal=rep(minmeansal17),temp.dur=rep(av.tempdur17),2), data=svwk2017)

mycols <- c("#a6611a","#dfc27d","#018571","#80cdc1")
fit <- list(Highsal = fitmeanhigh, Lowsal = fitmeanlow)
ggsurvplot(fit, data = svwk2017, combine = TRUE, # Combine curves
           #risk.table = TRUE,                  # Add risk table
           #conf.int = TRUE,                    # Add confidence interval
           #conf.int.style = "step",            # CI style, use "step" or "ribbon"
           censor = FALSE,                     # Remove censor points
           legend.labs=c("deep-high","shallow-high","deep-low", "shallow-low"),
           tables.theme = theme_cleantable(),  # Clean risk table
           xlim=c(2,11),
           break.x.by=1,
           palette = mycols)

#looking at marginal effect of temp.dur, holding everything else constant. 
maxtempdur <-max(svwk2017$temp.dur)#195, min will obviously be zero. 
fitmeanhigh <-survfit(cox_prev17, newdata=data.frame(depth=c("deep","shallow"),mean.temp=rep(av.meantemp17),mean.sal=rep(av.meansal17),temp.dur=rep(195),2), data=svwk2017)
fitmeanlow <-survfit(cox_prev17, newdata=data.frame(depth=c("deep","shallow"),mean.temp=rep(av.meantemp17),mean.sal=rep(av.meansal17),temp.dur=rep(0),2), data=svwk2017)

mycols <- c("#a6611a","#dfc27d","#018571","#80cdc1")

fit <- list(Highsal = fitmeanhigh, Lowsal = fitmeanlow)
ggsurvplot(fit, data = svwk2017, combine = TRUE, # Combine curves
           #risk.table = TRUE,                  # Add risk table
           #conf.int = TRUE,                    # Add confidence interval
           #conf.int.style = "step",            # CI style, use "step" or "ribbon"
           censor = FALSE,                     # Remove censor points
           legend.labs=c("deep-high","shallow-high","deep-low", "shallow-low"),
           tables.theme = theme_cleantable(),  # Clean risk table
           xlim=c(2,11),
           break.x.by=1,
           palette = mycols)

####################
############### Predictive temperature plots ###########
#we are probably not going to run this, but here's the code anyway. 

# NOTE you will still have to update this code to run this, it has not been updated. 

#let's do only temperature and look at some crazy values at only Deep.  
fitmean30 <-survfit(cox_prev17, newdata=data.frame(depth=c("deep"),mean.temp=rep(30),mean.sal=rep(30.78107),temp.dur=rep(0),2), data=svwk2017)
fitmean27 <-survfit(cox_prev17, newdata=data.frame(depth=c("deep"),mean.temp=rep(27.5),mean.sal=rep(30.78107),temp.dur=rep(0),2), data=svwk2017)
fitmean25 <-survfit(cox_prev17, newdata=data.frame(depth=c("deep"),mean.temp=rep(25),mean.sal=rep(30.78107),temp.dur=rep(0),2), data=svwk2017)
fitmean22 <-survfit(cox_prev17, newdata=data.frame(depth=c("deep"),mean.temp=rep(22),mean.sal=rep(30.78107),temp.dur=rep(0),2), data=svwk2017)
fitmean18.5 <-survfit(cox_prev17, newdata=data.frame(depth=c("deep"),mean.temp=rep(17),mean.sal=rep(30.78107),temp.dur=rep(0),2),data=svwk2017)

deeponly <-filter(nd, depth=="D")
mycolstemp <-c("#993404","#d95f0e", "#fe9929", "#cccccc", "#bdd7e7")

fit <- list(temp_30 = fitmean30, temp_27 =fitmean27,temp_25 =fitmean25,temp_22 =fitmean22,temp_18.5 =fitmean18.5)
ggsurvplot(fit, data = svwk2017, combine = TRUE, # Combine curves
           #risk.table = TRUE,                  # Add risk table
           #conf.int = TRUE,                    # Add confidence interval
           #conf.int.style = "step",            # CI style, use "step" or "ribbon"
           censor = FALSE,                     # Remove censor points
           legend.labs=c("30 C","27 C","25 C", "22 C", "18.5 C"),
           tables.theme = theme_cleantable(),  # Clean risk table
           xlim=c(2,11),
           #ylim=c(0.25,1.00),
           break.x.by=1,
           palette = mycolstemp)


#let's do only temperature and look at some crazy values at only shallow.  
fitmean30 <-survfit(cox_prev17, newdata=data.frame(depth=c("shallow"),mean.temp=rep(30),mean.sal=rep(29.35297),temp.dur=rep(15.45455),2), data=svwk2017)
fitmean27 <-survfit(cox_prev17, newdata=data.frame(depth=c("shallow"),mean.temp=rep(27.5),mean.sal=rep(29.35297),temp.dur=rep(15.45455),2), data=svwk2017)
fitmean25 <-survfit(cox_prev17, newdata=data.frame(depth=c("shallow"),mean.temp=rep(25),mean.sal=rep(29.35297),temp.dur=rep(15.45455),2), data=svwk2017)
fitmean22 <-survfit(cox_prev17, newdata=data.frame(depth=c("shallow"),mean.temp=rep(22),mean.sal=rep(29.35297),temp.dur=rep(15.45455),2), data=svwk2017)
fitmean18.5 <-survfit(cox_prev17, newdata=data.frame(depth=c("shallow"),mean.temp=rep(17),mean.sal=rep(29.35297),temp.dur=rep(15.45455),2),data=svwk2017)

shalonly <-filter(nd, depth=="S")
mycolstemp <-c("#993404","#d95f0e", "#fe9929", "#cccccc", "#bdd7e7")

fit <- list(temp_30 = fitmean30, temp_27 =fitmean27,temp_25 =fitmean25,temp_22 =fitmean22,temp_18.5 =fitmean18.5)
ggsurvplot(fit, data = svwk2017, combine = TRUE, # Combine curves
           #risk.table = TRUE,                  # Add risk table
           #conf.int = TRUE,                    # Add confidence interval
           #conf.int.style = "step",            # CI style, use "step" or "ribbon"
           censor = FALSE,                     # Remove censor points
           legend.labs=c("30 C","27 C","25 C", "22 C", "18.5 C"),
           tables.theme = theme_cleantable(),  # Clean risk table
           xlim=c(2,11),
           #ylim=c(0.25,1.00),
           break.x.by=1,
           palette = mycolstemp)

###################
#### Survival probability vs. time ####
#change out the model and dataset for year. 
cb <-prediction(cox_prev17, type="expected") %>% mutate(survprob=exp(-fitted)) %>% mutate(inst = 1-exp((1-survprob)/7))
cbgg <- ggplot(cb, aes(x=week, y=survprob, color=station))
cbgg + geom_line()+theme_bw() + ylim(0,1.0)

#######
######### Violin plots of covariates ###############
# are the model covariates actually different by site.

#2017
nd <-ddply(svwk2017,(station~week), summarise, mindo=mean(min.do),avdo=mean(mean.do),mxdo=mean(max.do),mintemp=mean(min.temp), mxtemp =mean(max.temp),avtemp=mean(mean.temp), dodur=mean(do.dur),tdur=mean(temp.dur), avsal=mean(mean.sal),avsat=mean(ssat))
nd <-mutate(nd, site=substr(station,1,1),depth=substr(station,2,2))
#growth covariates mean.temp
ggplot(nd,aes(x=station,y=avtemp))+ geom_violin(trim=FALSE, fill='#A4A4A4', color="black")+ geom_boxplot(width=0.1) + labs(x="station", y = "Mean Temperature (celcius)")+ theme_classic() #this is the real one because not based on how many fish.

#survival covariates ~depth*mean.temp+temp.dur+mean.sal 
ggplot(nd,aes(x=depth, y=avtemp))+ geom_violin(trim=FALSE, fill='#A4A4A4', color="black")+ geom_boxplot(width=0.1) + labs(x="depth", y = "Mean Temperature (celcius)")+ theme_classic() #this is the real one because not based on how many fish.
ggplot(nd,aes(x=depth,y=tdur))+ geom_violin(trim=FALSE, fill='#A4A4A4', color="black")+ geom_boxplot(width=0.1) + labs(x="depth", y = "Duration > 27 degrees C (minutes)")+ theme_classic() #this is the real one because not based on how many fish.
ggplot(nd,aes(x=depth,y=avsal))+ geom_violin(trim=FALSE, fill='#A4A4A4', color="black")+ geom_boxplot(width=0.1) + labs(x="depth", y = "Mean Salinity (ppt)")+ theme_classic() #this is the real one because not based on how many fish.



#2016
nd <-ddply(svwk2016,(station~week), summarise, mindo=mean(min.do),avdo=mean(mean.do),mxdo=mean(max.do),mintemp=mean(min.temp), mxtemp =mean(max.temp),avtemp=mean(mean.temp), dodur=mean(do.dur),tdur=mean(temp.dur), avsal=mean(mean.sal),avsat=mean(ssat))
nd <-mutate(nd, site=substr(station,1,2),depth=substr(station,3,3))

#growth covariates: site + max.do + min.temp + temp.dur
ggplot(nd,aes(x=site,y=mxdo))+ geom_violin(trim=FALSE, fill='#A4A4A4', color="black")+ geom_boxplot(width=0.1) + labs(x="site", y = "Maximum Dissolved Oxygen (mg/L)")+ theme_classic() #this is the real one because not based on how many fish.
ggplot(nd,aes(x=site,y=tdur))+ geom_violin(trim=FALSE, fill='#A4A4A4', color="black")+ geom_boxplot(width=0.1) + labs(x="site", y = "Duration > 27 degrees C (minutes)")+ theme_classic() #this is the real one because not based on how many fish.
ggplot(nd,aes(x=site,y=mintemp))+ geom_violin(trim=FALSE, fill='#A4A4A4', color="black")+ geom_boxplot(width=0.1) + labs(x="site", y = "Minimum Temperature (Celcius)")+ theme_classic() #this is the real one because not based on how many fish.

#survival covariates site + min.temp

############

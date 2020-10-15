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

#setwd("/Users//tdolan/Documents//WIP research/Caging paper/data/heather meeting/caging_rproj")
setwd("/Users/tdolan/Documents/R-Github/WFCageStudy")


svwk2016 <-read.csv("weeklysummarydata_2016_10-8-19.csv", header=TRUE) # the weekly summary data, covariates are unscaled
svwk2017 <-read.csv("weeklysummarydata_2017.csv", header=TRUE) # the weekly summary data, covariates are unscaled
cagetotals17 <-read.csv("cagetotals_17fixed.csv", header=TRUE)
cagetotals16 <-read.csv("cagetotals_16.csv", header=TRUE)

#2016 best models - dcor16 is best
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
min.mintemp.2016e <-10
min.mintemp.2016l <-18.9
meantemp17.early <-19.9
meantemp17.late <-20.36
meansal17.early <-31.14
meansal17.late <-31.08
tempdur17.early <-0
tempdur17.late <-0
tempdur16.early <-0
tempdur16.late <-2


#marginal effects with site.
fitearly16 <-survfit(dcor16, newdata=data.frame(site=c("CP"),min.temp=rep(min.mintemp.2016e),1), data=svwk2016)
fitlate16 <-survfit(dcor16, newdata=data.frame(site=c("CP"),min.temp=rep(min.mintemp.2016l),1), data=svwk2016)
fitearly17 <-survfit(cox_prev17, newdata=data.frame(depth=c("deep","shallow"),mean.temp=rep(meantemp17.early),mean.sal=rep(meansal17.early),temp.dur=rep(tempdur17.early),2), data=svwk2017)
fitlate17 <-survfit(cox_prev17, newdata=data.frame(depth=c("deep","shallow"),mean.temp=rep(meantemp17.late),mean.sal=rep(meansal17.late),temp.dur=rep(tempdur17.late),2), data=svwk2017)

#2016 - survival using real life data! 
mycols <- c("#a6bddb", "#1c9099")
fit <- list(earlytemp16 = fitearly16, latetemp16 = fitlate16)
ggsurvplot(fit, data = svwk2016, combine = TRUE, # Combine curves
           risk.table = FALSE,                  # Add risk table
           conf.int = TRUE,                   # Add confidence interval
           conf.int.style = "ribbon",            # CI style, use "step" or "ribbon"
           censor = FALSE,                     # Remove censor points
           legend.labs=c("2016_early","2016_late"),
           tables.theme = theme_cleantable(),  # Clean risk table
           xlim=c(2,11),
           break.x.by=1,
           palette = mycols)
#to me this plot suggests a significant site/min.temp interaction... 


#2017 - survival using real life data! 
mycols <- c("#748499","#a6bddb","#10565B","#1c9099")
fit <- list(earlytemp17 = fitearly17, latetemp17 = fitlate17)
ggsurvplot(fit, data = svwk2017, combine = TRUE, # Combine curves
           risk.table = FALSE,                  # Add risk table
           conf.int = FALSE,                   # Add confidence interval
           conf.int.style = "ribbon",            # CI style, use "step" or "ribbon"
           censor = FALSE,                     # Remove censor points
           legend.labs=c("2017 early deep","2017 early shallow","2017 late deep","2017 late shallow"),
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
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

#2017 
#cb17 <-prediction(cox_prev17, type="expected") %>% mutate(survprob=exp(-fitted), fittedup =fitted + se.fitted, fittedlo = fitted-se.fitted) %>% mutate(inst = 1-exp((1-survprob)/7), spup=exp(-fittedup),splo=exp(-fittedlo))
cb17 <-prediction(cox_prev17, type="expected", calculate_se = TRUE)

#resample to get rid of NANs
for(i in 1:length(cb17$se.fitted)){
  if(is.nan.data.frame(cb17$se.fitted[i])){
    cb17$se.fitted[i]<-sample(na.omit(cb17$se.fitted),size=1,replace=TRUE)}
  if(cb17$se.fitted[i] == 0){
    cb17$se.fitted[i]<-sample(na.omit(cb17$se.fitted),size=1,replace=TRUE)}
  if(cb17$fitted[i]== 0){
    cb17$fitted[i]<-sample(na.omit(cb17$fitted),size=1,replace=TRUE)}
  }
cb17$se.fitted <- as.numeric(as.character(cb17$se.fitted))
cb17$fitted <- as.numeric(as.character(cb17$fitted))

cagetotals17 <-read.csv("cagetotals_17fixed.csv", header=TRUE)
cb17 <- mutate(cb17, survprob=exp(-fitted),u=fitted+se.fitted, l=fitted-se.fitted)%>%
  mutate(upper.se=exp(-u), lower.se=exp(-l)) %>% 
  mutate(inst = 1-exp((1-survprob)/7), nd=10-survivors, ubar =survprob +u, lbar=survprob-l)
cagedeaths17 <-dplyr::select(cagetotals16, Site, depth, Week, CageID, deaths) %>% dplyr::rename(site=Site,cage=CageID, week=Week)
cb17 <-left_join(cb17, cagedeaths16, by=c("site","depth","week","cage"))
cb17 <-mutate(cb17, death.per=deaths/survivors, surv.per=(survivors-deaths)/survivors)



coeff <-1
cb17 %>%
  filter(!is.nan(se.fitted))%>%
  filter(week > 0)%>%
ggplot(aes(x=week))+
  #geom_point(aes(y=surv.per),color="black", alpha=0.2, shape =1)+
  geom_jitter(aes(y=surv.per),width = 0.8, height = 0.05, alpha=0.8, shape =19, size=1, color="dark gray")+
 # geom_point(aes(y=surv.per/coeff),color="black", alpha=0.5, shape =1)+
  #geom_smooth(aes(y=surv.per/coeff, method="loess"),color="blue",linetype="dashed", size=0.5, se=FALSE)+
  geom_line(aes(y=survprob), color="steelblue3")+
  #scale_y_continuous(limits=c(0,1),name = " ",)+
  geom_ribbon(aes(ymin=upper.se, ymax=lower.se), fill="steelblue3",alpha=0.4)+
    #sec.axis=sec_axis(trans=~.*1, name="percent survived"))+
    #+theme_bw() + ylim(0,1.0)+
  #scale_x_discrete(limits=c(3,6,9),labels=c("3"="29-Jun","6"="16-Jul","9"="3-Aug"), name="")+
  facet_grid(site~depth)+
  #theme(axis.line.y.right = element_line(color = "blue"), 
  #       axis.ticks.y.right = element_line(color = "blue"),
  #      axis.line.y.left = element_line(color = "red"), 
   #     axis.ticks.y.left = element_line(color = "red"),
   #     panel.background = element_rect(fill = "white", colour = "white"))
  theme_few()
#ggsave("survpred17jitter.png", path="/Users/tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()


#2016
cb16 <-prediction(dcor16, type="expected", calculate_se = TRUE) %>% 
  mutate(survprob=exp(-fitted),u=fitted+se.fitted, l=fitted-se.fitted)%>%
  mutate(upper.se=exp(-u), lower.se=exp(-l)) %>% 
  mutate(inst = 1-exp((1-survprob)/7), nd=10-survivors, ubar =survprob +u, lbar=survprob-l)
cagedeaths16 <-dplyr::select(cagetotals16, Site, depth, Week, CageID, deaths) %>% dplyr::rename(site=Site,cage=CageID, week=Week)
cb16 <-left_join(cb16, cagedeaths16, by=c("site","depth","week","cage"))
cb16 <-mutate(cb16, death.per=deaths/survivors, surv.per=(survivors-deaths)/survivors)


coeff <-1
cb16 %>%
  filter(!is.nan(se.fitted))%>%
  filter(week > 0)%>%
  ggplot(aes(x=week))+
  geom_jitter(aes(y=surv.per),width = 0.8, height = 0.05, alpha=0.8, shape =19, size=1, color="dark gray")+
 # geom_point(aes(y=surv.per/coeff),color="blue", alpha=0.02)+
 # geom_smooth(aes(y=surv.per/coeff, method="loess"),color="blue",linetype="dashed", size=0.5, se=FALSE)+
  geom_line(aes(y=survprob), color="steelblue3")+
  geom_ribbon(aes(ymin=lower.se, ymax=upper.se), fill="steelblue3",alpha=0.4)+
  scale_y_continuous(limits=c(0,1),name = " ")+
  #sec.axis=sec_axis(trans=~.*1, name="percent survived"))+
  #+theme_bw() + ylim(0,1.0)+
  #scale_x_discrete(limits=c(3,6,9),labels=c("3"="30-Jun","6"="20-Jul","9"="10-Aug"), name="")+
  facet_grid(site~depth)+
  theme_few()
#ggsave("survpred16_jitter.png", path="/Users/tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()





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

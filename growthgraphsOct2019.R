
#### Visualizations for lme growth models for caging study 
### October 19, 2019
### model selection done in the script "growthrevisited2019.R" in the folder caging_rproj

#packages you will need
library("MASS")
library("ggfortify")
library("MuMIn")
library("car")
library("multcomp")
library("nlme")
library("lattice")
library("car")
library("ggplot2")
library("reshape2")
library("plyr")
library("dplyr")
library("tidyr")
library("corrplot")
library("moments")
library("reshape")
library("lme4")
library("ggpubr")
library("plotrix")

#keeping everything in this folder
setwd("/Users//tdolan/Documents/Documents//R-Github/WFCageStudy")

#load the data
svrt16 <-read.csv("svrt2016forgrowth.csv", header=TRUE)
svwk2016 <-read.csv("weeklysummarydata_2016_10-8-19.csv", header=TRUE) # the weekly summary data, covariates are unscaled
svrt17 <-read.csv("svrt2017forgrowth.csv", header=TRUE)
svwk2017 <-read.csv("weeklysummarydata_2017.csv", header=TRUE) # the weekly summary data, covariates are unscaled

###Generate the datasets####
#2016
grow16 <-filter(svrt16, week > 1) # we have no environmental data for week 0 or 1 because that's before we started recording. 
grow16 <-mutate(grow16,av_length=ifelse(av_length==0,NA,av_length))#convert av_length of 0 to NA
grow16 <-filter(grow16, survivors > 0) #get rid of weeks where no fish are left.
#grow16 <-mutate_at(grow16,vars(week,site,depth,cage),factor) #factorize the variables. 
#create a dataset where variables are scaled. 
grow.scale16 <- mutate_at(grow16, vars(biom, min.do, max.do, mean.do, sd.do, min.temp, sd.do, min.temp, max.temp, sd.T, do.dur, temp.dur, mean.sal, sd.sal, ssat), scale)

#2017
grow17 <-filter(svrt17, week > 0) # we have no environmental data for week 0 or 1 because that's before we started recording. 
grow17 <-dplyr::rename(grow17, av_length=avlength)
grow17 <-mutate(grow17,av_length=ifelse(av_length==0,NA,av_length))#convert av_length of 0 to NA
grow17 <-filter(grow17, survivors > 0)#get rid of weeks where no fish are left.
#grow17 <-mutate_at(grow17,vars(week,site,depth,cage),factor) #factorize the variables. 
#create a dataset where variables are scaled. 
grow.scale17 <- mutate_at(grow17, vars(biom, min.do, max.do, mean.do, sd.do, min.temp, sd.do, min.temp, max.temp, sd.T, do.dur, temp.dur, mean.sal, sd.sal, ssat), scale)
##############
###### Visualization no model, just averages ######
#actual fishlength data with the average OVERLAID
ggplot(data=svwk2016,aes(x=week,y=fishlength, color = cage))+
  geom_point(alpha=1/10)+
  #scale_x_discrete(breaks=c("2","3","4","5","6","7","8","9","11"), labels=c("14-Jun","24-Jun","4-Jul","14-Jul","24-Jul","3-Aug","13-Aug","23-Aug","2-Sep"))+  # well that aggressively didn't work. 
  geom_line(aes(y=av_length, color=cage)) +
  facet_grid(site~depth)+
  theme_pubclean()

ggplot(data=svwk2017,aes(x=week,y=fishlength, color = cage))+
  geom_point(alpha=1/10)+
  #scale_x_discrete(breaks=c("2","3","4","5","6","7","8","9","11"), labels=c("14-Jun","24-Jun","4-Jul","14-Jul","24-Jul","3-Aug","13-Aug","23-Aug","2-Sep"))+  # well that aggressively didn't work. 
  geom_line(aes(y=avlength, color=cage)) +
  facet_grid(site~depth)+
  theme_pubclean()
########

##2016 top models
pwin16 <-lme(av_length ~ site + max.do + min.temp + temp.dur,random=~week|cage, data=grow16, method="ML") 

#2017 top models 
d73A <-lme(av_length ~ mean.temp ,random=~week|cage, data=grow17, method="ML")

##########
########### Residuals plots #########
#remember to change out the model and year. 
plot(d73A) #plot residuals vs. fitted values
grow17$predicted <-predict(d73A)
grow17$residuals <-residuals(d73A)

#residuals v. week
ggplot(grow17, aes(x=week, y=residuals))+
  geom_point(shape=21)+
  geom_smooth(method="lm",color="blue",linetype="dotted")+
  geom_hline(yintercept=0, linetype="dotted",col="red")+
  theme_bw()

#resids to fitted values. # I am actually not going to run this one. 
ggplot(grow16,aes(x=week,y=av_length)) +
  geom_smooth(method="lm",se=FALSE,color="lightgrey")+                    # this is the light grey smoothing function. 
  geom_segment(aes(xend=week,yend=predicted), alpha = .2) + #segment connecting predicted to residual. alpha fades the color
  #geom_point()+                           #the raw data, we can silence this to see the residuals better.
  geom_point(aes(color=abs(residuals)))+    #coloring the residuals
  scale_color_continuous(low="lightgrey",high="red") + # very bad residuals in dark blue 
  guides(color=FALSE)+                      #color legend removed. 
  ylab("Average length (mm)")+
  geom_point(aes(y=predicted),shape=3,alpha=1/5)+    #add the predicted values as crosses  
  theme_classic() 
##################
######## growth models over data, by cage #############

#over fishlength & cage average 2016
ggplot(data=svwk2016,aes(x=week,y=fishlength, color = cage))+
  geom_point(alpha=1/3, shape=1)+
  geom_point(data=grow16,aes(x=week,y=av_length, color = cage))+
  #scale_x_discrete(breaks=c("2","3","4","5","6","7","8","9","11"), labels=c("14-Jun","24-Jun","4-Jul","14-Jul","24-Jul","3-Aug","13-Aug","23-Aug","2-Sep"))+  # well that aggressively didn't work. 
  #geom_line(aes(y=av_length, color=cage)) + #actual average
  ylab("Average length (mm)")+
  geom_line(data = cbind(grow16, pred = predict(pwin16, level=1)), aes(y = pred)) + #predicted average
  facet_grid(site~depth)+
  theme_classic()

#over fishlength & cage average 2017
ggplot(data=svwk2017,aes(x=week,y=fishlength, color = cage))+
  geom_point(alpha=1/3, shape=1)+
  geom_point(data=grow17,aes(x=week,y=av_length, color = cage))+
  #scale_x_discrete(breaks=c("2","3","4","5","6","7","8","9","11"), labels=c("14-Jun","24-Jun","4-Jul","14-Jul","24-Jul","3-Aug","13-Aug","23-Aug","2-Sep"))+  # well that aggressively didn't work. 
  #geom_line(aes(y=av_length, color=cage)) + #actual average
  ylab("Average length (mm)")+
  geom_line(data = cbind(grow17, pred = predict(d73A, level=1)), aes(y = pred)) + #predicted average
  facet_grid(site~depth)+
  theme_classic()
############################
############ models by cage over data ######

#2016
newdf16 <-dplyr::select(grow16, cage,site,depth,station,week,av_length,min.temp,max.do,temp.dur)
newdf16$pred <-predict(pwin16, newdata=newdf16, level=1)

#now extract the confidence intervals. 
#create design matrix
Designmat <- model.matrix(eval(eval(pwin16$call$fixed)[-2]), newdf16[-ncol(newdf16)])
#compute standard error for predictions
predvar <- diag(Designmat %*% pwin16$varFix %*% t(Designmat))
newdf16$SE <- sqrt(predvar) 
newdf16$SE2 <- sqrt(predvar+pwin16$sigma^2)

p1 <- ggplot(newdf16,aes(x=week,y=pred)) + 
  geom_line() +
  geom_ribbon(aes(ymin=pred-2*SE2,ymax=pred+2*SE2),alpha=0.2,fill="red") +
  geom_ribbon(aes(ymin=pred-2*SE,ymax=pred+2*SE),alpha=0.2,fill="blue") +
  geom_point(data=svwk2016,aes(x=week,y=fishlength),alpha=0.2,fill="gray") +
  geom_point(data=newdf16,aes(x=week,y=av_length),fill="blue",shape=23) +
  ylab("Average Length (mm)")+
  #scale_y_continuous(name="average length (mm)", limits=c(40,100))+
  #facet_grid(site~depth)+
  #facet_wrap(~site)+
  facet_wrap(~cage,ncol=4)+
  theme_minimal()
p1

#2017
newdf17 <-dplyr::select(grow17, cage,site,depth,station,week,av_length,mean.temp)
newdf17$pred <-predict(d73A, newdata=newdf17, level=1)
#now extract the confidence intervals. 
#create design matrix
Designmat <- model.matrix(eval(eval(d73A$call$fixed)[-2]), newdf17[-ncol(newdf17)])
#compute standard error for predictions
predvar <- diag(Designmat %*% d73A$varFix %*% t(Designmat))
newdf17$SE <- sqrt(predvar) 
newdf17$SE2 <- sqrt(predvar+d73A$sigma^2)

p2 <- ggplot(newdf17,aes(x=week,y=pred)) + 
  geom_line() +
  geom_ribbon(aes(ymin=pred-2*SE2,ymax=pred+2*SE2),alpha=0.2,fill="red") +
  geom_ribbon(aes(ymin=pred-2*SE,ymax=pred+2*SE),alpha=0.2,fill="blue") +
  geom_point(data=svwk2017,aes(x=week,y=fishlength),alpha=0.2,fill="gray") +
  geom_point(data=newdf17,aes(x=week,y=av_length),fill="blue",shape=23) +
  ylab("Average Length (mm)")+
  #scale_y_continuous(name="average length (mm)", limits=c(40,100))+
  #facet_grid(site~depth)+
  #facet_wrap(~site)+
  facet_wrap(~cage,ncol=4)+
  theme_minimal()
p2
##############
########## marginal effects site 2016, ribbon & sawtooth #########
#there is no depth separation in 2017 because depth is not in the model. 

#across all sites
avmaxdo <-mean(grow16$max.do)#13.92222
avmintemp <-mean(grow16$min.temp)#19.38018
avtdur <-mean(grow16$temp.dur)#581.1594

#re-run from here
newdf5 <-dplyr::select(grow16, site, week, cage, depth,av_length)

#prediction level
lev <-1
fit.stat16 <- predict(pwin16, newdata=newdf5, data=grow16,level=lev)

#populate the new dataframe
newdf5$max.do=rep(avmaxdo)
newdf5$min.temp=rep(avmintemp)
newdf5$temp.dur=rep(avtdur)

#CP
newdf5_cp <-filter(newdf5,site=="CP")
fit.site16cp <- predict(pwin16, newdata=newdf5_cp, data=grow16,level=lev)
Designmatcp <- model.matrix(eval(eval(pwin16$call$fixed)[-2]), newdf5_cp[-ncol(newdf5_cp)])
newdf5_cp <- cbind(newdf5_cp,fit.site16cp)
#compute standard error for predictions
predvarcp <- diag(Designmatcp %*% pwin16$varFix %*% t(Designmatcp))
newdf5_cp$SE <- sqrt(predvarcp) 
newdf5_cp$SE2 <- sqrt(predvarcp+pwin16$sigma^2)

#MD
newdf5_md <-filter(newdf5,site=="MD")
fit.site16md <- predict(pwin16, newdata=newdf5_md, data=grow16,level=lev)
Designmatmd <- model.matrix(eval(eval(pwin16$call$fixed)[-2]), newdf5_md[-ncol(newdf5_md)])
newdf5_md <- cbind(newdf5_md,fit.site16md)
#compute standard error for predictions
predvarmd <- diag(Designmatmd %*% pwin16$varFix %*% t(Designmatmd))
newdf5_md$SE <- sqrt(predvarmd) 
newdf5_md$SE2 <- sqrt(predvarmd+pwin16$sigma^2)

#RS
newdf5_rs <-filter(newdf5,site=="RS")
fit.site16rs <- predict(pwin16, newdata=newdf5_rs, data=grow16,level=lev)
Designmatrs <- model.matrix(eval(eval(pwin16$call$fixed)[-2]), newdf5_rs[-ncol(newdf5_rs)])
newdf5_rs <- cbind(newdf5_rs,fit.site16rs)
#compute standard error for predictions
predvarrs <- diag(Designmatrs %*% pwin16$varFix %*% t(Designmatrs))
newdf5_rs$SE <- sqrt(predvarrs) 
newdf5_rs$SE2 <- sqrt(predvarrs+pwin16$sigma^2)

#p2 <- ggplot(newdf5,aes(x=week,y=av_length,color=cage)) + #if you want the cage level ribbons, and no sawtooth, must turn off geom_line
p2 <- ggplot(newdf5,aes(x=week)) + 
  #geom_line(data=newdf5_cp,aes(x=week,y=fit.site16cp),color="darkgreen") +
  geom_pointrange(data=newdf5_cp,aes(x=week,y=fit.site16cp, ymin=fit.site16cp-2*SE,ymax=fit.site16cp+2*SE),alpha=0.3,color="darkgreen")+
  #geom_point(data=newdf5_cp,aes(x=week,y=fit.site16cp),color="darkgreen") +
  #geom_ribbon(data=newdf5_cp,aes(ymin=fit.site16cp-2*SE2,ymax=fit.site16cp+2*SE2),alpha=0.3,fill="darkgreen") +
  #geom_ribbon(data=newdf5_cp,aes(ymin=fit.site16cp-2*SE,ymax=fit.site16cp+2*SE),alpha=0.1,fill="darkgreen") +
  
  #geom_line(data=newdf5_md,aes(x=week,y=fit.site16md),color="firebrick") +
  geom_pointrange(data=newdf5_md,aes(x=week,y=fit.site16md, ymin=fit.site16md-2*SE,ymax=fit.site16md+2*SE),alpha=0.3,color="firebrick")+
  #geom_point(data=newdf5_md,aes(x=week,y=fit.site16md),color="firebrick") +
 # geom_ribbon(data=newdf5_md,aes(ymin=fit.site16md-2*SE2,ymax=fit.site16md+2*SE2),alpha=0.3,fill="firebrick") +
  #geom_ribbon(data=newdf5_md,aes(ymin=fit.site16md-2*SE,ymax=fit.site16md+2*SE),alpha=0.1,fill="firebrick") +
  
  #geom_line(data=newdf5_rs,aes(x=week,y=fit.site16rs),color="blue") +
  geom_pointrange(data=newdf5_rs,aes(x=week,y=fit.site16rs, ymin=fit.site16rs-2*SE,ymax=fit.site16rs+2*SE),alpha=0.3,color="blue")+
  #geom_point(data=newdf5_rs,aes(x=week,y=fit.site16rs),color="blue") +
  #geom_ribbon(data=newdf5_rs,aes(ymin=fit.site16rs-2*SE2,ymax=fit.site16rs+2*SE2),alpha=0.3,fill="blue") +
  #geom_ribbon(data=newdf5_rs,aes(ymin=fit.site16rs-2*SE,ymax=fit.site16rs+2*SE),alpha=0.1,fill="blue") +
  
  #geom_point(data=svwk2016,aes(x=week,y=fishlength),alpha=0.2,fill="gray") +
  #geom_point(data=newdf5,aes(x=week,y=av_length),fill="blue",shape=23) +
  ylab("Average Length (mm)")+
  scale_y_continuous("y")+
  #facet_grid(site~depth)+
  facet_wrap(~site)+
  #facet_wrap(~cage,ncol=4)+
  theme_minimal()
p2

#this looks weird. What are my options to get rid of the sawtooth?
# 1. I could make a separate model for every cage, with bands.
# 2> I could fith a geom_smooth over the data
# I could combine all the data into one
# i could plot by cage and do a sneaky layers plot? 
# i could get rid of the SE data? & plot 1 line per cage. 
#############
########## further attempts at marginal effect of site  #########

#across all sites
avmaxdo <-mean(grow16$max.do)#13.92222
avmintemp <-mean(grow16$min.temp)#19.38018
avtdur <-mean(grow16$temp.dur)#581.1594

#re-run from here
newdf5 <-dplyr::select(grow16, site, week, cage, depth,av_length)

#prediction level
lev <-1
predall <- predict(pwin16, newdata=newdf5, data=grow16,level=lev)

#populate the new dataframe
newdf5$max.do=rep(avmaxdo)
newdf5$min.temp=rep(avmintemp)
newdf5$temp.dur=rep(avtdur)

#predcp
newdf5_cp <-filter(newdf5,site=="CP")
fit.site16cp <- predict(pwin16, newdata=newdf5_cp, data=grow16,level=lev)
Designmatcp <- model.matrix(eval(eval(pwin16$call$fixed)[-2]), newdf5_cp[-ncol(newdf5_cp)])
newdf5_cp <- cbind(newdf5_cp,fit.site16cp)
#compute standard error for predictions
predvarcp <- diag(Designmatcp %*% pwin16$varFix %*% t(Designmatcp))
newdf5_cp$SE <- sqrt(predvarcp) 
newdf5_cp$SE2 <- sqrt(predvarcp+pwin16$sigma^2) 
newdf5_cp <-rename(newdf5_cp,fit.site = fit.site16cp)
#MD
newdf5_md <-filter(newdf5,site=="MD")
fit.site16md <- predict(pwin16, newdata=newdf5_md, data=grow16,level=lev)
Designmatmd <- model.matrix(eval(eval(pwin16$call$fixed)[-2]), newdf5_md[-ncol(newdf5_md)])
newdf5_md <- cbind(newdf5_md,fit.site16md)
#compute standard error for predictions
predvarmd <- diag(Designmatmd %*% pwin16$varFix %*% t(Designmatmd))
newdf5_md$SE <- sqrt(predvarmd) 
newdf5_md$SE2 <- sqrt(predvarmd+pwin16$sigma^2)
newdf5_md <-rename(newdf5_md,fit.site = fit.site16md)
#RS
newdf5_rs <-filter(newdf5,site=="RS")
fit.site16rs <- predict(pwin16, newdata=newdf5_rs, data=grow16,level=lev)
Designmatrs <- model.matrix(eval(eval(pwin16$call$fixed)[-2]), newdf5_rs[-ncol(newdf5_rs)])
newdf5_rs <- cbind(newdf5_rs,fit.site16rs)
#compute standard error for predictions
predvarrs <- diag(Designmatrs %*% pwin16$varFix %*% t(Designmatrs))
newdf5_rs$SE <- sqrt(predvarrs) 
newdf5_rs$SE2 <- sqrt(predvarrs+pwin16$sigma^2)
newdf5_rs <-rename(newdf5_rs,fit.site = fit.site16rs)

newdf6 <-bind_rows(newdf5_cp,newdf5_md,newdf5_rs)

ddnewdf6 <-ddply(newdf6,week~site,summarise,av_length=mean(av_length), max.do=max(max.do),min.temp=min(min.temp),temp.dur=mean(temp.dur),fit.site=mean(fit.site),SE=mean(SE),SE2=mean(SE2))

#a line for each predicted cage trajectory, over actual average lengths. 
ggplot(data=grow16,aes(x=week, y=av_length,color=site))+
  #geom_point()+
  ylab("Average length (mm)")+
  #geom_line(data = cbind(grow16, pred = predict(pwin16, level=1)), aes(y = pred)) + #predicted average
  geom_line(data=ddnewdf6, aes(x=week,y=fit.site))+
  #geom_ribbon(data=newdf6,aes(ymin=fit.site-2*SE2,ymax=fit.site+2*SE2,x=week),alpha=0.1,fill="gray",inherit.aes=FALSE)+
  geom_point(data=newdf6,aes(y=fit.site,x=week),alpha=0.3)+
  #geom_smooth(data=newdf6,method=lm,se=FALSE)+
  geom_linerange(data=newdf6,aes(ymin=fit.site-2*SE2,ymax=fit.site+2*SE2,x=week),alpha=0.1,inherit.aes=FALSE)+
  facet_wrap(~site)+
  theme_classic()

#a line for the mean predicted cage, over actual average lengths.##
  ## a very hacky way to do it###

ggplot(data=ddnewdf6,aes(y=fit.site,x=week,color=site))+
  geom_line()+
  ylab("Average length (mm)")+
  #geom_point(data=newdf6,aes(y=fit.site,x=week),alpha=0.5)+
  geom_ribbon(data=ddnewdf6,aes(ymin=fit.site-2*SE2,ymax=fit.site+2*SE2,x=week),fill="gray",alpha=0.5,inherit.aes=FALSE)+
  #geom_smooth(data=newdf6,method=lm,se=FALSE)+
  #geom_linerange(data=newdf6,aes(ymin=fit.site-2*SE2,ymax=fit.site+2*SE2,x=week),alpha=0.1,inherit.aes=FALSE)+
  facet_wrap(~site)+
  theme_classic()
#####################
####### marginal effects of other variables 2016 ######

#across all sites
avmaxdo <-mean(grow16$max.do)#13.92222
avmintemp <-mean(grow16$min.temp)#19.38018
avtdur <-mean(grow16$temp.dur)#581.1594

#re-run from here
newdf <-dplyr::select(grow16, site, week, cage, depth,av_length)
newdf$max.do=rep(avmaxdo)
newdf$min.temp=rep(avmintemp)
newdf$temp.dur=rep(avtdur)

fit.site16 <- predict(pwin16, newdata=newdf, data=grow16)
newdf <- cbind(newdf,fit.site16)

#max.do - high, median and low.
#high
max(grow16$max.do)  
newdf$max.do=rep(23.8)
fit.maxdo_high <-predict(pwin16, newdata=newdf, data=grow16)
newdf <- cbind(newdf,fit.maxdo_high)
#mean
mean(grow16$max.do)
newdf$max.do=rep(13.92222)
fit.maxdo_mean <-predict(pwin16, newdata=newdf, data=grow16)
newdf <- cbind(newdf,fit.maxdo_mean)
#low
min(grow16$max.do)
newdf$max.do=rep(8.666918)
fit.maxdo_low <-predict(pwin16, newdata=newdf, data=grow16)
newdf <- cbind(newdf,fit.maxdo_low)

#min.temp - high, median and low.
newdf$max.do =rep(avmaxdo) #return max.do to it's mean before you start!
#high
max(grow16$min.temp)  
newdf$min.temp=rep(23.93)
fit.mintemp_high <-predict(pwin16, newdata=newdf, data=grow16)
newdf <- cbind(newdf,fit.mintemp_high)
#mean
mean(grow16$min.temp)
newdf$min.temp=rep(19.38018)
fit.mintemp_mean <-predict(pwin16, newdata=newdf, data=grow16)
newdf <- cbind(newdf,fit.mintemp_mean)
#low
min(grow16$min.temp)
newdf$min.temp=rep(12.31)
fit.mintemp_low <-predict(pwin16, newdata=newdf, data=grow16)
newdf <- cbind(newdf,fit.mintemp_low)

#temp.dur - high, median and low.
newdf$min.temp =rep(avmintemp) #return min.temp to it's mean before you start!
#high
max(grow16$temp.dur)  
newdf$temp.dur=rep(8910)
fit.tdur_high <-predict(pwin16, newdata=newdf, data=grow16)
newdf <- cbind(newdf,fit.tdur_high)
#mean
mean(grow16$temp.dur)
newdf$temp.dur=rep(581.1594)
fit.tdur_mean <-predict(pwin16, newdata=newdf, data=grow16)
newdf <- cbind(newdf,fit.tdur_mean)
#low
min(grow16$temp.dur)
newdf$temp.dur=rep(0)
fit.tdur_low <-predict(pwin16, newdata=newdf, data=grow16)
newdf <- cbind(newdf,fit.tdur_low)

#now we want to make basically a dataframe of all the fits and the actual average length, this will help with plots later. 
allfits <- dplyr::select(newdf,-max.do,-min.temp,-temp.dur,-av_length)
allfits <-full_join(grow16,allfits, by=c("cage","site","depth","week"))
allfits <-dplyr::select(allfits, -biom,-min.do,-max.do,-mean.do,-sd.do,-min.temp,-max.temp,-mean.temp,-sd.T,-do.dur,-temp.dur,-mean.sal,-sd.sal,-ssat,-residuals,-SD,-survivors, -X)

#melt version
af_melt <-gather(allfits,"fit", "pred_value", 6:17)
af_melt$fit <-as.factor(af_melt$fit)

#just another way to produce the original model over data, but quickly. 
allfits %>% 
  #filter(station=="CPS") %>%
  ggplot(aes(x=week,y=av_length, color=cage))+
  geom_point(alpha=1/2)+
  ylab("Average length (mm)")+
  geom_line(aes(y = predicted),linetype="dashed") + #predicted from original model
  #geom_line(aes(y = av_length)) + #actual model
  facet_wrap(~site)+
  theme_classic()

af_melt %>%
  #filter(station=="CPD") %>%
  filter(fit=="fit.maxdo_high"| fit=="fit.maxdo_low" |fit=="fit.maxdo_mean") %>%
  ggplot(aes(x=week,y=pred_value, color=fit))+
  ylab("Average length (mm)")+
  geom_line()+
  facet_wrap(~station)+
  #facet_grid(site~depth)+
  theme_bw()

#so you could plot every cage with it's own geom_line... that is an option. 
af_meltDO <- dplyr::filter(af_melt,fit=="fit.maxdo_high"|fit=="fit.maxdo_low"|fit=="fit.maxdo_mean")
ggplot(data=af_meltDO, aes(x=week,y=pred_value, color=fit))+
  ylab("Average length (mm)")+
  geom_line(data=subset(af_meltDO, cage=="CPD1" ),aes(x=week,y=pred_value))+
  geom_line(data=subset(af_meltDO, cage=="CPD2" ),aes(x=week,y=pred_value))+
  geom_line(data=subset(af_meltDO, cage=="CPD3" ),aes(x=week,y=pred_value))+
  geom_line(data=subset(af_meltDO, cage=="CPD4" ),aes(x=week,y=pred_value))+
  geom_line(data=subset(af_meltDO, cage=="CPS1" ),aes(x=week,y=pred_value))+
  geom_line(data=subset(af_meltDO, cage=="CPS2" ),aes(x=week,y=pred_value))+
  geom_line(data=subset(af_meltDO, cage=="CPS3" ),aes(x=week,y=pred_value))+
  geom_line(data=subset(af_meltDO, cage=="CPS4" ),aes(x=week,y=pred_value))+
  facet_wrap(~site)+
  theme_minimal()
#starts to get messy looking already

#You could create a mean value for ALL the cages in a site. 
ddaf_melt <-ddply(af_melt,week~site~fit,summarise,avpred=mean(pred_value),sepred = std.error(pred_value))

ggplot(data=ddaf_melt, aes(x=week,y=avpred, color=fit))+
  #scale_color_manual(values=wes_palette(n=3, name="Darjeeling1"))+
  ylab("Average length (mm)")+
  #geom_line(data=subset(ddaf_melt, fit=="fit.maxdo_high"|fit=="fit.maxdo_low"|fit=="fit.maxdo_mean"),aes(x=week,y=avpred,color=fit))+
  #geom_linerange(data=subset(ddaf_melt, fit=="fit.maxdo_high"|fit=="fit.maxdo_low"|fit=="fit.maxdo_mean"),aes(x=week,ymin=avpred-2*sepred,ymax=avpred+2*sepred),alpha=0.3,)+
  #geom_line(data=subset(ddaf_melt, fit=="fit.mintemp_high"|fit=="fit.mintemp_low"|fit=="fit.mintemp_mean"),aes(x=week,y=avpred,color=fit))+
  #geom_linerange(data=subset(ddaf_melt, fit=="fit.mintemp_high"|fit=="fit.mintemp_low"|fit=="fit.mintemp_mean"),aes(x=week,ymin=avpred-2*sepred,ymax=avpred+2*sepred),alpha=0.3,)+
  geom_line(data=subset(ddaf_melt, fit=="fit.tdur_high"|fit=="fit.tdur_low"|fit=="fit.tdur_mean"),aes(x=week,y=avpred,color=fit))+
  geom_linerange(data=subset(ddaf_melt, fit=="fit.tdur_high"|fit=="fit.tdur_low"|fit=="fit.tdur_mean"),aes(x=week,ymin=avpred-2*sepred,ymax=avpred+2*sepred),alpha=0.3,)+
  facet_wrap(~site)+
  theme_classic()

############# margin plots for 2017 ##########

#across all sites
avmeantemp <-mean(grow17$mean.temp)#21.27109

#re-run from here
newdf <-dplyr::select(grow17, site, week, cage, depth,av_length)
newdf$mean.temp=rep(avmeantemp)

fit.site17 <- predict(d73A, newdata=newdf, data=grow17)
#newdf <- cbind(newdf,fit.site16)

#mean.temp - high, median and low.
#high
max(grow17$mean.temp)  
newdf$mean.temp=rep(23.69685)
fit.meantemp_high <-predict(d73A, newdata=newdf, data=grow17)
newdf <- cbind(newdf,fit.meantemp_high)
#mean
mean(grow17$mean.temp)
newdf$mean.temp=rep(21.27109)
fit.meantemp_mean <-predict(d73A, newdata=newdf, data=grow17)
newdf <- cbind(newdf,fit.meantemp_mean)
#low
min(grow17$mean.temp)
newdf$mean.temp=rep(17.32684)
fit.meantemp_low <-predict(d73A, newdata=newdf, data=grow17)
newdf <-as.data.frame(cbind(newdf,fit.meantemp_low))


#now we want to make basically a dataframe of all the fits and the actual average length, this will help with plots later. 
allfits17 <- dplyr::select(newdf,-mean.temp,-av_length)
allfits17 <-full_join(grow17,allfits17, by=c("cage","site","depth","week"))
allfits17 <-dplyr::select(allfits17, -biom,-min.do,-max.do,-mean.do,-sd.do,-min.temp,-max.temp,-mean.temp,-sd.T,-do.dur,-temp.dur,-mean.sal,-sd.sal,-ssat,-residuals,-SD,-survivors, -X)

#melt version
af_melt17 <-gather(allfits17,"fit", "pred_value", 8:11)
af_melt17$fit <-as.factor(af_melt17$fit)

#just another way to produce the original model over data, but quickly. 
allfits17 %>% 
  ggplot(aes(x=week,y=av_length, color=cage))+
  geom_point(alpha=1/2)+
  ylab("Average length (mm)")+
  geom_line(aes(y = predicted),linetype="dashed") + #predicted from original model
  #geom_line(aes(y = av_length)) + #actual model
  facet_wrap(~site)+
  theme_classic()

#You could create a mean value for ALL the cages in a site. 
ddaf_melt17 <-ddply(af_melt17,week~fit,summarise,avpred=mean(pred_value),sepred = std.error(pred_value))

ggplot(data=ddaf_melt17, aes(x=week,y=avpred, color=fit))+
  #scale_color_manual(values=wes_palette(n=3, name="Darjeeling1"))+
  ylab("Average length (mm)")+
  geom_line(data=subset(ddaf_melt17, fit=="fit.meantemp_high"|fit=="fit.meantemp_low"|fit=="fit.meantemp_mean"),aes(x=week,y=avpred,color=fit))+
  geom_linerange(data=subset(ddaf_melt17, fit=="fit.meantemp_high"|fit=="fit.meantemp_low"|fit=="fit.meantemp_mean"),aes(x=week,ymin=avpred-2*sepred,ymax=avpred+2*sepred),alpha=0.3,)+
  #facet_grid(depth~site)+
  theme_classic()

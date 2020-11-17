
#### Visualizations for lme growth models for caging study 
### October 19, 2019
### model selection done in the script "growthrevisited2019.R" in the folder caging_rproj

#### setup ####
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
library("stringr")
library("ggthemes")
library("FSA")
library("purrr")

#keeping everything in this folder
setwd("/Users/tdolan/Documents/R-Github/WFCageStudy")

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
############## mortality estimates ####
# a different catch curve for each site

byweek16 <-ddply(svwk2016, station~week, summarize, catch = n_distinct(fishID))
edays16= c(0,8,14,23,28,36,42,49,58,72)
byweek16 <-byweek16 %>% base::split(.$station) 

#cpd
cpd <-mutate(byweek16$CPD, edays=edays16)
cpd_reg <-catchCurve(catch~edays,data=cpd, ages2use=0:72)
summary(cpd_reg, verbose=TRUE)
#so we could repeat this for every station, but instead, I think it's more meaningful to compare just to cormorant point. 

cp16 <-filter(svwk2016, site=="CP")%>%ddply(~week, summarize, catch = n_distinct(fishID)) %>% mutate(edays=edays16)
cp16_reg <-catchCurve(catch~edays,data=cp16, ages2use=0:72)
summary(cp16_reg, verbose=TRUE)
confint(cp16_reg)

#then in 2017, all the stations are at CP, so use all of them? 
edays17 <-c(0,6,14,22,27,31,35,41,49,55,63)
cp17 <-ddply(svwk2017, ~week, summarize, catch = n_distinct(fishID)) %>% mutate(edays=edays17)
cp17_reg <-catchCurve(catch~edays,data=cp17, ages2use=0:63)
summary(cp17_reg, verbose=TRUE)
confint(cp17_reg)

#get the other model parameters for the table
#function to extract model p.value
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)}

#change out the data
cp17mod<-lm(log(catch)~edays,data = cp17)
summary(cp17mod)

cp16mod<-lm(log(catch)~edays,data = cp16)
summary(cp16mod)

#cagenumber palette
sandycolors <-c("#a6611a", "#dfc27d","#80cdc1","#018571")
cagenumcols <-c("#E69F00","#56B4E9","#009E73","#CC79A7")
dark2 <- c("#1b9e77","#d95f02", "#7570b3", "#e7298a")
###### Visualization no model, just averages ######
svwk2016 <-mutate(svwk2016, cagenum = str_sub(cage,-1,-1))
##number of survivors ###
ggplot(data=svwk2016,aes(x=week,y=survivors, color = cagenum))+
  #geom_point(alpha=1/10)+
  geom_line() +
  scale_color_manual(values=sandycolors)+
  #scale_color_grey()+
  #ylim(0,10)+ 
  scale_y_discrete(limits=c(0,5,10))+
  scale_x_discrete(limits=c(3,6,9),labels=c("3"="30-Jun","6"="20-Jul","9"="10-Aug"))+
  facet_grid(site~depth)+
  theme_few()


#2017
svwk2017 <-mutate(svwk2017, cagenum = str_sub(cage,-1,-1))
##number of survivors ###
ggplot(data=svwk2017,aes(x=week,y=survivors, color = cagenum))+
  #geom_point(alpha=1/10)+
  geom_line() +
  #scale_color_grey()+
  scale_color_manual(values=sandycolors)+
  #ylim(0,10)+ 
  scale_y_discrete(limits=c(0,5,10))+
  scale_x_discrete(limits=c(3,6,9),labels=c("3"="29-Jun","6"="16-Jul","9"="3-Aug"))+
  facet_grid(site~depth)+
  theme_few()

## GROWTH ##
#rugs for number of fish?
#actual fishlength data with the average OVERLAID or not. 
ggplot(data=svwk2016,aes(x=week,y=fishlength))+
  #geom_point(alpha=1/10)+
  ylab("length (mm)")+xlab("")+
  ylim(40,110)+
  geom_rug(alpha = 1/2, position = "jitter", sides="b")+
  scale_x_discrete(limits=c(3,6,9),labels=c("3"="30-Jun","6"="20-Jul","9"="10-Aug"))+
  geom_line(aes(y=av_length, linetype=cagenum)) +
  facet_grid(site~depth)+
  theme_few()
#ggsave("actualgrowth16.png", path="/Users/tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()

ggplot(data=svwk2017,aes(x=week,y=fishlength))+
  #geom_point(alpha=1/10)+
  ylab("length (mm)")+xlab("")+
  ylim(40,110)+
  geom_rug(alpha = 1/2, position = "jitter", sides="b")+
  scale_x_discrete(limits=c(3,6,9),labels=c("3"="29-Jun","6"="16-Jul","9"="3-Aug"))+
  geom_line(aes(y=avlength, linetype=cagenum)) +
  facet_grid(site~depth)+
  theme_few()
#ggsave("actualgrowth17.png", path="/Users/tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()
######## TOP MODELS #####

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

drabcolors <-c("#969696", "#a6bddb","#3690c0","#016450")
boldcolors <-c("#ca0020","blue","#404040","grey")
sandycolors <-c("#a6611a", "#dfc27d","#80cdc1","#018571")


#New Formulation:: 2016
grow16 <-mutate(grow16, cagenum = str_sub(cage,-1,-1))
ggplot(data=svwk2016,aes(x=week,y=fishlength, color =cagenum))+
  geom_rug(alpha = 1/2, position = "jitter", sides="b")+
  scale_color_manual(values=sandycolors)+
  ylab("length (mm)")+xlab("")+
  ylim(40,110)+
  scale_x_discrete(limits=c(3,6,9),labels=c("3"="30-Jun","6"="20-Jul","9"="10-Aug"))+
  geom_line(aes(y=av_length), linetype="solid") + #actual average
  ylab("Average length (mm)")+
  geom_line(data = cbind(grow16, pred = predict(pwin16, level=1)), aes(y = pred),linetype="dashed") + #predicted average
  facet_grid(site~depth)+
  theme_few()
ggsave("actualgrowth16wmodel.png", path="/Users/tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()

#New Formulation:: 2017
grow17 <-mutate(grow17, cagenum = str_sub(cage,-1,-1))
ggplot(data=svwk2017,aes(x=week,y=fishlength, color =cagenum))+
  geom_rug(alpha = 1/2, position = "jitter", sides="b")+
  scale_color_manual(values=sandycolors)+
  ylab("length (mm)")+xlab("")+
  ylim(40,110)+
  scale_x_discrete(limits=c(3,6,9),labels=c("3"="29-Jun","6"="16-Jul","9"="3-Aug"))+
  #geom_point(alpha=1/3, shape=1)+
  geom_line(aes(y=av_length), linetype="dashed", data=grow17) + #actual average
  ylab("Average length (mm)")+
  geom_line(data = cbind(grow17, pred = predict(d73A, level=1)), aes(y = pred),linetype="solid") + #predicted average
  facet_grid(site~depth)+
  theme_few()
ggsave("actualgrowth17wmodel.png", path="/Users/tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()

#Previous formulation::: over fishlength & cage average 2017
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

##### number of fish vs. growth #####
## and growth variation by cage
## growth rate by cage, is it higher where more fish are recovered? Maybe a lrt. 
## and difference between observed and predicted vs. number of fish left. 

### 2016 ##
vargrowth <-ddply(svwk2016,week~cage~cagenum~station, summarize, surv=max(survivors), varlength=sd(fishlength), meanlength=mean(fishlength))

changelength <- function(df){
  deltalength <-c()
  for (i in 2:length(df$meanlength)){
    deltalength[i] <-df$meanlength[i]-df$meanlength[i-1]}
  deltalength[1] <-0
  df <-cbind(df,deltalength)
}

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

vargrowth <-vargrowth %>% base::split(.$cage)%>%
  purrr::map_dfr(changelength)

#weekly growth vs. vs. number of survivors
p2 <- ggplot(vargrowth) + 
  geom_point(aes(x=surv,y=deltalength), position = position_jitter(w = 0)) +
  scale_x_discrete(limits=c(0,2,4,6,8,10))+
  geom_smooth(aes(x=surv,y=deltalength), method="lm") +
  xlab("number of survivors")+ylab("weekly growth (mm)")+
  theme_few()
p2
#ggsave("grow16growthvsurvivors.png", path="/Users/tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()

#weekly growth vs. variance in length
p2 <- ggplot(vargrowth) + 
  geom_point(aes(x=surv,y=varlength), position = position_jitter(w = 0)) +
  scale_x_discrete(limits=c(0,2,4,6,8,10))+
  geom_smooth(aes(x=surv,y=varlength), method="lm") +
  xlab("number of survivors")+ylab("standard deviation of length")+
  theme_few()
p2
#ggsave("grow16SDgrowthvsurvivors.png", path="/Users/tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()

#find initial length and then get overall slope for each cage. 
library("purrr")

models <-svwk2016 %>% base::split(.$cage)%>%
  purrr::map(function(df) lm(fishlength~week,data = df))
slopes <-models %>% map(summary) %>% map_dbl(~.$coefficients[2]) #slope

initial_length <-svwk2016 %>% filter(week==2)%>% base::split(.$cage) %>% map_dfr(summarize, initlenght=mean(fishlength)) 

slopes<-cbind(slopes,initial_length)
slopes <-mutate(slopes, growth=slopes/7)
#weekly growth vs. initial length
p2 <- ggplot(slopes) + 
  geom_point(aes(x=initlenght,y=growth), position = position_jitter(w = 0)) +
  #scale_x_discrete(limits=c(0,2,4,6,8,10))+
  geom_smooth(aes(x=initlenght,y=growth), method="lm") +
  xlab("initial length (mm)")+ylab("average growth (mm/day)")+
  theme_few()
p2
#ggsave("initiallengthvsgrowth16.png", path="/Users/tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()

### 2017 ##
vargrowth <-ddply(svwk2017,week~cage~cagenum~station, summarize, surv=max(survivors), varlength=sd(fishlength), meanlength=mean(fishlength))
vargrowth <-vargrowth %>% base::split(.$cage)%>%
  purrr::map_dfr(changelength)

#weekly growth vs. vs. number of survivors
p2 <- ggplot(vargrowth) + 
  geom_point(aes(x=surv,y=deltalength), position = position_jitter(w = 0)) +
  scale_x_discrete(limits=c(0,2,4,6,8,10))+
  geom_smooth(aes(x=surv,y=deltalength), method="lm") +
  xlab("number of survivors")+ylab("weekly growth (mm)")+
  theme_few()
p2
#ggsave("grow17growthvsurvivors.png", path="/Users/tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()

#weekly growth vs. variance in length
p2 <- ggplot(vargrowth) + 
  geom_point(aes(x=surv,y=varlength), position = position_jitter(w = 0)) +
  scale_x_discrete(limits=c(0,2,4,6,8,10))+
  geom_smooth(aes(x=surv,y=varlength), method="lm") +
  xlab("number of survivors")+ylab("standard deviation of length")+
  theme_few()
p2
#ggsave("grow17SDgrowthvsurvivors.png", path="/Users/tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()


models <-svwk2017 %>% base::split(.$cage)%>%
  purrr::map(function(df) lm(fishlength~week,data = df))
slopes <-models %>% map(summary) %>% map_dbl(~.$coefficients[2]) #slope

initial_length <-svwk2017 %>% filter(week==2)%>% base::split(.$cage) %>% map_dfr(summarize, initlenght=mean(fishlength)) 

slopes<-cbind(slopes,initial_length)
slopes <-mutate(slopes, growth=slopes/7)
#weekly growth vs. initial length
p2 <- ggplot(slopes) + 
  geom_point(aes(x=initlenght,y=growth), position = position_jitter(w = 0)) +
  #scale_x_discrete(limits=c(0,2,4,6,8,10))+
  geom_smooth(aes(x=initlenght,y=growth), method="lm") +
  xlab("initial length (mm)")+ylab("average growth (mm/day)")+
  theme_few()
p2
#ggsave("initiallengthvsgrowth.png", path="/Users/tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()

########## marginal effect of site  #########

#across all sites
avmaxdo <-mean(grow16$max.do)#13.92222
avmintemp <-mean(grow16$min.temp)#19.38018
avtdur <-mean(grow16$temp.dur)#581.1594

#re-run from here
newdf5 <-dplyr::select(grow16, site, week, cage, depth,av_length)

#populate the new dataframe
newdf5$max.do=rep(avmaxdo)
newdf5$min.temp=rep(avmintemp)
newdf5$temp.dur=rep(avtdur)

#prediction level
lev <-1
predall <- predict(pwin16, newdata=newdf5, data=grow16,level=lev)

#predcp
newdf5_cp <-filter(newdf5,site=="CP")
fit.site16cp <- predict(pwin16, newdata=newdf5_cp, data=grow16,level=lev)
newdf5_cp <- cbind(newdf5_cp,fit.site16cp)
Designmatcp <- model.matrix(eval(eval(pwin16$call$fixed)[-2]), newdf5_cp[-ncol(newdf5_cp)])
#compute standard error for predictions
predvarcp <- diag(Designmatcp %*% pwin16$varFix %*% t(Designmatcp))
newdf5_cp$SE <- sqrt(predvarcp) 
newdf5_cp$SE2 <- sqrt(predvarcp+pwin16$sigma^2) 
newdf5_cp <-dplyr::rename(newdf5_cp,fit.site = fit.site16cp)
#MD
newdf5_md <-filter(newdf5,site=="MD")
fit.site16md <- predict(pwin16, newdata=newdf5_md, data=grow16,level=lev)
newdf5_md <- cbind(newdf5_md,fit.site16md)
Designmatmd <- model.matrix(eval(eval(pwin16$call$fixed)[-2]), newdf5_md[-ncol(newdf5_md)])
#compute standard error for predictions
predvarmd <- diag(Designmatmd %*% pwin16$varFix %*% t(Designmatmd))
newdf5_md$SE <- sqrt(predvarmd) 
newdf5_md$SE2 <- sqrt(predvarmd+pwin16$sigma^2)
newdf5_md <-dplyr::rename(newdf5_md,fit.site = fit.site16md)
#RS
newdf5_rs <-filter(newdf5,site=="RS")
fit.site16rs <- predict(pwin16, newdata=newdf5_rs, data=grow16,level=lev)
newdf5_rs <- cbind(newdf5_rs,fit.site16rs)
Designmatrs <- model.matrix(eval(eval(pwin16$call$fixed)[-2]), newdf5_rs[-ncol(newdf5_rs)])
#compute standard error for predictions
predvarrs <- diag(Designmatrs %*% pwin16$varFix %*% t(Designmatrs))
newdf5_rs$SE <- sqrt(predvarrs) 
newdf5_rs$SE2 <- sqrt(predvarrs+pwin16$sigma^2)
newdf5_rs <-dplyr::rename(newdf5_rs,fit.site = fit.site16rs)

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

ggplot(data=ddnewdf6,aes(y=fit.site,x=week))+
  geom_line()+
  ylab("length (mm)")+ xlab("")+
  ylim(35,75)+
  #geom_point(data=newdf6,aes(y=fit.site,x=week),alpha=0.5)+ #this is the predicted data if we held all to 
  #geom_point(data=svwk2016,aes(y=av_length, x=week))+ #this is the actual data
  geom_ribbon(data=ddnewdf6,aes(ymin=fit.site-2*SE2,ymax=fit.site+2*SE2,x=week),fill="grey",inherit.aes=FALSE,alpha=0.2,linetype="dotted",size=0.5)+
  #geom_smooth(data=newdf6,method=lm,se=FALSE)+
  #geom_linerange(data=newdf6,aes(ymin=fit.site-2*SE2,ymax=fit.site+2*SE2,x=week),alpha=0.1,inherit.aes=FALSE)+
  scale_x_discrete(limits=c(3,6,9),labels=c("3"="30-Jun","6"="20-Jul","9"="10-Aug"))+
  facet_wrap(~site)+
  theme_few()
#ggsave("grow16sitemodel.png", path="/Users/tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()
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
allfits <-dplyr::select(allfits, -biom,-min.do,-max.do,-mean.do,-sd.do,-min.temp,-max.temp,-mean.temp,-sd.T,-do.dur,-temp.dur,-mean.sal,-sd.sal,-ssat,-SD,-survivors, -X)

#melt version
af_melt <-pivot_longer(allfits,cols=c("fit.site16","fit.maxdo_high","fit.maxdo_mean","fit.maxdo_low","fit.mintemp_high","fit.mintemp_mean","fit.mintemp_low","fit.tdur_high","fit.tdur_mean","fit.tdur_low"), names_to = "fit")
af_melt$fit <-as.factor(af_melt$fit)

#You could create a mean value for ALL the cages in a site. 
ddaf_melt <-ddply(af_melt,week~site~fit,summarise,avpred=mean(value),sepred = std.error(value))


#TDUR
#with ribbon instead of linerange
ggplot(data=ddaf_melt, aes(x=week,y=avpred, color=fit))+
  ylab("length (mm)")+xlab("")+
  geom_line(data=subset(ddaf_melt, fit=="fit.tdur_high"|fit=="fit.tdur_low"|fit=="fit.tdur_mean"),aes(x=week,y=avpred))+
  geom_ribbon(data=subset(ddaf_melt, fit=="fit.tdur_high"|fit=="fit.tdur_low"|fit=="fit.tdur_mean"),
              aes(x=week,ymin=avpred-2*sepred,ymax=avpred+2*sepred,fill=fit),alpha=0.2,linetype="blank")+
  scale_x_discrete(limits=c(3,6,9),labels=c("3"="30-Jun","6"="20-Jul","9"="10-Aug"))+
  facet_wrap(~site)+
  theme_few()
#ggsave("growmargin16TDUR.png", path="/Users/tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()

#MINTEMP
ggplot(data=ddaf_melt, aes(x=week,y=avpred, color=fit))+
  ylab("length (mm)")+xlab("")+
  geom_line(data=subset(ddaf_melt, fit=="fit.mintemp_high"|fit=="fit.mintemp_low"|fit=="fit.mintemp_mean"),aes(x=week,y=avpred))+
  geom_ribbon(data=subset(ddaf_melt, fit=="fit.mintemp_high"|fit=="fit.mintemp_low"|fit=="fit.mintemp_mean"),
              aes(x=week,ymin=avpred-2*sepred,ymax=avpred+2*sepred,fill=fit),alpha=0.2,linetype="blank")+
  scale_x_discrete(limits=c(3,6,9),labels=c("3"="30-Jun","6"="20-Jul","9"="10-Aug"))+
  facet_wrap(~site)+
  theme_few()
#ggsave("growmargin16MINTEMP.png", path="/Users/tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()

#MAXDO
ggplot(data=ddaf_melt, aes(x=week,y=avpred, color=fit))+
  ylab("length (mm)")+xlab("")+
  geom_line(data=subset(ddaf_melt, fit=="fit.maxdo_high"|fit=="fit.maxdo_low"|fit=="fit.maxdo_mean"),aes(x=week,y=avpred))+
  geom_ribbon(data=subset(ddaf_melt, fit=="fit.maxdo_high"|fit=="fit.maxdo_low"|fit=="fit.maxdo_mean"),
              aes(x=week,ymin=avpred-2*sepred,ymax=avpred+2*sepred,fill=fit),alpha=0.2,linetype="blank")+
  scale_x_discrete(limits=c(3,6,9),labels=c("3"="30-Jun","6"="20-Jul","9"="10-Aug"))+
  facet_wrap(~site)+
  theme_few()
#ggsave("growmargin16MAXDO.png", path="/Users/tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()

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
allfits17 <-dplyr::select(allfits17, -biom,-min.do,-max.do,-mean.do,-sd.do,-min.temp,-max.temp,-mean.temp,-sd.T,-do.dur,-temp.dur,-mean.sal,-sd.sal,-ssat,-SD,-survivors, -X)

#melt version
af_melt17 <-pivot_longer(allfits17,cols=c("fit.meantemp_high","fit.meantemp_mean","fit.meantemp_low"), names_to = "fit")
af_melt17$fit <-as.factor(af_melt17$fit)

#You could create a mean value for ALL the cages in a site. 
ddaf_melt17 <-ddply(af_melt17,week~fit,summarise,avpred=mean(value),sepred = std.error(value))

ggplot(data=ddaf_melt17, aes(x=week,y=avpred, color=fit))+
  #scale_color_manual(values=wes_palette(n=3, name="Darjeeling1"))+
  ylab("Average length (mm)")+
  geom_line(data=subset(ddaf_melt17, fit=="fit.meantemp_high"|fit=="fit.meantemp_low"|fit=="fit.meantemp_mean"),aes(x=week,y=avpred,color=fit))+
  geom_linerange(data=subset(ddaf_melt17, fit=="fit.meantemp_high"|fit=="fit.meantemp_low"|fit=="fit.meantemp_mean"),aes(x=week,ymin=avpred-2*sepred,ymax=avpred+2*sepred),alpha=0.3,)+
  #facet_grid(depth~site)+
  theme_classic()

#MEANTEMP17
ggplot(data=ddaf_melt17, aes(x=week,y=avpred, color=fit))+
  ylab("length (mm)")+xlab("")+
  geom_line(data=subset(ddaf_melt17, fit=="fit.meantemp_high"|fit=="fit.meantemp_low"|fit=="fit.meantemp_mean"),aes(x=week,y=avpred,color=fit))+
  geom_ribbon(data=subset(ddaf_melt17, fit=="fit.meantemp_high"|fit=="fit.meantemp_low"|fit=="fit.meantemp_mean"),
              aes(x=week,ymin=avpred-2*sepred,ymax=avpred+2*sepred,fill=fit),alpha=0.2,linetype="blank")+
  scale_x_discrete(limits=c(3,6,9),labels=c("3"="29-Jun","6"="16-Jul","9"="3-Aug"))+
  theme_few()
#ggsave("growmargin17MEANTEMP.png", path="/Users/tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()

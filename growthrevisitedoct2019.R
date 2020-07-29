#Growth revisited October 9, 2019
# taking the rmd version of this and making it into a script

#this is only for model selection. Graphing will be done in growthgraphsOct2019.R
#see also saved calculations_1_15_19 in caging_rproj

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

#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github/WFCageStudy")

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

## Visualization no model, just averages ##
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

####VIF function code####
#https://stackoverflow.com/questions/26633483/collinearity-after-accounting-for-random-mixed-effects

#2016
vif.lme16 <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)] }
  #exclude week
  wk <- sum(1 * (nam == "week" | nam == "week"))
  if (wk > 0) {
    v <- v[-(1:wk), -(1:wk), drop = FALSE]
    nam <- nam[-(1:wk)] }
  #exclude categorical variable siteMD
  md <- sum(1 * (nam == "siteMD" | nam == "siteMD"))
  if (md > 0) {
    v <- v[-(1:md), -(1:md), drop = FALSE]
    nam <- nam[-(1:md)] }
  #exclude categorical variable siteRS
  rs <- sum(1 * (nam == "siteRS" | nam == "siteRS"))
  if (rs > 0) {
    v <- v[-(1:rs), -(1:rs), drop = FALSE]
    nam <- nam[-(1:rs)] }
  #exclude categorical variable depthShallow
  ds <- sum(1 * (nam == "depthshallow" | nam == "depthshallow"))
  if (ds > 0) {
    v <- v[-(1:ds), -(1:ds), drop = FALSE]
    nam <- nam[-(1:ds)] }
  #exclude categorical variable week:siteMD
  wmd <- sum(1 * (nam == "siteMD:week" | nam == "siteMD:week"))
  if (wmd > 0) {
    v <- v[-(1:wmd), -(1:wmd), drop = FALSE]
    nam <- nam[-(1:wmd)] }
  #exclude categorical variable week:siteRS
  wrs <- sum(1 * (nam == "siteRS:week" | nam == "siteRS:week"))
  if (wrs > 0) {
    v <- v[-(1:wrs), -(1:wrs), drop = FALSE]
    nam <- nam[-(1:wrs)] }
  #exclude categorical variable week:depthshallow
  wds <- sum(1 * (nam == "depthshallow:week" | nam == "depthshallow:week"))
  if (wds > 0) {
    v <- v[-(1:wds), -(1:wds), drop = FALSE]
    nam <- nam[-(1:wds)] }
  #exclude categorical variable siteMD:depthshallow
  smdds <- sum(1 * (nam == "siteMD:depthshallow" | nam == "siteMD:depthshallow"))
  if (smdds > 0) {
    v <- v[-(1:smdds), -(1:smdds), drop = FALSE]
    nam <- nam[-(1:smdds)] }
  #exclude categorical variable siteRS:depthshallow
  srsds <- sum(1 * (nam == "siteRS:depthshallow" | nam == "siteRS:depthshallow"))
  if (srsds > 0) {
    v <- v[-(1:srsds), -(1:srsds), drop = FALSE]
    nam <- nam[-(1:srsds)] }
  #exclude categorical variable week:siteMD:depthshallow
  wksmdds <- sum(1 * (nam == "siteMD:depthshallow:week" | nam == "siteMD:depthshallow:week"))
  if (wksmdds > 0) {
    v <- v[-(1:wksmdds), -(1:wksmdds), drop = FALSE]
    nam <- nam[-(1:wksmdds)] }
  #exclude categorical variable week:siteRS:depthshallow
  wksrsds <- sum(1 * (nam == "siteRS:depthshallow:week" | nam == "siteRS:depthshallow:week"))
  if (wksrsds > 0) {
    v <- v[-(1:wksrsds), -(1:wksrsds), drop = FALSE]
    nam <- nam[-(1:wksrsds)] } 
  #I think i got all the categorical variables
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v }


#2017
#### VIF.lme function
vif.lme17 <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)] }
  #exclude week
  wk <- sum(1 * (nam == "week" | nam == "week"))
  if (wk > 0) {
    v <- v[-(1:wk), -(1:wk), drop = FALSE]
    nam <- nam[-(1:wk)] }
  #exclude categorical variable siteMD
  md <- sum(1 * (nam == "siteN" | nam == "siteN"))
  if (md > 0) {
    v <- v[-(1:md), -(1:md), drop = FALSE]
    nam <- nam[-(1:md)] }
  #exclude categorical variable siteRS
  rs <- sum(1 * (nam == "siteS" | nam == "siteS"))
  if (rs > 0) {
    v <- v[-(1:rs), -(1:rs), drop = FALSE]
    nam <- nam[-(1:rs)] }
  #exclude categorical variable depthShallow
  ds <- sum(1 * (nam == "depthshallow" | nam == "depthshallow"))
  if (ds > 0) {
    v <- v[-(1:ds), -(1:ds), drop = FALSE]
    nam <- nam[-(1:ds)] }
  #exclude categorical variable week:siteMD
  wmd <- sum(1 * (nam == "siteN:week" | nam == "siteN:week"))
  if (wmd > 0) {
    v <- v[-(1:wmd), -(1:wmd), drop = FALSE]
    nam <- nam[-(1:wmd)] }
  #exclude categorical variable week:siteRS
  wrs <- sum(1 * (nam == "siteS:week" | nam == "siteS:week"))
  if (wrs > 0) {
    v <- v[-(1:wrs), -(1:wrs), drop = FALSE]
    nam <- nam[-(1:wrs)] }
  #exclude categorical variable week:depthshallow
  wds <- sum(1 * (nam == "depthshallow:week" | nam == "depthshallow:week"))
  if (wds > 0) {
    v <- v[-(1:wds), -(1:wds), drop = FALSE]
    nam <- nam[-(1:wds)] }
  #exclude categorical variable siteMD:depthshallow
  smdds <- sum(1 * (nam == "siteN:depthshallow" | nam == "siteN:depthshallow"))
  if (smdds > 0) {
    v <- v[-(1:smdds), -(1:smdds), drop = FALSE]
    nam <- nam[-(1:smdds)] }
  #exclude categorical variable siteRS:depthshallow
  srsds <- sum(1 * (nam == "siteS:depthshallow" | nam == "siteS:depthshallow"))
  if (srsds > 0) {
    v <- v[-(1:srsds), -(1:srsds), drop = FALSE]
    nam <- nam[-(1:srsds)] }
  #exclude categorical variable week:siteMD:depthshallow
  wksmdds <- sum(1 * (nam == "siteN:depthshallow:week" | nam == "siteN:depthshallow:week"))
  if (wksmdds > 0) {
    v <- v[-(1:wksmdds), -(1:wksmdds), drop = FALSE]
    nam <- nam[-(1:wksmdds)] }
  #exclude categorical variable week:siteRS:depthshallow
  wksrsds <- sum(1 * (nam == "siteS:depthshallow:week" | nam == "siteS:depthshallow:week"))
  if (wksrsds > 0) {
    v <- v[-(1:wksrsds), -(1:wksrsds), drop = FALSE]
    nam <- nam[-(1:wksrsds)] } 
  #I think i got all the categorical variables
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v }
##########

#####Create the saturated models #####
##VIF table creation##
#First, create a saturated model for both years and do VIF on all of them. Then look at the average vif across years. 
# R won't run the full models, so conduct vif using lm versions, see if that works. 

#ok try this:
sat16 <-lme(av_length ~ site + depth + min.do + mean.do + max.do + min.temp + max.temp + mean.temp + do.dur + temp.dur + mean.sal + ssat ,random=~week|cage, data=grow16, method="ML")




sat17 <-lme(av_length ~ site + depth + min.do + mean.do + max.do + min.temp + max.temp + mean.temp + do.dur + temp.dur + mean.sal + ssat ,random=~week|cage, data=grow17, method="ML")

d16 <-dredge(sat16)
fullwin16 <-lme(av_length ~ site + max.do + min.temp + temp.dur,random=~week|cage, data=grow16, method="ML")
fw16.2 <-lme(av_length ~ site + max.do + min.temp,random=~week|cage, data=grow16, method="ML")
fw16.3 <-lme(av_length ~ site + max.do + min.temp + temp.dur + min.do ,random=~week|cage, data=grow16, method="ML")
d17 <- dredge(sat17)
fullwin17 <-lme(av_length ~  max.do + mean.do + mean.temp + min.do,random=~week|cage, data=grow17, method="ML")
fw17.2 <-lme(av_length ~ site + max.do + mean.do + mean.temp + min.do,random=~week|cage, data=grow17, method="ML")
fw17.3 <-lme(av_length ~ mean.temp,random=~week|cage, data=grow17, method="ML")

vif16 <- as.data.frame(vif.lme16(sat16)) 
vif17 <- as.data.frame(vif.lme17(sat17)) %>% add_rownames(var="var")
vifs <-bind_cols(vif16, vif17) 
vifs <- mutate(vifs, avgvif=((vif.lme16(sat16)+vif.lme17(sat17))/2)) 
vifs

#mean temp has the highest vif in all which means that it's correlated with the most things, probably as the summer progresses, more fish die,
## therefore it's a good predictor, but not a good explanatory variable. so let's take that out first. 
### Although it's vif is not very high for summer 2016! So I think theres a good argument for leaving it in for 2016
sat16a <-update(sat16, .~.-mean.temp)
sat17a <-update(sat17, .~.-mean.temp)
vif16 <- as.data.frame(vif.lme16(sat16a)) 
vif17 <- as.data.frame(vif.lme17(sat17a)) %>% add_rownames(var="var")
vifs <-bind_cols(vif16, vif17) 
vifs <- mutate(vifs, avgvif=((vif.lme16(sat16a)+vif.lme17(sat17a))/2)) 
vifs

#next remove mean.do - vifs are still high. We are aiming for 5. 
sat16b <-update(sat16a, .~.-mean.do)
sat17b <-update(sat17a, .~.-mean.do)
vif16 <- as.data.frame(vif.lme16(sat16b)) 
vif17 <- as.data.frame(vif.lme17(sat17b)) %>% add_rownames(var="var")
vifs <-bind_cols(vif16, vif17) 
vifs <- mutate(vifs, avgvif=((vif.lme16(sat16b)+vif.lme17(sat17b))/2)) 
vifs

#All are well under 5! I think start here for actual model selection. 

#2016
d6 <-dredge(sat16b) #dredge it.
d6 <-as.data.frame(d6)
#2016 dredge top 3 models: 
dg16A <-lme(av_length ~ site + max.do + min.temp + max.temp + temp.dur ,random=~week|cage, data=grow16, method="ML") #1134.1
dg16B <-lme(av_length ~ site + max.do + min.temp + max.temp ,random=~week|cage, data=grow16, method="ML") #1135.1
dg16C <-lme(av_length ~ site + max.do + min.temp + max.temp + min.do ,random=~week|cage, data=grow16, method="ML") #1135.4
#write.csv(d6,file="growthdredgewin16.csv")
fits16 <-stepAIC(sat16b, direction="both", k=2) #step it
summary(fits16)
# 2016 top model from Step
step16A <-lme(av_length ~ site + max.do + min.temp + temp.dur ,random=~week|cage, data=grow16, method="ML") #1134.068
#this is the same model as pwin16

#Now perform VIF process from the full model separately. 
vif.lme16(sat16) #remove mean.do first. 
sat16c <-update(sat16, .~.-mean.do)
vif.lme16(sat16c) #fully under 5... 
ds16 <-dredge(sat16c)
#this generally does not result in better models than the VIF together. so keep vif together. 

#Look at winning models from previous attempts before max.do was corrected. Are they better?
pwin16 <-lme(av_length ~ site + max.do + min.temp + temp.dur,random=~week|cage, data=grow16, method="ML") #1134.068
pwin16.2 <-lme(av_length ~ site + max.do + min.temp ,random=~week|cage, data=grow16, method="ML") #1135.092
pwin16.3 <- lme(av_length ~ site + max.do + min.temp + min.do + temp.dur,random=~week|cage, data=grow16, method="ML") #1135.386
AICc(pwin16,pwin16.2,pwin16.3)

#try eliminating by categories
temp16 <- lme(av_length ~ site + depth + max.temp + mean.temp  + min.temp + temp.dur,random=~week|cage, data=grow16, method="ML")
summary(temp16)
dt16 <-dredge(temp16)
# i think it makes sense to have both max temp and min temp in the model... that seems fine to me. mean temp seems less important. 
do16 <-lme(av_length ~ site + depth + max.do + min.do + mean.do + do.dur + ssat,random=~week|cage, data=grow16, method="ML")
summary(do16)
# to me this argues that min.do is not a good model term, and mean.do should be used. 

hyp16 <-lme(av_length ~ site + depth + mean.do + min.temp + max.temp + mean.sal ,random=~week|cage, data=grow16, method="ML")
dh16 <-dredge(hyp16)
hwin16 <-lme(av_length ~ site + max.temp + min.temp ,random=~week|cage, data=grow16, method="ML")

#Another way is to look at the coefficients of the saturated model and remove ones that have a sign that doesn't make sense. 




summary(sat17)
wd17 <-lme(av_length ~ site + depth + min.do + max.do + max.temp + mean.temp + mean.sal + ssat ,random=~week|cage, data=grow17, method="ML")
dwd17 <-dredge(wd17, m.max =3) #here we limit the model to three terms. 
dwd17A <-lme(av_length ~ mean.temp ,random=~week|cage, data=grow17, method="ML")
dwd17B <-lme(av_length ~ mean.temp + min.do ,random=~week|cage, data=grow17, method="ML")
dwd17C <-lme(av_length ~ mean.temp + depth ,random=~week|cage, data=grow17, method="ML")
dwd17.nolimits <-dredge(wd17) #here we do not limit the model to 3 terms
dwdnl17A <-lme(av_length ~ mean.temp + max.do + min.do + ssat ,random=~week|cage, data=grow17, method="ML")
dwdnl17B <-lme(av_length ~ mean.temp + max.do + min.do + ssat + site ,random=~week|cage, data=grow17, method="ML")
dwdnl17C <-lme(av_length ~ mean.temp  + min.do  ,random=~week|cage, data=grow17, method="ML")


#Going by AICc alone, compare the models to find the top 3: 
faceofflme16 <- as.data.frame(model.sel(dg16A,dg16B,dg16C,step16A,pwin16,pwin16.2,pwin16.3,hwin16,fullwin16,fw16.2, fw16.3, b16,b16.2,b16.3))
#So all the pwins are still the top 3 models! But their coefficents might be slightly different than last time so you still have to re-run the things. 
##write.csv(faceofflme16,file="faceofflme16.csv")

#2017
d7 <-dredge(sat17b) #dredge it.
dg17A <-lme(av_length ~  max.temp + min.temp,random=~week|cage, data=grow17, method="ML")
dg17B <-lme(av_length ~  min.do + min.temp + temp.dur ,random=~week|cage, data=grow17, method="ML")
dg17c <-lme(av_length ~  max.temp + min.do + min.temp,random=~week|cage, data=grow17, method="ML")

d7.limits3 <-dredge(sat17,m.max=3)
d73A <-lme(av_length ~ mean.temp ,random=~week|cage, data=grow17, method="ML")
d73B <-lme(av_length ~ mean.temp + min.do ,random=~week|cage, data=grow17, method="ML")
d73C <-lme(av_length ~ mean.temp + depth ,random=~week|cage, data=grow17, method="ML")

#Try some interactions
d73BInt <-lme(av_length ~ mean.temp*depth + min.do ,random=~week|cage, data=grow17, method="ML")

#write.csv(d6,file="growthdredgewin16.csv")
fits17 <-stepAIC(sat17b, direction="both", k=2) #step it
summary(fits17)

#previous winners from 2017
pwin17 <- lme(av_length ~  min.do + mean.do + max.do + mean.temp,random=~week|cage, data=grow17, method="ML")
pwin17.2 <- lme(av_length ~ site + min.do + mean.do + max.do + mean.temp,random=~week|cage, data=grow17, method="ML")
pwin17.3 <- lme(av_length ~ depth +  min.do + mean.do + max.do + mean.temp,random=~week|cage, data=grow17, method="ML")

faceofflme17 <-model.sel(d73A,d73B,d73C,dg17A,dg17B,dg17c, pwin17,pwin17.2, pwin17.3, fits17, fullwin17, fw17.2, fw17.3, dwd17A,dwd17B,dwd17C, dwdnl17A,dwdnl17B,dwdnl17C, b17, b17.2, b17.3)

##write.csv(faceofflme17, file="faceofflme17.csv")
#there are many different model selection procedures to try, but ultimately there are a few rules:
# no models with more than 2 terms in the same "category"
# no models where the coefficient is in the wrong direction
# no models with more than 4 terms total! too complicated. 

#if we follow our rules, and dredge from the full model (m.max=4), we get:
b16 <- lme(av_length ~ site + max.do + min.temp + temp.dur,random=~week|cage, data=grow16, method="ML")
b16.2 <- lme(av_length ~ site + max.do + min.temp,random=~week|cage, data=grow16, method="ML")
b16.3 <- lme(av_length ~ site + max.do + mean.temp + temp.dur,random=~week|cage, data=grow16, method="ML")

#2017
b17 <-lme(av_length ~ mean.temp,random=~week|cage, data=grow17, method="ML")
b17.2 <- lme(av_length ~  min.do + mean.temp,random=~week|cage, data=grow17, method="ML")
b17.3 <-lme(av_length ~  depth + mean.temp, random=~week|cage, data=grow17, method="ML")

#make a function that plots the model and prints AIC
mod16 <-function(m){
  newdf <-dplyr::select(grow16, cage,av_length,depth,week,min.do, max.do, mean.do, min.temp, max.temp, mean.temp, do.dur, temp.dur, mean.sal, ssat)
  newdf$pred <-predict(r16, newdata=newdf, level=1)
  p <- ggplot(newdf,aes(x=week,y=pred)) + 
    geom_point(data=newdf,aes(x=week,y=pred),fill="red",shape=23) +
    geom_point(data=svwk2016,aes(x=week,y=av_length),alpha=0.2,fill="gray") +
    facet_wrap(~cage)+
    theme_pubr()
  return(list(AICc(m),p))}


#Following the tutorial. 
#the AICc to beat is 1134.068

#first fit a random model on cages
#https://rpsychologist.com/r-guide-longitudinal-lme-lmer
r16 <- lmer(av_length~ 1 + (1|cage), data=grow16)
mod16(r16) #1490.375, fit is not great. 
#add the element of time
t16 <- lmer(av_length~ week + (week|cage), data=grow16)  # has already beaten our previous AICc (1128.441) I would call this the "fully random model" 
mod16(t16)
# conditional growth (site only)
cg16 <- lmer(av_length~ week*site + (week|cage), data=grow16) # even better. 1105.768
mod16(cg16)
# conditional growth random intercept only
ri16 <- lmer(av_length~ week*site + (1|cage), data=grow16) # i'm not sure the plot function is working. AICc 1205.408
mod16(ri16) 
# conditional growth random random slope only
rs16 <- lmer(av_length~ week*site + (0 + week|cage), data=grow16) # 1131.004
mod16(rs16)
# conditional growth, dropping intercept-slope covariance
nocovar16 <- lmer(av_length~ week*site + (week||cage), data=grow16) #1124.193
mod16(nocovar16)

# do we have to add an interaction effect of week on every term? Does it understand that week is time??? I think we went over this. 

#Is it even possible to add additional terms to the model using this format? In this format week is in the model interacting with site. 

#lets try with our best model from before. 
pwin16 <-lme(av_length ~ site + max.do + min.temp + temp.dur,random=~week|cage, data=grow16, method="ML")
pwin16.ws <-lme(av_length ~ week*site + max.do + min.temp + temp.dur,random=~week|cage, data=grow16, method="ML")
pwin16.wall <-lme(av_length ~ week*(site + max.do + min.temp + temp.dur),random=~week|cage, data=grow16, method="ML")
AICc(pwin16,pwin16.ws, pwin16.wall) #wall is best!

#try changing the models and changing the predition level. 
#plot the model over the data
# in your weird sawtooth graphs, you solve that by making a separate line for each cage, or by taking some kind of average of each prediction and plotting that. 
newdf <-dplyr::select(grow16, site, cage,av_length,depth,week,min.do, max.do, mean.do, min.temp, max.temp, mean.temp, do.dur, temp.dur, mean.sal, ssat)
newdf$pred <-predict(pwin16.wall, newdata=newdf, level=1)
  p <- ggplot(newdf,aes(x=week,y=pred)) + 
    geom_line(color="red") +
    geom_point(data=svwk2016,aes(x=week,y=av_length),alpha=0.2,fill="gray") +
    facet_wrap(~cage)+
    #facet_grid(site~depth)+
    theme_pubr()
  p







#####################################################

#To find the top models go to the worksheet "saved calculations_1_15_19" in the caging_rproj folder  
#there were a lot of good models and we didn't just choose the best by AICc see the xcel sheet. 
  
#2016 top models
  pwin16 <-lme(av_length ~ site + max.do + min.temp + temp.dur,random=~week|cage, data=grow16, method="ML") 
  pwin16.2 <-lme(av_length ~ site + max.do + min.temp ,random=~week|cage, data=grow16, method="ML") 
  b16.3 <- lme(av_length ~ site + max.do + mean.temp + temp.dur,random=~week|cage, data=grow16, method="ML")
  win16 <-model.sel(pwin16,pwin16.2,b16.3)
  
#2017 top models 
  d73A <-lme(av_length ~ mean.temp ,random=~week|cage, data=grow17, method="ML")
  dwdnl17A <-lme(av_length ~ mean.temp + max.do + min.do + ssat ,random=~week|cage, data=grow17, method="ML")
  d73B <-lme(av_length ~ mean.temp + min.do ,random=~week|cage, data=grow17, method="ML")
  win17 <-model.sel(d73A,dwdnl17A,d73B)
#the graphing will be done in a new script called "growthgraphsOct2019.R"
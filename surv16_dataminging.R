# Survival analysis 2016 making the data:

library("survival")
library("car")
library("ggplot2")
library("reshape2")
library("plyr")
library("moments")
library("reshape")
library("tidyr")
library("dplyr")




#keeping everything in this folder
setwd("/Users//tdolan/Documents/R-Github/WFCageStudy")

six <-read.csv("Indexed_envirodata16_v2.csv")

###################### sal data minging ######
#First get missing data 
new.sal <-dplyr::select(six,week, index, contains("sal"))
new.sal <-slice(new.sal, 4:7585) # 6/13/2016 15:15 on 6/15 to 8/31/2016 14:30.
new.sal <- as.data.frame(new.sal)
new.sal[new.sal=="#N/A"] <-NA#replace #N/A with actual NA
new.sal[new.sal==".."] <-NA

###################### First actually you want to clean up the data. Impute data where it's an na.######
#CPD
for(i in 1:length(new.sal$CPD.sal)){
  if(is.na(new.sal$CPD.sal[i])){
    new.sal$CPD.sal[i]<-sample(na.omit(new.sal$CPD.sal),size=1,replace=TRUE)}}
#CPS
for(i in 1:length(new.sal$CPS.sal)){
  if(is.na(new.sal$CPS.sal[i])){
    new.sal$CPS.sal[i]<-sample(na.omit(new.sal$CPS.sal),size=1,replace=TRUE)}}
#RSD
for(i in 1:length(new.sal$RSD.sal)){
  if(is.na(new.sal$RSD.sal[i])){
    new.sal$RSD.sal[i]<-sample(na.omit(new.sal$RSD.sal),size=1,replace=TRUE) }}

#CPD
for(i in 1:length(new.sal$CPD.sal)){
  new.sal$CPD.sal[i]<-10.7592+0.7130*new.sal$CPD.sal[i]}
#plot(new.sal$index, new.sal$CPD.sal, main="CPD Salinity",ylab="ppt")

###################### CPS##########
#first correct the drops and zeros. 
for(i in 1:length(new.sal$CPS.sal)){
  if(new.sal$CPS.sal[i]< 23){
    new.sal$CPS.sal[i]<-10.157950+0.531815*new.sal$CPD.sal[i]}}

#now de trend with the correction factor - which is the relationship between the cage and field at CPD, not the relationship between CPS and corrected CPD. 
for(i in 1:length(new.sal$CPS.sal)){
  new.sal$CPS.sal[i]<-10.7592+0.7130*new.sal$CPS.sal[i]}
#plot(new.sal$index,new.sal$CPS.sal,main="CPS salinity", ylab="ppt")
#this looks reasonable. 

###################### RSD##################
#for RSD, find out where it started to go wrong, only impute from CPD after that point. 
for(i in 1:length(new.sal$RSD.sal)){
  new.sal$RSD.sal[i]<-10.7592+0.7130*new.sal$RSD.sal[i]}
for(i in 1:length(new.sal$RSD.sal)){
  if(new.sal$RSD.sal[i] < 15){
    new.sal$RSD.sal[i]<-sample(na.omit(new.sal$RSD.sal),size=1,replace=TRUE) }}

for(i in 1:length(new.sal$RSD.sal)){
  if(new.sal$RSD.sal[i] < 25){
    new.sal$RSD.sal[i]<-1.43597+0.85076*new.sal$CPD.sal[i]}}

###################### RSS ########
#impute the rest of the summer past week 2 from RSD
for(i in 1:length(new.sal$RSS.sal)){
  if(is.na(new.sal$RSS.sal[i])){
    new.sal$RSS.sal[i]<-19.88741+0.37552*new.sal$RSD.sal[i]}}

for(i in 1:length(new.sal$RSS.sal)){
  if(new.sal$RSS.sal[i] < 25){
    new.sal$RSS.sal[i]<-19.88741+0.37552*new.sal$RSD.sal[i]}}

###################### MDD#################################
#first generally apply the correction factor, then re-plot
for(i in 1:length(new.sal$MDD.sal)){
  new.sal$MDD.sal[i]<-10.7592+0.7130*new.sal$MDD.sal[i]}

#the whole thing is a piece of shit and has to come from CPD.
for(i in 1:length(new.sal$MDD.sal)){ #first get rid of the goddamn NA
  if(is.na(new.sal$MDD.sal[i])){
    new.sal$MDD.sal[i]<-sample(na.omit(new.sal$MDD.sal),size=1,replace=TRUE)}} #this one was a mess. 
for(i in 1:length(new.sal$MDD.sal)){  # then impute from CPD
  new.sal$MDD.sal[i]<--17.61644 + 1.41582*new.sal$CPD.sal[i]}
#Then do we re-apply the correction factor? I think not.

###################### finally, MDS #############################
new.sal$MDS.sal <- as.numeric(as.character(new.sal$MDS.sal))
for(i in 1:length(new.sal$MDS.sal)){
  if(is.na(new.sal$MDS.sal[i])){
    new.sal$MDS.sal[i]<-sample(na.omit(new.sal$MDS.sal),size=1,replace=TRUE)}}
#apply correction factor
for(i in 1:length(new.sal$MDS.sal)){
  new.sal$MDS.sal[i]<-10.7592+0.7130*new.sal$MDS.sal[i]}
new.sal$MDS.sal <- as.numeric(as.character(new.sal$MDS.sal))
#correct the zeros & weird numbers
for(i in 1:length(new.sal$MDS.sal)){
  if(new.sal$MDS.sal[i] < 24){
    new.sal$MDS.sal[i]<-4.0425671+0.823900*new.sal$MDD.sal[i]}}

###################### sal gather ##### 
sal.gather <-melt(new.sal, id.vars="index",measure.vars=c("CPD.sal","CPS.sal","RSD.sal","RSS.sal", "MDD.sal","MDS.sal"))
levels(sal.gather$variable) <- c("CPD","CPS","RSD","RSS","MDD","MDS") #renaming levels will help us combine data later
sal.gather <-filter(sal.gather, value > 20, value) #subsetting between 0 and 23, because outside of this range isn't plausible. 

################################### TEMPERATURE ####################################
new.T <-dplyr::select(six,week, index, contains(".T"))
new.T <-as.data.frame(slice(new.T, 4:7585)) # 6/13/2016 15:15 on 6/15 to 8/31/2016 14:30.
new.T[new.T=="#N/A"] <-NA
new.T[new.T==".."] <-NA

#fix NA
#CPD
for(i in 1:length(new.T$CPD.T)){
  if(is.na(new.T$CPD.T[i])){
    new.T$CPD.T[i]<-sample(na.omit(new.T$CPD.T),size=1,replace=TRUE)}}
new.T$CPD.T <- as.numeric(as.character(new.T$CPD.T))

#CPS
for(i in 1:length(new.T$CPS.T)){
  if(is.na(new.T$CPS.T[i])){
    new.T$CPS.T[i]<-sample(na.omit(new.T$CPS.T),size=1,replace=TRUE)}}
new.T$CPS.T <- as.numeric(as.character(new.T$CPS.T))

#RSD
for(i in 1:length(new.T$RSD.T)){
  if(is.na(new.T$RSD.T[i])){
    new.T$RSD.T[i]<-sample(na.omit(new.T$RSD.T),size=1,replace=TRUE) }}
new.T$RSD.T <- as.numeric(as.character(new.T$RSD.T))

#RSS
for(i in 1:length(new.T$RSS.T)){
  if(is.na(new.T$RSS.T[i])){
    new.T$RSS.T[i]<-sample(na.omit(new.T$RSS.T),size=1,replace=TRUE) }}
new.T$RSS.T <- as.numeric(as.character(new.T$RSS.T))

#MDD
for(i in 1:length(new.T$MDD.T)){
  if(is.na(new.T$MDD.T[i])){
    new.T$MDD.T[i]<-sample(na.omit(new.T$MDD.T),size=1,replace=TRUE) }}
new.T$MDD.T <- as.numeric(as.character(new.T$MDD.T))

#MDS
for(i in 1:length(new.T$MDS.T)){
  if(is.na(new.T$MDS.T[i])){
    new.T$MDS.T[i]<-sample(na.omit(new.T$MDS.T),size=1,replace=TRUE) }}
new.T$MDS.T <- as.numeric(as.character(new.T$MDS.T))

#Temperature timeseries
windows(30,20) #change to quartz
par(mfrow=c(2,3))
par(mar=c(4,4,1,1))
#MDD.T
plot(new.T$index,new.T$MDD.T, ylim=c(10,30),ylab="(C)",main="MDD.T", cex=2)
abline(h=18.5)
abline(h=27)
#MDS.T
plot(new.T$index,new.T$MDS.T, ylim=c(10,30),ylab="(C)",main="MDS.T", cex=2)
abline(h=18.5)
abline(h=27)
#CPD.T
plot(new.T$index,new.T$CPD.T, ylim=c(10,30),ylab="(C)",main="CPD.T", cex=2)
abline(h=18.5)
abline(h=27)
#CPS.T
plot(new.T$index,new.T$CPS.T, ylim=c(10,30),ylab="(C)",main="CPS.T", cex=2)
abline(h=18.5)
abline(h=27)
#RSD.T
plot(new.T$index,new.T$RSD.T, ylim=c(10,30),ylab="(C)",main="RSD.T", cex=2)
abline(h=18.5)
abline(h=27)
#RSS.T
plot(new.T$index,new.T$RSS.T, ylim=c(10,30),ylab="(C)",main="RSS.T", cex=2)
abline(h=18.5)
abline(h=27)




######## Temperature Violin Plots#####

#gather
T.gather <-melt(new.T, id.vars="index",measure.vars=c("CPD.T","CPS.T","RSD.T","RSS.T", "MDD.T","MDS.T"))
levels(T.gather$variable) <- c("CPD","CPS","RSD","RSS","MDD","MDS") #renaming levels will help us combine data later
T.gather <-filter(T.gather, value > 10, value) #subsetting between 0 and 10, because outside of this range isn't plausible. 

################################### DISSOLVED OXYGEN ############

new.DO <-dplyr::select(six,week, index, contains("DO"))
new.DO <-as.data.frame(slice(new.DO, 4:7585)) # 6/13/2016 15:15 on 6/15 to 8/31/2016 14:30.
new.DO[new.DO=="#N/A"] <-NA
new.DO[new.DO==".."] <-NA


#fix NA
#CPD
for(i in 1:length(new.DO$CPDDO)){
  if(is.na(new.DO$CPDDO[i])){
    new.DO$CPDDO[i]<-sample(na.omit(new.DO$CPDDO),size=1,replace=TRUE)}}
new.DO$CPDDO <- as.numeric(as.character(new.DO$CPDDO))

#CPS
for(i in 1:length(new.DO$CPSDO)){
  if(is.na(new.DO$CPSDO[i])){
    new.DO$CPSDO[i]<-sample(na.omit(new.DO$CPSDO),size=1,replace=TRUE)}}
new.DO$CPSDO <- as.numeric(as.character(new.DO$CPSDO))

#RSD
for(i in 1:length(new.DO$RSDDO)){
  if(is.na(new.DO$RSDDO[i])){
    new.DO$RSDDO[i]<-sample(na.omit(new.DO$RSDDO),size=1,replace=TRUE) }}
new.DO$RSDDO <- as.numeric(as.character(new.DO$RSDDO))

#RSS
for(i in 1:length(new.DO$RSSDO)){
  if(is.na(new.DO$RSSDO[i])){
    new.DO$RSSDO[i]<-sample(na.omit(new.DO$RSSDO),size=1,replace=TRUE) }}
new.DO$RSSDO <- as.numeric(as.character(new.DO$RSSDO))

#MDD
for(i in 1:length(new.DO$MDDDO)){
  if(is.na(new.DO$MDDDO[i])){
    new.DO$MDDDO[i]<-sample(na.omit(new.DO$MDDDO),size=1,replace=TRUE) }}
new.DO$MDDDO <- as.numeric(as.character(new.DO$MDDDO))

#MDS
for(i in 1:length(new.DO$MDSDO)){
  if(is.na(new.DO$MDSDO[i])){
    new.DO$MDSDO[i]<-sample(na.omit(new.DO$MDSDO),size=1,replace=TRUE) }}
new.DO$MDSDO <- as.numeric(as.character(new.DO$MDSDO))

######## DO gather ######
DO.gather <-melt(new.DO, id.vars="index",measure.vars=c("CPDDO","CPSDO","RSDDO","RSSDO", "MDDDO","MDSDO"))
levels(DO.gather$variable) <- c("CPD","CPS","RSD","RSS","MDD","MDS") #renaming levels will help us combine data later
DO.gather <-filter(DO.gather, value >= 0, value) #changed the subset filters on this! 




# combine the data


vert <- full_join(T.gather,DO.gather,by=c("index","variable")) #not sure this is in the order we want, but i think it's ok
vert <-full_join(vert,sal.gather, by= c("index","variable"))#add in salinity
names(vert) <- c("index","station","temp","DO","sal")

#add a column that is percent DO 
#we got the pressure from july 30, 2016 at Gabreski
#documentation https://www.rdocumentation.org/packages/rMR/versions/1.1.0/topics/DO.unit.convert 
library("rMR")

vert <-mutate(vert, perDO = DO.unit.convert(DO, DO.units.in = "mg/L", DO.units.out = "pct", bar.units.in = "atm", bar.press =1.0056, temp.C = vert$temp, bar.units.out = "atm", salinity =vert$sal, salinity.units="pp.thou"))

#add in cage week. 
sixset <- dplyr::select(six,"index","week")
vert2 <-full_join(vert,sixset, by="index")
#vert2 <-as.data.frame(filter(vert2,week!=0))#get rid of cage week 0
vert2 <-na.omit(vert2) #filters the na. out. 

#separate station into site vs. depth. 
vert2 <-mutate(vert2,site=substr(station,1,2),depth=substr(station,3,3))

#########

################################## Weekly aggregation of data########################################

##we need to figure out a way to make it clear that the position of these values is relative, so a real index.
vert2_2016 <-group_by(vert2,index)
vert2_2016 <-as.data.frame(arrange(vert2,station)) #i think that worked?
write.csv(vert2_2016, file="vert2_2016.csv")


#we need now to do the cumulative columns. 

#we have to compute summary stats by week for each station. Then we can append it to the cage data

vert2_2016 <-as.data.frame(group_by(vert2_2016,week))

#### IMPORTANT CHANGES ### 

# 1. we made max do a function of percent DO instead of absolute because it depends on the saturation ability. 
# 2. we created a duration of time above 115% DO 

# 4. we took out quantiles because they're redundant. I think. 

vert4_2016 <-ddply(vert2,week~station,summarise,min.do=min(DO),max.do=max(DO),mean.do=mean(DO),sd.do=sd(DO),min.temp=min(temp), max.temp=max(temp), mean.temp=mean(temp),sd.T=sd(temp),do.dur=(sum(DO < 2)*15),temp.dur=(sum(temp >= 27)*15),mean.sal=mean(sal),sd.sal=sd(sal),ssat =(sum(perDO >= 115)*15))
#we multiplied it by 15 so it's the actual duration in minutes. This is going to be suuuuper zero inflated so I am not sure it's informative. 
allsummer2016 <-ddply(vert2,~station,summarise,min.do=min(DO),max.do=max(DO),mean.do=mean(DO),sd.do=sd(DO),min.temp=min(temp), max.temp=max(temp), mean.temp=mean(temp),sd.T=sd(temp),do.dur=(sum(DO < 2)*15),temp.dur=(sum(temp >= 27)*15),mean.sal=mean(sal),sd.sal=sd(sal),ssat =(sum(perDO >= 115)*15))



#use vert4 if you want conditions that week. 
write.csv(vert4_2016, file ="vert4_2016.csv")

#now the cage data. 

ct<-read.csv("cagetotals_16.csv") 
names(ct) <-c("cageID", "Site","depth","sitedepth","week","Fish1","Fish2","Fish3","Fish4","Fish5","Fish6","Fish7", "Fish8","Fish9","Fish10","numdead","survivors","av_length","SD", "Survival1", "Survival2")

ct <-replace(ct,is.na(ct),0)
tl <-na.pass(mutate(ct, tl= as.numeric(ct$Fish1)+as.numeric(ct$Fish2)+as.numeric(ct$Fish3)+as.numeric(ct$Fish4)+as.numeric(ct$Fish5)+as.numeric(ct$Fish6)+as.numeric(ct$Fish7)+as.numeric(ct$Fish8)+as.numeric(ct$Fish9)+as.numeric(ct$Fish10)))

# new way: from Ken rose paper. recall this gives weights in milligrams, not grams like above.  
tl <-na.pass(mutate(tl,w1=((as.numeric(Fish1)/10.723)^(1/0.28)), w2=((as.numeric(Fish2)/10.723)^(1/0.28)),w3=((as.numeric(Fish3)/10.723)^(1/0.28)),w4=((as.numeric(Fish4)/10.723)^(1/0.28)),w5=((as.numeric(Fish5)/10.723)^(1/0.28)),w6=((as.numeric(Fish6)/10.723)^(1/0.28)),w7=((as.numeric(Fish7)/10.723)^(1/0.28)),w8=((as.numeric(Fish8)/10.723)^(1/0.28)),w9=((as.numeric(Fish9)/10.723)^(1/0.28)),w10=((as.numeric(Fish10)/10.723)^(1/0.28))))

#add total biomass column (weight of all fish in that cage that week in grams)
tl <-na.pass(mutate(tl,biom=w1+w2+w3+w4+w5+w6+w7+w8+w9+w10))# it works i think. 

#first, plot biomass by week. gotta summarize the data by week and sitedepth. 
tl<-as.data.frame(group_by(tl,sitedepth))
biol <-ddply(tl,week~sitedepth,summarise,mean=mean(biom),sd=sd(biom))

###### Length by cage by week. barplots to see progression of each fish... #######

tl.gather <-melt(tl, id.vars=c("cageID","Site", "depth", "sitedepth", "week"),measure.vars=c("Fish1","Fish2","Fish3","Fish4", "Fish5","Fish6","Fish7","Fish8","Fish9","Fish10"))

######################## Non-cumulative, Non-Scaled Dataset #########################################


# when we use the non-cumulative data: We actually don't want to set up the dataset using the svrt dataset. 
# instead it's going to be set up where every fish is observed every week, and just gets a status indicator as to whether it's alive or dead. 


#### START HERE ####


#### START HERE ####
cagetotal <-read.csv("cagetotals_16.csv", header=TRUE) 

cgt1 <-melt(cagetotal, id.vars=c("CageID","Site", "depth", "sitedepth", "Week", "Survival1"),measure.vars=c("X1","X2","X3","X4", "X5","X6","X7","X8","X9","X10")) %>% group_by(Week) %>% group_by(CageID) %>% group_by(variable)

#create an individual fish id
cgt2 <-unite_(cgt1, "index", c("CageID","Week","variable"),remove=FALSE)
cgt2 <-unite_(cgt2, "fishID", c("CageID","variable"),remove=FALSE)
cgt2 <-dplyr::rename(cgt2,cage=CageID,site=Site,station=sitedepth,week=Week, Fish=variable,fishlength=value)

#So what we have here is an observation for every fish for every week. and what we need is environmental & summary data to go along with that. 
#Here "week zero" means the day the fish was caught. 


#we need to append biomass, numfish and avlength to cgt2
tl2 <-select(tl,cageID,Site,depth,sitedepth,week,survivors,av_length,SD,biom) #week x cageID
tl2 <- dplyr::rename(tl2,station=sitedepth)

#gather up all the environmental data where impacts are non-cumulative.
#incase you need to reset the creation of vert4. 
vert4 <-ddply(vert2,week~station,summarise,min.do=min(DO),max.do=max(DO),mean.do=mean(DO),sd.do=sd(DO),min.temp=min(temp), max.temp=max(temp), mean.temp=mean(temp),sd.T=sd(temp),do.dur=(sum(DO < 2)*15),temp.dur=(sum(temp >= 27)*15),mean.sal=mean(sal),sd.sal=sd(sal),ssat =(sum(perDO >= 115)*15))

vert4 <- mutate(vert4, week=as.integer(week)) #use vert4 because you're not doing cumulative data. 
# week 1 here really means the end of week 1. Whereas week 1 in cage data means the 2nd measurement of the fish. 
#So IT IS WEEK 2! Environmental data starts on week 2! 
vert4 <-mutate(vert4, week=week+1)# now the environmental data is properly synched to the survival data. 

svrt <- left_join(tl2,vert4,by=c("week","station")) 
svrt <- dplyr::rename(svrt, cage=cageID,site=Site, week = week) 

###############At this point 1/15/19, we're going to write svrt to a csv to use in the growth models. ###############################
write.csv(svrt, file="svrt2016forgrowth_regenerated.csv")

# Attach the survival data for  
tl4 <-right_join(svrt,cgt2, by=c("cage","week","depth","site","station"))

#delete the entries where fishlen = NA - because we don't take observations for a fish after that fish has died. 
tl4 <-filter(tl4,!is.na(fishlength))

#now we somehow need to incorporate the death.week. 
# for each unique value in fishID, count how many instances of that value there are, then give the variable "deathweek" that value. 
cgt3 <- dplyr::count(tl4,fishID)  # so n now is number of observations of that organism. the maximum should be 12. 

tl4 <-right_join(tl4,cgt3, by="fishID")  # cool so now an "n" of 12 means that fish was observed 12 times. the first was when it was put in. 
# so if n= 1 that means the fish died in the first week.
# if n = 11 that means the fish died between week 10 and week 11, on week 11 the fish was not there.
# if n = 12 the fish survived the experiment because it was observed when the cage was dismantled. 

#fish that were obeserved 12 times survived the experiment. 
#fish that were observed 11 times died during interval 11 (week 10)


tl4 <- mutate(tl4, death.week=ifelse(n==12, 12, n-1)) #because we don't count the first measurement from the day the fish was put in as it having survived a week. 
tl4 <-mutate(tl4, event = ifelse(death.week==12,3,NA)) %>% mutate(event2=ifelse(death.week==12,0,NA))#status = 0 if fish had 12 observations, meaning it survived the experiment, otherwise leave it blank. 

tl5 <-mutate(tl4, start=week-1, end=week)
tl5 <-mutate(tl5, event=ifelse(is.na(event) & death.week==end,1,3)) %>% mutate(event2=ifelse(is.na(event2) & death.week == end,1,0)) #if the end of the week of the of the observation = number of times animal was observed, the animal died in the preceding interval.
#tl5 <-filter(tl5, week > 0)
tl5 <-filter(tl5, week > 1)

View(tl5)

write.csv(tl5,file="weeklysummarydata_2016_10-8-19_regenerated.csv")


################################################# Fake week 12 version ############################

#start with the creation of tl4. 

tl4 <-right_join(svrt,cgt2, by=c("cage","week","depth","site","station"))
tl4 <-filter(tl4,!is.na(fishlength))
cgt3 <- dplyr::count(tl4,fishID)
tl4 <-right_join(tl4,cgt3, by="fishID") 

#the easiest way to do this would be to impute fake environmental data for weeks 0 and 1. 
#min.do
for(i in 1:length(tl4$min.do)){
  if(is.na(tl4$min.do[i])){
    tl4$min.do[i]<-sample(na.omit(tl4$min.do),size=1,replace=TRUE)}}
#max.do
for(i in 1:length(tl4$max.do)){
  if(is.na(tl4$max.do[i])){
    tl4$max.do[i]<-sample(na.omit(tl4$max.do),size=1,replace=TRUE)}}
#mean.do
for(i in 1:length(tl4$mean.do)){
  if(is.na(tl4$mean.do[i])){
    tl4$mean.do[i]<-sample(na.omit(tl4$mean.do),size=1,replace=TRUE)}}
#sd.do
for(i in 1:length(tl4$sd.do)){
  if(is.na(tl4$sd.do[i])){
    tl4$sd.do[i]<-sample(na.omit(tl4$sd.do),size=1,replace=TRUE)}}
#min.temp
for(i in 1:length(tl4$min.temp)){
  if(is.na(tl4$min.temp[i])){
    tl4$min.temp[i]<-sample(na.omit(tl4$min.temp),size=1,replace=TRUE)}}
#sd.T
for(i in 1:length(tl4$sd.T)){
  if(is.na(tl4$sd.T[i])){
    tl4$sd.T[i]<-sample(na.omit(tl4$sd.T),size=1,replace=TRUE)}}
#do.dur
for(i in 1:length(tl4$do.dur)){
  if(is.na(tl4$do.dur[i])){
    tl4$do.dur[i]<-sample(na.omit(tl4$do.dur),size=1,replace=TRUE)}}
#temp.dur
for(i in 1:length(tl4$temp.dur)){
  if(is.na(tl4$temp.dur[i])){
    tl4$temp.dur[i]<-sample(na.omit(tl4$temp.dur),size=1,replace=TRUE)}}
#mean.sal
for(i in 1:length(tl4$mean.sal)){
  if(is.na(tl4$mean.sal[i])){
    tl4$mean.sal[i]<-sample(na.omit(tl4$mean.sal),size=1,replace=TRUE)}}
#sd.sal
for(i in 1:length(tl4$sd.sal)){
  if(is.na(tl4$sd.sal[i])){
    tl4$sd.sal[i]<-sample(na.omit(tl4$sd.sal),size=1,replace=TRUE)}}
#ssat
for(i in 1:length(tl4$ssat)){
  if(is.na(tl4$ssat[i])){
    tl4$ssat[i]<-sample(na.omit(tl4$ssat),size=1,replace=TRUE)}}
#max.temp
for(i in 1:length(tl4$max.temp)){
  if(is.na(tl4$max.temp[i])){
    tl4$max.temp[i]<-sample(na.omit(tl4$max.temp),size=1,replace=TRUE)}}
#mean.temp
for(i in 1:length(tl4$mean.temp)){
  if(is.na(tl4$mean.temp[i])){
    tl4$mean.temp[i]<-sample(na.omit(tl4$mean.temp),size=1,replace=TRUE)}}

tw4 <- mutate(tl4, death.week=ifelse(n==12, 12, n-1)) #because we don't count the first measurement from the day the fish was put in as it having survived a week. 
tw4 <-mutate(tw4, event = ifelse(death.week==12,3,NA)) %>% mutate(event2=ifelse(death.week==12,0,NA))#status = 0 if fish had 12 observations, meaning it survived the experiment, otherwise leave it blank. 

tw5 <-mutate(tw4, start=week-1, end=week)
tw5 <-mutate(tw5, event=ifelse(is.na(event) & death.week==end,1,3)) %>% mutate(event2=ifelse(is.na(event2) & death.week == end,1,0)) #if the end of the week of the of the observation = number of times animal was observed, the animal died in the preceding interval.
tw5 <-filter(tw5, week >= 1)

View(tw5)

#write.csv(tw5,file="weeklysummarydata_2016_extraweek.csv")








###################################################################NOW MAKE THE SCALED VERSION ############################################################

tl6 <- mutate_at(tl5, vars(biom, min.do, max.do, mean.do, sd.do, min.temp, sd.do, min.temp, max.temp, sd.T, do.dur, temp.dur, mean.sal, sd.sal, ssat), scale)

#write.csv(tl6,file="weeklysummScaledata_2016.csv")



#################################################################NOW MAKE THE CUMULATIVE VERSION ###############################################################


#Use vert5 if you want cumulative conditions. 
wk1 <- filter(vert2, week==1) %>% dplyr::group_by(station) %>% summarise(min.do=min(DO),max.do=max(DO),mean.do=mean(DO),sd.do=sd(DO),min.temp=min(temp), max.temp=max(temp), mean.temp=mean(temp),sd.T=sd(temp),do.dur=(sum(DO < 2)*15),temp.dur=(sum(temp >= 27)*15),mean.sal=mean(sal),sd.sal=sd(sal),ssat =(sum(perDO >= 115)*15)) %>% mutate(week="1")
wk2 <- filter(vert2, week <= 2) %>% group_by(station) %>% summarise(min.do=min(DO),max.do=max(DO),mean.do=mean(DO),sd.do=sd(DO),min.temp=min(temp), max.temp=max(temp), mean.temp=mean(temp),sd.T=sd(temp),do.dur=(sum(DO < 2)*15),temp.dur=(sum(temp >= 27)*15),mean.sal=mean(sal),sd.sal=sd(sal),ssat =(sum(perDO >= 115)*15)) %>% mutate(week="2")
wk3 <- filter(vert2, week <= 3) %>% group_by(station) %>% summarise(min.do=min(DO),max.do=max(DO),mean.do=mean(DO),sd.do=sd(DO),min.temp=min(temp), max.temp=max(temp), mean.temp=mean(temp),sd.T=sd(temp),do.dur=(sum(DO < 2)*15),temp.dur=(sum(temp >= 27)*15),mean.sal=mean(sal),sd.sal=sd(sal),ssat =(sum(perDO >= 115)*15)) %>% mutate(week="3")
wk4 <- filter(vert2, week <= 4) %>% group_by(station) %>% summarise(min.do=min(DO),max.do=max(DO),mean.do=mean(DO),sd.do=sd(DO),min.temp=min(temp), max.temp=max(temp), mean.temp=mean(temp),sd.T=sd(temp),do.dur=(sum(DO < 2)*15),temp.dur=(sum(temp >= 27)*15),mean.sal=mean(sal),sd.sal=sd(sal),ssat =(sum(perDO >= 115)*15)) %>% mutate(week="4")
wk5 <- filter(vert2, week <= 5) %>% group_by(station) %>% summarise(min.do=min(DO),max.do=max(DO),mean.do=mean(DO),sd.do=sd(DO),min.temp=min(temp), max.temp=max(temp), mean.temp=mean(temp),sd.T=sd(temp),do.dur=(sum(DO < 2)*15),temp.dur=(sum(temp >= 27)*15),mean.sal=mean(sal),sd.sal=sd(sal),ssat =(sum(perDO >= 115)*15)) %>% mutate(week="5")
wk6 <- filter(vert2, week <= 6) %>% group_by(station) %>% summarise(min.do=min(DO),max.do=max(DO),mean.do=mean(DO),sd.do=sd(DO),min.temp=min(temp), max.temp=max(temp), mean.temp=mean(temp),sd.T=sd(temp),do.dur=(sum(DO < 2)*15),temp.dur=(sum(temp >= 27)*15),mean.sal=mean(sal),sd.sal=sd(sal),ssat =(sum(perDO >= 115)*15)) %>% mutate(week="6")
wk7 <- filter(vert2, week <= 7) %>% group_by(station) %>% summarise(min.do=min(DO),max.do=max(DO),mean.do=mean(DO),sd.do=sd(DO),min.temp=min(temp), max.temp=max(temp), mean.temp=mean(temp),sd.T=sd(temp),do.dur=(sum(DO < 2)*15),temp.dur=(sum(temp >= 27)*15),mean.sal=mean(sal),sd.sal=sd(sal),ssat =(sum(perDO >= 115)*15)) %>% mutate(week="7")
wk8 <- filter(vert2, week <= 8) %>% group_by(station) %>% summarise(min.do=min(DO),max.do=max(DO),mean.do=mean(DO),sd.do=sd(DO),min.temp=min(temp), max.temp=max(temp), mean.temp=mean(temp),sd.T=sd(temp),do.dur=(sum(DO < 2)*15),temp.dur=(sum(temp >= 27)*15),mean.sal=mean(sal),sd.sal=sd(sal),ssat =(sum(perDO >= 115)*15)) %>% mutate(week="8")
wk9 <- filter(vert2, week <= 9) %>% group_by(station) %>% summarise(min.do=min(DO),max.do=max(DO),mean.do=mean(DO),sd.do=sd(DO),min.temp=min(temp), max.temp=max(temp), mean.temp=mean(temp),sd.T=sd(temp),do.dur=(sum(DO < 2)*15),temp.dur=(sum(temp >= 27)*15),mean.sal=mean(sal),sd.sal=sd(sal),ssat =(sum(perDO >= 115)*15)) %>% mutate(week="9")
wk10 <- filter(vert2, week <= 10) %>% group_by(station) %>% summarise(min.do=min(DO),max.do=max(DO),mean.do=mean(DO),sd.do=sd(DO),min.temp=min(temp), max.temp=max(temp), mean.temp=mean(temp),sd.T=sd(temp),do.dur=(sum(DO < 2)*15),temp.dur=(sum(temp >= 27)*15),mean.sal=mean(sal),sd.sal=sd(sal),ssat =(sum(perDO >= 115)*15)) %>% mutate(week="10")
wk11 <- filter(vert2, week <= 11) %>% group_by(station) %>% summarise(min.do=min(DO),max.do=max(DO),mean.do=mean(DO),sd.do=sd(DO),min.temp=min(temp), max.temp=max(temp), mean.temp=mean(temp),sd.T=sd(temp),do.dur=(sum(DO < 2)*15),temp.dur=(sum(temp >= 27)*15),mean.sal=mean(sal),sd.sal=sd(sal),ssat =(sum(perDO >= 115)*15)) %>% mutate(week="11")

vert5 <-bind_rows(wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11)

vert5 <- mutate(vert5, week=as.integer(week)) #use vert4 because you're not doing cumulative data. 
# week 1 here really means the end of week 1. Whereas week 1 in cage data means the 2nd measurement of the fish. 
#So IT IS WEEK 2! Environmental data starts on week 2! 
vert5 <-mutate(vert5, week=week+1)# now the environmental data is properly synched to the survival data. 

svrt5 <- left_join(tl2,vert5,by=c("week","station")) 
svrt5 <- dplyr::rename(svrt5, cage=cageID,site=Site, week = week) 


# Attach the survival data for  
tc4 <-right_join(svrt5,cgt2, by=c("cage","week","depth","site","station"))

#delete the entries where fishlen = NA - because we don't take observations for a fish after that fish has died. 
tc4 <-filter(tc4,!is.na(fishlength))

#now we somehow need to incorporate the death.week. 
# for each unique value in fishID, count how many instances of that value there are, then give the variable "deathweek" that value. 
cgt3 <- dplyr::count(tc4,fishID)  # so n now is number of observations of that organism. the maximum should be 12. 

tc4 <-right_join(tc4,cgt3, by="fishID")  # cool so now an "n" of 12 means that fish was observed 12 times. the first was when it was put in. 
# so if n= 1 that means the fish died in the first week.
# if n = 11 that means the fish died between week 10 and week 11, on week 11 the fish was not there.
# if n = 12 the fish survived the experiment because it was observed when the cage was dismantled. 

#fish that were obeserved 12 times survived the experiment. 
#fish that were observed 11 times died during interval 11 (week 10)


tc4 <- mutate(tc4, death.week=ifelse(n==12, 12, n-1)) #because we don't count the first measurement from the day the fish was put in as it having survived a week. 
tc4 <-mutate(tc4, event = ifelse(death.week==12,3,NA)) %>% mutate(event2=ifelse(death.week==12,0,NA))#status = 0 if fish had 12 observations, meaning it survived the experiment, otherwise leave it blank. 

tc5 <-mutate(tc4, start=week-1, end=week)
tc5 <-mutate(tc5, event=ifelse(is.na(event) & death.week==end,1,3)) %>% mutate(event2=ifelse(is.na(event2) & death.week == end,1,0))#if the end of the week of the of the observation = number of times animal was observed, the animal died in the preceding interval.
tc5 <-filter(tc5, week > 1)

View(tc5)

#write.csv(tc5,file="cumsummarydata_2016.csv")


#create the dataset for death week only. 

tc12 <-filter(tc5, death.week==12 & week==11) #these are all the fish that survived the experiment
dim(tc12)
tc6 <-filter(tc5, event==1) #these are all the fish that died during the experiment
dim(tc6)
140+24 #### TOO FEW FISH! i think because a lot of them died in that first week.... 
tc7 <-bind_rows(tc6,tc12)

#write.csv(tc7,file="cumsum2016_deathwkonly.csv")


################### Create a special version of the weekly summary data with an additional week. Impute data for that week ######

#the purpose of this is to test the assumptions with degrees of freedom. 


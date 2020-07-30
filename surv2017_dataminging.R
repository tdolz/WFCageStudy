#packages you will need

library("survival")
library("car")
library("ggplot2")
library("reshape2")
library("plyr")
library("moments")
library("reshape")
library("tidyr")
library("dplyr")
library("margins")
library("rMR")


#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github/WFCageStudy")

#inspect the data.
sev<-read.csv("indexed_all17_10_30_17.csv")

##### GITHUB VERSION ######

################################## Dissolved Oxygen ############################
new.do <-dplyr::select(sev,cage.week, index, contains("DO"))
new.do <-slice(new.do, 578:7910) # from 15:00 on 6/15 to midnight on august 31.

### NSDO is fine as is. 
################ MDDO #############
#need to impute data to cover week of 15:00 6/15 to 8:00 6/21 (1126)
#do this by fitting a line to the data between the MDD0 and MSDO and using the equatio of the line to predict missing values for MDDO
#plot(new.do$MSDO,new.do$MDDO)
#abline(lm(new.do$MDDO~new.do$MSDO)) # y ~ x
#summary(lm(new.do$MDDO~new.do$MSDO))
#imputing the way we did in problem set 13 is inappropriate because there are temporal and diel trends. better to base it off of MSDO

#imputing the MDDO data by using the slope from MSDO to predict empty MDDO
for(i in 1:length(new.do$MDDO)){
  if(new.do$MDDO[i]==0){
    new.do$MDDO[i]<-4.925872+0.326993*new.do$MSDO[i]
  }
} 

################ MSDO ######### 
for(i in 1:length(new.do$MSDO)){  #just fix the zeros
  if(new.do$MSDO[i]< 1){
    new.do$MSDO[i]<-sample(na.omit(new.do$MSDO),size=1,replace=TRUE)
  }
}
################ SDDO ############## 
for(i in 1:length(new.do$SDDO)){
  if(new.do$SDDO[i]< 1){
    new.do$SDDO[i]<-sample(na.omit(new.do$SDDO),size=1,replace=TRUE)
  }
}
################ SSDO ############
for(i in 1:length(new.do$SSDO)){
  if(new.do$SSDO[i]< 1){
    new.do$SSDO[i]<-sample(na.omit(new.do$SSDO),size=1,replace=TRUE)
  }
}
################ NDDO ############## 
#a hot mess. Impute from NSDO
#plot(new.do$index,new.do$NDDO)
#plot(new.do$NSDO,new.do$NDDO)
#abline(lm(new.do$NDDO~new.do$NSDO)) # y ~ x
#summary(lm(new.do$NDDO~new.do$NSDO))
# y= mx + b -> NDDO = 0.141495532x + 5.018696125
for(i in 1:length(new.do$NDDO)){
  if(new.do$NDDO[i]<3){
    new.do$NDDO[i]<-5.018696125+0.141495532*new.do$NSDO[i]
  }
} 

################ DO Time Series ##########
#DO timeseries
#windows(30,20)
#par(mfrow=c(2,3))
#par(mar=c(4,4,1,1))
#MDDO
#plot(new.do$index,new.do$MDDO, ylim=c(0,20),ylab="(mg/L)",main="MDDO", cex=2)
#abline(h=2)
#abline(h=5)
#MSDO
#plot(new.do$index,new.do$MSDO, ylim=c(0,20),ylab="(mg/L)",main="MSDO", cex=2)
#abline(h=2)
#abline(h=5)
#NDDO
#plot(new.do$index,new.do$NDDO, ylim=c(0,20),ylab="(mg/L)",main="NDDO", cex=2)
#abline(h=2)
#abline(h=5)
#NSDO
#plot(new.do$index,new.do$NSDO, ylim=c(0,20),ylab="(mg/L)",main="NSDO", cex=2)
#abline(h=2)
#abline(h=5)
#SDDO
#plot(new.do$index,new.do$SDDO, ylim=c(0,20),ylab="(mg/L)",main="SDDO", cex=2)
#abline(h=2)
#abline(h=5)
#SSDO
#plot(new.do$index,new.do$SSDO, ylim=c(0,20),ylab="(mg/L)",main="SSDO", cex=2)
#abline(h=2)
#abline(h=5)

################ DO Boxplots ####################################

names(new.do) <- c("week","index", "MDDO","MSDO","NDDO","NSDO", "SDDO","SSDO")
head(new.do)
new.do <-as.data.frame(new.do)

do.gather <-melt(new.do, id.vars="index",measure.vars=c("MDDO","MSDO","NDDO","NSDO", "SDDO","SSDO"))
levels(do.gather$variable) <- c("MD","MS","ND","NS","SD","SS") #renaming levels will help us combine data later
do.gather <-filter(do.gather, value > 0 & value < 20) #subsetting between 0 and 20, because outside of this range isn't plausible. 

##make a nicer boxplot
s <-ggplot(do.gather,aes(factor(variable),value))
s + geom_violin(scale="count", aes(fill=factor(variable)), kernel="gaussian", trim=FALSE)+
  geom_boxplot(width=0.2,notch=TRUE)+
  labs(title="Dissolved Oxygen", x="station", y="(mg/L)")+
  stat_summary(fun.y = "mean", geom="point", shape =8, size=1, color="purple")+
  geom_hline(yintercept = 2)+
  geom_hline(yintercept = 5)





#################################### SALINITY #####################################################################

new.sal <-dplyr::select(sev,cage.week, index, contains("sal"))
new.sal <-slice(new.sal, 578:7910) # from 15:00 on 6/15 to midnight on august 31.

#plot to see where we need to impute
#windows(10,10)
#ggplot(new.sal, aes(index, y = value, color = variable))+
#  geom_point(aes(y=MS.sal,col="MS.sal"))+
#  geom_point(aes(y=MD.sal,col="MD.sal"))+
#  geom_point(aes(y=ND.sal,col="ND.sal"))+
#  geom_point(aes(y=NS.sal,col="NS.sal"))+
#  geom_point(aes(y=SD.sal,col="SD.sal"))+
#  geom_point(aes(y=ss.sal.new,col="ss.sal.new"))

#we have to determine a relationship that is based on the parts of each time series that is not compromised. 
#ggplot(new.sal, aes(index, y = value, color = variable))+
#  geom_point(aes(y=MS.sal,col="MS.sal"))+
#  geom_point(aes(y=MD.sal,col="MD.sal"))

################ MD.sal ################ 
#develope a relationship that excludes the missing data by filtering it out temporarily. 
#rel.sal <- filter(new.sal, MS.sal > 23, MD.sal > 23)
#plot(rel.sal$MS.sal, rel.sal$MD.sal)
#abline(lm(rel.sal$MD.sal~rel.sal$MS.sal)) # y ~ x we are predicting missing values for MD.sal here (the orange)
#summary(lm(rel.sal$MD.sal~rel.sal$MS.sal))
#y = mx + b: MD.sal = -0.096570x+33.166188


#imputing the MD.sal data by using the slope from MS.sal to predict MD.sal
for(i in 1:length(new.sal$MD.sal)){
  if(new.sal$MD.sal[i]< 23){
    new.sal$MD.sal[i] <- -0.096570*new.sal$MS.sal[i] + 33.166188
  }
}

################ MS.sal ##############
#now to predict MS from MD (predict the blue)
#plot(rel.sal$MD.sal, rel.sal$MS.sal)
#abline(lm(rel.sal$MS.sal~rel.sal$MD.sal))
#summary(lm(rel.sal$MS.sal~rel.sal$MD.sal))
#y=mx + b: MS.sal = -0.33153x + 38.04454

#imputing the MS.sal data by using the slope from MD.sal to predict MS.sal
for(i in 1:length(new.sal$MS.sal)){
  if(new.sal$MS.sal[i]< 23){
    new.sal$MS.sal[i] <- -0.33153*new.sal$MD.sal[i] + 38.04454
  }
}

################ Impute the rest  ###########
#now impute the rest by drawing randomly from the rest of the trend. 
for(i in 1:length(new.sal$ND.sal)){
  if(new.sal$ND.sal[i]< 26){
    new.sal$ND.sal[i]<-sample(na.omit(new.sal$ND.sal),size=1,replace=TRUE)
  }
}
#NS.sal 
for(i in 1:length(new.sal$NS.sal)){
  if(new.sal$NS.sal[i]< 24){
    new.sal$NS.sal[i]<-sample(na.omit(new.sal$NS.sal),size=1,replace=TRUE)
  }
}
#SD.sal
for(i in 1:length(new.sal$SD.sal)){
  if(new.sal$SD.sal[i]< 24){
    new.sal$SD.sal[i]<-sample(na.omit(new.sal$SD.sal),size=1,replace=TRUE)
  }
}
#SS.sal
for(i in 1:length(new.sal$ss.sal.new)){
  if(new.sal$ss.sal.new[i]< 26){
    new.sal$ss.sal.new[i]<-sample(na.omit(new.sal$ss.sal.new),size=1,replace=TRUE)
  }
}

##### plots #####
#salnity timeseries
windows(30,20)
par(mfrow=c(2,3))
par(mar=c(4,4,1,1))
#MD.sal
plot(new.sal$index,new.sal$MD.sal, ylim=c(20,35),ylab="(ppt)",main="MD.sal", cex=2)
abline(h=25)
abline(h=30)
#MS.sal
plot(new.sal$index,new.sal$MS.sal, ylim=c(20,35),ylab="(ppt)",main="MS.sal", cex=2)
abline(h=25)
abline(h=30)
#ND.sal
plot(new.sal$index,new.sal$ND.sal, ylim=c(20,35),ylab="(ppt)",main="ND.sal", cex=2)
abline(h=25)
abline(h=30)
#NS.sal
plot(new.sal$index,new.sal$NS.sal, ylim=c(20,35),ylab="(ppt)",main="NS.sal", cex=2)
abline(h=25)
abline(h=30)
#SD.sal
plot(new.sal$index,new.sal$SD.sal, ylim=c(20,35),ylab="(ppt)",main="SD.sal", cex=2)
abline(h=25)
abline(h=30)
#SS.sal
plot(new.sal$index,new.sal$ss.sal.new, ylim=c(20,35),ylab="(ppt)",main="SS.sal", cex=2)
abline(h=25)
abline(h=30)


############## Salinity Violin plots ##################
sal.gather <-melt(new.sal, id.vars="index",measure.vars=c("MD.sal","MS.sal","ND.sal","NS.sal", "SD.sal","ss.sal.new"))
levels(sal.gather$variable) <- c("MD","MS","ND","NS","SD","SS") #renaming levels will help us combine data later
sal.gather <-filter(sal.gather, value > 22, value) #subsetting between 0 and 23, because outside of this range isn't plausible. 

#sa <-ggplot(sal.gather,aes(factor(variable),value))
#sa + geom_violin(scale="count", aes(fill=factor(variable)), kernel="gaussian", trim=FALSE)+
#  geom_boxplot(width=0.2,notch=TRUE)+
#  labs(title="Salinity", x="station", y="(ppt)")+
#  stat_summary(fun.y = "mean", geom="point", shape =8, size=1, color="purple")



################################### Temperature ###############################################

new.T <-select(sev,cage.week, index, contains("Temp"))
new.T <-slice(new.T, 578:7910)
sevset <- select(sev,"index","cage.week")
new.T <-full_join(new.T,sevset, by="index")
new.T <-as.data.frame(new.T)
names(new.T) <- c("week","index", "MDTemp","MSTemp","NDTemp","NSTemp", "SDTemp","SSTemp","cage.week")
#earlier we had decided we were going to subset the data from  15:00 on 6/15 to midnight on august 31.
#MSTemp.d was basically zeros, so we went to the cagedata folder to hproj17 and found more up to date data, which I inserted into Sev.

#imputing temperature data:DECIDED NOT TO 
#SS, NS, MS, is fine.
#SD needs some fixing for the drops. 

#for(i in new.T$SDTemp[2016]:new.T$SDTemp[7333]){
 # if(new.T$SDTemp[i]< 18){
   # new.T$SDTemp[i]<-sample(na.omit(new.T$SDTemp),size=1,replace=TRUE)}}

################ Temperature time series ####################################################################################################################
#let's make time series plots first, that way if something is terrible, it will be obvious. 
#temp timeseries
#windows(30,20)
#par(mfrow=c(2,3))
#par(mar=c(4,4,1,1))
#MDTemp
#plot(new.T$index,new.T$MDTemp.c, ylim=c(10,30),ylab="(C)",main="MD Temp", cex=2)
#abline(h=18.5)
#abline(h=27)
#MSTemp
#plot(new.T$index,new.T$MSTemp.d, ylim=c(10,30),ylab="(C)",main="MS Temp", cex=2)
#abline(h=18.5)
#abline(h=27)
#NDTemp
#plot(new.T$index,new.T$NDTemp, ylim=c(10,30),ylab="(C)",main="ND Temp", cex=2)
#abline(h=18.5)
#abline(h=27)
#NSTemp
#plot(new.T$index,new.T$NSTemp.d, ylim=c(10,30),ylab="(C)",main="NS Temp", cex=2)
#abline(h=18.5)
#abline(h=27)
#SDTemp
#plot(new.T$index,new.T$SDTemp, ylim=c(10,30),ylab="(C)",main="SD Temp", cex=2)
#abline(h=18.5)
#abline(h=27)
#SSTemp
#plot(new.T$index,new.T$SSTemp, ylim=c(10,30),ylab="(C)",main="SS Temp", cex=2)
#abline(h=18.5)
#abline(h=27)

################ Temperature Violin Plots######
#melt the data frame
T.gather <-melt(new.T, id.vars="index",measure.vars=c("MDTemp","MSTemp","NDTemp","NSTemp", "SDTemp","SSTemp"))
levels(T.gather$variable) <- c("MD","MS","ND","NS","SD","SS")
T.gather <- filter(T.gather,T.gather$value > 5) #values at zero are likely instrument malfunctioning. 

#make nicer boxplot
#p <-ggplot(T.gather,aes(factor(variable),value))
#p + geom_violin(scale="count", aes(fill=factor(variable)), kernel="gaussian", trim=FALSE)+
#  geom_boxplot(width=0.2,notch=TRUE)+
#  labs(title="Temperature", x="station", y="C")+
#  stat_summary(fun.y = "mean", geom="point", shape =8, size=1, color="purple")+
#  geom_hline(yintercept = 18.5)+
#  geom_hline(yintercept = 27)


################################### SUMMARY STATS #######################


################ data minging: create vert, vert2 #####

#combine do.gather and T.gather into one column
vert <- full_join(T.gather,do.gather,by=c("index","variable")) #not sure this is in the order we want, but i think it's ok
vert <-full_join(vert,sal.gather, by= c("index","variable"))#add in salinity
names(vert) <- c("index","station","temp","DO","sal")

#add a column that is percent DO 
#we got the pressure from july 30, 2016 at Gabreski
#documentation https://www.rdocumentation.org/packages/rMR/versions/1.1.0/topics/DO.unit.convert 
library("rMR")
vert <-mutate(vert, perDO = DO.unit.convert(DO, DO.units.in = "mg/L", DO.units.out = "pct", bar.units.in = "atm", bar.press =1.0056, temp.C = vert$temp, bar.units.out = "atm", salinity =vert$sal, salinity.units="pp.thou"))

#add in cage week. 
sevset <- select(sev,"index","cage.week")
vert2 <-full_join(vert,sevset, by="index") %>% mutate(week=cage.week) %>% select(-cage.week)
#vert2 <-as.data.frame(filter(vert2,week!=0))#get rid of cage week 0
vert2 <-na.omit(vert2) #filters the na. out. 

#separate station into site vs. depth. 
vert2 <-mutate(vert2,site=substr(station,1,1),depth=substr(station,2,2))

################################## Weekly aggregation of data########################################

##we need to figure out a way to make it clear that the position of these values is relative, so a real index.
vert2 <-group_by(vert2,index)
vert2 <-as.data.frame(arrange(vert2,station)) #i think that worked?
#vert2 <- tibble::rowid_to_column(vert2, "ID") #adds an index, which may be useful. 

write.csv(vert2, "vert2_2017.csv")

#### IMPORTANT CHANGES ### 

# 1. we made max do a function of percent DO instead of absolute because it depends on the saturation ability. 
# 2. we created a duration of time above 115% DO 
# 3. we created some additional metrics relateing to quantiles of the environmental data. 


vert4 <-ddply(vert2,week~station,summarise,min.do=min(DO),max.do=max(DO),mean.do=mean(DO),sd.do=sd(DO),min.temp=min(temp), max.temp=max(temp), mean.temp=mean(temp),sd.T=sd(temp),do.dur=(sum(DO < 2)*15),temp.dur=(sum(temp >= 27)*15),mean.sal=mean(sal),sd.sal=sd(sal),ssat =(sum(perDO >= 115)*15))
write.csv(vert4, file="vert4_2017.csv")
allsummer2017 <-ddply(vert2,~station,summarise,min.do=min(DO),max.do=max(DO),mean.do=mean(DO),sd.do=sd(DO),min.temp=min(temp), max.temp=max(temp), mean.temp=mean(temp),sd.T=sd(temp),do.dur=(sum(DO < 2)*15),temp.dur=(sum(temp >= 27)*15),mean.sal=mean(sal),sd.sal=sd(sal),ssat =(sum(perDO >= 115)*15))


################################## now add the cage data.#########################################
#ct<-read.csv("cagetotals_17fixed.csv") #in this version, NOPE 
#ct <- ct[,1:15]

ct <- dplyr::rename(ct, cageID = CageID, Fish1=X1, Fish2=X2, Fish3=X3, Fish4=X4, Fish5=X5, Fish6=X6, Fish7=X7, Fish8=X8, Fish9=X9, Fish10=X10)

ct <-mutate(ct,Site=substr(cageID,1,1),depth=substr(cageID,2,2),sitedepth=substr(cageID,1,2))
ct$depth <-ifelse(ct$depth=="D","deep","shallow")


ct <-replace(ct,is.na(ct),0) #you need to do this because otherwise the NA handleing is off later, although there must be a better way/ 

###### adding biomass & length #####
#Janet said that biomass in the cage is an important part of growth & survivorship.... So I computed total length for the cage
tl <-na.pass(mutate(ct, tl= as.numeric(ct$Fish1)+as.numeric(ct$Fish2)+as.numeric(ct$Fish3)+as.numeric(ct$Fish4)+as.numeric(ct$Fish5)+as.numeric(ct$Fish6)+as.numeric(ct$Fish7)+as.numeric(ct$Fish8)+as.numeric(ct$Fish9)+as.numeric(ct$Fish10)))
#but what we would we ideally want is an actual measure of biomass. 
# so we can use this equation W = cL^b, where W = weight in grams, c and b are constants and L is length
#aka log(W) = log(c)+ B log(L)
# for winter flounder according to Lux 1969, wf b= 3.070 and logc = -5.440

#so add new weight columns --> check if you did this right. 
# old way
#tl <-na.pass(mutate(tl,w1=(10^(-5.440))*as.numeric(Fish1)^(3.070), w2=(10^(-5.440))*as.numeric(Fish2)^(3.070),w3=(10^(-5.440))*as.numeric(Fish3)^(3.070), w4=(10^(-5.440))*as.numeric(Fish4)^(3.070),w5=(10^(-5.440))*as.numeric(Fish5)^(3.070),w6=(10^(-5.440))*as.numeric(Fish6)^(3.070),w7=(10^(-5.440))*as.numeric(Fish7)^(3.070),w8=(10^(-5.440))*as.numeric(Fish8)^(3.070),w9=(10^(-5.440))*as.numeric(Fish9)^(3.070),w10=(10^(-5.440))*as.numeric(Fish10)^(3.070)))
# new way: from Ken rose paper. recall this gives weights in milligrams, not grams like above.  
tl <-na.pass(mutate(tl,w1=((as.numeric(Fish1)/10.723)^(1/0.28)), w2=((as.numeric(Fish2)/10.723)^(1/0.28)),w3=((as.numeric(Fish3)/10.723)^(1/0.28)),w4=((as.numeric(Fish4)/10.723)^(1/0.28)),w5=((as.numeric(Fish5)/10.723)^(1/0.28)),w6=((as.numeric(Fish6)/10.723)^(1/0.28)),w7=((as.numeric(Fish7)/10.723)^(1/0.28)),w8=((as.numeric(Fish8)/10.723)^(1/0.28)),w9=((as.numeric(Fish9)/10.723)^(1/0.28)),w10=((as.numeric(Fish10)/10.723)^(1/0.28))))

#add total biomass column (weight of all fish in that cage that week in grams)
tl <-na.pass(mutate(tl,biom=w1+w2+w3+w4+w5+w6+w7+w8+w9+w10))# it works i think. 

r <-ggplot(tl, aes(Week,biom,color=sitedepth)) + 
  geom_point(size=6)+
  labs(title="total biomass", x="week", y="milligrams")
r 

###### plot biomass, length, number of fish by week for each station ######

#first, plot biomass by week. gotta summarize the data by week and sitedepth. 
tl<-as.data.frame(group_by(tl,sitedepth))
biol <-ddply(tl,Week~sitedepth,summarise,mean=mean(biom),sd=sd(biom))
#ggplot(biol, aes(week,mean,color=sitedepth))+
#  geom_line(size=1)+
#labs(title="mean total biomass", x="week", y="milligrams")

#av_length
biolen <-ddply(tl,Week~sitedepth,summarise,mean=mean(as.numeric(avlength)),sd=sd(as.numeric(avlength)))
#ggplot(biolen, aes(week,mean,color=sitedepth))+
#  geom_line(size=1)+
#  labs(title="mean of cage average lengths", x="week", y="mm")

#mortality plot (num fish)
mor <-ddply(tl,Week~sitedepth,summarise,mean=mean(count),sd=sd(count))
#ggplot(mor, aes(week,mean,color=sitedepth))+
#  geom_line(size=1)+
#  labs(title="mean number of fish alive", x="week", y="average count")

###### Length by cage by week. barplots to see progression of each fish... #######

tl.gather <-melt(tl, id.vars=c("cageID","Site", "depth", "sitedepth", "Week"),measure.vars=c("Fish1","Fish2","Fish3","Fish4", "Fish5","Fish6","Fish7","Fish8","Fish9","Fish10"))

md <-filter(tl.gather,sitedepth=="MD")
ms <-filter(tl.gather,sitedepth=="MS")
nd <-filter(tl.gather,sitedepth=="ND")
ns <-filter(tl.gather,sitedepth=="NS")
sd <-filter(tl.gather,sitedepth=="SD")
ss <-filter(tl.gather,sitedepth=="SS")

#change the data source to create new plot
windows(30,20)
plot.md <- ggplot(data=md,aes(x=variable,y=value))+
  geom_bar(stat="identity")+ylim(0,100)+xlim("Fish1","Fish2","Fish3","Fish4","Fish5","Fish6","Fish7", "Fish8","Fish9","Fish10")+
  geom_text(aes(label=value),vjust=-0.3,size=2)+
facet_grid(Week~cageID)

#try to do the same but for histogram?
#this isn't working for skew/kurtosis because not enough data.
#it's basically just like the barplot. 
plot.ss <- ggplot(data=ss,aes(value))+
  geom_histogram(stat="bin", binwidth = 1)+
  facet_grid(Week~cageID)


###### VIRIDIS PLOTS for cages ################################################


###### calculate skewness & kurtosis by week for each SITEDEPTH######
#first get rid of zeros on tl.gather - important!
tl.gath <-filter(tl.gather, value != 0)
skewcalc <-ddply(tl.gath,Week~sitedepth,summarise,mean=mean(value),sd=sd(value),skew=skewness(value),kurt=kurtosis(value))

psc <- ggplot(data=skewcalc,aes(x=Week,y=skew))+
  geom_line(size=1)+
  geom_hline(yintercept=0,linetype="dashed",color="gray")+
  facet_wrap(~sitedepth,ncol=2)
psk <-ggplot(data=skewcalc,aes(x=Week,y=kurt))+
  geom_line(size=1)+
  geom_hline(yintercept=0,linetype="dashed",color="gray")+
  facet_wrap(~sitedepth,ncol=2)


#######################################NON CUMULATIVE NON SCALED DATASET ############################

# when we use the non-cumulative data: We actually don't want to set up the dataset using the svrt dataset. 
# instead it's going to be set up where every fish is observed every week, and just gets a status indicator as to whether it's alive or dead. 

cagetotal <-read.csv("cagetotals_17fixed.csv", header=TRUE) 

cgt1 <-melt(cagetotal, id.vars=c("CageID","Site", "depth", "sitedepth", "Week"),measure.vars=c("X1","X2","X3","X4", "X5","X6","X7","X8","X9","X10")) %>% group_by(Week) %>% group_by(CageID) %>% group_by(variable)

#create an individual fish id
cgt2 <-unite_(cgt1, "index", c("CageID","Week","variable"),remove=FALSE)
cgt2 <-unite_(cgt2, "fishID", c("CageID","variable"),remove=FALSE)
cgt2 <-dplyr::rename(cgt2,cage=CageID,site=Site,station=sitedepth,week=Week, Fish=variable,fishlength=value)

#So what we have here is an observation for every fish for every week. and what we need is environmental & summary data to go along with that. 
#Here "week zero" means the day the fish was caught. 

#we need to append biomass, numfish and avlength to cgt2
tl2 <-dplyr::select(tl,-Fish1,-Fish2,-Fish3,-Fish4,-Fish5,-Fish6,-Fish7,-Fish8,-Fish9,-Fish10, -w1,-w2,-w3,-w4,-w5,-w6,-w7,-w8,-w9,-w10) #week x cageID
tl2 <- dplyr::rename(tl2,station=sitedepth, week=Week)

#gather up all the environmental data where impacts are non-cumulative.
#incase you need to reset the creation of vert4. 
vert4 <-ddply(vert2,week~station,summarise,min.do=min(DO),max.do=max(DO),mean.do=mean(DO),sd.do=sd(DO),min.temp=min(temp), max.temp=max(temp), mean.temp=mean(temp),sd.T=sd(temp),do.dur=(sum(DO < 2)*15),temp.dur=(sum(temp >= 27)*15),mean.sal=mean(sal),sd.sal=sd(sal),ssat =(sum(perDO >= 115)*15))

vert4 <- mutate(vert4, week=as.integer(week)) #use vert4 because you're not doing cumulative data. 
# for 2017 environmental data starts on week 1. 
#vert4 <-mutate(vert4, week=week+1)# now the environmental data is properly synched to the survival data. 

svrt <- left_join(tl2,vert4,by=c("week","station")) 
svrt <- dplyr::rename(svrt, cage=cageID,site=Site, week = week) 

####################################### write svrt to csv for growth ################################################################
write.csv(svrt, file="svrt2017forgrowth_regenerated.csv")



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

tl4 <- mutate(tl4, death.week=n) #this i think is very different from the  
tl4 <-mutate(tl4, event = ifelse(death.week==12,3,NA)) %>% mutate(event2=ifelse(death.week==12,0,NA))#status = 0 if fish had 12 observations, meaning it survived the experiment, otherwise leave it blank. 

tl5 <-mutate(tl4, start=week, end=week+1)

tl5 <-mutate(tl5, event=ifelse(is.na(event) & death.week==end,1,3)) %>% mutate(event2=ifelse(is.na(event2) & death.week == end,1,0)) #if the end of the week of the of the observation = number of times animal was observed, the animal died in the preceding interval.
tl5 <-filter(tl5, week > 0)

View(tl5)

write.csv(tl5,file="weeklysummarydata_2017_regenerated.csv")


############################################MAKE THE SCALED DATASET ################################################

tl6 <- mutate_at(tl5, vars(biom, min.do, max.do, mean.do, sd.do, min.temp, sd.do, min.temp, max.temp, sd.T, do.dur, temp.dur, mean.sal, sd.sal, ssat), scale)

#write.csv(tl6,file="weeklysummScaledata_2017.csv")

########################################### MAKE THE CUMLATIVE VERSION #########################################################################

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

#vert5 <-mutate(vert5, week=week+1)# now the environmental data is properly synched to the survival data. NOT FOR 2017

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


tc4 <- mutate(tc4, death.week=n)# i think.
tc4 <-mutate(tc4, event = ifelse(death.week==12,3,NA)) %>% mutate(event2=ifelse(death.week==12,0,NA))#status = 0 if fish had 12 observations, meaning it survived the experiment, otherwise leave it blank. 

tc5 <-mutate(tc4, start=week, end=week+1)
tc5 <-mutate(tc5, event=ifelse(is.na(event) & death.week==end,1,3)) %>% mutate(event2=ifelse(is.na(event2) & death.week == end,1,0))#if the end of the week of the of the observation = number of times animal was observed, the animal died in the preceding interval.
tc5 <-filter(tc5, week > 0)

View(tc5)

#write.csv(tc5,file="cumsummarydata_2017.csv")




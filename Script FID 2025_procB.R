##############
# Workspace #
#############

setwd("C:/Users/jeffreym/OneDrive - Tartu ülikool/Dokumendid/Donnees chevreuils")

###########
# library #
###########

library(MASS)
library(MuMIn)
library(nlme)
library(lme4)
library(VGAM)
library(AICcmodavg)
library(rptR)
library(visreg)
library(lubridate)
library(car)
library(lmerTest)
library(dplyr)
library(interactions)

################
#  1.1. data   #
################
# Complete dataset
rm(list = ls())
DataManip=read.csv("Data-FID-2025.csv",header=T, dec=".")
summary(DataManip)
dim(DataManip) # 895 obs. of 82

# Transformations
DataManip$Nest=as.factor(DataManip$nest)
DataManip$Start=as.numeric(DataManip$SD_m)
DataManip$Flight=as.numeric(DataManip$FID_m)
DataManip$DateA=as.factor(DataManip$Date)
DataManip$Year=as.factor(DataManip$Year)
DataManip$Repet=as.numeric(DataManip$Repetition)
DataManip$Height=as.numeric(DataManip$height_m)
DataManip$DensityGull=as.numeric(DataManip$Count_Larcan_nests_Hannah)
DataManip$Area=as.factor(DataManip$Zone)
DataManip$Egg1=as.numeric(DataManip$Egg1)
DataManip$Egg2=as.numeric(DataManip$Egg2)
DataManip$Egg3=as.numeric(DataManip$Egg3)
DataManip$Egg4=as.numeric(DataManip$Egg4)
DataManip$Chick1=as.numeric(DataManip$Chick1)
DataManip$Chick2=as.numeric(DataManip$Chick2)
DataManip$Chick3=as.numeric(DataManip$Chick3)
DataManip$Chick4=as.numeric(DataManip$Chick4)
DataManip$EggHatched=as.numeric(DataManip$EggHatched)
DataManip$Clutch=as.numeric(DataManip$ClutchMass)
DataManip$Brood=as.numeric(DataManip$ChickNmass)
DataManip$Success=as.numeric(DataManip$HatchSucc)
DataManip$ChickP=as.numeric(DataManip$ProbaChick)
DataManip$EggNest=as.factor(DataManip$TotEgg)
DataManip$FID_Blup=as.numeric(DataManip$BLUP_25)
DataManip$FirstEgg=as.factor(DataManip$first_egg_date)
attach(DataManip)

# Creattion of "Julian date"
DataManip$DateJ=strptime(DataManip$DateA,"%d/%m/%Y")$yday+1
DataManip$DateJ
DataManip$DateJ=as.numeric(DataManip$DateJ)
attach(DataManip)

# Creation of Laying date of 1st egg since 1st of January
DataManip$DateEgg1=strptime(DataManip$FirstEgg,"%d/%m/%Y")$yday+1
DataManip$DateEgg1
DataManip$DateEgg1=as.numeric(DataManip$DateEgg1)
attach(DataManip)

# Splitting areas between North (more exposed to humans) and other parts
DataManip$ProxiH[DataManip$Area == "Head"] <- "Less used"
DataManip$ProxiH[DataManip$Area == "North"] <- "More used"
DataManip$ProxiH[DataManip$Area == "Leg"] <- "Less used"
DataManip$ProxiH[DataManip$Area == "Middle"] <- "Less used"
DataManip$ProxiH=as.factor(DataManip$ProxiH)

#########################
### 1. Repeatability ###
########################
DataRep=DataManip[,c("Nest","Flight","Repet","Start","DateJ","Year")]
DataRepsNA=na.omit(DataRep) # 879 observations

# Manual calculation of corrected repeatability
lmerRep=lmer(Flight ~ Repet + Start + DateJ + (1|Nest), data = DataRepsNA)

plot(lmerRep)
qqnorm(residuals(lmerRep))
hist(residuals(lmerRep))
vif(lmerRep)
summary(lmerRep)

repeatability = 45.42/(45.42+45.52)
repeatability # 0.50

# Extraction of BLUPs for FIDs
df_BLUPS_B <- data_frame(Nest = row.names(ranef(lmerRep)$Nest),
                         BLUP_B = ranef(lmerRep)$Nest[,"(Intercept)"])

# Rptr package for repeatability
rep=rptGaussian(Flight ~ (1|Nest), grname = "Nest",
                data = DataRepsNA, nboot = 1000, npermut = 1000)

print(rep) # R = 0.52 without correction; CI [0.45, 0.58]
summary(rep)
plot(rep, main="Repeatability of Flight initiation distance")

## corrected for starting distance, repetition number, Julian date ##
rep=rptGaussian(Flight ~ Start+Repet+DateJ+(1|Nest), grname = "Nest",
                data = DataRepsNA, nboot = 1000, npermut = 1000)

print(rep) # R = 0.50 with correction; CI [0.43, 0.57]
summary(rep)
plot(rep, main="Repeatability of Flight initiation distance")

#######################################
### 2. FID - Nesting site features ###
######################################

########################
### 2.1 Basic model ###
#######################

DataTest=DataManip[,c("Nest","Flight","Height","Repet","Start","DensityGull","ProxiH","Year")]
DataTestsNA=na.omit(DataTest) # 879 observations

# Data check
summary(DataTestsNA)
str(DataTestsNA)

# Model
lmtest= lmer(Flight~ Height + Repet + Start + DensityGull + ProxiH
             + DensityGull*Height + DensityGull*ProxiH
             + (1|Nest) + (1|Year),
             data = DataTestsNA, REML=F,na.action=na.fail)

vif(lmtest)
Anova(lmtest, type = 2)
summary(lmtest)

# Model selection
options(na.action = na.fail)
dd <- dredge(lmtest)
d=subset(dd, delta < 2)
d

nst=nested(dd)
d=subset(dd, delta < 2 & !nested(.)) # table with only non-nested models
d

# Results
avgm <- model.avg(dd,subset=c("28","32"))
summary(avgm)
confint(avgm)

# Final model
lmerfin= lmer(Flight~ Height + Repet + Start + DensityGull + ProxiH
              + (1|Nest) + (1|Year),
              data = DataTestsNA, REML=F,na.action=na.fail)

r.squaredGLMM(lmerfin) # 0.26 & 0.64
plot(lmerfin,ask=TRUE) ## ok
hist(residuals(lmerfin))

## Plot for height
newdata=expand.grid(Height=seq(min(DataTestsNA$Height),
                                 max(DataTestsNA$Height),0.01),
                    Repet=mean(DataTestsNA$Repet),
                    Start=mean(DataTestsNA$Start),
                    DensityGull=mean(DataTestsNA$DensityGull),
                    ProxiH=levels(DataTestsNA$ProxiH)[1])

pred=predictSE(lmerfin,newdata=newdata,se.fit=TRUE)
newdata=cbind(newdata,pred)
newdata$low=newdata$fit-1.96*newdata$se.fit
newdata$upp=newdata$fit+1.96*newdata$se.fit

plot(Flight~Height,data=DataTestsNA,type="n",
     #xlim=c(0.3,1.8),ylim=c(-5,70),
     xaxt='n', yaxt='n',
     xlab="Height of the nesting site (m)",
     ylab="Flight initiation distance (m)",
     cex.lab=1.5,
     cex.axis=1.5)
axis(side=1, at=seq(0, 1.8, by=0.2))
axis(side=2, at=seq(-5, 70, by=10))
polygon(c(rev(newdata$Height),newdata$Height),c(rev(newdata$low),newdata$upp),
        col = adjustcolor("grey90", alpha.f = 0.7),
        border=NA)
points(Flight~Height,data=DataTestsNA,pch=20,lwd=4,col="dodgerblue")
lines(newdata$fit~newdata$Height,col="black",lwd=4)
lines(newdata$low~newdata$Height,col="black",lwd=3,lty=4)
lines(newdata$upp~newdata$Height,col="black",lwd=3,lty=4)

## Plot for density
newdata=expand.grid(DensityGull=seq(min(DataTestsNA$DensityGull),
                                    max(DataTestsNA$DensityGull),0.1),
                    Repet=mean(DataTestsNA$Repet),
                    Start=mean(DataTestsNA$Start),
                    Height=mean(DataTestsNA$Height),
                    ProxiH=levels(DataTestsNA$ProxiH)[1])

pred=predictSE(lmerfin,newdata=newdata,se.fit=TRUE)
newdata=cbind(newdata,pred)
newdata$low=newdata$fit-1.96*newdata$se.fit
newdata$upp=newdata$fit+1.96*newdata$se.fit

plot(Flight~DensityGull,data=DataTestsNA,type="n",
     #xlim=c(0,20),ylim=c(-10,70),
     xlab="DensityGull of Common gull in a 5-meter radius",
     ylab="Flight initiation distance (m)",
     cex.lab=1.5,
     cex.axis=1.5)
polygon(c(rev(newdata$DensityGull),newdata$DensityGull),c(rev(newdata$low),newdata$upp),
        col = adjustcolor("grey90", alpha.f = 0.7),
        border=NA)
points(Flight~DensityGull,data=DataTestsNA,pch=20,lwd=4,col="dodgerblue")
lines(newdata$fit~newdata$DensityGull,col="black",lwd=4)
lines(newdata$low~newdata$DensityGull,col="black",lwd=3,lty=4)
lines(newdata$upp~newdata$DensityGull,col="black",lwd=3,lty=4)

#########################
### 3. FID - Fitness ###
########################

###################################
### 3.1 Number of eggs hatched ###
##################################

DataTest=DataManip[,c("Nest","Flight","Height","DensityGull","Repet","Start","Year",
                      "EggHatched","FID_Blup","ProxiH")]

DataTestsNA=na.omit(DataTest) #239 observations
DataTestsNA=subset(DataTestsNA,Nest!="24_24" & Nest!="42_24") #eggs non-incubated on these 2 nests due to external factors

glm1=glm(EggHatched ~ FID_Blup + Height + DensityGull + Year
         + FID_Blup*Height + FID_Blup*DensityGull,
         data = DataTestsNA, family = poisson(link="log"), na.action=na.fail)

vif(glm1)
Anova(glm1, type = 2)
summary(glm1)

# Model selection
options(na.action = na.fail)
dd <- dredge(glm1)
d=subset(dd, delta < 2)
d

# Results
glmfin=glm(EggHatched ~ 1,data = DataTestsNA, family = poisson(link="log"), na.action=na.fail)

summary(glmfin)

#################################################
### 3.2 Probability to have at least 1 chick ###
################################################

DataTest=DataManip[,c("Nest","Flight","Height","DensityGull","Repet","Start","Year",
                      "ChickP","EggNest","FID_Blup","ProxiH")]

DataTestsNA=na.omit(DataTest) #239 observations
DataTestsNA=subset(DataTestsNA,Nest!="24_24" & Nest!="42_24") #eggs non-incubated on these 2 nests due to external factors

glm1=glm(ChickP ~ FID_Blup + Height + DensityGull + Year
         + FID_Blup*Height + FID_Blup*DensityGull,
         data = DataTestsNA, family=binomial(link="logit"))

vif(glm1)
Anova(glm1, type = 2)
summary(glm1)

# Model selection
options(na.action = na.fail)
dd <- dredge(glm1)
d=subset(dd, delta < 2)
d

# Results
glmfin=glm(ChickP ~ Year,
           data = DataTestsNA, family=binomial(link="logit"))

Anova(glmfin)
summary(glmfin)
confint(glmfin)

########################
### 3.3 Clutch mass ###
#######################

DataTest=DataManip[,c("Nest","Flight","Height","DensityGull","Repet","Start","Year",
                      "Clutch","EggNest","FID_Blup","ProxiH")]

DataTestsNA=na.omit(DataTest) #239 observations
DataTestsNA=subset(DataTestsNA,EggNest!="4" & EggNest!="5" & EggNest!="6") #without second clutches / 226 obs.

# Model 1
lmtest= lm(Clutch ~ FID_Blup + Height + DensityGull + Year
           + FID_Blup*Height + FID_Blup*DensityGull,
           data = DataTestsNA)

vif(lmtest)
Anova(lmtest, type = 2)
summary(lmtest)

# Model selection
options(na.action = na.fail)
dd <- dredge(lmtest)
d=subset(dd, delta < 2)
d

nst=nested(dd)
d=subset(dd, delta < 2 & !nested(.)) # table with only non-nested models
d

# Results
avgm <- model.avg(dd,subset=c("47","39"))
summary(avgm)
confint(avgm)

# Model 1
lmfin= lm(Clutch ~ FID_Blup + Height + Year
          +FID_Blup*Height,
          data = DataTestsNA)

Anova(lmfin, type = 3)
summary(lmfin)
hist(residuals(lmfin))
plot(lmfin)

## Clutch-FID-Height graph ##
POvJa<-subset(DataTestsNA,Height<0.7)
POvJc<-subset(DataTestsNA,Height>0.7)

## POvJa ##
lmer2 = lm(Clutch ~ FID_Blup + Height + Year
           +FID_Blup*Height,
           data = DataTestsNA)

newdata=expand.grid(FID_Blup=seq(min(POvJa$FID_Blup),
                               max(POvJa$FID_Blup),3),
                    Year=levels(DataTestsNA$Year)[2],
                    Height=mean(POvJa$Height))

pred=predict(lmer2,newdata=newdata,se.fit=TRUE)
newdata=cbind(newdata,pred)
newdata$low=newdata$fit-1.96*newdata$se.fit
newdata$upp=newdata$fit+1.96*newdata$se.fit

plot(Clutch~FID_Blup,data=DataTestsNA,type="n",
     xlim=c(-15,25),ylim=c(10,250),
     xlab="FID BLUPs (Bold to shy)",
     ylab="Total clutch mass (m)",
     cex.lab=1.5,
     cex.axis=1.5)

legend(x=12,y=240, inset=.01, legend=c("Height > 0.7", "Height < 0.7"),
       col=c("darkred", "dodgerblue"), bty="n", lty=1, lwd=3, cex=0.9, text.font=4)

polygon(c(rev(newdata$FID_Blup),newdata$FID_Blup),c(rev(newdata$low),newdata$upp),
        col = adjustcolor("lightblue", alpha.f = 0.7),
        border=NA)

points(Clutch~FID_Blup,data=POvJa,pch=20,lwd=4,col="dodgerblue")
lines(newdata$fit~newdata$FID_Blup,col="dodgerblue",lwd=4)
lines(newdata$low~newdata$FID_Blup,col="dodgerblue",lwd=3,lty=4)
lines(newdata$upp~newdata$FID_Blup,col="dodgerblue",lwd=3,lty=4)

## POvJc ##
newdata=expand.grid(FID_Blup=seq(min(POvJc$FID_Blup),
                                 max(POvJc$FID_Blup),2),
                    Year=levels(DataTestsNA$Year)[2],
                    Height=mean(POvJc$Height))

pred=predict(lmer2,newdata=newdata,se.fit=TRUE)
newdata=cbind(newdata,pred)
newdata$low=newdata$fit-1.96*newdata$se.fit
newdata$upp=newdata$fit+1.96*newdata$se.fit

polygon(c(rev(newdata$FID_Blup),newdata$FID_Blup),c(rev(newdata$low),newdata$upp),
        col = adjustcolor("pink", alpha.f = 0.7),
        border=NA)

points(Clutch~FID_Blup,data=POvJc,pch=20,lwd=4,col="darkred")
lines(newdata$fit~newdata$FID_Blup,col="darkred",lwd=4)
lines(newdata$low~newdata$FID_Blup,col="darkred",lwd=3,lty=4)
lines(newdata$upp~newdata$FID_Blup,col="darkred",lwd=3,lty=4)

#######################
### 3.4 Brood mass ###
######################

DataTest=DataManip[,c("Nest","Flight","Height","DensityGull","Repet","Start","Year",
                      "Brood","EggHatched","FID_Blup","ProxiH","EggNest")]

DataTestsNA=na.omit(DataTest) # 239 observations
DataTestsNA=subset(DataTestsNA,Nest!="24_24" & Nest!="42_24" & EggNest!="4" & EggNest!="5" & EggNest!="6")
#External factors + second clutches removed - 224 obs.

# Model
lmtest= lm(Brood ~ FID_Blup + Height + DensityGull + Year
           + FID_Blup*Height + FID_Blup*DensityGull,
           data = DataTestsNA)

vif(lmtest)
Anova(lmtest, type = 2)
summary(lmtest)

# Model selection
options(na.action = na.fail)
dd <- dredge(lmtest)
d=subset(dd, delta < 2)
d

lmfin= lm(Brood ~ Year, data = DataTestsNA)

Anova(lmfin)
summary(lmfin)
confint(lmfin)
plot(lmfin)

#################################
### 3.5 First egg laying date ###
#################################

DataTest=DataManip[,c("Nest","Height","DensityGull","Repet","Start","Year",
                      "DateEgg1","FID_Blup","ProxiH")]

DataTestsNA=na.omit(DataTest) #239 observations

# Model 1
lmtest= lm(DateEgg1 ~ FID_Blup + Height + DensityGull + Year +
           + FID_Blup*Height + FID_Blup*DensityGull,
           data = DataTestsNA)

vif(lmtest)
Anova(lmtest, type = 3)
summary(lmtest)

# Model selection
options(na.action = na.fail)
dd <- dredge(lmtest)
d=subset(dd, delta < 2)
d

# Model final
lmfin= lm(DateEgg1 ~ FID_Blup + DensityGull + Year
          + FID_Blup*DensityGull,
          data = DataTestsNA)

Anova(lmfin, type = 3)
summary(lmfin)
hist(residuals(lmfin))
plot(lmfin)

## Laying date-FID-Density graph ##
POvJa<-subset(DataTestsNA,DensityGull<5)
POvJb<-subset(DataTestsNA,DensityGull>5 & DensityGull<10)
POvJc<-subset(DataTestsNA,DensityGull>10)

## POvJa ##
lmer2 = lm(DateEgg1 ~ FID_Blup + DensityGull + Year
           + FID_Blup*DensityGull,
           data = DataTestsNA)

newdata=expand.grid(FID_Blup=seq(min(POvJa$FID_Blup),
                                 max(POvJa$FID_Blup),3),
                    Year=levels(DataTestsNA$Year)[2],
                    DensityGull=mean(POvJa$DensityGull))

pred=predict(lmer2,newdata=newdata,se.fit=TRUE)
newdata=cbind(newdata,pred)
newdata$low=newdata$fit-1.96*newdata$se.fit
newdata$upp=newdata$fit+1.96*newdata$se.fit

plot(DateEgg1~FID_Blup,data=DataTestsNA,type="n",
     xlim=c(-15,30),ylim=c(119,131),
     xlab="FID BLUPs (Bold to shy)",
     ylab="First egg laying date (since 1st of January)",
     cex.lab=1.5,
     cex.axis=1.5)

legend(x=12,y=131, inset=.01, legend=c("...", "...", "..."),
       col=c("black", "darkred", "dodgerblue"), bty="n", lty=1, lwd=3, cex=0.9, text.font=4)

polygon(c(rev(newdata$FID_Blup),newdata$FID_Blup),c(rev(newdata$low),newdata$upp),
        col = adjustcolor("lightblue", alpha.f = 0.7),
        border=NA)

points(DateEgg1~FID_Blup,data=POvJa,pch=20,lwd=4,col="dodgerblue")
lines(newdata$fit~newdata$FID_Blup,col="dodgerblue",lwd=4)
lines(newdata$low~newdata$FID_Blup,col="dodgerblue",lwd=3,lty=4)
lines(newdata$upp~newdata$FID_Blup,col="dodgerblue",lwd=3,lty=4)

## POvJb ##
newdata=expand.grid(FID_Blup=seq(min(POvJb$FID_Blup),
                                 max(POvJb$FID_Blup),3),
                    Year=levels(DataTestsNA$Year)[2],
                    DensityGull=mean(POvJb$DensityGull))

pred=predict(lmer2,newdata=newdata,se.fit=TRUE)
newdata=cbind(newdata,pred)
newdata$low=newdata$fit-1.96*newdata$se.fit
newdata$upp=newdata$fit+1.96*newdata$se.fit

polygon(c(rev(newdata$FID_Blup),newdata$FID_Blup),c(rev(newdata$low),newdata$upp),
        col = adjustcolor("pink", alpha.f = 0.7),
        border=NA)

points(DateEgg1~FID_Blup,data=POvJb,pch=20,lwd=4,col="darkred")
lines(newdata$fit~newdata$FID_Blup,col="darkred",lwd=4)
lines(newdata$low~newdata$FID_Blup,col="darkred",lwd=3,lty=4)
lines(newdata$upp~newdata$FID_Blup,col="darkred",lwd=3,lty=4)

## POvJc ##
newdata=expand.grid(FID_Blup=seq(min(POvJc$FID_Blup),
                                 max(POvJc$FID_Blup),3),
                    Year=levels(DataTestsNA$Year)[2],
                    DensityGull=mean(POvJc$DensityGull))

pred=predict(lmer2,newdata=newdata,se.fit=TRUE)
newdata=cbind(newdata,pred)
newdata$low=newdata$fit-1.96*newdata$se.fit
newdata$upp=newdata$fit+1.96*newdata$se.fit

polygon(c(rev(newdata$FID_Blup),newdata$FID_Blup),c(rev(newdata$low),newdata$upp),
        col = adjustcolor("grey", alpha.f = 0.7),
        border=NA)

points(DateEgg1~FID_Blup,data=POvJc,pch=20,lwd=4,col="black")
lines(newdata$fit~newdata$FID_Blup,col="black",lwd=4)
lines(newdata$low~newdata$FID_Blup,col="black",lwd=3,lty=4)
lines(newdata$upp~newdata$FID_Blup,col="black",lwd=3,lty=4)

#####################
### 3.6 Egg mass ###
####################

# Dataset
rm(list = ls())
DataManip=read.csv("FID-Egg-Chick-23-24-25.csv",header=T, dec=".")
summary(DataManip)
dim(DataManip) # 2090 obs. of 20 variables
attach(DataManip)

# Transformations
DataManip$Nest=as.factor(DataManip$nest)
DataManip$Flight=as.numeric(DataManip$FID_m)
DataManip$Year=as.factor(DataManip$Year)
DataManip$DensityGull=as.numeric(DataManip$Count_Larcan_nests_Hannah)
DataManip$Dist=as.numeric(DataManip$Distance_coast)
DataManip$Height=as.numeric(DataManip$height_m)
DataManip$Area=as.factor(DataManip$Zone)
DataManip$Egg=as.numeric(DataManip$Egg)
DataManip$EggNb=as.factor(DataManip$EggNb)
DataManip$FID_Blup=as.numeric(DataManip$BLUP)
DataManip$Growth=as.numeric(DataManip$Growth)
DataManip$ChickMass=as.numeric(DataManip$Chick)

attach(DataManip)

DataTest=DataManip[,c("Nest","Height","DensityGull","Year",
                      "Egg","EggNb","FID_Blup","Area","ProxiH")]
DataTestsNA=na.omit(DataTest) #657 observations

# Model 1
DataTestsNA=subset(DataTestsNA,EggNb!="4") #without 4th egg

lmertest= lmer(Egg ~ FID_Blup + EggNb + Height + DensityGull
               + FID_Blup*Height + FID_Blup*DensityGull
               + (1|Nest) + (1|Year),
               data = DataTestsNA, REML=F,na.action=na.fail)

vif(lmertest)
Anova(lmertest)
summary(lmertest)

# Model selection
options(na.action = na.fail)
dd <- dredge(lmertest)
d=subset(dd, delta < 2)
d

# Model 1
lmerfin= lmer(Egg ~ EggNb
              + (1|Nest) + (1|Year),
              data = DataTestsNA, REML=F,na.action=na.fail)

summary(lmerfin)
r.squaredGLMM(lmerfin) # 0.09 & 0.78

plot(lmerfin)
visreg(lmerfin)

##########################
### 3.6 Mass at birth ###
#########################

DataTest=DataManip[,c("Nest","Height","DensityGull","Year",
                      "ChickMass","EggNb","FID_Blup","Area","ProxiH")]

DataTestsNA=na.omit(DataTest) #454 observations

# Model 1
DataTestsNA=subset(DataTestsNA,EggNb!="4") #without 4th egg

lmtest= lmer(ChickMass ~ FID_Blup + EggNb + Height + DensityGull
             + FID_Blup*Height + FID_Blup*DensityGull
             + (1|Nest) + (1|Year),
             data = DataTestsNA, REML=F,na.action=na.fail)

vif(lmtest)
Anova(lmtest)
summary(lmtest)

# Model selection
options(na.action = na.fail)
dd <- dredge(lmtest)
d=subset(dd, delta < 2)
d

# Results
lmfin= lmer(ChickMass ~ EggNb + (1|Nest) + (1|Year),
            data = DataTestsNA, REML=F,na.action=na.fail)

Anova(lmfin)
summary(lmfin)
confint(lmfin)
r.squaredGLMM(lmerfin) # 0.10 & 0.78

plot(lmfin)
visreg(lmfin)
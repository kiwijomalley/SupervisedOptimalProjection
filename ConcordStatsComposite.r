## Developed by James O'Malley June 29, 2024 ##

## Step 2 of data analysis process 
##  ConcordStatsComposite.r estimates hierarchical p2 and p2 models for relational data
##  to establish whether the direction of patient-sharing that has the greater weight is the more 
##   likely to correspond to a true nomination.
##  Also evaluates descriptive statistics

## Goal: refine concordance analysis of NPO survey nominations and patient-sharing in Medicare ##
## Runs on the 40 different ways of forming weighted physician networks ##
## Can explore impact of using different thresholds to classified an edge as existing, and evaluate ##
##  whether optimal thresholded smart network constructions can do better than optimal thresholded ##
##  basic methods ##
## The data sets for each type include the survey respondents as npi1 and have varying definitions 
##  for who they were potentially eligible to nominate

#install.packages("lme4",repos=NULL,contriburl="file:/opt/software/cran/src/contrib")
#install.packages("nloptr",repos=NULL,contriburl="file:/opt/software/cran/src/contrib")
#install.packages("statnet",repos=NULL,contriburl="file:/opt/software/cran/src/contrib")
#install.packages("amen",repos=NULL,contriburl="file:/opt/software/cran/src/contrib")
#install.packages("latentnet",repos=NULL,contriburl="file:/opt/software/cran/src/contrib")
#install.packages("eigenmodel",repos=NULL,contriburl="file:/opt/software/cran/src/contrib")
#install.packages("coda",repos=NULL,contriburl="file:/opt/software/cran/src/contrib")
#install.packages("geepack",repos=NULL,contriburl="file:/opt/software/cran/src/contrib")
#install.packages("igraph",repos=NULL,contriburl="file:/opt/software/cran/src/contrib")
#install.packages("dyads",repos=NULL,contriburl="file:/opt/software/cran/src/contrib")
#install.packages("plyr",repos=NULL,contriburl="file:/opt/software/cran/src/contrib")

library(lme4)
library(nloptr)
library(foreign)
library(statnet)
library(amen)
library(latentnet)
library(eigenmodel)
library(coda)
library(geepack)
library(survival)
library(igraph)
library(dyads)
library(plyr)

rsource <- "/drives/drive1/home/f001580/Desktop/drives/54054-Linux/54054dua/programs/jomalley/"
setwd(rsource)
datdir="../../idata/jomalley/"
datdirpy="../../idata/jomalley/PHN-subnetwork/2018/batch-1/"
outdir="../../idata/jomalley/"

## Data-Dictionary of NPO data merged with CMS data ##
#npi1: cms npi of physician 1
#npi2: cms npi of physician 2
#dyadid: unique id of dyad
#upat: number of shared patients
#gvisit: Geometric mean of number of patient encounters with physicians of shared patients
#w1,...,w40: weights based on the 40 ways of building weighted physician networks from CMS data
#tieid: tieid using npis (can drop)
#InfOwn: Information wanted at own hospital
#RefOwn: Referral at own hospital
#InfOth: Information wanted other hospital
#RefOth: Referral other hospital      
#named: Any survey nomination (any of the above four)
#HospRespName: Name of hospital of survey respondent
#HospNomID: Numerical ID of hospital of nominated physician
#SpecialtyNom: Specialty of nominated physician
#id1: ID of physician 1 (starting from 1)      
#id2: ID of physician 2 (starting from 1)
#dyadstatus: 2 = mutual, 1 = directed, 0 = null over all nominations
#mInfOwn: dyadstatus for information wanted own hospital (2, 1, 0)
#mRefOwn: dyadstatus for referral own hospital (2, 1, 0)
#mInfOth: dyadstatus for information wanted other hospital (2, 1, 0)
#mRefOth: dyadstatus for referral other hospital (2, 1, 0)
## The next four variables originate from the old way of forming networks; they do not seem to be used at all! ##
#nominee1: whether physician 1 was a nominee (1 = yes, 0 = no). 0 = responding physician takes precedence
#sound1: whether physician employed by Sound 
#nominee2: whether 2 was a nominee (1 = yes, 0 = no). 0 = responding physician takes precedence     
#sound2: whether physician 2 was employed by Sound

#The following variables are from the NPO sampling frame and so are missing for non-NPO nominees
#  The suffix "k" takes the values 1 for physician 1 and 2 for physician 2
# hospidk = Hospital ID of Sound physician 
# EmploymentTypek = Regular (W2) or various types of Consultant (1099)
# Statusk = Active or not (inactive physicians already excluded)
# ColleagueTypek = Role (most are physician, only physicians were in survey)
# Suffixk = MD or DO
# in_encounterk = Flag for whether record is in encounter data?
# in_mipsk = Flag for whether record is in the MIPS?
# nAffilsk = Number of hospital affiliations
# Regionk = Sound Region of Country
# TeachingStatusk = Teaching status of hospital
# Localityk = U (urban) verus R (rural)
# Bedsizek = Size of hospital as reflected by number of beds
# nSoundPhysk = Number of Sound physicians at hospital
# nRegionk = Numerical code for region of country
# Teachingk = 1 if teaching, 0 otherwise
# Ruralk = 1 if rural, 0 otherwise
# Bed_Sizek = Numerical version of Bedsize variable
# Sampledk = 1 if among 600 sampled physicians, 0 otherwise

diagprops <- function(gold,test,N=1e7) {
 th <- sort(unique(test)) #smallest to biggest
 print(length(th))
 nth <- length(th)
 sens <- rep(0,nth); spec <- rep(0,nth)
 for (i in 1:nth) {
  sens[i] <- mean(test[gold==1]>th[i])
  spec[i] <- mean(test[gold==0]<=th[i])
 }
 mis <- sens+spec-1
 Jindex <- max(mis) #Youden Index or maximum improvement over chance
 pos <- order(-mis)[1]
 optth <- th[pos] #Optimal threshold on original scale (e.g., number of shared patients)
 sensopt <- sens[pos]
 specopt <- spec[pos]
 ppvopt <- mean(gold[test>optth])
 npvopt <- 1-mean(gold[test<=optth])
 aucprob <- auc_probability(as.logical(gold),test,N)
 aucarea <- auc_area(sens,spec)
 return(c(optth,sensopt,specopt,ppvopt,npvopt,nth,Jindex,aucprob,aucarea))
}

auc_probability <- function(gold,test,N=1e7) {
 pos <- sample(test[gold], N, replace=TRUE)
 neg <- sample(test[!gold], N, replace=TRUE)
 # sum( (1 + sign(pos - neg))/2)/N # does the same thing
 return((sum(pos > neg) + sum(pos == neg)/2) / N) # give partial credit for ties
}

auc_area <- function(sens,spec) {
 n <- length(sens)
 fpr <- 1-spec
 sens <- c(1,sens); fpr <- c(1,fpr)
 return(sum((sens[1:n]+sens[2:(n+1)])*(fpr[1:n]-fpr[2:(n+1)]))/2)
}

#Build models that predict directed dyads that use directional shared-patient edge weights
# That is, the composite prediction in which one edge exists and one doesn't 

#Make or load tie-level data set containing CMS and survey information for concordance analysis
# This has an edgelist structure in that only ties that are present are listed?
makedata <- 0 #Avoid reforming data set
type <- 1 #Which data set to use to evaluate diagnostic properties. Use 1 and 3 for p2ML and p2, respectively
if (makedata==1) {
 source('ConcordDataWrangle.r')
}
if (type==1) {
 data <- read.csv(paste(datdir,"ConcordSameHospPristine_19.csv",sep=""),header=TRUE)
} else if (type==2) {
 data <- read.csv(paste(datdir,"ConcordSameHospRough_19.csv",sep=""),header=TRUE)
} else if (type==3) {
 data <- read.csv(paste(datdir,"ConcordAllSoundBilling_19.csv",sep=""),header=TRUE)
} else if (type==4) {
 #data <- read.csv(paste(datdir,"ConcordanceTieData19.csv",sep=""),header=TRUE)
 data <- read.csv(paste(datdir,"ConcordAllSound_19.csv",sep=""),header=TRUE) 
}
sensdata <- read.csv(paste(datdir,"ConcordSensitivityAll_19.csv",sep=""),header=TRUE) #Not currently using

#Add the number of patients seen on the same day by the pair of providers as an additional measure
data$nsameday <- data$w6+data$w6r-data$w30
data$wtsameday <- data$w7+data$w7r-data$w28 #Not as meaningful as same-day as mixing concepts

#Decide whether to subset ties for which both survey respondent and nominee are in CMS 
table(data$nocmsresp,data$nocmsnom)
dropextr <- 1
if (dropextr==1) {
 data <- data[(data$nocmsresp==0 & data$nocmsnom==0),]
}

#Form the 4 data sets for internal nominations
data$intsurvnom <- ifelse(data$InfOwn+data$RefOwn>0,1,0)
data$extsurvnom <- ifelse(data$InfOth+data$RefOth>0,1,0) #Ideally, there wouldn't be any of these!

#Form label of a mutual eligible dyad (both physicians have to respond)
freqdyad <- table(data$dyadid)
mutdyadid <- as.numeric(names(freqdyad))
mutdyadid <- mutdyadid[(freqdyad==2)]
data$mutelig <- 1*(data$dyadid %in% mutdyadid) #Potential mutual dyads (reciprocated ties): both physicians responded

#Form indicator of a mutual dyad
freqdyadnom <- tapply(data$named,data$dyadid,'sum') #Shows if dyad had 0, 1, or 2 nominations
mutdyadid <- as.numeric(names(freqdyadnom))
freqdyadnom <- as.vector(freqdyadnom)
mutdata <- data.frame(dyadid=mutdyadid,dyadtype=freqdyadnom)
data <- merge(data,mutdata,by="dyadid",all.x=TRUE)

table(data$mutelig)
table(data$dyadtype)
table(data$named)
table(data$mutelig,data$dyadtype)
table(data$dyadtype,data$named) #The 25+ ties with named=0 and mutual=1 are reverse edges 
# (respondent didn't name but someone named them and this shows up in another nomination)

#Restrict to dyads eligible to the dyadic (both physicians responded) for types 1 and 2
mutdata <- data[data$mutelig==1,]
table(mutdata$dyadtype)
table(mutdata$named)
table(mutdata$dyadtype,mutdata$named) 
table(mutdata$shospid1==mutdata$shospid2) #Checks that only consider within hospital nominations for type 1 and type 2

## Evaluation of descriptive statistics ##

DiagStats <- read.csv(paste(outdir,"DiagStatsInPristine_19.csv",sep=""),header=TRUE) #Data to get optima theta by design
Joptthts <- DiagStats$optth_s; Joptthtr <- DiagStats$optth_r
tht <- Joptthts[1:20]; thtr <- Joptthtr[1:20]
Jopt <- 0
if (Jopt==0) {
 tht <- rep(0,20) #Threshold for patient-sharing; unfair to some designs to use fixed-threshold (using threshold based on J would be fairer)
}

#Strong sensitivity and specificity
mutdyads <- mutdata[mutdata$dyadtype==2,] #Any name generator of mutual dyad status
#Strong specificity
nulldyads <- mutdata[mutdata$dyadtype==0,] #Any name generator of mutual dyad status
#Mutual sensitivity or ddauc: want directed dyads
dirdyads <- mutdata[mutdata$dyadtype==1,] #Any name generator of survey-based directional dyad status

#y <- dirdyads$named #Any name generator (within and external to hospital)
y <- dirdyads$intsurvnom #Within-hospital name generator of survey-based edge status
strsens <- rep(0,nrow=20); strspec <- rep(0,nrow=20); infotest <- rep(0,20); edgePPV <- rep(0,20); dyadPPV <- rep(0,20)
ddauc <- matrix(0,nrow=20,ncol=3)
for (i in 1:nrow(ddauc)) {
 x <- paste("w",i,sep="")
 strsens[i] <- mean(mutdyads[,x]>tht[i])
 strspec[i] <- mean(nulldyads[,x]<=tht[i])
 xr <- paste(x,"r",sep="")
 up1 <- mean(dirdyads[,x])+sqrt(var(dirdyads[,x])) #Threshold at which excess difference occurs
 ddauc[i,1] <- mean(dirdyads[y==1,x]>dirdyads[y==1,xr])+0.5*mean(dirdyads[y==1,x]==dirdyads[y==1,xr]) #directional AUC
 ddauc[i,2] <- (mean(dirdyads[y==1,x])-mean(dirdyads[y==1,xr]))/sqrt(var(dirdyads[,x])) #Mean difference
 ddauc[i,3] <- mean(dirdyads[y==1,x]-dirdyads[y==1,xr]>up1)-mean(dirdyads[y==1,xr]-dirdyads[y==1,x]>up1) #excess diff (1 sd)
 infotest[i] <- strsens[i]-mean(dirdyads[y==1,x]>tht[i])
 #Edge PPV
 patsharemut <- mutdata[(dirdyads[y==1,x]>tht[i] & dirdyads[y==1,xr]>tht[i]),]
 edgePPV[i] <- mean(patsharemut$intsurvnom)
 #Dyad PPV
 dyadPPV[i] <- mean(patsharemut$dyadtype==2) #Technically, this is for any name generator, not within-hospital
}

#Output descriptive statistics results
outdesc <- rbind(strsens,strspec,t(ddauc),infotest,edgePPV,dyadPPV)
dirname <- paste(outdir,"DescriptiveTests1.csv",sep="") #1 for type==1, 3 for type==3
#write.csv(outdesc,file=dirname,row.names=FALSE)

#Use dyads package to estimate a hierarchical P2 model (p2ML procedure)
# Define X and R (D1 and D2 in the package) as matrices whose elements are the directed and undirected 
# patient-sharing networks. Ignore actor and receiver covariates but allow random-effects
# Use list to combine matrices of possibly different dimension (plyr function handles this)!

#Check if dyads span hospitals. they shouldn't for types 1 and 2 but will for 3 and 4!
ndyad1 <- tapply(mutdata$shospid1,mutdata$dyadid,'mean')
#print(ndyad1) #This is invariant

#Form network objects based on shospid1 (=shospid2)
# First sort by shospid 1 and form overall survey-based network (Y) and shared-patient network (X) for w7 (used to use w8)

#Network datasets for specific hospitals and for NPO
if (type %in% c(1,2)) {
 netdata <- mutdata[order(mutdata$npi1),] #To be consistent with order of physician-level covariates when using tapply!
} else {
 netdata <- data[data$npi2 %in% unique(data$npi1),] #Restricts data to respondents (assume they each could have named each other but does not require that they actually did)
 netdata <- netdata[order(netdata$npi1),]
}

source('p2est.r')
wantall <- 1
if (wantall==1) {
 dsgns <- 1:20
 nd <- length(dsgns)
 outfit <- matrix(0,nrow=nd,ncol=3*4)
 outtest <- matrix(0,nrow=nd,ncol=3*4)
 for (i in 1:nd) {
  p2mod <- p2est(netdata,type,dsgns[i])
  outfit[i,] <- c(p2mod$Mp[9,c(1:3,8)],p2mod$Mp[10,c(1:3,8)],p2mod$Mp[12,c(1:3,8)]) #7, 8, 9 for unadjusted; 
  outtest[i,] <- c(p2mod$Mptest[9,c(1:3,8)],p2mod$Mptest[10,c(1:3,8)],p2mod$Mptest[12,c(1:3,8)])
 }
 outfit <- data.frame(outfit); outtest <- data.frame(outtest)
 names(outfit) <- paste(rep(c("ddir","dsym","msym"),each=4),rep(c("Est","SE","L2.5","U97.5"),times=3),sep="")
 names(outtest) <- paste(rep(c("ddir","dbas","mbas"),each=4),rep(c("Est","SE","L2.5","U97.5"),times=3),sep="")
} else {
 optdsgn <- 7 #7, 17
 p2mod <- p2est(netdata,type,optdsgn)
 optfit <- p2mod$Mp
 opttest <- p2mod$Mptest
}

if (wantall==1) {
 if (type %in% c(1,2)) {
  #Output from fitted models
  sink(paste(rsource,"AllFitHierP2Models1.txt",sep="")) #Use 1 for type=1 and 2 for type=2
  print("Top 10 Designs: Hierarchical P2 Models of Within-Hospital NPO Networks")
 } else {
  #Output from fitted models
  sink(paste(rsource,"AllFitP2Models3.txt",sep="")) #Add U if unadjusted, H if account for homophily, 1 for type = 1 or 2
  print("Top 10 Designs: P2 Models of NPO Network")
 }
 print(outfit)
 print(outtest)
} else {
 if (type %in% c(1,2)) {
  #Output from fitted models
  sink(paste(rsource,"FittedHierP2Models1.txt",sep="")) #Use 1 for type=1 and 2 for type=2
  print("Hierarchical P2 Models of Within-Hospital NPO Networks")
 } else {
  #Output from fitted models
  sink(paste(rsource,"FittedP2Models3_7.txt",sep="")) #Add U if unadjusted, H if account for homophily, 3 for type = 3 or 4, _ModNum
  print("Optimal P2 Models of NPO Network")
 }
 print(optfit)
 print(opttest)
}
sink()

## Descriptive statistics for paper ##
data <- read.csv(paste(datdir,"SoundPhysAttrData.csv",sep=""),header=TRUE) #Physician attribute data
table(data$nominee)
table(data$sound_provider)
table(data$nominee,data$sound_provider)

data1 <- read.csv(paste(datdir,"ConcordSameHospPristine_19.csv",sep=""),header=TRUE)
table(data1$named,data1$samehosp)
length(unique(data1$npi1))
length(unique(data1$npi2))
length(unique(c(data1$npi1,data1$npi2)))

ind <- (data1$npi2 %in% npiresp)
data1both <- data1[ind,] #Data set used in the dyadic analysis (both physicians responded to the survey)
length(unique(data1both$npi1))
length(unique(data1both$npi2))

data2 <- read.csv(paste(datdir,"ConcordAllSoundBilling_19.csv",sep=""),header=TRUE)
table(data2$named,data2$samehosp)

data3 <- read.csv(paste(datdir,"ConcordAllSound_19.csv",sep=""),header=TRUE) 
data3$samehosp <- ifelse(data3$hospid1==data3$hospid2,1,0)
data3$samehosp[is.na(data3$samehosp)] <- 0 #If can't match hospitals then assume they're different
table(data3$named,data3$samehosp) #The lower counts obtained here compared to 
table(data3$InfOwn,data3$samehosp) 



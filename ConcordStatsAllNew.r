## Developed by James O'Malley June 29, 2024 ##

## Step 1 of data analysis process
## Performs concordance analysis of NPO survey nominations and patient-sharing in Medicare ##
## Also makes plots of networks formed using different projections ##
## Runs on the 40 different ways of forming weighted physician networks ##
## We assume that all Sound physicians at a hospital know each other ##

#Need to install packages using the following due to working insie Granite
#install.packages("Matrix",repos=NULL,contriburl="file:/opt/software/cran/src/contrib")
#install.packages("lme4",repos=NULL,contriburl="file:/opt/software/cran/src/contrib")
#install.packages("nloptr",repos=NULL,contriburl="file:/opt/software/cran/src/contrib")
#install.packages("DescTools",repos=NULL,contriburl="file:/opt/software/cran/src/contrib") #ROC analysis

library(Matrix)
library(lme4)
library(nloptr)
library(foreign)
library(DescTools)
library(igraph)

source <- "/drives/drive1/home/f001580/Desktop/drives/54054-Linux/54054dua/programs/jomalley/"
setwd(source)
datdir="../../idata/jomalley/"
datdirpy="../../idata/jomalley/PHN-subnetwork/2019/batch-1/"
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
#nominee1: whether physician 1 was a nominee (1 = yes, 0 = no). 0 = responding physician takes precedance
#sound1: whether physician employed by NPO
#nominee2: whether 2 was a nominee (1 = yes, 0 = no). 0 = responding physician takes precedence     
#sound2: whether physician 2 was employed by NPO

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
# nSoundPhysk = Number of NPO physicians at hospital
# nRegionk = Numerical code for region of country
# Teachingk = 1 if teaching, 0 otherwise
# Ruralk = 1 if rural, 0 otherwise
# Bed_Sizek = Numerical version of Bedsize variable
# Sampledk = 1 if among 600 sampled physicians, 0 otherwise

#Compute diagnostic properties of test based on shared-patients
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

#Make or load tie-level data set containing CMS and survey information for concordance analysis
makedata <- 0
type <- 1 #Specifications for data set to use for evaluation of diagnostic properties
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
 data <- read.csv(paste(datdir,"ConcordanceTieData19.csv",sep=""),header=TRUE)
}
sensdata <- read.csv(paste(datdir,"ConcordSensitivityAll_19.csv",sep=""),header=TRUE) #Not currently using

#Add the number of patients seen on the same day by the pair of providers as an additional measure
data$nsameday <- data$w6+data$w6r-data$w30
data$wtsameday <- data$w7+data$w7r-data$w28 #Not as meaningful as same-day as mixing concepts

#Decide whether to subset
table(data$nocmsresp,data$nocmsnom)
dropextr <- 1
if (dropextr==1) {
 data <- data[(data$nocmsresp==0 & data$nocmsnom==0),] 
}

#Form the 4 data sets for internal nominations
data$intsurvnom <- ifelse(data$InfOwn+data$RefOwn>0,1,0)
data$extsurvnom <- ifelse(data$InfOth+data$RefOth>0,1,0) #Ideally, there wouldn't be any of these!

#Reduce data set when type=3 to just the sociocentric dyads involving only survey responders
onlyresp <- 0
if (type==3 & onlyresp==1) {
 data <- data[data$npi2 %in% unique(data$npi1),] #A total of 166 physicians
}

#Plot network
covdata <- tapply(data$shospid1,data$id1,'mean') #When analyze pristine can use shosp1
id1 <- as.numeric(names(covdata)); hospid <- as.numeric(covdata)
covdata <- data.frame(id1=id1,hospital=hospid)

edgelist <- data.frame(id1=data$id1,id2=data$id2,w=data$w17,innom=data$intsurvnom,named=data$named,same=data$samehosp) #w7 and w17 are optimal for symmetric and directed decomposition, resp
if (type==1) {
 edgelist <- edgelist[((edgelist$w*(edgelist$same==1)>0 | edgelist$innom==1) & edgelist$id2<=mxid),]
} else {
 edgelist <- edgelist[((edgelist$w>0 | edgelist$named==1) & (edgelist$id2 %in% edgelist$id1)),]
}
nodes <- unique(c(edgelist$id1,edgelist$id2))
covdata <- covdata[covdata$id1 %in% nodes,]
names(covdata) <- c("nodes","hospital")

gnet <- graph_from_data_frame(d=edgelist,vertices=covdata,directed=TRUE)
V(gnet)$color <- V(gnet)$hospital
if (type==1) {
 E(gnet)$width <- edgelist$innom+1
 E(gnet)$color[edgelist$w>0 & edgelist$innom==0] <- 'black' #False positives
 E(gnet)$color[edgelist$w>0 & edgelist$innom==1] <- 'green' #True positives
 E(gnet)$color[edgelist$w==0 & edgelist$innom==1] <- 'red' #False negatives (colorless = true negatives)
} else {
 E(gnet)$width <- edgelist$named+1
 E(gnet)$color[edgelist$w>0 & edgelist$named==0] <- 'black' #False positives
 E(gnet)$color[edgelist$w>0 & edgelist$named==1] <- 'green' #True positives
 E(gnet)$color[edgelist$w==0 & edgelist$named==1] <- 'red' #False negatives (colorless = true negatives)
}
par(mar=c(0,0,0,0))
layout <- layout_with_kk(gnet)
#layout <- layout_with_fr(gnet)
#layout <- layout_nicely(gnet)
plot(gnet,
     layout=layout,
     vertex.color=V(gnet)$color,
     vertex.label.color="black",
     vertex.label.cex=0.75,
     vertex.size=5,
     edge.curved=0.25,
     edge.color=E(gnet)$color, #"grey20"
     edge.width=1, #E(gnet)$width,
     edge.arrow.size=0.3)
if (type==1) {
 dev.copy2pdf(file=paste(outdir,"NPOHospw17kk.pdf",sep=""),width=6,height=6) #If type=1, use kk for w17 to show more interesting plot
} else if (type==3) {
 dev.copy2pdf(file=paste(outdir,"NPOAllRespw7.pdf",sep=""),width=6,height=6) #If type=3 and onlyresp=1
}

#The following ordering of the designs appears to confirm that the ABA designs only favored reverse edges
Interm <- rep(rep(c(0,1),each=5),times=4)
Cont <- 1-Interm #1 = Constraint enforced (no intermediate visits allowed)
Revisit <- rep(rep(c(0,1),each=10),times=2) #1 = Constraint enforced (only count if ABA occurs)
Mult <- rep(seq(1,5),times=8)
#The most simple directed measures are no Cont and no Revisit using mult=1 or 2 or 4; these are w6, w7, and w9

#Patient-sharing variables
wtcols <- paste("w",1:40,sep="")
wtcolsr <- paste(paste("w",1:40,sep=""),"r",sep="")
pats <- c("upat","upatsame","gvisit","upatdir","upatdirr",wtcols,"nsameday")
Directed=c(0,0,0,1,1,rep(rep(c(1,0),each=20),times=1),0) #1 = Directed, 0 = Undirected
npats <- length(pats)
testdata <- data[,c(pats)]
patsr <- c(pats[1:5],wtcolsr,pats[46])
testdatar <- data[,patsr]

diagprop1 <- matrix(0,nrow=npats,ncol=9)
diagprop2 <- matrix(0,nrow=npats,ncol=9)
for (i in 1:npats) {
 diagprop1[i,] <- diagprops(gold=data$named,testdata[,i],N=1e7)
 diagprop2[i,] <- diagprops(gold=data$intsurvnom,testdata[,i],N=1e7)
}
diagprop1 <- data.frame(diagprop1); diagprop2 <- data.frame(diagprop2)
names(diagprop1) <- c("optth","sensopt","specopt","ppvopt","npvopt","distval","Youden","aucprob","aucarea")
names(diagprop2) <- c("optth","sensopt","specopt","ppvopt","npvopt","distval","Youden","aucprob","aucarea")
outallnom <- data.frame(Method=pats,Directed=Directed,diagprop1)
outinnom <- data.frame(Method=pats,Directed=Directed,diagprop2)

#Evaluate diagnostic accuracy of reversed patient-sharing values for the 20 directed designs
nr <- sum(Directed)
diagprop1 <- matrix(0,nrow=npats,ncol=9)
diagprop2 <- matrix(0,nrow=npats,ncol=9)
for (i in 4:(nr+3)) { #Don't redo undirected edge-weights
 diagprop1[i,] <- diagprops(gold=data$named,testdatar[,i],N=1e7)
 diagprop2[i,] <- diagprops(gold=data$intsurvnom,testdatar[,i],N=1e7)
}
diagprop1 <- data.frame(diagprop1); diagprop2 <- data.frame(diagprop2)
names(diagprop1) <- c("optth","sensopt","specopt","ppvopt","npvopt","distval","Youden","aucprob","aucarea")
names(diagprop2) <- c("optth","sensopt","specopt","ppvopt","npvopt","distval","Youden","aucprob","aucarea")
outallnomr <- data.frame(Method=pats,Directed=Directed,diagprop1)
outinnomr <- data.frame(Method=pats,Directed=Directed,diagprop2)

#Analyses to demonstrate utility of retaining directed information
outdiagall <- matrix(0,nrow=npats,ncol=5); outdiagin <- matrix(0,nrow=npats,ncol=5)
outfitall <- matrix(0,nrow=npats,ncol=6); outfitin <- matrix(0,nrow=npats,ncol=6)
outfit1all <- matrix(0,nrow=npats,ncol=6); outfit1in <- matrix(0,nrow=npats,ncol=6)
outfit2all <- matrix(0,nrow=npats,ncol=9); outfit2in <- matrix(0,nrow=npats,ncol=9)
for (i in 6:(nr+3)) {
 xsame <- testdata[,i]
 xrev <- testdatar[,i]
 xsum <- xsame+xrev #New predictor that focuses on symmetric information
 xdirdiff <- (xsame-xrev)/2 #New predictor that focuses on directional information
 fit <- summary(glm(data$named~xsame+xrev+npi,family=binomial(link="logit")))$coef
 fit0 <- glm(data$named~data$upat+npi,family=binomial(link="logit"))
 fit1 <- glm(data$named~xsum+data$upat+npi,family=binomial(link="logit"))
 fit1d <- glm(data$named~xdirdiff+xsum+npi,family=binomial(link="logit"))
 fit2 <- glm(data$named~xdirdiff+xsum+data$upat+npi,family=binomial(link="logit"))
 anovas <- anova(fit0,fit1,fit2)$Deviance[2:3]
 cstats <- c(Cstat(fit0),Cstat(fit1),Cstat(fit2))
 outdiagall[i,] <- c(anovas,cstats)
 outfitall[i,] <- c(fit[2:3,1],fit[2:3,3],fit[2:3,4])
 ofit1 <- summary(fit1d)$coef
 outfit1all[i,] <- c(ofit1[2:3,1],ofit1[2:3,3],ofit1[2:3,4])
 ofit2 <- summary(fit2)$coef
 outfit2all[i,] <- c(ofit2[2:4,1],ofit2[2:4,3],ofit2[2:4,4])
 fit <- summary(glm(data$intsurvnom~xsame+xrev+npi,family=binomial(link="logit")))$coef
 fit0 <- glm(data$intsurvnom~data$upat+npi,family=binomial(link="logit"))
 fit1 <- glm(data$intsurvnom~xsum+data$upat+npi,family=binomial(link="logit"))
 fit1d <- glm(data$intsurvnom~xdirdiff+xsum+npi,family=binomial(link="logit"))
 fit2 <- glm(data$intsurvnom~xdirdiff+xsum+data$upat+npi,family=binomial(link="logit"))
 anovas <- anova(fit0,fit1,fit2)$Deviance[2:3]
 cstats <- c(Cstat(fit0),Cstat(fit1),Cstat(fit2))
 outdiagin[i,] <- c(anovas,cstats)
 outfitin[i,] <- c(fit[2:3,1],fit[2:3,3],fit[2:3,4])
 ofit1 <- summary(fit1d)$coef
 outfit1in[i,] <- c(ofit1[2:3,1],ofit1[2:3,3],ofit1[2:3,4])
 ofit2 <- summary(fit2)$coef
 outfit2in[i,] <- c(ofit2[2:4,1],ofit2[2:4,3],ofit2[2:4,4])
 #Retain optimal values
 if (i==12) { #Optimizes the predictiveness of upat+xsum
  xsumopt <- xsum 
  xdirdiffn <- xdirdiff #Retain this just to confirm that non-predictive
 } else if (i==22) {
  xdirdiffopt <- xdirdiff #Optimizes the predictiveness of upat+xsum+xdirdiffopt
 }
}
outfit1all <- data.frame(outfit1all)
names(outfit1all) <- c("diffest","sumest","diffz","sumz","diffp","sump")
outfit1in <- data.frame(outfit1in)
names(outfit1in) <- c("diffest","sumest","diffz","sumz","diffp","sump")
outfit2all <- data.frame(outfit2all)
names(outfit2all) <- c("diffest","sumest","uest","diffz","sumz","uz","diffp","sump","up")
outfit2in <- data.frame(outfit2in)
names(outfit2in) <- c("diffest","sumest","uest","diffz","sumz","uz","diffp","sump","up")

#Super model: Use upat, xsum when no constraints and multiplicity is 2 or 3, xdirdiff when revisit = 1 and mult is 2
print('Evaluation of Optimal Network Construction: the regression coefficients and preds can be used to build network!')
fitopt <- glm(data$intsurvnom~data$upat+xsumopt+xdirdiffn+xdirdiffopt,family=binomial(link="logit"))
fitopt1 <- glm(data$intsurvnom~xsumopt+xdirdiffn+xdirdiffopt,family=binomial(link="logit"))
fitopt2 <- glm(data$intsurvnom~xsumopt+xdirdiffopt,family=binomial(link="logit"))
print(anova(fitopt2,fitopt1,fitopt)) #Threshold for stat sig at 95% level is difference of 3.8415 or greater, so upat model just sig better!
print(summary(fitopt2))
print(Cstat(fitopt2)) #But fitopt2 might be the best way to build the network? Either this or fitopt!

#Form optimal combinations
corsamerev <- rep(0,npats)
rhoall <- rep(0,npats)
rhoin <- rep(0,npats)
diagprop1 <- matrix(0,nrow=npats,ncol=9)
diagprop2 <- matrix(0,nrow=npats,ncol=9)
for (i in 6:(nr+3)) {
 corsamerev[i] <- cor(testdata[,i],testdatar[,i])
 rhoall[i] <- outfitall[i,1]/(outfitall[i,1]+outfitall[i,2])
 testoptall <- rhoall[i]*testdata[,i]+(1-rhoall[i])*testdatar[,i]
 rhoin[i] <- outfitin[i,1]/(outfitin[i,1]+outfitin[i,2])
 testoptin <- rhoin[i]*testdata[,i]+(1-rhoin[i])*testdatar[,i]
 diagprop1[i,] <- diagprops(gold=data$named,testoptall,N=1e7)
 diagprop2[i,] <- diagprops(gold=data$intsurvnom,testoptin,N=1e7)
}
diagprop1 <- data.frame(diagprop1); diagprop2 <- data.frame(diagprop2)
names(diagprop1) <- c("optth","sensopt","specopt","ppvopt","npvopt","distval","Youden","aucprob","aucarea")
names(diagprop2) <- c("optth","sensopt","specopt","ppvopt","npvopt","distval","Youden","aucprob","aucarea")
outallnomopt <- data.frame(Method=pats,Directed=Directed,diagprop1)
outinnomopt <- data.frame(Method=pats,Directed=Directed,diagprop2)

#Output data
nc <- ncol(outallnom)
allnom <- data.frame(cbind(outallnom,outallnomr[,3:nc],outallnomopt[,3:nc],corsamerev,rhoall))
innom <- data.frame(cbind(outinnom,outinnomr[,3:nc],outinnomopt[,3:nc],corsamerev,rhoin))
named <- c("Method","Directed","optth","sensopt","specopt","ppvopt","npvopt","distval","Youden","aucprob","aucarea")
names(allnom) <- c(paste(named,rep("s",times=nc),sep="_"),paste(named[3:nc],rep("r",times=(nc-2)),sep="_"),paste(named[3:nc],rep("opt",times=(nc-2)),sep="_"),"corsamerev","rho")
names(innom) <- c(paste(named,rep("s",times=nc),sep="_"),paste(named[3:nc],rep("r",times=(nc-2)),sep="_"),paste(named[3:nc],rep("opt",times=(nc-2)),sep="_"),"corsamerev","rho")

#Add new diagnostic output to output matrices
outdiagin <- data.frame(outdiagin)
names(outdiagin) <- c("anova01","anova12","cstat0","cstat1","cstat2")
outdiagall <- data.frame(outdiagall)
names(outdiagall) <- c("anova01","anova12","cstat0","cstat1","cstat2")
innom <- data.frame(innom,outdiagin,outfit1in,outfit2in) #outfit1in
allnom <- data.frame(allnom,outdiagall,outfit1all,outfit2all) #outfit1all

#Add factors to output matrices
Revisit <- c(rep(NA,5),Revisit,NA)
Cont <- c(rep(NA,5),Cont,NA)
Mult <- c(rep(NA,5),Mult,NA)
innom <- data.frame(innom,Revisit=Revisit,Cont=Cont,Mult=Mult)
allnom <- data.frame(allnom,Revisit=Revisit,Cont=Cont,Mult=Mult)

namecomp <- c("Pristine","Rough","Billing","AllTies")
dirname <- paste(outdir,paste(paste("DiagStatsNamedFE",namecomp[type],sep=""),"19.csv",sep="_"),sep="") #r if add results for responders only
write.csv(allnom,file=dirname,row.names=FALSE)
dirname <- paste(outdir,paste(paste("DiagStatsInFE",namecomp[type],sep=""),"19.csv",sep="_"),sep="")
write.csv(innom,file=dirname,row.names=FALSE)

#Summary analysis of overall factors (can also do this outside of Granite) to obtain overall impact of Revisit and Continuity
innoma <- innom[c(6:25),c(3:(ncol(innom)-3))]
innomby <- innom$Mult[c(6:25)]
mnsmult <- matrix(0,nrow=5,ncol=ncol(innoma))
for (i in 1:ncol(innoma)) {
 mnsmult[,i] <- tapply(innoma[,i],innomby,'mean') #Average value for each level of Multiplicity averaging over Revisit and Continuity
}

innomby <- innom$Cont[c(6:25)]
mnscont <- matrix(0,nrow=2,ncol=ncol(innoma))
for (i in 1:ncol(innoma)) {
  mnscont[,i] <- tapply(innoma[,i],innomby,'mean')
}

innomby <- innom$Revisit[c(6:25)]
mnsrevisit <- matrix(0,nrow=2,ncol=ncol(innoma))
for (i in 1:ncol(innoma)) {
  mnsrevisit[,i] <- tapply(innoma[,i],innomby,'mean')
}

innomby <- 10*innom$Revisit[c(6:25)]+innom$Cont[c(6:25)]
mnsrevcont <- matrix(0,nrow=4,ncol=ncol(innoma))
for (i in 1:ncol(innoma)) {
  mnsrevcont[,i] <- tapply(innoma[,i],innomby,'mean') #Average value for each combination of Revisit and Continuity averaging over Multiplicity
}

mns <- data.frame(rbind(mnsmult,mnscont,mnsrevisit,mnsrevcont))
names(mns) <- names(innom)[3:36]
mns$Cont <- c(rep(0,5),0,1,0,0,0,1,0,1)
mns$Revisit <- c(rep(0,5),0,0,0,1,0,0,1,1)
mns$Mult <- c(1:5,rep(0,8))
dirname <- paste(outdir,paste(paste("SumStatsFactorsFE",namecomp[type],sep=""),"19.csv",sep="_"),sep="")
write.csv(mns,file=dirname,row.names=FALSE)


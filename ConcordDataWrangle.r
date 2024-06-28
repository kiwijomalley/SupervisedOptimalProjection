## Developed by James O'Malley June 29, 2024 ##

## Step 2 of data wrangling process
##  Avoids needing to redevelop data from raw each time
##  Compute the number of visits and number unique benes for each Part B physician dyad
##  We can assume that all NPO physicians at a hospital know each other

library(lme4)
library(foreign)

source <- "/drives/drive1/home/f001580/Desktop/drives/54054-Linux/54054dua/programs/jomalley/"
setwd(source)
datdir="../../idata/jomalley/"
datdirpy="../../idata/jomalley/PHN-subnetwork/2019/batch-1/"
outdir="../../idata/jomalley/"

find_mode <- function(x) {
 u <- unique(x)
 tab <- tabulate(match(x,u))
 u[tab==max(tab)]
}


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
#sound1: whether physician employed by NPO
#nominee2: whether 2 was a nominee (1 = yes, 0 = no). 0 = responding physician takes precedence     
#sound2: whether physician 2 was employed by NPO

#The following variables are from the NPO sampling frame and so are missing for non-NPO nominees
#  The suffix "k" takes the values 1 for physician 1 and 2 for physician 2
# hospidk = Hospital ID of NPO physician 
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

#Make or load tie-level data set containing CMS and survey information for concordance analysis
makedata <- 0
if (makedata==1) {
 source('ConcordDataBuildNew.r')
} else {
 data <- read.csv(paste(datdir,"ConcordanceTieData19.csv",sep=""),header=TRUE) 
}

#Checking data
ntie <- nrow(data)
ndyad <- length(unique(data$dyadid))
freqdyadid <- table(data$dyadid)
print(table(freqdyadid))

insoundonly <- (data$named==1 & data$nocmsnom==1)
#The number of ties only in Sound equals the number of dyads only appearing once in data set
print(c(ntie-ndyad,sum(insoundonly)))

#Load nomination data to restrict analytic data sets to samples for which know true nomination status
# Do this just for internal nominations, at least initially, as only have valid nomination status for them
# There are a number of ways of doing this!

## Identify physicians estimated to be at risk of being nominated ##
# First, load data set that contains hospital affiliations and form a variable that indicates whether the 
#  physicians are at the same hospital. These ought to be the physicians with a fair chance of being nominated.

#Load accurate hospital data
hospdataall <- read.csv(paste(datdir,"NPI_sound_billing.csv",sep=""),header=TRUE)
hospdata <- hospdataall[!duplicated(hospdataall[,1]),] #Unique hospital IDs; 5 duplicates but first hospital seems likeliest
hospdata <- hospdata[hospdata$numEncounters>=1,] #Try 1, 10, 100

#Check that these NPIs are in the CMS data
checknpi <- (hospdata$NPI %in% unique(c(data$npi1,data$npi2)))
table(checknpi)
#Compute number of potential edges contained within these hospitals
ind <- (data$npi1 %in% hospdata$NPI & data$npi2 %in% hospdata$NPI)
table(ind)

#Checks that edges are in AllNPIs19
allnpi19 <- read.csv(paste(datdir,"AllNPIs19.csv",sep=""),header=TRUE)
table(hospdata$NPI %in% allnpi19$npi)
table(allnpi19$npi %in% data$npi1)
table(allnpi19$npi %in% data$npi2)

data <- data[,1:136] #In case need to rebuild dataset (136 if just weighted, 216 if binary too)
#Merge Carly's hospital ID with the data and evaluate if overlap with hospid is high
data <- merge(data,hospdata[,1:2],by.x='npi1',by.y='NPI',all.x=TRUE)
names(data)[ncol(data)] <- 'shospid1'
data <- merge(data,hospdata[,1:2],by.x='npi2',by.y='NPI',all.x=TRUE)
names(data)[ncol(data)] <- 'shospid2'
table(data$shospid1==data$hospid1) #The agreement is no higher among survey respondents so don't enlarge shospid1, shospid2
table(data$shospid2==data$hospid2)

#Form unique hospital ID per physicians using the most common one for all observations of patient; might generate earlier?
mhosp1 <- unlist(tapply(data$shospid1,data$npi1,'find_mode'))
mhosp2 <- unlist(tapply(data$shospid2,data$npi2,'find_mode'))
hospdatau <- matrix(as.numeric(cbind(names(mhosp1),as.vector(mhosp1))),ncol=2,byrow=FALSE)
hospdatau <- data.frame(npi=hospdatau[,1],shospidu=hospdatau[,2]) #Hospital attribute dataset
data <- merge(data,hospdatau,by.x="npi1",by.y="npi",all.x=TRUE)
data <- merge(data,hospdatau,by.x="npi2",by.y="npi",all.x=TRUE,suffixes=c(1,2))

#Pristine analysis: Reduce data set to the ties for which both physicians are at the same hospital using strict definition
# of same hospital and one of the physicians is the responder
data$samehosp <- ifelse(data$shospid1==data$shospid2,1,0)
data$samehosp[is.na(data$samehosp)] <- 0 #If can't match hospitals then assume they're different

dataown <- data[data$samehosp==1,]
dataown <- dataown[(dataown$npi1!=dataown$npi2),]
#print(anyDuplicated(dataown))
table(dataown$nocmsresp,dataown$nocmsnom)
write.csv(dataown,paste(outdir,"ConcordSameHospPristine_19.csv",sep=""),row.names=FALSE)

#Less pristine internal hospital nomination data (does not use Carly's assignments)
# This seems to yield data set with non-reciprocated edges in some hospitals
data$samehosp <- ifelse(data$shospid1==data$shospid2 | data$hospid1==data$hospid2 | data$shospid1==data$hospid2 | data$hospid1==data$shospid2,1,0)
data$samehosp[is.na(data$samehosp)] <- 0 #If can't match hospitals then assume they're different

#Re-assign hospital affiliation (can get 3 IDs the same so have the second group of 3 conditions)
expandshosp <- 1 #This seems to help get hospital affiliation to be more accurate
if (expandshosp==1) {
 h11 <- (data$shospid1!=data$shospid2 & data$hospid1==data$hospid2 & data$shospid1!=data$hospid2 & data$hospid1!=data$shospid2)
 hs1 <- (data$shospid1!=data$shospid2 & data$hospid1!=data$hospid2 & data$shospid1==data$hospid2 & data$hospid1!=data$shospid2)
 h1s <- (data$shospid1!=data$shospid2 & data$hospid1!=data$hospid2 & data$shospid1!=data$hospid2 & data$hospid1==data$shospid2)
 h11n <- (data$shospid1!=data$shospid2 & data$hospid1==data$hospid2 & data$shospid1==data$hospid2 & data$hospid1!=data$shospid2)
 hs1n <- (data$shospid1!=data$shospid2 & data$hospid1==data$hospid2 & data$shospid1!=data$hospid2 & data$hospid1==data$shospid2)
 h1sn <- (data$shospid1!=data$shospid2 & data$hospid1!=data$hospid2 & data$shospid1==data$hospid2 & data$hospid1==data$shospid2)
 h11[is.na(h11)] <- 0; hs1[is.na(hs1)] <- 0; h1s[is.na(h1s)] <- 0
 h11n[is.na(h11n)] <- 0; hs1n[is.na(hs1n)] <- 0; h1sn[is.na(h1sn)] <- 0
 data$shospid1[h11==1] <- data$hospid1[h11==1]; data$shospid2[h11==1] <- data$hospid2[h11==1]
 data$shospid2[hs1==1] <- data$hospid2[hs1==1]
 data$shospid1[h1s==1] <- data$hospid1[h1s==1]
 data$shospid2[h11n==1] <- data$hospid1[h11n==1]
 data$shospid1[hs1n==1] <- data$hospid2[hs1n==1]
 data$shospid2[h1sn==1] <- data$shospid1[h1sn==1]
}

dataown <- data[data$samehosp==1,]
dataown <- dataown[(dataown$npi1!=dataown$npi2),]
#print(anyDuplicated(dataown))
table(dataown$nocmsresp,dataown$nocmsnom)
write.csv(dataown,paste(outdir,"ConcordSameHospRough_19.csv",sep=""),row.names=FALSE)

#Nominations allowing all employees at Carly's Sound hospitals that are actively billing to be the eligible set
dataown <- data[data$npi1 %in% hospdata$NPI & data$npi2 %in% hospdata$NPI,]
#dataown <- dataown[(dataown$nocmsresp==0 & dataown$nocmsnom==0),]
#print(anyDuplicated(dataown))
table(dataown$nocmsresp,dataown$nocmsnom)
write.csv(dataown,paste(outdir,"ConcordAllSoundBilling_19.csv",sep=""),row.names=FALSE)

#Note: If keep observations with nocmsresp==0 then the following calculation yields a different data set
# The reason is that nominations not in CMS are replicated across all respondents as opposed to only
# being included on the edge of the physician who nominated them. Retaining them like this means that
# only sensitivity and NPV are impacted by such physicians. If replicate them to include the nulls
# across the 166 respondents (i.e., 165 more null ties) then specificity would also be impacted.
npiresp <- unique(data$npi1[data$npi1 %in% hospdata$NPI])
npihosp <- hospdata$NPI
outdata <- c()
for (i in 1:length(npiresp)) {
 npipair <- cbind(rep(npiresp[i],times=length(npihosp)),npihosp)
 outdata <- rbind(outdata,npipair)
}
outdata <- outdata[!(outdata[,1]==outdata[,2]),]
outdata <- data.frame(npi1=outdata[,1],npi2=outdata[,2])

#All sources sensitivity data 
sensdata <- data[data$named==1,]
write.csv(sensdata,paste(outdir,"ConcordSensitivityAll_19.csv",sep=""),row.names=FALSE)

#All sources data
write.csv(data,paste(outdir,"ConcordAllSound_19.csv",sep=""),row.names=FALSE)
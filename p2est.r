## Developed by James O'Malley June 29, 2024 ##

#p2est: R function for setting up and estimating p2 and hierarchical p2 models, respectively.
# Use type=1 or 2 for the hierarchical p2 model
# Use type=3 or 4 for the p2 model 
# The difference between type=1 vs type=2 and type=3 vs type=4 is immaterial once inside this function
#  The difference reflects whether a single model is being estimated or multiple are (i.e., this procedure
#  is called multiple times by an outside function)

p2est <- function(netdata=netdata,type=type,ndsgn=design) {
 pid <- unique(netdata$npi1)
 np <- length(pid) #Number of respondents
 Y <- matrix(0,np,np)
 Xbase <- matrix(0,np,np)
 X <- matrix(0,np,np)
 Xr <- matrix(0,np,np)
 Xdir <- matrix(0,np,np)
 Xsym <- matrix(0,np,np)
 PhysHosp1 <- tapply(netdata$shospidu1,netdata$npi1,'mean') #Hospital affiliation of responding physician
 wtk <- paste("w",ndsgn,sep="")
 wtrk <- paste(paste("w",ndsgn,sep=""),"r",sep="")
 for (i in 1:nrow(netdata)) { #Uses a calculation that avoids double loop
  pos1 <- 1*(pid==netdata$npi1[i])
  pos2 <- 1*(pid==netdata$npi2[i])*netdata$named[i] #Wipes out lone cell that equals 1 if not named
  Y <- Y + matrix(pos1,ncol=1) %*% matrix(pos2,nrow=1)
  wbase <- 1*(pid==netdata$npi2[i])*netdata$upat[i] #Undirected base weight (upat=w30)
  w <- 1*(pid==netdata$npi2[i])*netdata[i,wtk] #Directed weight (multiply by 0 if want total directed weight)
  wr <- 1*(pid==netdata$npi2[i])*netdata[i,wtrk] #Reverse direction weight
  wdir <- 1*(pid==netdata$npi2[i])*(netdata[i,wtk]-netdata[i,wtrk]) #Directed weight (multiply by 0 if want total directed weight)
  wsym <- 1*(pid==netdata$npi2[i])*(netdata[i,wtk]+netdata[i,wtrk]) #Undirected weight
  Xbase <- Xbase + matrix(pos1,ncol=1) %*% matrix(w,nrow=1)
  X <- X + matrix(pos1,ncol=1) %*% matrix(w,nrow=1)
  Xr <- Xr + matrix(pos1,ncol=1) %*% matrix(wr,nrow=1)
  Xdir <- Xdir + matrix(pos1,ncol=1) %*% matrix(wdir,nrow=1)
  Xsym <- Xsym + matrix(pos1,ncol=1) %*% matrix(wsym,nrow=1) 
 }

 #Be careful to order covariates in the same order as physicians are ordered in the network
 #Form hospital ID factor variables
 u <- unique(netdata$shospidu1)
 nu <- length(u) #Number of hospitals
 PhysHosp <- as.numeric(tapply(netdata$shospidu1,netdata$npi1,'mean'))
 PhysHospAttr <- matrix(0,nrow=np,ncol=nu)
 PhysHospSize <- rep(0,nu)
 for (i in 1:nu) {
  PhysHospAttr[i,] <- 1*(PhysHosp[i]==u)
  PhysHospSize[i] <- sum(netdata$shospidu1==u[i]) #This is a count of dyads (=NumPhys*(NumPhys-1))
 }
 PhysHospAttr1 <- PhysHospAttr[,1] #Make physician 1 the excluded category, can use as sender and receiver covariate?
 
 #Can adjust for individual columns of hospital attribution but not all levels
 #A lot of physicians were the only respondent at a hospital so perhaps better to just adjust for hospital network size
 if (type>2) {
  np2 <- length(unique(netdata$npi2)) #Number of unique physicians who could have been named
  Physbyhosp <- table(netdata$npi1,netdata$shospidu1)/(np-1) #Gets correct number of physicians per hospital
 } else {
  Physbyhosp <- table(netdata$npi1,netdata$shospidu1) 
  Physbyhosp <- 1*(Physbyhosp>0)
 }
 HospSize <- apply(Physbyhosp,2,'sum')
 HospSizeData <- data.frame(ID=names(HospSize),NPhys=as.numeric(HospSize))
 netdata <- merge(netdata,HospSizeData,by.x='shospidu1',by.y='ID')
 PhysHospSize <- as.numeric(tapply(netdata$NPhys,netdata$npi1,'mean')) 
 #PhysHospSize <- (1+sqrt(1+4*PhysHospSize))/2 #As the number of edges=NumPhys(NumPhys-1)

 #Form hospital homophily variable
 if (type %in% c(3,4)) {
  PhysHosp1 <- matrix(netdata$shospidu1,nrow=np,byrow=FALSE)
  PhysHosp2 <- matrix(netdata$shospidu2,nrow=np,byrow=TRUE)
  HomoHosp <- 1*matrix(netdata$shospidu1==netdata$shospidu2,ncol=np,byrow=FALSE)
  AugHomoHosp <- matrix(0,np,np)
  AugHomoHosp[c(2:np),1] <- HomoHosp[,1]
  for (i in 2:(np-1)) {
   AugHomoHosp[c(1:(i-1),(i+1):np),i] <- HomoHosp[,i]
  }
  AugHomoHosp[1:(np-1),np] <- HomoHosp[,np]
 }

 #Form subnetworks by extracting set of physicians with each hospital and then getting all their ties
 Ysub <- c()
 Xbsub <- c()
 Xsub <- c()
 Xrsub <- c()
 Xdirsub <- c()
 Xsymsub <- c()
 PhysHospSizeMat <- c()
 for (i in 1:nu) {
  pidh <- unique(netdata$npi1[netdata$shospid1==u[i]]) #Physicians at a given hospital
  rcphys <- seq(1,np)[pid %in% pidh]
  Ysub[[i]] <- Y[rcphys,rcphys]
  Xbsub[[i]] <- Xbase[rcphys,rcphys]
  Xsub[[i]] <- X[rcphys,rcphys]
  Xrsub[[i]] <- Xr[rcphys,rcphys]
  Xdirsub[[i]] <- Xdir[rcphys,rcphys]
  Xsymsub[[i]] <- Xsym[rcphys,rcphys]
  PhysHospSizeMat[[i]] <- PhysHospSize[rcphys]*diag(length(rcphys)) #All physicians at the same hospital get same value of covariate
 }
 Xbt <- do.call(plyr::rbind.fill.matrix, Xbsub) #NA's denote ties that do not exist as hospital smaller than biggest hospital
 Xt <- do.call(plyr::rbind.fill.matrix, Xsub) #NA's denote ties that do not exist as hospital smaller than biggest hospital
 Xrt <- do.call(plyr::rbind.fill.matrix, Xrsub) #NA's denote ties that do not exist as hospital smaller than biggest hospital
 Xdirt <- do.call(plyr::rbind.fill.matrix, Xdirsub) #NA's denote ties that do not exist as hospital smaller than biggest hospital
 Xsymt <- do.call(plyr::rbind.fill.matrix, Xsymsub) #NA's denote ties that do not exist as hospital smaller than biggest hospital
 HospSizePhys <- do.call(plyr::rbind.fill.matrix, PhysHospSizeMat) #NA

 nobs <- as.vector(table(netdata$shospid1))
 nphys <- (1+sqrt(1+4*nobs))/2 #nobs=nphys*(nphys-1) #Only sure to work if perfect partition of physicians to hospitals
 print(sum(nphys)-ncol(Y)) #Check that total number of physicians sums to dimension of Y

 #Scale shared patient relationships by shared patient counts
 Xbases <- Xbase/100
 Xs <- X/100
 Xrs <- Xr/100
 Xdirs <- Xdir/100
 Xsyms <- Xsym/100
 #PhysHospSizes <- PhysHospSize/100
 Xbts <- Xt/100
 Xts <- Xt/100
 Xrts <- Xrt/100
 Xdirts <- Xdirt/100
 Xsymts <- Xsymt/100
 #HospSizePhys <- HospSizePhys/100

 if (type %in% c(1,2)) {
  #Hierarchical P2 model - with no dependence on shared-patient network. This reveals evidence of heterogeneity in reciprocity
  #Xt and Rt coefficients give association of shared patient-value with likelihood of directed and mutual shared-patient relationships
  M2 <- p2ML(Ysub,sender=~PhysHospSize,receiver=~PhysHospSize,density=~Xdirts+Xsymts,reciprocity=~Xsymts,adapt=10,burnin=1000)
  M2test <- p2ML(Ysub,sender=~PhysHospSize,receiver=~PhysHospSize,density=~Xdirts+Xbts,reciprocity=~Xsymts+Xbts,adapt=10,burnin=1000)
  return(list(Mp=summary(M2),Mptest=summary(M2test)))
 } else {
  #P2 model. Important that exclude reciprocity if want model with no reciprocity covariates!
  Mp2 <- p2(Y,sender=~PhysHospSize,receiver=~PhysHospSize,density=~Xdirs+Xsyms+AugHomoPhys,reciprocity=~Xsyms,adapt=10,burnin=1000,sample=4000)
  Mp2test <- p2(Y,sender=~PhysHospSize,receiver=~PhysHospSize,density=~Xdirs+Xbases+AugHomoPhys,reciprocity=~Xsyms+Xbases,adapt=10,burnin=1000,sample=4000)
  return(list(Mp=Mp2,Mptest=Mp2test)) 
 }
}
  
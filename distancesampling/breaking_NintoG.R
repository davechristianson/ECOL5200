

rm(list=ls())
library(extraDistr)
library(nimble)
library(mcmcplots)
library(sf)
library(jpeg)
groups<-st_cast(st_read("ecol5200distance.gpkg",layer="groups"),to="POINT")

# negative hypergeometric provide the number of samplings, without replacement befor
# r white balls are pulled from n black balls and m white balls
# negative hypergeometric will max out at n + (r-1), because once all the black black balls have been
# pulled, only white balls remain 

# simulated population and data
# n+ m is population size?
set.seed(123456)
N=150
n=N #black balls, maximum group size
m = 10 # white balls 
r=1 
#set nn to number of groups-1
cuts<-(rnhyper(1000000,n,m,r))
max(cuts)
tcuts<-table(cuts)
pcuts<-(tcuts/1000000)
rcuts<-list()
groupsizes<-list()
nn<-10000
for(i in 1:nn){
  rcuts[[i]]<-sort(rnhyper(1,n,m,r))
  groupsizes[[i]]<-diff(c(0,rcuts[[i]]))
}
ng=10
groups<- sample.int(N,ng,replace=F,prob=c(pcuts,rep(0,N-length(pcuts))))
hist(diff(c(0,sort(groups),N)))
grousizes<-diff(c(0,sort(groups),N))



## start simulation informed by data
set.seed(123456)
# sorting simultaneously detected groups by size at each site, to determine the rank of indivdivals
# here we are going to use the observed data to reconstruct the 'cutting' of population size N into G groups
B<-400 # truncation distance
groups<-groups[groups$Species=="pronghorn"&groups$dist_to_transect<B&groups$Total.Group.Size<400,]


# first confirm no strong bias against small groups at large distances
plot(log(groups$Total.Group.Size)~groups$dist_to_transect)
abline(lm(log(groups$Total.Group.Size)~groups$dist_to_transect))

#extract relevant data
obsgs<-groups$Total.Group.Size
obsdx<-groups$dist_to_transect
obssite<-groups$site_id
obsrep<-groups$replicate



# extract the groups at each site in each replicate and sort by size
exclusivegroups<-list()
j=0
for(i in 1:length(groups$Total.Group.Size)){
  for(k in 1:6){
    if(length(groups$Total.Group.Size[groups$site_id==i&groups$replicate==k])==0){
      next
    }
    else {  
    j=j+1
    exclusivegroups[[j]]<-sort(as.numeric(groups$Total.Group.Size[groups$site_id==i&groups$replicate==k]))
    }
  }
}

# now identify the 'ranking' of the first individual in each group 
# ranking is the '1st' individual in a groups position in the larger population
breaks<-lapply(exclusivegroups,cumsum)
hist(unlist(breaks)[unlist(breaks)<400])

# simulate population
set.seed(123456)
nsites<-100
habitat<-rnorm(nsites)
mean.lam<-50
beta0<-log(mean.lam)
beta1<-0.5
lambda<-exp(beta0 + beta1*habitat)
N<-rnbinom(nsites,size=2, mu=lambda)
jpeg("N_hist_sim.jpg",width=6.5,height=6.5,units="in",res=600)
hist(N, main="",mar=c(5,5,1,1))
dev.off()
# determine number of groups
mu<-22 #(mean group size minus 1)
G<-apply(cbind(N,rpois(nsites, N/mu)+1),1,min) #realized group sizes, number of groups increase proportionally to population size, Number of groups must be <= population size


# now determine the groups at each site, using the 'breakpoints' from observed data
cuts<-table(factor(unlist(breaks),levels=1:max(N)))/sum(unlist(breaks))

# create a list of groupsizes at each site
breaklist<-list()
GSlist<-list()
DXlist<-list()
for(i in 1:nsites){
  if(N[i]==0){
    breaklist[[i]]<-0
    GSlist[[i]]<-NULL
    DXlist[[i]]<-NULL
  }
  else{
    breaklist[[i]]<-c(0,sort(sample.int(N[i]-1,size=G[i]-1,replace=F,prob=cuts[1:(N[i]-1)])))
    GSlist[[i]]<-diff(c(breaklist[[i]],N[i]))
    DXlist[[i]]<-runif(sum(unlist(GSlist[[i]])>0),1,B)
  }
}
hist(unlist(GSlist)) # group size distribution
hist(unlist(DXlist)) # distance distrubtion (should be uniform)
table(unlist(GSlist))
unlist(lapply(GSlist,sum))==N
# determine the distances and detections of each group at each site
nD<-10
delta<-B/nD
DClist<-list()
for(i in 1:nsites){
  if(is.null(unlist(DXlist[[i]]))){
    next
  }
  else{
    DClist[[i]]<-as.numeric(cut(unlist(DXlist[[i]]),breaks=seq(0,B,by=delta)))  
  }
}
# DClist[unlist(lapply(DXlist,function (x) !is.null(x)))]<-lapply(DXlist[unlist(lapply(DXlist,function (x) !is.null(x)))],function (x) as.numeric(cut(x,breaks=seq(0,B,by=delta))))
midpt<-(delta/2)+delta*((1:nD)-1)

# number of groups OCCURRING in each distance class
NGDClist<-lapply(DClist,function (x) table(factor(x,levels=1:nD)))

# which groups are detected
pnaught<-0.8 # overall detection probability

sigma0<-B/2
alpha0<-log(sigma0)
alpha1<-0.2 #  groups >5 will show almost no decline in detection probability when = 1

Ylist<-list()
Slist<-list()
Plist<-list()
for(i in 1:nsites){
  if(is.null(GSlist[[i]])){
    next
  }
  else{
    Slist[[i]]<-exp(alpha0 + alpha1*log(GSlist[[i]])) # shape parameter dependent on group size
    Plist[[i]]<-pnaught*exp(-(midpt[DClist[[i]]])^2/(2*(Slist[[i]])^2)) # p depenent on sigma and dclass
    Ylist[[i]]<-rbinom(G[i],1,Plist[[i]])
  }
}
plot(unlist(Plist)~unlist(DXlist))

# depending on the model fit, detections as vectors of unique observations may be desirable
# or as a matrix of detections by site and distance class.
# first create a matrix of detections by distance class (becomes detached from the groupsize)
ymat<-matrix(0,nsites,nD)
for(i in 1:nsites){
    ymat[i,]<-as.numeric(table(factor(Ylist[[i]]*DClist[[i]],levels=1:nD)))
}
# next a vector of observations (retains connection to groupsize)
d<-(unlist(DXlist)*unlist(Ylist))[unlist(Ylist)!=0]
dclass<-(unlist(DClist)*unlist(Ylist))[unlist(Ylist)!=0]
gs<-(unlist(GSlist)*unlist(Ylist))[unlist(Ylist)!=0]
gsite<-rep(c(1:nsites)[rowSums(ymat)>0],rowSums(ymat)[rowSums(ymat)>0])

M <- 400                        # Size of augmented data set is M
nz <- M-length(dclass)              # Number of "pseudo-groups" added
y <- c(rep(1,length(dclass)),rep(0,nz))      # Indicator of capture (== 1 for all obs. groups)
ngroup <- length(dclass)             # Number of observed groups
site <- c(gsite, rep(NA,nz)) # Site they belong to is unknown 
d <- c(d, rep(NA,nz))    # Their distance data are missing ...
groupsize <- c(gs, rep(NA,nz)) # .... as is their size
zst <- y                        # Starting values for data augmentation variable

## Bundle data and produce summary
str(win.data <- list (y=y, B=B, ngroup=ngroup, nsites=nsites, d=d, habitat=habitat,site=site, nz=nz, groupsize=groupsize-1))


## Define model in BUGS langauge as nimbleCode object
# "model1.txt" in book code
Section9p2p2_code <- nimbleCode({ 
  
  ## Prior distributions for model parameters
  alpha0 ~ dunif(-10,10)
  alpha1 ~ dunif(-10,10)
  beta0 ~ dunif(-10,10)
  beta1 ~ dunif(-10,10)
  lambda.group ~ dgamma(0.1, 0.1)
  r ~ dunif(0,1000)
  r.group ~ dunif(0,10)
  probs.lambda.group<-r.group/(r.group+lambda.group)
  ## psi is a derived parameter
  psi <- sum(lambda[1:nsites])/(ngroup+nz)
  
  ## Individual level model: observations and process
  for(i in 1:(ngroup+nz)){
    z[i] ~ dbern(psi)                   # Data augmentation variables
    d[i] ~ dunif(0, B)                  # Distance is uniformly distributed
    groupsize[i] ~ dnegbin(probs.lambda.group,r.group)  # Group size is Poisson
    
    log(sigma[i]) <- alpha0 +  alpha1*log(groupsize[i]+1)
    mu[i] <- z[i]*exp(-d[i]*d[i]/(2*sigma[i]*sigma[i])) #p dep on dist class 
    
    ## here using the half normal detection function (Buckland et al. 2001)
    y[i] ~ dbern(mu[i])
    
    site[i] ~ dcat(site.probs[1:nsites]) # Population distribution among sites
    zg[i]<- z[i]*(groupsize[i] + 1)      # Number of individuals in that group
  }
  
  for(s in 1:nsites){
    ## Model for population size of groups
    N[s] ~ dpois(lambda[s])
    probs.lambda[s] <- r/(r+lambda[s])
    log(lambda[s])<- beta0 + beta1*habitat[s]
    site.probs[s]<- lambda[s]/sum(lambda[1:nsites])
  }
  
  # Derived quantities
  G <- sum(z[1:(ngroup+nz)])        # Total number of groups
  Ntotal <- sum(zg[1:(ngroup+nz)])  # Total population size (all groups combined)
}
)

# define MCMC settings, inits function and parameters to save
ni <- 60000   ;   nb <- 20000   ;   nt <- 100   ;   nc <- 5 ## MCMC settings
inits <- function(){list(alpha0=0,     ## Inits as a function
                         alpha1=0.5, 
                         beta0=0, 
                         beta1=0, 
                         z=zst)}

params <- c("alpha0", "alpha1", "beta0", "beta1", "psi", "Ntotal", "G", 
            "lambda.group","r","r.group")                     ## define params for monitors

# Call NIMBLE, plot posterior distributions
out1 <- nimbleMCMC(code = Section9p2p2_code, 
                   constants = win.data,
                   inits = inits,
                   monitors = params,
                   niter = ni,
                   nburnin = nb,
                   samplesAsCodaMCMC = TRUE)
# this is the classic DA approach where lambda estimates group abundance

library(mcmcplots)
colnames(out1)
mcmcplot(out1)
print(out1,dig=3)
samplesSummary(out1,round=3)




plot(log(unlist(GSlist)),unlist(Plist))

####################################
#3333333333333333333333333333333333#
####################################
#Traditional mutinomial formulation of HDS, N = G*gamma

gobs <- apply(ymat, 1, sum)  # Total detections per site and occasion
str( data3 <- list(y=ymat, nsites=nsites, nD=nD, midpt=midpt, delta=delta, habitat=habitat, B=B, gobs = gobs,ngroup=ngroup,site=gsite,groupsize=gs-1) )


## Define model in BUGS

Section9p5p4p1_code <- nimbleCode({
  # Prior distributions
  beta0 ~ dnorm(0, 0.01)  # Intercept for log(lambda)
  mean.lam <- exp(beta0)
  beta1 ~ dnorm(0, 0.01)  # Coefficient of lambda on habitat
  #phi ~ dunif(0,1)        # Probability of availability
  sigma ~ dunif(0.01,800)   # Distance function parameter
  gamma0 ~ dnorm(0,0.01)
  r.group ~ dunif(0,10)
  # Detection probs for each distance interval and related things
  for(b in 1:nD){
    log(g[b]) <- -midpt[b]*midpt[b]/(2*sigma*sigma) # half-normal 
    f[b] <- 1/nD    # radial density function
    cellprobs[b] <- g[b]*f[b]
    cellprobs.cond[b] <- cellprobs[b]/sum(cellprobs[1:nD])
  }
  cellprobs[nD+1]<- 1-sum(cellprobs[1:nD])
  
  for (s in 1:nsites) {
    y[s,1:nD] ~ dmulti(cellprobs.cond[1:nD], gobs[s])
    gobs[s] ~ dbin(pdet[s], G[s])
    pdet[s] <- sum(cellprobs[1:nD])   # Distance class probabilities
    #for (k in 1:K) {
      #
      #pmarg[s,k] <- pdet[s,k]*phi         # Marginal probability
      
      # Model part 4: distance class frequencies
      #
      # Model part 3: total number of detections:
      #
      # nobs[s,k] ~ dbin(pdet[s,k], Navail[s,k]) # Alternative formulation
      # Model part 2: Availability. Not used in this model but simulated.
      #Navail[s,k] ~ dbin(phi, M[s]) 
    #}  # end k loop
    # Model part 1: Abundance model
    G[s] ~ dpois(lambda[s])    
    log(lambda[s]) <- beta0 + beta1*habitat[s]
  }  # End s loop
  
  
  for(i in 1:ngroup){
    groupsize[i] ~ dnegbin(probs[site[i]],r.group)
  }
  for(s in 1:nsites){
    log(gamma[s]) <- gamma0 
    probs[s]<-r.group/(r.group+gamma[s])
  }
  
  # Derived quantities
  Gtot <- sum(G[1:nsites])
  Ntot<- sum(G[1:nsites]*(1+gamma[1:nsites]))

} # End model
)


# Assemble the initial values and parameters to save for JAGS

Gst <- gobs  +2
inits3 <- function(){
  list(G=Gst, sigma = 327, beta0=log(2), beta1=0.25)
}
params3 <- c("sigma", "beta0", "mean.lam", "gamma0","beta1", "Gtot", "Ntot","cellprobs.cond")

# MCMC settings
ni <- 60000   ;   nb <- 10000   ;   nt <- 100   ;   nc <- 5

# Run nimble:
# Some additional init may be needed for nimble.
out3 <- nimbleMCMC(
  code = Section9p5p4p1_code,
  constants = data3,
  inits = inits3,
  monitors = params3,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)


samplesSummary(out3,round=3)
mcmcplot(out3)




####################################
#4444444444444444444444444444444444#
####################################
#using pdet to estimate N from total nobs at a site

nobs<-as.numeric(tapply(gs,factor(gsite,levels=1:nsites),sum,na.rm=T))
nobs[is.na(nobs)]<-0

str( data4 <- list(y=ymat, nsites=nsites, nD=nD, midpt=midpt, delta=delta, habitat=habitat, B=B, gobs = gobs,ngroup=ngroup,site=gsite,groupsize=gs-1), nobs=nobs)


## Define model in BUGS

PdetToNobsCode <- nimbleCode({
  # Prior distributions
  beta0 ~ dnorm(0, 0.01)  # Intercept for log(lambda)
  mean.lam <- exp(beta0)
  beta1 ~ dnorm(0, 0.01)  # Coefficient of lambda on habitat
  #phi ~ dunif(0,1)        # Probability of availability
  sigma ~ dunif(0.01,800)   # Distance function parameter
  gamma0 ~ dnorm(0,0.01)
  gamma1 ~ dnorm(0,0.01)
  r ~ dunif(0,10)
  r.group ~ dunif(0,10)
  # Detection probs for each distance interval and related things
  for(b in 1:nD){
    log(g[b]) <- -midpt[b]*midpt[b]/(2*sigma*sigma) # half-normal 
    f[b] <- 1/nD    # radial density function
    cellprobs[b] <- g[b]*f[b]
    cellprobs.cond[b] <- cellprobs[b]/sum(cellprobs[1:nD])
  }
  cellprobs[nD+1]<- 1-sum(cellprobs[1:nD])
  
  for (s in 1:nsites) {
    y[s,1:nD] ~ dmulti(cellprobs.cond[1:nD], gobs[s])
    gobs[s] ~ dbin(pdet[s], G[s])
    G[s] ~ dnegbin(probs.lambda[s],r)
    
    
    N[s] ~ dpois(lambda[s])
    probs.lambda[s]<-  r/(r+(lambda[s]/gamma[s]))
    pdet[s] <- sum(cellprobs[1:nD])   # Distance class probabilities
    #for (k in 1:K) {
    #
    #pmarg[s,k] <- pdet[s,k]*phi         # Marginal probability
    
    # Model part 4: distance class frequencies
    #
    # Model part 3: total number of detections:
    #
    # nobs[s,k] ~ dbin(pdet[s,k], Navail[s,k]) # Alternative formulation
    # Model part 2: Availability. Not used in this model but simulated.
    #Navail[s,k] ~ dbin(phi, M[s]) 
    #}  # end k loop
    # Model part 1: Abundance model
    log(lambda[s]) <- beta0 + beta1*habitat[s]
  }  # End s loop
  
  
  for(i in 1:ngroup){
    groupsize[i] ~ dnegbin(probs[site[i]],r.group)
  }
  for(s in 1:nsites){
    log(gamma[s]) <- gamma0 + gamma1*habitat[s]
    probs[s]<-r.group/(r.group+gamma[s])
  }
  
  # Derived quantities
  mean.gamma<-mean(gamma[1:nsites])
  Ntot<- sum(N[1:nsites])
  Gtot<- sum(G[1:nsites])
} # End model
)


# Assemble the initial values and parameters to save for JAGS

Nst <- nobs  + 10
Gst <- gobs +2
inits4 <- function(){
  list(N=Nst, G=Gst, sigma = 327, gamma0 = log(10), beta0=log(2), beta1=0.25,r.group=0.5)
}
params4 <- c("sigma", "beta0", "mean.lam", "beta1", "Ntot","Gtot","gamma0","r","r.group","mean.gamma")

# MCMC settings
ni <- 60000   ;   nb <- 10000   ;   nt <- 100   ;   nc <- 5

# Run nimble:
# Some additional init may be needed for nimble.
out4 <- nimbleMCMC(
  code = PdetToNobsCode,
  constants = data4,
  inits = inits4,
  monitors = params4,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)


samplesSummary(out4,round=3)
mcmcplot(out4)



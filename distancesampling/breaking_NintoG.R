

rm(list=ls())
library(extraDistr)
library(nimble)
library(mcmcplots)
library(sf)
library(jpeg)
library(BiasedUrn)

# negative hypergeometric provide the number of samplings, without replacement before
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

# Wallis Noncentral (Multivariate) Hypergeometric Distribution
#rMWNCHypergeo(nran=1,m=rep(1,482),n=5,odds=c)


groups<-st_cast(st_read("ecol5200distance.gpkg",layer="groups"),to="POINT")
B<-300 # truncation distance
groups<-groups[groups$Species=="pronghorn"&groups$dist_to_transect<B,]



## start simulation informed by data
set.seed(123456)
# sorting simultaneously detected groups by size at each site, to determine the rank of indivdivals
# here we are going to use the observed data to reconstruct the 'cutting' of population size N into G groups


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
hist(unlist(breaks))
dev.off()
jpeg("pronghorn_breakpoints.jpeg",height=6.5, width=6.5, units="in",res=600)
hist(unlist(breaks),main="",xlab="Group Breakpoints (Cumulative Sums of Indviduals)",breaks=seq(0,500,1))
dev.off()
# same graph as above but cumulative sum of individuals as frequency
hist(rep(unlist(breaks),unlist(breaks)),breaks=seq(0,500,1))

# simulated group ranks (cumulative sums) using pronghorn group ranks
# Wallis Noncentral (Multivariate) Hypergeometric Distribution only allows 32 bins by default
rMWNCHypergeo(nran=1,m=rep(1,488),n=5,odds=as.numeric(table(factor(unlist(breaks),levels=1:488))),precision=0.1)

# compare the function sample.int(replace = F) with Wallenius' Multivariate Noncentral Hypergometric
# they are identical
probs<-c(100,50,25,rep(1, 20),100,200)
simsamp<-list()
for(i in 1:10000){
  simsamp[[i]]<-as.numeric(table(factor(sample.int(n=25,size=5,replace=F,prob=probs),levels=1:25)))
}
simsampmat<-do.call("rbind",simsamp)

par(mfrow=c(3,1))
plot(rowSums(rMWNCHypergeo(nran=10000,m=rep(1,25),n=5,odds=probs)))
plot(colSums(simsampmat))
plot(rowSums(rMFNCHypergeo(nran=10000,m=rep(1,25),n=5,odds=probs)))

# now construct a new population of pronghorn using the observered cumultive sums
simsamp<-list()
for(i in 1:10000){
  simsamp[[i]]<-as.numeric(table(factor(sample.int(n=488,size=5,replace=F,prob=as.numeric(table(factor(unlist(breaks),levels=1:488)))),levels=1:488)))
}

par(mfrow=c(2,1))
plot(colSums(do.call("rbind",simsamp)))

#Compare populations of 1 10 and 100 individuals
#Group sizes of mean of 1, 2 and 10 individuals


# now simulate population
set.seed(1234567)
nsites<-1000
habitat<-rnorm(nsites)
mean.lam<-1
beta0<-log(mean.lam)
beta1<-0.33
lambda<-exp(beta0 + beta1*habitat)
hist(lambda)
N<-rnbinom(nsites,size=5, mu=lambda)
hist(N)
dev.off()
jpeg("N_hist_sim.jpg",width=6.5,height=6.5,units="in",res=600)
hist(N, main="",mar=c(5,5,1,1))
dev.off()

# determine number of groups
# there are several approaches that could be used.  
# one appproach is to start with the idealized group size and then extract that from N
# here, gammatr and gammatr0 refer to  group size - 1 (0 truncated) 
gamma0_tr<-1.5


# sim1: Group size held constant: G increase with N.  simplest formulation, mean group size is completely independent of N.  G must increase as N increases
mu_GS_tr<-exp(gamma0_tr)
G<-apply(cbind(N,rpois(nsites, N/(mu_GS_tr+1))+1),1,min) #realized group sizes, number of groups increase proportionally to population size, Number of groups must be <= population size, first +1 is to account for -1 in group size model, second +1 is to truncate the poisson at 1.
real_mu_site_GS<-N[N!=0]/G[N!=0]
real_mu_GS<-mean(real_mu_site_GS) # realized mean group size at every site, not = gamma because of the schoastic random poisson draw of group number and the marginal constrations imposed by N
plot(real_mu_site_GS~N[N!=0])
plot(G[N!=0]~N[N!=0])
plot(mu_GS_tr~N[N!=0])

summary(glm(real_mu_site_GS~N[N!=0],family=quasipoisson))



# Sim 2: Number of groups held constant. group size increases proportionately by (1/G or 1/5 here) with increasing density (group abundance is random poisson draw at constant rate), 
mu_G<-5
G<-apply(cbind(N,rpois(nsites, mu_G)+1),1,min) #realized group sizes, number of groups increase proportionally to population size, Number of groups must be <= population size
real_mu_site_GS<-N[N!=0]/G[N!=0]
real_mu_GS<-mean(real_mu_site_GS)
plot(real_mu_site_GS~N[N!=0])
plot(G[N!=0]~N[N!=0])
coef(glm(real_mu_site_GS~N[N!=0],family="quasipoisson"))
coef(lm(real_mu_site_GS~N[N!=0]))


# sim 3: inverse response: mean group size covaries inversely with N via the same site covariate that positively influences N
# this is indirect negative covariation in density and group sizes. number of groups must increase as N increases
gamma1_tr<--0.5
mu_GS_tr<-exp(gamma0 + gamma1_tr*habitat)
G<-apply(cbind(N,rpois(nsites, N/(mu_GS_tr+1))+1),1,min) #realized group sizes, number of groups increase proportionally to population size, Number of groups must be <= population size
real_mu_site_GS<-N[N!=0]/G[N!=0]
real_mu_GS<-mean(real_mu_site_GS)
plot(real_mu_site_GS~N[N!=0])
plot(G[N!=0]~N[N!=0])
plot(mu_GS_tr~N[N!=0])

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
  else if(N[i]==1){ # if only one individual in the population, only 1 group is possible
    breaklist[[i]]<-0
    GSlist[[i]]<-1
    DXlist[[i]]<-runif(1,1,B)
  }
  else{
    breaklist[[i]]<-c(0,sort(sample.int(N[i]-1,size=G[i]-1,replace=F,prob=cuts[1:(N[i]-1)])))
    GSlist[[i]]<-diff(c(breaklist[[i]],N[i]))
    DXlist[[i]]<-runif(sum(unlist(GSlist[[i]])>0),1,B)
  }
}

dev.off()
jpeg("sim_vs_real_gs_dist.jpg",width=6.5,height=6.5,res=600,units="in")
par(mfrow=c(2,1))
hist(unlist(GSlist),breaks=seq(0,500,by=5)) # group size distribution in simulation
hist(groups$Total.Group.Size,breaks=seq(0,500,by=5))
dev.off()

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
dev.off()
jpeg("detp_vs_dist_and_gs_sim.jpeg",width=6.5,height=6.5,units="in",res=600)
plot(unlist(Plist)~unlist(DXlist),cex=log(unlist(GSlist)+0.5),xlab="distance (m)",ylab="detection probability")
dev.off()

min(unlist(GSlist)) # no zeros, correct?
mean(unlist(GSlist))

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

# other pieces of data needed in some models
# number of groups observed at each site
gobs <- apply(ymat, 1, sum)  # Total detections per site and occasion

# number of individuals observed at each site
nobs<-as.numeric(tapply(gs,factor(gsite,levels=1:nsites),sum,na.rm=T))
nobs[is.na(nobs)]<-0

###########################################
## MODEL 1 conventional DA approach #######
###########################################


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
  gamma0 ~ dgamma(0.1, 0.1) 
  gamma1 ~ dunif(-10,10)
  #r ~ dunif(0,50)
  #r.gamma ~ dunif(0,20)
  
  ## psi is a derived parameter
  psi <- sum(lambda[1:nsites])/(ngroup+nz)
  
  ## Individual level model: observations and process
  for(i in 1:(ngroup+nz)){
    z[i] ~ dbern(psi)                   # Data augmentation variables
    d[i] ~ dunif(0, B)                  # Distance is uniformly distributed
    #groupsize[i] ~ dnegbin(probs.gamma[site[i]],r.gamma)  # Group size is Poisson
     
    
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
    #probs.lambda[s] <- r/(r+lambda[s])
    log(lambda[s]) <- beta0 + beta1*habitat[s]
    site.probs[s] <- lambda[s]/sum(lambda[1:nsites])
    log(gamma[s]) <- gamma0 + gamma1*habitat[s]
    #probs.gamma[s]<-r.gamma/(r.gamma+gamma[s])
  }
  for(i in 1:(ngroup+nz)){
    groupsize[i] ~ dpois(gamma[site[i]])
  }
  
  # Derived quantities
  G <- sum(z[1:(ngroup+nz)])        # Total number of groups
  Ntotal <- sum(zg[1:(ngroup+nz)])  # Total population size (all groups combined)
  mean.sig<-mean(sigma[1:ngroup])
}
)

# define MCMC settings, inits function and parameters to save
ni <- 60000   ;   nb <- 20000   ;   nt <- 100   ;   nc <- 5 ## MCMC settings
inits <- function(){list(alpha0=0,     ## Inits as a function
                         alpha1=0.5, 
                         beta0=0, 
                         beta1=0, 
                         z=zst)}

params <- c("mean.sig","alpha0", "alpha1", "beta0", "beta1", "psi", "Ntotal", "G", 
            "gamma0","gamma1")                     ## define params for monitors

# Call NIMBLE, plot posterior distributions
out1 <- nimbleMCMC(code = Section9p2p2_code, 
                   constants = win.data,
                   inits = inits,
                   monitors = params,
                   niter = ni,
                   nburnin = nb,
                   samplesAsCodaMCMC = TRUE)
# this is the classic DA approach where lambda estimates group abundance


mcmcplot(out1)
print(out1,dig=3)
samplesSummary(out1,round=3)




plot(log(unlist(GSlist)),unlist(Plist))

####################################
#3333333333333333333333333333333333#
####################################
#Traditional mutinomial formulation of HDS, N = G*gamma

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

Gst <- gobs+2
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


str( data4 <- list(y=ymat, nsites=nsites, nD=nD, midpt=midpt, delta=delta, habitat=habitat, B=B, gobs = gobs,site=gsite, nobs=nobs))


## Define model in BUGS

PdetToNobsCode <- nimbleCode({
  # Prior distributions
  beta0 ~ dnorm(0, 0.01)  # Intercept for log(lambda)
  mean.lam <- exp(beta0)
  beta1 ~ dnorm(0, 0.01)  # Coefficient of lambda on habitat
  #phi ~ dunif(0,1)        # Probability of availability
  sigma ~ dunif(0.01,800)   # Distance function parameter
  #gamma0 ~ dnorm(0,0.01)
  #gamma1 ~ dnorm(0,0.01)
  r ~ dunif(0.9,7)
  #r.group ~ dunif(0,10)
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
    #gobs[s] ~ dbin(pdet[s], G[s])
    #G[s] ~ dnegbin(probs.lambda[s],r)
    
    
    N[s] ~ dnegbin(probs.lambda[s],r)
    probs.lambda[s]<-  r/(r+(lambda[s]))
    pdet[s] <- sum(cellprobs[1:nD])   # Distance class probabilities
    #for (k in 1:K) {
    #
    #pmarg[s,k] <- pdet[s,k]*phi         # Marginal probability
    
    # Model part 4: distance class frequencies
    #
    # Model part 3: total number of detections:
    #
    
    nobs[s] ~ dbin(pdet[s], N[s]) # Alternative formulation
    # Model part 2: Availability. Not used in this model but simulated.
    #Navail[s,k] ~ dbin(phi, M[s]) 
    #}  # end k loop
    # Model part 1: Abundance model
    log(lambda[s]) <- beta0 + beta1*habitat[s]
  }  # End s loop
  
  
  #for(i in 1:ngroup){
  #  groupsize[i] ~ dnegbin(probs[site[i]],r.group)
  #}
  #for(s in 1:nsites){
  #  log(gamma[s]) <- gamma0 + gamma1*habitat[s]
  #  probs[s]<-r.group/(r.group+gamma[s])
  #}
  
  # Derived quantities
  #mean.gamma<-mean(gamma[1:nsites])
  #Ntot<- sum(N[1:nsites])
} # End model
)


# Assemble the initial values and parameters to save for JAGS

Nst <- nobs  + 10
#Gst <- gobs +2
inits4 <-  list(N=Nst, sigma = 250, beta0=log(50), beta1=0.25, r=3)

params4 <- c("sigma", "pdet","mean.lam","beta1")

# MCMC settings
ni <- 60000   ;   nb <- 10000   ;   nt <- 100   ;   nc <- 5

# Run nimble:
# Some additional init may be needed for nimble.

nmodel4<-nimbleModel(
  code = PdetToNobsCode,
  constants = data4,
  inits = inits4)


comp4<-configureMCMC(nmodel4)

comp4$removeSamplers(c("r","beta0","beta1"))

comp4$addSampler(target=c("r","beta0","beta1"),type="AF_slice")

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






####################################
#55555555555555555555555555555555555#
####################################
# mutinomial formulation of HDS, N = G*gamma but with nobs

str( data5 <- list(y=ymat, nsites=nsites, nD=nD, midpt=midpt, delta=delta, habitat=habitat, B=B, gobs = gobs,ngroup=ngroup,site=gsite,groupsize=gs-1,nobs=nobs) )


## Define model in BUGS

PdetNobsGammaCode <- nimbleCode({
  # Prior distributions
  beta0 ~ dnorm(0, 0.01)  # Intercept for log(lambda)
  mean.lam <- exp(beta0)
  beta1 ~ dnorm(0, 0.01)  # Coefficient of lambda on habitat
  #phi ~ dunif(0,1)        # Probability of availability
  sigma ~ dunif(0.01,800)   # Distance function parameter
  gamma0 ~ dnorm(0,0.01)
  r.group ~ dunif(0.5,5)

  # Detection probs for each distance interval and related things
  for(b in 1:nD){
    log(g[b]) <- -midpt[b]*midpt[b]/(2*sigma*sigma) # half-normal 
    f[b] <- 1/nD    # radial density function
    cellprobs[b] <- g[b]*f[b]
    cellprobs.cond[b] <- cellprobs[b]/sum(cellprobs[1:nD])
  }
  cellprobs[nD+1]<- 1-sum(cellprobs[1:nD])
  pdet <- sum(cellprobs[1:nD])
  
  for (s in 1:nsites) {
    y[s,1:nD] ~ dmulti(cellprobs.cond[1:nD], gobs[s]) # detection probability of individuals resolved
    #pdet[s] <- sum(cellprobs[1:nD])   # Distance class probabilities
    #for (k in 1:K) {
    #
    #pmarg[s,k] <- pdet[s,k]*phi         # Marginal probability
    
    # Model part 4: distance class frequencies
    #
    # Model part 3: total number of detections:
    #
    nobs[s] ~ dbin(pdet, G[s]*gamma[s]) # Alternative formulation
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
inits5 <- function(){
  list(G=Gst, sigma = 211, beta0=log(3), beta1=0.25,r.group=1)
}
params5 <- c("sigma", "beta0", "mean.lam", "gamma0","beta1","Ntot","Gtot")
# MCMC settings
ni <- 60000   ;   nb <- 10000   ;   nt <- 100   ;   nc <- 5

# Run nimble:
# Some additional init may be needed for nimble.
out5 <- nimbleMCMC(
  code = PdetNobsGammaCode,
  constants = data5,
  inits = inits5,
  monitors = params5,
  nburnin = nb,
  niter = ni,
  samplesAsCodaMCMC = TRUE
)


samplesSummary(out5,round=3)
mcmcplot(out3)



rm(list=ls())
library(sf)
library(rjags)
sites<-st_cast(st_read("ecol5200distance.gpkg",layer="sites"),to="LINESTRING")
groups<-st_cast(st_read("ecol5200distance.gpkg",layer="groups"),to="POINT")

# distance sampling parameters
nD<-6
B<-0.600
delta<-0.100
midpt<-0.5*(delta)+seq(0,(B-delta),by=delta)


# site replicate info

repcols<-paste0("rep.",seq(1,max(sites$totalreps)))
rep_effort<-1*st_drop_geometry(sites[,repcols])
K<-max(sites$totalreps)
nreps<-sites$totalreps

# spatial covariates (by site)
nsites<-as.numeric(max(sites$site_id))
roads<-sites$road_kmspersqkm
land<-sites$prop_publand

# organize detections in 3-D array 
g3d<-array(data=0, dim=c(nsites,nD,K))


groups$dclass<-as.numeric(cut(groups$dist_to_transect,seq(0,B*1000,delta*1000)))         
#exclude observations that can't be binned into a dclass (farther out than B)
groups<-st_drop_geometry(groups[!is.na(groups$dclass)&groups$Species=="pronghorn",])

# get rid of the one outlier...
groups<-groups[groups$Total.Group.Size!=426,]


testarr<-table(factor(groups$site_id,levels=levels(as.factor(sites$site_id))),
                groups$dclass,groups$replicate,useNA="ifany")
# names are now out of sequential order
for (i in 1:K){
  g3d[,,i]<-testarr[order(as.numeric(row.names(testarr[,,i]))),,i]
}



           
# Bundle and summarize the data set
gobs <- apply(g3d, c(1,3), sum) # Total detections per site and occasion
nobs.max<-tapply(groups$Total.Group.Size,factor(groups$site_id,levels=1:nsites),sum,na.rm=T)
nobs.max[is.na(nobs.max)]<-0
nobs <- tapply(groups$Total.Group.Size,list(factor(groups$site_id,levels=1:nsites),groups$replicate),sum,na.rm=T)
nobs[is.na(nobs)]<-0

  
str( data <- list(g3d=g3d, nsites=nsites, K=K, nD=nD, midpt=midpt, delta=delta,
                  roads=roads,land=land, B=B, gobs = gobs, nobs=nobs,
                  nreps=nreps,
                  area=sites$site_length_km*2*B,
                  groupsize=groups$Total.Group.Size,
                  site=groups$site_id,
                  ngroup=length(groups$site_id)))


#Define model. From 9.5.4.1, page 501.
cat("
model {
# Prior distributions
beta0 ~ dnorm(0, 0.01) # Intercept for log(lambda.group)
mean.lam <- exp(beta0)
beta1 ~ dnorm(0, 0.01) # Coefficient of roads on lambda.group
beta2 ~ dnorm(0,0.01) # Coefficient of public land
phi0 ~ dunif(0,1) # Probability of availability
alpha0 ~ dnorm(0,0.01)
alpha1 ~ dnorm(0,0.01)
alpha2 ~ dnorm(0,0.01)
#sigma ~ dunif(0.01,5) # Distance function parameter
gamma0 ~ dnorm(0,0.01)
gamma1 ~ dnorm(0,0.01)
gamma2 ~ dnorm(0,0.01)
sd~dunif(0,10)
tau<-1/pow(sd,2)
r ~ dunif(0,10)
r.lam ~ dunif(0,10)


  for (s in 1:nsites) {
    log(sigma[s])<- alpha0+alpha1*roads[s]+alpha2*land[s]
  
    for(b in 1:nD){
      log(g[b,s]) <- -midpt[b]*midpt[b]/(2*sigma[s]*sigma[s]) # half-normal
      f[b,s] <- delta/B # density function
      cellprobs[b,s] <- g[b,s]*f[b,s]
      cellprobs.cond[b,s] <- cellprobs[b,s]/sum(cellprobs[1:nD,s])
    }
    cellprobs[nD+1,s] <- 1-sum(cellprobs[1:nD,s])

    #logit(phi[s]) <- alpha0 + alpha1*roads[s] + alpha2*land[s]
    
    for (k in 1:nreps[s]) {
      pdet[s,k] <- sum(cellprobs[1:nD,s]) # Distance class probabilities
      g3d[s,1:nD,k] ~ dmulti(cellprobs.cond[1:nD,s], gobs[s,k])
      
      # Model  1: using standard model of group abundance multiplied by mean group size
      gobs[s,k] ~ dbin(pdet[s,k], G[s])
      
      #Model 2:  lambda=individuals, Group abundance is Individual abundance divided by mean group size
      #gobs[s,k] ~ dbin(pdet[s,k], G[s])
      
      # Model 3: lambda=indivduals, indivdiual det prob inherited from group det prob
      #nobs[s,k] ~ dbin(pdet[s,k], N[s])
      
    } 
    eps.site[s] ~ dnorm(0,tau)
    lambda[s] <- exp(llambda.lim[s])
    llambda.lim[s] <- min(250, max(-250, llambda[s])) # 'Stabilize' log
    llambda[s] <- log(area[s]) + beta0 + beta1*roads[s] +beta2*(exp(land[s])/(1+exp(land[s]))) + eps.site[s]
    
    #Model 1: lambda=groups, conventional model of group abundance multiplied by mean groups size
    G[s] ~ dpois(lambda[s]) # super-population of group sizes, is superpopulation/mean groupsize

    #Model 2:  lambda=individuals, Group abundance is Individual abundance divided by mean group size
    #N[s] ~ dnegbin(probs.lam[s], r.lam)
    #probs.lam[s] <- r.lam/(r.lam+lambda[s])
    #G[s] ~ dpois(N[s]/gamma[s]) # super-population of group sizes, is superpopulation/mean groupsize

    #Model 3:  lambda=indivduals, indivdiual det prob inherited from group det prob
    #N[s] ~ dnegbin(probs.lam[s], r.lam)
    #probs.lam[s] <- r.lam/(r.lam+lambda[s])

    #for (k in 1:K){
      #Gavail[s,k] ~ dbin(phi[s], G[s])
      #Gavail[s,k] ~ dbin(phi0, G[s])
    #}

  } 

  for(i in 1:ngroup){
    groupsize[i] ~ dnegbin(probs[site[i]],r)
  }

  for(s in 1:nsites){
    log(gamma[s]) <- gamma0 + gamma1*roads[s] + gamma2*(exp(land[s])/(1+exp(land[s])))
    #log(gamma[s]) <- gamma0 
    probs[s]<-r/(r+gamma[s])
  }
  
  #Model 1
  Ntot<- sum(lambda[]*gamma[])
  
  #Model 2
  #Ntot<-sum(N[])
  
  #Model 3
  #Ntot<-sum(N[])

} 
",file="model.txt")
# Assemble the initial values and parameters to save for JAGS
Gavail.st <- apply(g3d, c(1,3),sum)+1
#Nst <- apply(Gavail.st, c( 1), max)+ 2
Nst <- nobs.max+10
inits <- function(){
  list(N=Nst, r=1.1,r.lam=5, phi0=0.2, beta0=1, beta1=0, Gavail=Gavail.st)
}
params <- c("r", "r.lam","mean.pdet","sigma","phi0","alpha1","alpha2","mean.lam", "beta1","beta2","mean.gamma","gamma0","gamma1","gamma2","Ntot","lambda","gamma","Gavail","G","N")
# MCMC settings
ni <- 60000 ; nb <- 10000 ; nt <- 100 ; nc <- 5
# Run  JAGS

library("jagsUI") # JAGS works but WinBUGS does not!

# Run JAGS: This fails quite often ('invalid parent node'), just keep trying
#Model 1
outTE1 <- jags(data, inits, params, "model.txt", n.thin=nt,n.chains=nc,
               n.burnin=nb,n.iter=ni, parallel = TRUE)
print(outTE1,3)

#Model 2
outTE2 <- jags(data, inits, params, "model.txt", n.thin=nt,n.chains=nc,
               n.burnin=nb,n.iter=ni, parallel = TRUE)
print(outTE2,3)

#Model 3
outTE3 <- jags(data, inits, params, "model.txt", n.thin=nt,n.chains=nc,
               n.burnin=nb,n.iter=ni, parallel = TRUE)
print(outTE3,3)

outTE3$mean$lambda



outTE1$mean$gamma[8]
outTE1$mean$Gavail[8]
outTE1$mean$lambda[8]
outTE1$mean$G[8]
outTE1$mean$N[8]
outTE1$mean$phi[8]


traceplot(outTE1)  
sum(data$groupsize)

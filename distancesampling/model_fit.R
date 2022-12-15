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
y3d<-array(data=0, dim=c(nsites,nD,K))
groups$dclass<-as.numeric(cut(groups$dist_to_transect,seq(0,B*1000,delta*1000)))         
#exclude observations that can't be binned into a dclass (farther out than B)
groups<-st_drop_geometry(groups[!is.na(groups$dclass)&groups$Species=="pronghorn",])
testarr<-table(factor(groups$site_id,levels=levels(as.factor(sites$site_id))),
                groups$dclass,groups$replicate,useNA="ifany")
# names are now out of sequential order
for (i in 1:K){
  y3d[,,i]<-testarr[order(as.numeric(row.names(testarr[,,i]))),,i]
}

           
# Bundle and summarize the data set
nobs <- apply(y3d, c(1,3), sum) # Total detections per site and occasion
str( data <- list(y3d=y3d, nsites=nsites, K=K, nD=nD, midpt=midpt, delta=delta,
                  roads=roads,land=land, B=B, nobs = nobs, nreps=nreps,
                  area=sites$site_length_km*2*B,
                  groupsize=groups$Total.Group.Size,
                  site=groups$site_id,
                  ngroup=length(groups$site_id)))


#Define model. From 9.5.4.1, page 501.
cat("
model {
# Prior distributions
beta0 ~ dnorm(0, 0.01) # Intercept for log(lambda)
mean.lam <- exp(beta0)
beta1 ~ dnorm(0, 0.01) # Coefficient of roads on lambda
beta2 ~ dnorm(0,0.01) # Coefficient of public land
phi0 ~ dunif(0,1) # Probability of availability
alpha0 <- exp(phi0)/(1+exp(phi0))
alpha1 ~ dnorm(0,0.01)
alpha2 ~ dnorm(0,0.01)
sigma ~ dunif(0.01,5) # Distance function parameter
gamma0 ~ dnorm(0,0.01)
gamma1 ~ dnorm(0,0.01)
gamma2 ~ dnorm(0,0.01)

r ~ dunif(0,20)


  for(b in 1:nD){
    log(g[b]) <- -midpt[b]*midpt[b]/(2*sigma*sigma) # half-normal
    f[b] <- delta/B # density function
    cellprobs[b] <- g[b]*f[b]
    cellprobs.cond[b] <- cellprobs[b]/sum(cellprobs[1:nD])
  }
  cellprobs[nD+1] <- 1-sum(cellprobs[1:nD])

  for (s in 1:nsites) {
    
    logit(phi[s]) <- alpha0 + alpha1*roads[s] + alpha2*land[s]

    for (k in 1:nreps[s]) {
      pdet[s,k] <- sum(cellprobs[1:nD]) # Distance class probabilities
      pmarg[s,k] <- pdet[s,k]*phi[s] # Marginal probability
      y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[1:nD], nobs[s,k])
      nobs[s,k] ~ dbin(pmarg[s,k], M[s])
    } 

    for (k in 1:K){
      Navail[s,k] ~ dbin(phi[s], M[s])
    }
  M[s] ~ dpois(lambda[s])
  log(lambda[s]) <- log(area[s]) + beta0 + beta1*roads[s] +beta2*land[s]
  } 

  for(i in 1:ngroup){
    groupsize[i] ~ dnegbin(probs[site[i]],r)
  }

  for(s in 1:nsites){
    log(gamma[s]) <- gamma0 + gamma1*roads[s] + gamma2*land[s]
    probs[s]<-r/(r+gamma[s])
  }
  
# Derived quantities
  Mtot <- sum(M[])

  for(k in 1:K){
    Ntot[k]<- sum(Navail[,k])
    D[k]<-sum(Navail[,k]*gamma[])/(2*B*434.8) # total length of all transects is 434.8 km
  }
} 
",file="model.txt")
# Assemble the initial values and parameters to save for JAGS
Navail.st <- apply(y3d, c(1,3),sum)
Mst <- apply(Navail.st, c( 1), max) + 2
inits <- function(){
  list(M=Mst, sigma = 0.3, phi0=0.2, beta0=-0.5, beta1=0)
}
params <- c("sigma", "phi0","alpha1","alpha2","mean.lam", "beta1","beta2","gamma0","gamma1","gamma2","D")
# MCMC settings
ni <- 60000 ; nb <- 10000 ; nt <- 100 ; nc <- 5
# Run  JAGS

library("jagsUI") # JAGS works but WinBUGS does not!

# Run JAGS: This fails quite often ('invalid parent node'), just keep trying
outTE1 <- jags(data, inits, params, "model.txt", n.thin=nt,n.chains=nc,
               n.burnin=nb,n.iter=ni, parallel = TRUE)

print(outTE1,3)
traceplot(outTE1)  
sum(data$groupsize)

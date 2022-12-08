
library(rjags)

# Bundle and summarize the data set
nobs <- apply(y3d, c(1,3), sum)  # Total detections per site and occasion

# DC add vectors of road density and prop. public land for all sites date as spatial covariates for lambda
# DC add vector of the total number of reps at each site, nreps 
# DC add constant for total number of groups, ngroup
# DC add vector of group sizes for all groups, groupsize
# DC add vector of site ID's for all groups, site
str( data <- list(y3d=y3d, nsites=nsites, K=K, nD=nD, midpt=midpt, delta=delta,
                  roads=roads, land=land, B=B, nobs = nobs, nreps=nreps,
                  ngroup=sum(y3d), groupsize=groups$Total.Group.Size,
                  site=groups$site_id, area=sites$site_length_km*2*B) )



cat("
model {

  beta0 ~ dnorm(0, 0.01)  
  mean.lam <- exp(beta0)
  beta1 ~ dnorm(0, 0.01) # DC this is now the road density coefficeint  
  beta2 ~ dnorm(0, 0.01) # DC add another prior for land coeffient
  phi ~ dunif(0,1)        
  sigma ~ dunif(0.01,5)   
  gamma0 ~ dnorm(0, 0.01) # DC add a prior for groupsize intercept
  gamma1 ~ dnorm(0, 0.01) # DC add a prior road density coefficient on groupsize
  gamma2 ~ dnorm(0, 0.01) # DC add another prior for land coefficent on groupsize

  for(b in 1:nD){
    log(g[b]) <- -midpt[b]*midpt[b]/(2*sigma*sigma) 
    f[b] <- delta/B     #DC change to uniform density function for line transects
    cellprobs[b] <- g[b]*f[b]
    cellprobs.cond[b] <- cellprobs[b]/sum(cellprobs[1:nD])
  }
  cellprobs[nD+1]<- 1-sum(cellprobs[1:nD])
  for (s in 1:nsites) {
    for (k in 1:nreps[s]) {  #DC change constant K to s indexed nreps
      pdet[s,k] <- sum(cellprobs[1:nD])   
      pmarg[s,k] <- pdet[s,k]*phi        

      y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[1:nD], nobs[s,k])

      nobs[s,k] ~ dbin(pmarg[s,k], M[s])

    }  
    for (k in 1:K){
        Navail[s,k] ~ dbin(phi, M[s])  #move Navail from the limited k (above) to full k loop over K=6 replicates at each site
    }

    M[s] ~ dpois(lambda[s])
    log(lambda[s]) <- log(area[s]) + beta0 + beta1*roads[s] + beta2*land[s]  #DC add roads and land and offset term for search area size
    log(gamma[s]) <- gamma0 + gamma1*roads[s] + gamma2*land[s] #DC add GLM for expected group size at each site
  }  
  
  for(i in 1:ngroup){
    groupsize[i] ~ dpois(gamma[site[i]]) # need to add in a new loop over all groups
  }
  
  Mtot <- sum(M[])
  for(k in 1:K){
    Ntot[k]<- sum(Navail[,k])
    D[k]<-sum(Navail[,k]*gamma[])/(2*B*434.8) # total length of all transects is 434.8 km
  }
} 
",file="model_orig.txt")


Navail.st <- apply(y3d, c(1,3),sum)
Mst <- apply(Navail.st, c( 1), max)  + 2 # change +2 to +10
inits <- function(){
  list(M=Mst, sigma = 0.2, phi=0.3, beta0=log(2), beta1=0.5) # DC change intial sigma and phi
}
params <- c("sigma", "phi", "beta0", "mean.lam", "beta1","beta2","gamma0","gamma1","gamma2","D", "Ntot")

# MCMC settings
ni <- 60000   ;   nb <- 10000   ;   nt <- 100   ;   nc <- 5 # change thinning rate and chains


# Run JAGS: This fails quite often ('invalid parent node'), just keep trying
outTE1 <- jags(data, inits, params, "model_orig.txt",
               n.thin=nt,n.chains=nc, n.burnin=nb,n.iter=ni, parallel = TRUE)
traceplot(outTE1)   
print(outTE1, 3)            # ART 4 min

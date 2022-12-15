
# this script will successfully estimate pronghorn density as INDIVIDUALS per sq km. only when you replace the ################################ sections with appropriate code, will it run correctly

rm(list=ls())
library(sf)
library(jagsUI)

# the data are stored in a geopackage as spatial vector objects. you need site information (which also holds replicate information at each site). you also need group information (the observations of large mammals)
sites<-st_cast(st_read("ecol5200distance.gpkg",layer="sites"),to="LINESTRING")
groups<-st_cast(st_read("ecol5200distance.gpkg",layer="groups"),to="POINT")

# visual inspection (your plot window needs to be pretty big to see this), remember, each segment on a transect route iss a unique 'site'
plot(sites)
plot(groups)

# distance sampling parameters that describe our study and model design. you will want to make surey everything is in km's
nD<-6 # number of distance class bins
B<-0.6 # strip half-width in km                
delta<-B/nD # distance between bins in km
midpt<-seq(delta/2,B,delta) # midpoint of each distance class bin in km

# site replicate info: this code will determine how many times each site was survyed replicates are only identified sequentially, not strictly chronologically therefore, replicate 1 simply refers to the first time a site was surveyed and replicate 'nrep[s]' refers to the last time a site was surveyed the sequence 1:nrep[s] identifies all the replicates in which a site was surveyed
repcols<-paste0("rep.",seq(1,max(sites$totalreps)))
rep_effort<-1*st_drop_geometry(sites[,repcols])
K<-max(sites$totalreps)
nreps<-sites$totalreps

# extract spatial covariates (by site) to explain variation in abundance
roads<-sites$road_kmspersqkm
land<-sites$prop_publand


# convert distances to site to distance classes. distances were saved as meters but we want to convert to km's
groups$dclass<-as.numeric(cut(groups$dist_to_transect,seq(0,B*1000,delta*1000)))     

# exclude observations of species other than pronghorn and exclude pronghorn greater than B from a site
groups<-st_drop_geometry(groups[!is.na(groups$dclass)&groups$Species=="pronghorn",])


# extract pronghorn group covariates
groupsize=groups$Total.Group.Size
site=groups$site_id

# now organize detections in 3-D array.  create empty array to hold the number of gorups counted in each distance class bin at each site in each replicate up to the maximum number of replicates (6)
y3d<-array(data=NA, dim=c(length(sites$site_id),nD,K)) 

# create a table of counts of groups in each distance class at each site in each replicate, this takes the group data and restructures it as site by replicate data. Groups lose their identity. This table will also produce the 'zeroes' needed to fit our data correclty, filling in a 0 count where and when no groups are found.
testarr<-table(factor(groups$site_id,levels=levels(as.factor(sites$site_id))),
               groups$dclass,groups$replicate)

# populate the y3d array with the site counts in each distance bin in correct site order (the table function puts them in alphabetical order and this loop puts them back in numerical order)
for (i in 1:K){
  y3d[,,i]<-testarr[order(as.numeric(row.names(testarr[,,i]))),,i]
}

# The three-part mutinomial uses our group observations twice, once in the y3d array (broken down by distance class) and once in 'nobs' a total count across all distance bins of groups by site and replicate.  This code also introduces 'zeroes' into our data whereever and whenever we went through a site and did not see any pronghorn. 
nobs <- apply(y3d, c(1,3), sum) # Total detections per site and occasion

# a few other tidbits of data that make the model easier to write
ngroup=length(groups$site_id) # Total nubmer of pronghorn groups observed (needed for loop through each group for groupsize)
nsites<-as.numeric(max(sites$site_id)) # Total number of sites (needed for loop through each site)
area<-sites$site_length_km*2*B # surveyed area at each site (needed for site-specific density)

# bundle data for jags

str( data <- list(y3d=y3d, 
                  nsites=nsites, 
                  K=K, 
                  nD=nD, 
                  midpt=midpt, 
                  delta=delta,
                  roads=roads,
                  land=land, 
                  B=B, 
                  nobs = nobs, 
                  nreps=nreps,
                  area=area,
                  groupsize=groupsize,
                  site=site,
                  ngroup=ngroup))


# this code is adapted from chapter 9, section 9.5.4.1, page 501
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
    for (k in 1:nreps[s]) { 
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
    D[k]<-sum(Navail[,k]*gamma[])/(2*B*434.8) # total LENGTH of all transects is 434.8 km
  }
  Davg<-mean(D[])  # overall average density of pronghorn across all replicates
  
} 
",file="model_orig.txt")


Navail.st <- apply(y3d, c(1,3),sum)
Mst <- apply(Navail.st, c( 1), max)  + 2 # change +2 to +10
inits <- function(){
  list(M=Mst, sigma = 0.2, phi=0.3, beta0=log(2), beta1=0.5) # DC change intial sigma and phi
}
params <- c("Davg","sigma", "phi", "beta0", "mean.lam", "beta1","beta2","gamma0","gamma1","gamma2","D", "Ntot")

# MCMC settings
ni <- 60000   ;   nb <- 10000   ;   nt <- 100   ;   nc <- 5 # change thinning rate and chains


# Run JAGS: This fails quite often ('invalid parent node'), just keep trying
outTE1 <- jags(data, inits, params, "model_orig.txt",
               n.thin=nt,n.chains=nc, n.burnin=nb,n.iter=ni, parallel = TRUE)
#traceplot(outTE1)   
print(outTE1, 3)            # ART 4 min

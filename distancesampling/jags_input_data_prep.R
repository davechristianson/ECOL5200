rm(list=ls())
library(sf)
library(rjags)

# the data are stored in a geopackage as spatial vector objects
# you need site information (which also holds replicate information at each site)
# you also need group information (the observations of large mammals)
sites<-st_cast(st_read("ecol5200distance.gpkg",layer="sites"),to="LINESTRING")
groups<-st_cast(st_read("ecol5200distance.gpkg",layer="groups"),to="POINT")

plot(sites)
points(groups)

# distance sampling parameters
nD<-#####################
B<-####################                 
delta<-#########################
midpt<-###########################


# site replicate info
# this will collect the necessary information about replicates from the surveys
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

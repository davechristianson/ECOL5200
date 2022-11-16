
## start all scripts with emptying the environment
rm(list=ls())

library(googlesheets4)
library(sp)
library(sf)
library(geosphere)
library(dplyr)

# first, we will clean up the detections data from the Google worksheet and extract the information we need for each observation.
###########################################################
# detections

# here is a backup line of code if trouble connecting to Google sheet.
#detections<-read.csv("Detections - fall22.csv",skip=1)

# read data from Google sheet (select '1:Yes' when prompted and sign into googl and agree to give tidyverse access)
detections<-read_sheet("https://docs.google.com/spreadsheets/d/1_TOxrZ5MtqTASe3ZSz2RnX2sM2wivyXYvDPUEKaim6g/edit#gid=0",skip=1)
detections<-as.data.frame(detections)

################################################################
# some species were entered in more than one way. this will make them consistent
detections$Species[detections$Species %in% c("pronghorn","Pronghorn")]<-"pronghorn"
detections$Species[detections$Species %in% c("Mule Deer","mule deer")]<-"muledeer"
detections$Species[detections$Species %in% c("End","END")]<-"end"
detections$Species[detections$Species %in% c("Start","START")]<-"start"

# some longitudes were not entered as negative, but should be.
detections$`Detection Longitude`<- ((detections$`Detection Longitude`>0)*-1*detections$`Detection Longitude`) + ((detections$`Detection Longitude`<0)*1*detections$`Detection Longitude`) 

# now the data should be ready to be convert to an sf (spatial) object in R
detections<-st_as_sf(detections,coords=c("Detection Longitude","Detection Latitude"),crs = st_crs(4326))

# some distances are not numbers - they are surrounded by text
detections$`Detection Distance`<-as.numeric(gsub(".*?([0-9]+).*", "\\1",  detections$`Detection Distance`))

# now use destPoint() function to get the lat and long of the detected animal groups (not the observer)
# destPoint() uses observer's lat and long as well as the bearing and distance to project the animal groups location
df_sp <- as_Spatial(detections)
destPoint(as.matrix(coordinates(df_sp))[!is.na(as.numeric(df_sp$Detection.Distance)),],d=as.numeric(df_sp$Detection.Distance)[!is.na(as.numeric(df_sp$Detection.Distance))], b=as.numeric(df_sp$Detection.Angle)[!is.na(as.numeric(df_sp$Detection.Distance))]) %>%
  as.data.frame() %>%
  st_as_sf(coords = c('lon', 'lat')) -> df_sf
st_crs(df_sf)<-st_crs(4326)

# now create a new spatial object that contains only groups and the groups' lat and long
groups<-detections[!is.na(as.numeric(df_sp$Detection.Distance)),] # only detections (remove START and END)
groups[,c("obs_longitude","obs_latitude")]<-st_coordinates(groups) # save the observer lat and long as data columns
st_crs(groups)<-st_crs(4326)
st_geometry(groups)<-st_geometry(df_sf) # new coordinates from destPoint() above

# two pieces of information are needed for hierarchical modelling of the data - shortest distance to transect line and ID of the nearest transect line
sites<-st_cast(st_read("ecol5200distance.gpkg",layer="sites"),to="LINESTRING") # read in the vector layer of transects from the geopackage.
sites$site_id<-row.names(sites) # created a site_id column to identify each site as an integer
groups$dist_to_transect<-apply(st_distance(groups,st_transform(sites,crs=st_crs(groups)), by_element=F),1,min) # measure distance from group to nearest transect
groups$site_id<-apply(st_distance(groups,st_transform(sites,crs=st_crs(groups)), by_element=F),1,which.min)    # identify the transect that is closest to each site

#  make a neat (but unnecessary) vector of lines from the observer point to the detected group just for display purposes
lineofsight<-st_sfc(mapply(function(a,b){st_cast(st_union(a,b),"LINESTRING")}, groups$geometry, detections$geometry[!is.na(as.numeric(df_sp$Detection.Distance))], SIMPLIFY=FALSE))
st_crs(lineofsight)<-st_crs(4326)
lineofsight<-st_sf(lineofsight)
lineofsight[,c(names(st_drop_geometry(groups)))]<-as.matrix(st_drop_geometry(groups))
lineofsight$dist_to_transect<-apply(st_distance(groups,st_transform(sites,crs=st_crs(groups)), by_element=F),1,min)

st_write(detections,dsn="ecol5200distance.gpkg",layer="detections",append=F)
st_write(lineofsight,dsn="ecol5200distance.gpkg",layer="line_of_sight",append=F)
st_write(groups,dsn="ecol5200distance.gpkg",layer="groups",append=F)



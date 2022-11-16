##################################
#Site Covariates
## read in the sites data and add a new site_id
sites<-st_cast(st_read("ecol5200distance.gpkg",layer="sites"),to="LINESTRING")
sites$site_id<-row.names(sites)

# read in all roads vector file
allroads<-st_cast(st_read("ecol5200distance.gpkg",layer="wyohub_sa_roads"),to="MULTILINESTRING")

# buffer all site transect segmetns by 2000 meters
sites_buffer<-st_buffer(sites,dist=2000)

# read in public lands vector file
lands<-st_make_valid(st_cast(st_read("ecol5200distance.gpkg",layer="wyohub_sa_landuse"),to="MULTIPOLYGON"))

# get new road lines that are labelled by the site buffer polygon they are in
inters<-st_intersection(allroads,sites_buffer)

# calculate the lenghth of road lines in each site buffer polygon
road_lengths<-tapply(st_length(inters), inters$site_id,sum)

# convert road lengths to km's and divide by the area of the site buffered polygon, and append to vector
sites$road_kmspersqkm[order(sites_buffer$site_id)]<-(road_lengths/1000)/(as.numeric(st_area(st_transform(sites_buffer,crs=st_crs(lands))[sort(sites_buffer$site_id),]))/(1000*1000))
sites_buffer$road_kmspersqkm[order(sites_buffer$site_id)]<-(road_lengths/1000)/(as.numeric(st_area(st_transform(sites_buffer,crs=st_crs(lands))[sort(sites_buffer$site_id),]))/(1000*1000))

# get new land use polygons that are labelled by the site buffer polygon they are in 
land_inters<-st_intersection(lands,st_transform(sites_buffer,crs=st_crs(lands)))

# calculate the total area of public land in each site buffered polygon
land_areas<-tapply(st_area(land_inters[!land_inters$Name %in% c("Water","Private"),]), land_inters$site_id[!land_inters$Name %in% c("Water","Private")],sum)

# convert areas to square kms and divide by the area of the site buffer polygon to get proportion of site covered in public land, append to vector
sites$prop_publand[order(sites_buffer$site_id)]<-(land_areas/(1000*1000)) /(as.numeric(st_area(st_transform(sites_buffer,crs=st_crs(lands))[order(sites_buffer$site_id),]))/(1000*1000))
sites_buffer$prop_publand[order(sites_buffer$site_id)]<-(land_areas/(1000*1000)) /(as.numeric(st_area(st_transform(sites_buffer,crs=st_crs(lands))[order(sites_buffer$site_id),]))/(1000*1000))

st_write(sites,dsn="ecol5200distance.gpkg",layer="sites",append=F)
st_write(sites_buffer,dsn="ecol5200distance.gpkg",layer="sites_buffer",append=F)

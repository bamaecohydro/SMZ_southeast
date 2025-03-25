#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Title: SMZ demo script
# Date: 6/21/2024; edited 11/24/24
# Coder: Delaney Peterson and Nate Jones (natejones@ua.edu)
# Purpose: explore topographic metrics associated with smz shapefiles
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# NOTE: will have to run on all 6 of the data chunks, but probably will run
# East and West Gulf on the lab computer

#Piedmont.shp : DONE (1853)             {01:11 hr:mm}
#CPlainMisc.shp : DONE (1604)           {00:55 hr:mm}
#AppalachianPlateaus.shp : DONE (3037)  {01:54 hr:mm}
#Ouachita.shp : DONE [labcomp] (6250)   {01:15 hr:mm}
#EGCPlain.shp : DONE [labcomp] (25532)  {02:45 hr:mm}
#WGCPlain.shp : DONE [labcomp] (20842)  {~2:45 hr:mm}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setup workspace ---------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clear workspace
remove(list=ls())

#load libraries
library(tidyverse)
library(whitebox)
library(sf)
library(raster)
library(elevatr)
library(mapview)
library(parallel)

#set folders
data_dir <- "/Users/Delaney/R/data_dir/WaterSci/Final_Project/data/SMZ_physio/" 
temp_dir <- "/Users/Delaney/R/temp_dir/" 

#Load SMZ shape
SMZs <- st_read(paste0(data_dir, "AppalachianPlateaus.shp"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create function to find "channel head" ----------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
channel_head_fun<-function(n){

#load libraries~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(whitebox)
library(sf)
library(raster)
library(elevatr)

#Isolate SMZ of interest ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#subset to single SMZ
SMZ <- SMZs[n,] 

#identify unique identifier
OBJECTID <- SMZ$OBJECTID

#Define master crs
crs <- crs(SMZ)

#Download DEM
dem <- get_elev_raster(SMZ, z = 14)

# Delineate flow network ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Write raster to temp file
writeRaster(dem, paste0(temp_dir, 'dem_', OBJECTID,'.tif'), overwrite=T)

#Fill single cell pits
wbt_fill_single_cell_pits(
  dem = paste0(temp_dir, 'dem_', OBJECTID,'.tif'),
  output = paste0(temp_dir, 'dem_fill_', OBJECTID,'.tif'),  
  wd = temp_dir
)

#breach depressions 
wbt_breach_depressions(
  dem = paste0(temp_dir, 'dem_fill_', OBJECTID,'.tif'),
  output = paste0(temp_dir, 'dem_breach_', OBJECTID,'.tif'), 
  wd= temp_dir
)

#flow direction
wbt_d8_pointer(
  dem = paste0(temp_dir, 'dem_breach_', OBJECTID,'.tif'),
  output = paste0(temp_dir, 'fdr_', OBJECTID,'.tif'),
  wd = temp_dir)

#flow accumulation
wbt_d8_flow_accumulation(
  input = paste0(temp_dir, 'dem_breach_', OBJECTID,'.tif'),
  output = paste0(temp_dir, 'fac_', OBJECTID,'.tif'),
  pntr = F,
  wd = temp_dir
)
fac <- raster(paste0(temp_dir, 'fac_', OBJECTID,'.tif'))

#slope
wbt_slope(
  dem = paste0(temp_dir, 'dem_', OBJECTID,'.tif'),
  output = paste0(temp_dir, 'slope_', OBJECTID,'.tif'),
  wd = temp_dir
)
slope <- raster(paste0(temp_dir, 'slope_', OBJECTID,'.tif'))

#try curvature
wbt_total_curvature(
  dem = paste0(temp_dir, 'dem_', OBJECTID,'.tif'),
  output = paste0(temp_dir, 'curvature_', OBJECTID,'.tif'),
  wd = temp_dir
)
curv <- raster(paste0(temp_dir, 'curvature_', OBJECTID,'.tif'))

#Extract streams
wbt_extract_streams(
  flow_accum = paste0(temp_dir, 'fac_', OBJECTID,'.tif'),
  output = paste0(temp_dir, 'flow_net_', OBJECTID,'.tif'),
  threshold = 1000, 
  wd = temp_dir)

#convert to vector
wbt_raster_streams_to_vector(
  streams = paste0(temp_dir, 'flow_net_', OBJECTID,'.tif'), 
  d8_pntr = paste0(temp_dir, 'fdr_', OBJECTID,'.tif'), 
  output = paste0(temp_dir, 'flow_net_', OBJECTID,'.shp'), 
  wd = temp_dir)

#Read into R environement
flow_net<-st_read(
  paste0(temp_dir, 'flow_net_', OBJECTID,'.shp'), 
  crs = st_crs(crs))

#clip flownet to SMZ
flow_net <- flow_net[SMZ,]

#Identify channel head ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #For now, identify one channel head per SMZ. In the future, we could identify 
  #a channel head for each "reach" of the flownet flowing through the SMZ

#Convert SMZ to a multilinestring and then identify intersecting points
inlet_outlet_pnts <- st_cast(SMZ, "MULTILINESTRING") %>% st_intersection(., flow_net) %>% st_cast("POINT")

#Identify channel head based on fac value
channel_head <- inlet_outlet_pnts %>% 
  mutate(fac = extract(fac, inlet_outlet_pnts)) %>% 
  filter(fac == min(fac, na.rm =T)) %>% 
  slice(1)

#Estimate slope and watershed area of channel head ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output <- channel_head %>% 
  dplyr::select(OBJECTID) %>% 
  mutate(
    x          = st_coordinates(channel_head)[1], 
    y          = st_coordinates(channel_head)[2], 
    ws_area_m2 = extract(fac, channel_head)*res(fac)[1]*res(fac)[2],
    slope      = extract(slope, channel_head),
    curvature  = extract(curv, channel_head)) %>% 
  st_drop_geometry()

#Cleanup workspace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create list of temp files
files <- list.files(temp_dir) %>% as_tibble() %>% filter(str_detect(value, paste(OBJECTID)))

#Create cleanup function
cleanup_fun<-function(n){
  file <- files[n,]
  file.remove(paste0(temp_dir, file))
}
#apply cleanup function
lapply(seq(1,nrow(files)), cleanup_fun)

#print output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output
}

#Test
# t0<-Sys.time()
# channel_head_fun(5)
# tf<-Sys.time()
# tf-t0

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create errors handling wrapper ------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create wrapper function 
error_fun<-function(n){
  tryCatch(
    expr = channel_head_fun(n), 
    error = function(e)
      tibble(
        OBJECTID = SMZs$OBJECTID[n], 
        x        =   NA, 
        y        =   NA, 
        ws_area_m2 = NA, 
        slope      = NA,
        curvature  = NA)
  )
}  

#test function
# t0<-Sys.time()
# lapply(
#    X=seq(1, 5),
#    FUN=error_fun) %>%
#  bind_rows()
#tf<-Sys.time()
#tf-t0

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Run function in parallel ------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Determine number of processing cores available on your machine
n.cores<-detectCores()

#Create clusters
cl<-makeCluster(n.cores)

#Send libraries to cluster
clusterEvalQ(cl, {
  library(tidyverse)
  library(whitebox)
  library(sf)
  library(raster)
  library(elevatr)
})

#Export data to cluter environments
clusterExport(cl, c("channel_head_fun", "SMZs", "data_dir", "temp_dir"))

#Now run function
output<-parLapply(
  cl=cl,
  seq(1, nrow(SMZs)), 
  error_fun)

#Now, bind rows from list output
df<-output %>% bind_rows()
df

#Stop the clusters
stopCluster(cl)

write_csv(df, "/Users/Delaney/R/data_dir/WaterSci/Final_Project/data/output_aplateaus.csv")






mapview(SMZ)
#is there a difference between states????




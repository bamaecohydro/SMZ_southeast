#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Title: Filiter shape files
# Date: 4/1/2025
# Coder: Delaney Peterson and Nate Jones (natejones@ua.edu)
# Purpose: Filter SMZs where "inlet" occurs at percieved property boundary
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#It turns out we need to do a network analysis to find the most upstream reaches. 
#   Here...we will define reach as stream lengths between conluences.  We will merge all SMZs along 
#   each reach.  Then, we define the most upstream reaches using some sort of connectivity magic.
#   after that, then we will proced as planned.  

#   Yes -- I still need to figure out how do deal with property boundaries...

#   Also, we need to deal with threshold for each region...


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
data_dir <- "data/SMZ_physio/" 

#Load SMZ shape
SMZs <- st_read(paste0(data_dir, "WGCPlain.shp"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create filtering function to identify property lines --------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#channel_head_fun<-function(n){
  
  #For testing
  n <- 16

  #load libraries~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  library(tidyverse)
  library(whitebox)
  library(sf)
  library(raster)
  library(elevatr)

  #Create temp directory
  temp_dir <- tempdir()
  dir.create(temp_dir)
  
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
  
  mapview(flow_net) + mapview(SMZ) + mapview(inlet_outlet_pnts) + mapview(channel_head, col.regions="red") + mapview(SMZs, col.regions="green")
  
#   #Estimate slope and watershed area of channel head ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   output <- channel_head %>% 
#     dplyr::select(OBJECTID) %>% 
#     mutate(
#       x          = st_coordinates(channel_head)[1], 
#       y          = st_coordinates(channel_head)[2], 
#       ws_area_m2 = extract(fac, channel_head)*res(fac)[1]*res(fac)[2],
#       slope      = extract(slope, channel_head),
#       curvature  = extract(curv, channel_head)) %>% 
#     st_drop_geometry()
#   
#   #Cleanup workspace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   #Create list of temp files
#   files <- list.files(temp_dir) %>% as_tibble() %>% filter(str_detect(value, paste(OBJECTID)))
#   
#   #Create cleanup function
#   cleanup_fun<-function(n){
#     file <- files[n,]
#     file.remove(paste0(temp_dir, file))
#   }
#   #apply cleanup function
#   lapply(seq(1,nrow(files)), cleanup_fun)
#   
#   #print output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   output
# }
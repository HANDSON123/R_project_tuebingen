## Script for downloadind usefull datasets

# Load necessary library

library(malariaAtlas)
library(tidyverse)
library(dplyr)
library(terra)
library(sf)
library(dplyr)
library(chirps)
library(ggplot2)
library(ggmap)
library(curl)
library(readr)
library(readxl)
library(openxlsx)
library(wpgpDownloadR)

# Set the working directory

setwd("C:/Users/USER/Documents/Research project tubingen/R_project_tuebingen/")

### Downloading the rainfall dataset (for Cameroon, Congo, Benin, Gabon, Tanzania)


# URL
url <- "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/netcdf/chirps-v3.0.monthly.nc"


# Download

curl_download(url, "chirps_monthly.nc")

# load CHIRPS data
chirps_cmr <- rast("chirps_monthly.nc")

# Case 1: Cameroon

# load the shape file

shp_cmr <- st_read("cmr_data/CMR_shape_file/cameroun_districts_2024.shp")

district_vect <- vect(shp_cmr)


#chirps <- readRast("chirps_monthly.nc")
chirps_cmr <- crop(chirps_cmr, district_vect)
dates <- time(chirps_cmr)

idx <- which(format(dates, "%Y") %in% c("2020", "2021", "2022", "2023"))
chirps_2020_2023_cmr <- chirps_cmr[[idx]]


rain_df <- terra::extract(chirps_2020_2023_cmr, district_vect, mean, na.rm=TRUE)
rain_df$district <- district_vect$Admin2[rain_df$ID]

# Assign dates to column names
colnames(rain_df)[2:(length(idx)+1)] <- as.character(dates[idx])

rain_long <- pivot_longer(
  rain_df,
  cols = -c(ID, district),
  names_to = "date",
  values_to = "rainfall"
) %>% select(-c(ID))


write_csv(rain_long, "cmr_data/cmr_rainfall_data.csv")




#loading the population data

## 1: Benin

ben_pop <- rast("ben_data/ben_pd_2020_1km_UNadj.tif")

ben_pop_2020 <- extract(
  ben_pop,
  ben_vect,
  fun = mean,
  na.rm = TRUE)

ben_district_pop_2020 <- shp_ben %>% 
  st_drop_geometry() %>% 
  select(NAME_3) %>% 
  rename(District=NAME_3) %>% 
  mutate(population_density_2021 = ben_pop_2020[,2])

write_csv(ben_district_pop_2020, file = "ben_data/ben_pop_density_2020.csv")

## 2: Cameroon

#2011
cmr_pop <- rast("cmr_data/cmr_pd_2011_1km_UNadj.tif")
cmr_pop_2011 <- extract(
  cmr_pop,
  cmr_vect,
  fun = mean,
  na.rm = TRUE)


cmr_district_pop_2011 <- shp_cmr %>% 
  st_drop_geometry() %>% 
  select(Admin2) %>% 
  mutate(pop_density_2011 = cmr_pop_2011[,2])%>% 
  rename(District=Admin2)


write_csv(cmr_district_pop_2011, file = "cmr_data/cmr_pop_density_2011.csv")


#2020
cmr_pop <- rast("cmr_data/cmr_pd_2020_1km_UNadj.tif")
cmr_pop_2020 <- extract(
  cmr_pop,
  cmr_vect,
  fun = mean,
  na.rm = TRUE)

cmr_district_pop_2020 <- shp_cmr %>% 
  st_drop_geometry() %>% 
  select(Admin2) %>% 
  rename(District=Admin2) %>% 
  mutate(population_density_2021 = cmr_pop_2020[,2])

write_csv(cmr_district_pop_2020, file = "cmr_data/cmr_pop_density_2020.csv")

# 3: Congo

cong_pop <- rast("cong_data/cog_pd_2020_1km_UNadj.tif")
cong_pop_2020 <- extract(
  cong_pop,
  con_vect,
  fun = mean,
  na.rm = TRUE)

cong_district_pop_2020 <- shp_con %>% 
  st_drop_geometry() %>% 
  select(adm2_name) %>% 
  rename(District=adm2_name) %>% 
  mutate(population_density_2021 = cong_pop_2020[,2])

write_csv(cong_district_pop_2020, file = "cong_data/cong_pop_density_2020.csv")


#Rainfall data

##1: Benin

path1 <- "CHIRPS_PENTAD_precipitation_2021-01-01_2021-12-31_0p9311to13p9448N_8p2516to17p0602E.precipitation.tif"

ben_rainfall_2021 <- rast(path1)

rain_extracted <- extract(
  ben_rainfall_2021,
  ben_vect,
  fun   = "mean",
  na.rm = TRUE
)

ben_rainfall_final_2021 <- shp_ben %>% 
  st_drop_geometry() %>% 
  select(NAME_3) %>% 
  mutate(mean_rainfall_2021 = rain_extracted[, 2])%>% 
  rename(District=NAME_3)

write_csv(ben_rainfall_final_2021, "ben_data/ben_rainfall_final_2021.csv")



##2: Cameroon (2011 and 2021)

#2011

path <- "cmr_data/chirps_2011/CHIRPS_PENTAD_precipitation_2011-01-01_2011-12-31_-5p8264to13p9448N_0p6930to19p1696E.precipitation.tif"

cmr_rainfall_2011 <- rast(path)

rain_extracted <- extract(
  cmr_rainfall_2011,
  cmr_vect,
  fun   = "mean",
  na.rm = TRUE
)

cmr_rainfall_final_2011 <- shp_cmr %>% 
  st_drop_geometry() %>% 
  select(Admin2) %>% 
  mutate(mean_rainfall_2021 = rain_extracted[, 2])%>% 
  rename(District=Admin2)

write_csv(cmr_rainfall_final_2011, "cmr_data/cmr_rainfall_2011.csv")

#2021

cmr_rainfall_2021 <- rast(path1)

rain_extracted <- extract(
  cmr_rainfall_2021,
  cmr_vect,
  fun   = "mean",
  na.rm = TRUE
)

cmr_rainfall_final_2021 <- shp_cmr %>% 
  st_drop_geometry() %>% 
  select(Admin2) %>% 
  mutate(mean_rainfall_2021 = rain_extracted[, 2])%>% 
  rename(District=Admin2)

write_csv(cmr_rainfall_final_2021, "cmr_data/cmr_rainfall_final_2021.csv")

# 3: Congo 

cong_rainfall_2021 <- rast(path1)

rain_extracted <- extract(
  cong_rainfall_2021,
  con_vect,
  fun   = "mean",
  na.rm = TRUE
)

cong_rainfall_final_2021 <- shp_con %>% 
  st_drop_geometry() %>% 
  select(adm2_name) %>% 
  mutate(mean_rainfall_2021 = rain_extracted[, 2])%>% 
  rename(District=adm2_name)

write_csv(cong_rainfall_final_2021, "cong_data/cong_rainfall_final_2021.csv")

#temperature data

path2 <- "ERA5_LAND_DAILY_temperature_2m_2021-01-01_2021-12-31_1p1068to13p6033N_8p1637to16p8844E.temperature_2m.tif"

##1: Benin

ben_temperature_2021 <- rast(path2)

temperature_extracted <- extract(
  ben_temperature_2021,
  ben_vect,
  fun   = "mean",
  na.rm = TRUE)


ben_temperature_final_2021 <- shp_ben %>% 
  st_drop_geometry() %>% 
  select(NAME_3) %>% 
  mutate(mean_temperature_2021 = temperature_extracted[, 2])%>% 
  rename(District=NAME_3)

write_csv(ben_temperature_final_2021, "ben_data/ben_temperature_2021.csv")


## 2: Cameroon

# 2011

path <- "cmr_data/ERA5_temperature/ERA5_LAND_DAILY_temperature_2m_2011-01-01_2011-12-31_-5p4765to13p9448N_0p6930to19p3454E.temperature_2m.tif"

cmr_temperature_2011 <- rast(path)

temperature_extracted <- extract(
  cmr_temperature_2011,
  cmr_vect,
  fun   = "mean",
  na.rm = TRUE)


cmr_temperature_final_2011 <- shp_cmr %>% 
  st_drop_geometry() %>% 
  select(Admin2) %>% 
  mutate(mean_temperature_2011 = temperature_extracted[, 2])%>% 
  rename(District=Admin2)

write_csv(cmr_temperature_final_2011, "cmr_data/cmr_temperature_2011.csv")


# 2021

cmr_temperature_2021 <- rast(path2)

temperature_extracted <- extract(
  cmr_temperature_2021,
  cmr_vect,
  fun   = "mean",
  na.rm = TRUE)


cmr_temperature_final_2021 <- shp_cmr %>% 
  st_drop_geometry() %>% 
  select(Admin2) %>% 
  mutate(mean_temperature_2021 = temperature_extracted[, 2])%>% 
  rename(District=Admin2)

write_csv(cmr_temperature_final_2021, "cmr_data/cmr_temperature_2021.csv")

# 3: Congo

cong_temperature_2021 <- rast(path2)

temperature_extracted <- extract(
  cong_temperature_2021,
  con_vect,
  fun   = "mean",
  na.rm = TRUE)


cong_temperature_final_2021 <- shp_con %>% 
  st_drop_geometry() %>% 
  select(adm2_name) %>% 
  mutate(mean_temperature_2021 = temperature_extracted[, 2])%>% 
  rename(District=adm2_name)

write_csv(cong_temperature_final_2021, "cong_data/cong_temperature_2021.csv")


# Elevation data

## 1: Benin

elevation_ben <- elevation_30s(
  country = "BEN",       
  path    = "C:/Users/USER/Documents/Research project tubingen/R_project_tuebingen/ben_data/"    
)

elev_extracted_ben <- extract(
  elevation_ben,
  ben_vect,
  fun   = "mean",
  na.rm = TRUE
)

shp_ben$mean_elevation <- elev_extracted_ben[, 2]

elevation_df_ben <- shp_ben %>%
  st_drop_geometry() %>%          
  dplyr::select(                       
    NAME_3,                       
    mean_elevation                
  ) %>%
  rename(
    District = NAME_3)

write_csv(elevation_df_ben, "ben_data/ben_elevation_data.csv")

## 2: Cameroon

elevation_cmr <- elevation_30s(
  country = "CMR",       
  path    = "C:/Users/USER/Documents/Research project tubingen/R_project_tuebingen/cmr_data/"    
)
#plot(elevation_cmr, main = "Elevation Benin (metres)")

elev_extracted_cmr <- extract(
  elevation_cmr,
  cmr_vect,
  fun   = "mean",
  na.rm = TRUE
)

shp_cmr$mean_elevation <- elev_extracted_cmr[, 2]

elevation_df_cmr <- shp_cmr %>%
  st_drop_geometry() %>%          
  dplyr::select(                       
    Admin2,                       
    mean_elevation                
  ) %>%
  rename(
    District = Admin2)

write_csv(elevation_df_cmr, "cmr_data/cmr_elevation_data.csv")

# elevation data for Africa

elev_global <- elevation_global(res = 2.5, path ="C:/Users/USER/Documents/Research project tubingen/R_project_tuebingen/Africa_data/")

## 3: Congo

elevation_cong <- elevation_30s(
  country = "COG",       
  path    = "C:/Users/USER/Documents/Research project tubingen/R_project_tuebingen/cong_data/"    
)


elev_extracted_cong <- extract(
  elevation_cong,
  con_vect,
  fun   = "mean",
  na.rm = TRUE
)

elevation_df_cong <- shp_con %>%
  mutate(mean_elevation = elev_extracted_cong[, 2]) %>% 
  st_drop_geometry() %>%          
  dplyr::select(adm2_name, mean_elevation) %>%
  rename(District = adm2_name)

write_csv(elevation_df_cong, "cong_data/cong_elevation_data.csv")

## 4: Gabon

elevation_gab <- elevation_30s(
  country = "GAB",       
  path    = "C:/Users/USER/Documents/Research project tubingen/R_project_tuebingen/gab_data/"    
)

# Land cover

## 1: Benin
path = "C:/Users/USER/Documents/Research project tubingen/R_project_tuebingen/ben_data/"
lc_water   <- landcover(var = "water",    path = path)  
lc_wetland <- landcover(var = "wetland",  path = path)  
lc_trees   <- landcover(var = "trees",    path = path)  
lc_crops   <- landcover(var = "cropland",    path = path)  
lc_built   <- landcover(var = "built",    path = path)  

# Helper function — crop, mask and extract mean per district
extract_mean <- function(raster, districts_vect){
  r_crop <- crop(raster, districts_vect)
  r_mask <- mask(r_crop, districts_vect)
  extract(r_mask, districts_vect, fun = "mean", na.rm = TRUE)[, 2]
}

# Extract all classes
shp_ben$prop_water    <- extract_mean(lc_water,    ben_vect)
shp_ben$prop_wetland  <- extract_mean(lc_wetland,  ben_vect)
shp_ben$prop_trees    <- extract_mean(lc_trees,    ben_vect)
shp_ben$prop_cropland <- extract_mean(lc_crops, ben_vect)
shp_ben$prop_built    <- extract_mean(lc_built,    ben_vect)

# Save to CSV
landcover_df_ben <- shp_ben %>%
  st_drop_geometry() %>%
  dplyr::select(NAME_3,
                prop_water, prop_wetland, prop_trees,
                prop_cropland, prop_built) %>%
  rename(District = NAME_3)

write_csv(landcover_df_ben, "C:/Users/USER/Documents/Research project tubingen/R_project_tuebingen/ben_data/ben_landcover_per_district.csv")

lc_vars_ben <- landcover_df_ben %>%
  dplyr::select(prop_water, prop_wetland, prop_trees, 
                prop_cropland, prop_built)
cor_matrix <- cor(lc_vars_ben, use = "complete.obs")
print(round(cor_matrix, 2))

## 2: Cameroon

# Extract all classes
shp_cmr$prop_water    <- extract_mean(lc_water,    cmr_vect)
shp_cmr$prop_wetland  <- extract_mean(lc_wetland,  cmr_vect)
shp_cmr$prop_trees    <- extract_mean(lc_trees,    cmr_vect)
shp_cmr$prop_cropland <- extract_mean(lc_crops, cmr_vect)
shp_cmr$prop_built    <- extract_mean(lc_built,    cmr_vect)

# Save to CSV
landcover_df_cmr <- shp_cmr %>%
  st_drop_geometry() %>%
  dplyr::select(Admin2,
                prop_water, prop_wetland, prop_trees,
                prop_cropland, prop_built) %>%
  rename(District = Admin2)

write_csv(landcover_df_cmr, "cmr_data/cmr_landcover_per_district.csv")

# 3: Congo

# Extract all classes
shp_con$prop_water    <- extract_mean(lc_water,    con_vect)
shp_con$prop_wetland  <- extract_mean(lc_wetland,  con_vect)
shp_con$prop_trees    <- extract_mean(lc_trees,    con_vect)
shp_con$prop_cropland <- extract_mean(lc_crops, con_vect)
shp_con$prop_built    <- extract_mean(lc_built,    con_vect)

# Save to CSV
landcover_df_cong <- shp_con %>%
  st_drop_geometry() %>%
  dplyr::select(adm2_name,
                prop_water, prop_wetland, prop_trees,
                prop_cropland, prop_built) %>%
  rename(District = adm2_name)

write_csv(landcover_df_cong, "cong_data/cong_landcover_per_district.csv")


# Distance to water


# Load HydroRIVERS Africa
rivers <- st_read("HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp")

# Load HydroLAKES Africa
lakes  <- st_read("HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp")

# 1: Benin

shp_ben_country <- shp_ben %>%
  group_by(COUNTRY) %>% 
  summarise(geometry=st_union(geometry))

# Extract coastline from country boundary
# The outer boundary of the country includes the coastline
#coastline <- st_cast(shp_ben_country, "MULTILINESTRING")


# Clip both to your country
rivers_ben <- st_intersection(rivers, shp_ben_country)
lakes_ben  <- st_intersection(lakes,  shp_ben_country)


# Quick check
ggplot() +
  geom_sf(data = shp_ben_country, fill = NA, color = "black") +
  geom_sf(data = rivers_ben, color = "blue", size = 0.2)
# Convert lakes (polygons) to lines for distance calculation
lakes_lines <- st_cast(lakes_ben, "MULTILINESTRING")

# Make sure CRS matches
#coastline   <- st_transform(coastline,   crs = 4326)
rivers_ben  <- st_transform(rivers_ben,  crs = 4326)
lakes_lines <- st_transform(lakes_lines, crs = 4326)
shp_ben   <- st_transform(shp_ben,   crs = 4326)

# Combine rivers and lake boundaries
all_water <- bind_rows(
  rivers_ben  %>% dplyr::select(geometry),
  lakes_lines %>% dplyr::select(geometry)
)


# Get district centroids
district_centroids <- st_centroid(shp_ben)

# Compute distance from each centroid to nearest water feature
# Note: this may take a few minutes depending on how many features
dist_matrix <- st_distance(district_centroids, all_water)

# Minimum distance per district in metres
shp_ben$dist_to_water_m  <- apply(dist_matrix, 1, min)

# Convert to kilometres
shp_ben$dist_to_water_km <- shp_ben$dist_to_water_m / 1000

# Check
head(shp_ben[, c("NAME_2", "dist_to_water_km")])
summary(shp_ben$dist_to_water_km)

dist_water_df <- shp_ben %>%
  st_drop_geometry() %>%
  dplyr::select(NAME_3, dist_to_water_m, dist_to_water_km) %>%
  rename(district = NAME_3)

write_csv(dist_water_df, "ben_data/ben_distance_to_water.csv")


# 2: Cameroon

shp_cmr_country <- shp_cmr %>%
  group_by(Admin1) %>% 
  summarise(geometry=st_union(geometry))

# Clip both to your country
rivers_cmr <- st_intersection(rivers, shp_cmr_country)
lakes_cmr  <- st_intersection(lakes,  shp_cmr_country)

lakes_cmr <- st_cast(lakes_cmr, "MULTILINESTRING")

# Make sure CRS matches
rivers_cmr  <- st_transform(rivers_cmr,  crs = 4326)
lakes_cmr <- st_transform(lakes_cmr, crs = 4326)
shp_cmr   <- st_transform(shp_cmr,   crs = 4326)

# Combine rivers and lake boundaries
all_water <- bind_rows(
  rivers_cmr  %>% dplyr::select(geometry),
  lakes_cmr %>% dplyr::select(geometry)
)


# Get district centroids
district_centroids <- st_centroid(shp_cmr)

# Compute distance from each centroid to nearest water feature
# Note: this may take a few minutes depending on how many features
dist_matrix <- st_distance(district_centroids, all_water)

# Minimum distance per district in metres
shp_cmr$dist_to_water_m  <- apply(dist_matrix, 1, min)

# Convert to kilometres
shp_cmr$dist_to_water_km <- shp_cmr$dist_to_water_m / 1000

# Check
head(shp_cmr[, c("Admin2", "dist_to_water_km")])
summary(shp_cmr$dist_to_water_km)

dist_water_df <- shp_cmr %>%
  st_drop_geometry() %>%
  dplyr::select(Admin1, Admin2,
                dist_to_water_m,
                dist_to_water_km) %>%
  rename(region   = Admin1,
         District = Admin2)

write.csv(
  dist_water_df,
  "cmr_data/cmr_distance_to_water.csv")

# 3: Congo

shp_con_country <- shp_con %>%
  group_by(adm0_name) %>% 
  summarise(geometry=st_union(geometry))

# Clip both to your country
rivers_cong <- st_intersection(rivers, shp_con_country)S
lakes <- st_make_valid(lakes)
lakes_cong  <- st_intersection(lakes,  shp_con_country)

lakes_cong <- st_cast(lakes_cong, "MULTILINESTRING")

# Make sure CRS matches
rivers_cong  <- st_transform(rivers_cong,  crs = 4326)
lakes_cong <- st_transform(lakes_cong, crs = 4326)
shp_con   <- st_transform(shp_con,   crs = 4326)

# Combine rivers and lake boundaries
all_water <- bind_rows(
  rivers_cong  %>% dplyr::select(geometry),
  lakes_cong %>% dplyr::select(geometry)
)


# Get district centroids
district_centroids <- st_centroid(shp_con)

# Compute distance from each centroid to nearest water feature
# Note: this may take a few minutes depending on how many features
dist_matrix <- st_distance(district_centroids, all_water)

# Minimum distance per district in metres
shp_con$dist_to_water_m  <- apply(dist_matrix, 1, min)

# Convert to kilometres
shp_con$dist_to_water_km <- shp_con$dist_to_water_m / 1000

# Check
head(shp_con[, c("adm2_name", "dist_to_water_km")])
summary(shp_con$dist_to_water_km)

dist_water_df <- shp_con %>%
  st_drop_geometry() %>%
  dplyr::select(adm2_name,
                dist_to_water_m,
                dist_to_water_km) %>%
  rename(District = adm2_name)

write.csv(
  dist_water_df,
  "cong_data/cong_distance_to_water.csv")


# Vegetation index

## 1: Benin 

ndvi_path <- "ben_data/NDVI2021/"

# Load all 12 monthly NDVI tif files
ndvi_files <- list.files(
  ndvi_path,
  pattern    = "_1_km_monthly_NDVI.*\\.tif$",
  full.names = TRUE
)

# Stack all 12 months into one raster
ndvi_stack <- rast(ndvi_files)

# # MODIS NDVI scale factor = 0.0001
# # Values below -2000 before scaling are fill/invalid values
# 
# # Apply scale factor
# ndvi_stack_scaled <- ndvi_stack * 0.0001
# 
# # Remove invalid values
# ndvi_stack_scaled[ndvi_stack < -0.2] <- NA
# ndvi_stack_scaled[ndvi_stack >  1.0] <- NA

ndvi_annual <- mean(ndvi_stack, na.rm = TRUE)

# Reproject raster to match districts if needed
ndvi_annual <- project(ndvi_annual, "EPSG:4326")

ndvi_extracted <- extract(
  ndvi_annual,
  ben_vect,
  fun   = "mean",
  na.rm = TRUE
)

shp_ben$NDVI <- ndvi_extracted[, 2]

vegetation_index_ben <- shp_ben %>%
  st_drop_geometry() %>%          
  dplyr::select(                       
    NAME_3,                       
    NDVI                
  ) %>%
  rename(
    District = NAME_3)

write_csv(vegetation_index_ben, "ben_data/ben_vegetation_index_2021.csv")

## 2: Cameroon (2011 and 2021)

# 2011
ndvi_path <- "cmr_data/NDVI2011/"

# Load all 12 monthly NDVI tif files
ndvi_files <- list.files(
  ndvi_path,
  pattern    = "_1_km_monthly_NDVI.*\\.tif$",
  full.names = TRUE
)

# Stack all 12 months into one raster
ndvi_stack <- rast(ndvi_files)

ndvi_annual <- mean(ndvi_stack, na.rm = TRUE)

# Reproject raster to match districts if needed
ndvi_annual <- project(ndvi_annual, "EPSG:4326")

ndvi_extracted <- extract(
  ndvi_annual,
  cmr_vect,
  fun   = "mean",
  na.rm = TRUE
)

shp_cmr$NDVI <- ndvi_extracted[, 2]

vegetation_index_cmr <- shp_cmr %>%
  st_drop_geometry() %>%          
  dplyr::select(                       
    Admin2,                       
    NDVI                
  ) %>%
  rename(
    District = Admin2)

write_csv(vegetation_index_cmr, "cmr_data/cmr_vegetation_index_2011.csv")


# 2021
ndvi_path <- "cmr_data/NDVI2021/"

# Load all 12 monthly NDVI tif files
ndvi_files <- list.files(
  ndvi_path,
  pattern    = "_1_km_monthly_NDVI.*\\.tif$",
  full.names = TRUE
)

# Stack all 12 months into one raster
ndvi_stack <- rast(ndvi_files)

ndvi_annual <- mean(ndvi_stack, na.rm = TRUE)

# Reproject raster to match districts if needed
ndvi_annual <- project(ndvi_annual, "EPSG:4326")

ndvi_extracted <- extract(
  ndvi_annual,
  cmr_vect,
  fun   = "mean",
  na.rm = TRUE
)

shp_cmr$NDVI_1 <- ndvi_extracted[, 2]

vegetation_index_cmr <- shp_cmr %>%
  st_drop_geometry() %>%          
  dplyr::select(                       
    Admin2,                       
    NDVI_1                
  ) %>%
  rename(
    District = Admin2)

write_csv(vegetation_index_cmr, "cmr_data/cmr_vegetation_index_2021.csv")

# 3: Congo

ndvi_path <- "cong_data/NDVI2021/"

# Load all 12 monthly NDVI tif files
ndvi_files <- list.files(
  ndvi_path,
  pattern    = "_1_km_monthly_NDVI.*\\.tif$",
  full.names = TRUE
)

# Stack all 12 months into one raster
ndvi_stack <- rast(ndvi_files)

ndvi_annual <- mean(ndvi_stack, na.rm = TRUE)

# Reproject raster to match districts if needed
ndvi_annual <- project(ndvi_annual, "EPSG:4326")

ndvi_extracted <- extract(
  ndvi_annual,
  con_vect,
  fun   = "mean",
  na.rm = TRUE
)

vegetation_index_cong <- shp_con %>%
  mutate(NDVI_1 = ndvi_extracted[, 2]) %>% 
  st_drop_geometry() %>%          
  dplyr::select(adm2_name, NDVI_1) %>%
  rename(District = adm2_name)

write_csv(vegetation_index_cong, "cong_data/cong_vegetation_index_2021.csv")


temp_df <- as.data.frame(ben_temperature_2021, xy = TRUE, na.rm = FALSE)

p <- ggplot() +
  geom_raster(
    data = temp_df,
    aes(x = x, y = y, fill = `ERA5_LAND_DAILY_temperature_2m_2021-01-01_2021-12-31_1p1068to13p6033N_8p1637to16p8844E.temperature_2m`)
  ) +
  geom_sf(
    data = shp_ben,
    fill = NA,
    color = "red",
    size = 0.6
  ) +
  scale_fill_viridis_c(option = "magma", na.value = "white") +
  coord_sf() +
  theme_minimal() +
  labs(
    title = "ERA5-Land Temperature (2021) with Benin District Boundaries",
    fill  = "Temp (°C)"
  )

ggsave("ben_data/verification_ben.png", plot = p, width = 15, height = 10)


temp_df$na_flag <- is.na(temp_df[[3]])

library(stars)
temp_stars <- st_as_stars(ben_temperature_2021)

# Plot interactive map
mapview(temp_stars, col.regions = viridis::viridis(100), na.color = "white") +
  mapview(shp_ben, col.regions = NA, color = "red", lwd = 2)



################################################################################
#   Geocodes
################################################################################


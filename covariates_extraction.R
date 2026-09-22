################################################################################
#        Covariates extraction at sampled and prediction locations             #
################################################################################

# Load necessary libraries
library(tidyverse)
library(sf)
library(terra)

################################################################################
#                              Useful functions                                #
################################################################################

## Function: fill NA values in one column using nearest non‑NA neighbor
fill_by_nearest <- function(sf_obj, colname) {
  
  # rows with non‑missing values
  non_na <- sf_obj %>% filter(!is.na(.data[[colname]]))
  
  # rows with missing values
  na_rows <- sf_obj %>% filter(is.na(.data[[colname]]))
  
  # if no NA: return unchanged
  if (nrow(na_rows) == 0) return(sf_obj)
  
  # find nearest non‑NA feature for each NA row
  nearest_id <- st_nearest_feature(na_rows, non_na)
  
  # extract the corresponding values
  replacement_values <- non_na[[colname]][nearest_id]
  
  # fill the NA values
  sf_obj[[colname]][is.na(sf_obj[[colname]])] <- replacement_values
  
  sf_obj
}

################################################################################
#                                Loading the data                              #
################################################################################

cluster_data_cmr_sf  <- st_read("data/cmr_data/spatial_data_sample_cmr_cluster.gpkg")
CMR_df_1_sf          <- st_read("data/cmr_data/spatial_data_sample_cmr.gpkg")

GAB_df_1_sf          <- st_read("data/gab_data/spatial_data_sample_gab.gpkg")

spatial_data_ben     <- st_read("data/ben_data/spatial_data_sample_ben.gpkg")

spatial_data_cong    <- st_read("data/cong_data/spatial_data_sample_cong.gpkg")


########################## Loading shape files #################################


shp_ben <- read_sf("data/ben_data/BEN_shape_file/gadm41_BEN_3.shp")
shp_con <- read_sf("data/cong_data/CON_shape_file/cog_admin2.shp")
shp_cmr <- read_sf("data/cmr_data/CMR_shape_file/cameroun_districts_2024.shp")
shp_gab <- read_sf("data/gab_data/gab_shape/gadm41_GAB_2.shp")


######################### renaming the districts ###############################

shp_ben$NAME_3[shp_ben$NAME_3=="Gakpè"]="Gakpe"
shp_ben$NAME_3[shp_ben$NAME_3=="Tokpa Domé"]="Tokpa-dome"
shp_ben$NAME_3[shp_ben$NAME_3=="Ouakpé-Daho"]="Houakpe-daho"
shp_ben$NAME_3[shp_ben$NAME_3=="Aganmalomé"]="Aganmanlome"
shp_ben$NAME_3[shp_ben$NAME_3=="Agonkanmè"]="Agonkamey"

ben_vect <- vect(shp_ben)
cmr_vect <- vect(shp_cmr)
con_vect <- vect(shp_con)
gab_vect <- vect(shp_gab)

########################### Loading raster files ###############################

land_cover_water <- rast("data/ben_data/landuse/WorldCover_water_30s.tif")
land_cover_trees <- rast("data/ben_data/landuse/WorldCover_trees_30s.tif")
land_cover_builts <- rast("data/ben_data/landuse/WorldCover_built_30s.tif")
land_cover_cropland <- rast("data/ben_data/landuse/WorldCover_cropland_30s.tif")
land_cover_wetland <- rast("data/ben_data/landuse/WorldCover_wetland_30s.tif")

global_spatial_limit <- rast("data/Africa_data/spatial_limit_of_transmission/2010_Pf_Limits_Decompressed.geotiff")

precip_ratser_2021 <- list.files("data/CHIPRS 2021/", pattern = "chirps", full.names = TRUE)
annual_temperature_climate <- rast("data/Africa_data/ERA5_LAND_DAILY_temperature_2m_2021-01-01_2021-12-31_-39p5213to41p4698N_-22p2966to54p0315E.temperature_2m.tif")

path1 <- "data/CHIRPS_PENTAD_precipitation_2021-01-01_2021-12-31_0p9311to13p9448N_8p2516to17p0602E.precipitation.tif"

rainfall_2021 <- rast(path1)


##################################### Benin case ###############################

ndvi_path <- "data/ben_data/NDVI2021/"
# Load all 12 monthly NDVI tif files
ndvi_files <- list.files(
  ndvi_path,
  pattern    = "_1_km_monthly_NDVI.*\\.tif$",
  full.names = TRUE
)



ITN_intervention_2021 <- rast("data/ITN_2000/2025_GBD2024_Africa_ITN_2021.tif")

ben_pop <- rast("data/ben_data/ben_pd_2020_1km_UNadj.tif")

elevation_ben <- rast("data/ben_data/elevation/BEN_elv_msk.tif")
ben_ndvi_stack<- rast(ndvi_files)

ben_land_cover_water <- crop(land_cover_water, ben_vect)
ben_land_cover_water <- mask(ben_land_cover_water, ben_vect)

ben_land_cover_built <- crop(land_cover_builts, ben_vect)
ben_land_cover_built <- mask(ben_land_cover_built, ben_vect)

ben_land_cover_cropland <- crop(land_cover_cropland, ben_vect)
ben_land_cover_cropland <- mask(ben_land_cover_cropland, ben_vect)

ben_land_cover_trees <- crop(land_cover_trees, ben_vect)
ben_land_cover_trees <- mask(ben_land_cover_trees, ben_vect)

ben_land_cover_wetland <- crop(land_cover_wetland, ben_vect)
ben_land_cover_wetland <- mask(ben_land_cover_wetland, ben_vect)

ben_rainfall_2021 <- crop(rainfall_2021, ben_vect)
ben_rainfall_2021 <- mask(ben_rainfall_2021, ben_vect)

ben_temperature_2021 <- crop(annual_temperature_climate , ben_vect)
ben_temperature_2021 <- mask(ben_temperature_2021, ben_vect)

ben_ITN_intervention_2021  <- crop(ITN_intervention_2021, ben_vect)
ben_ITN_intervention_2021  <- mask(ben_ITN_intervention_2021, ben_vect)

# # create a circle around the GPS coordinates extracted manually to ensure of 
# # reasonable extracted value for the population density
# 
# village_coords_ben_1 <- village_coords_ben %>% select(Village, lat, lng)
# 
# village_coords_ben_1 <- st_as_sf(village_coords_ben_1, coords = c("lng", "lat"),
#                                  crs = 4326)
# # project to a metric system
# village_coords_ben_1_projected <- st_transform(village_coords_ben_1, crs = 32631)
# 
# # create 1000 buffer around each village point
# # the distance can be adjusted based on the result obtained
# village_coords_ben_1_projected_buffer <- st_buffer(village_coords_ben_1_projected,
#                                                    dist = 1000)
# 
# # project back to match the raster's crs if they differ
# 
# village_coords_ben_1_projected_buffer_crs <- st_transform(village_coords_ben_1_projected_buffer,
#                                                           crs = crs(ben_pop))
# 
# ben_pop_density_2020_df_buffer <- terra::extract(ben_pop, 
#                                                  village_coords_ben_1_projected_buffer_crs,
#                                                  fun = mean,
#                                                  na.rm = TRUE)
# 

ben_ITN_intervention_2021_df <- terra::extract(ben_ITN_intervention_2021,  spatial_data_ben)
ben_pop_density_2020_df <- terra::extract(ben_pop,  spatial_data_ben)
ben_rainfall_2021_df    <- terra::extract(ben_rainfall_2021, spatial_data_ben)
ben_temperature_2021_df <- terra::extract(ben_temperature_2021, spatial_data_ben)
ndvi_vals_ben_monthly <- terra::extract(ben_ndvi_stack, spatial_data_ben)
ndvi_mean_ben <- apply(ndvi_vals_ben_monthly[,-1], 1, mean, na.rm = TRUE)
elevation_ben_df <- terra::extract(elevation_ben, spatial_data_ben)
ben_land_cover_water_df <- terra::extract(ben_land_cover_water, spatial_data_ben)
ben_land_cover_built_df <- terra::extract(ben_land_cover_built, spatial_data_ben)
ben_land_cover_cropland_df <- terra::extract(ben_land_cover_cropland, spatial_data_ben)
ben_land_cover_trees_df <- terra::extract(ben_land_cover_trees, spatial_data_ben)
ben_land_cover_wetland_df <- terra::extract(ben_land_cover_wetland, spatial_data_ben)



spatial_data_ben_sf <- spatial_data_ben %>% 
  mutate(pop_density_2020 = ben_pop_density_2020_df[,2],
         annual_rainfall_2021 = ben_rainfall_2021_df[,2],
         annual_temperature_2021 = ben_temperature_2021_df[,2],
         vegetation_index = ndvi_mean_ben,
         elevation = elevation_ben_df[,2],
         lc_water = ben_land_cover_water_df[,2],
         lc_built = ben_land_cover_built_df[,2],
         lc_trees = ben_land_cover_trees_df[,2],
         lc_wetland = ben_land_cover_wetland_df[,2],
         lc_cropland = ben_land_cover_cropland_df[,2],
         annual_ITN_usage_2021 = ben_ITN_intervention_2021_df[,2])


# distance to water

bbox <- st_bbox(shp_ben)


rivers <- st_read("data/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp",
                  wkt_filter = st_as_text(st_as_sfc(bbox)))

lakes  <- st_read("data/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp",
                  wkt_filter = st_as_text(st_as_sfc(bbox)))

rivers_valid <- st_make_valid(rivers)
lakes_valid  <- st_make_valid(lakes)

rivers_valid <- ms_simplify(rivers_valid, keep = 0.3, keep_shapes = TRUE)
lakes_valid  <- ms_simplify(lakes_valid,  keep = 0.3, keep_shapes = TRUE)

rivers_valid <- st_transform(rivers_valid, st_crs(shp_ben))
lakes_valid  <- st_transform(lakes_valid,  st_crs(shp_ben))

# nearest river segment for each point
nn_riv <- st_nn(spatial_data_ben_sf, rivers_valid, k = 1, returnDist = TRUE)
dist_riv_min <- sapply(nn_riv$dist, `[`, 1)

# nearest lake polygon for each point
nn_lake <- st_nn(spatial_data_ben_sf, lakes_valid, k = 1, returnDist = TRUE)
dist_lake_min <- sapply(nn_lake$dist, `[`, 1)

spatial_data_ben_sf$dist_to_water <- pmin(dist_riv_min, dist_lake_min)

numeric_cols <- names(spatial_data_ben_sf)[sapply(spatial_data_ben_sf, is.numeric)]

for (col in numeric_cols) {
  spatial_data_ben_sf <- fill_by_nearest(spatial_data_ben_sf, col)
}


st_write(spatial_data_ben_sf, "data/ben_data/spatial_data_ben_sf.gpkg", append = FALSE)


############################# Congo case #######################################

Cong_pop <- rast("data/cong_data/cog_pd_2020_1km_UNadj.tif")

elevation_Cong <- rast("data/cong_data/elevation/COG_elv_msk.tif")
ndvi_path <- "data/cong_data/NDVI2021/"
# Load all 12 monthly NDVI tif files
ndvi_files <- list.files(
  ndvi_path,
  pattern    = "_1_km_monthly_NDVI.*\\.tif$",
  full.names = TRUE
)

Cong_ndvi_stack <- rast(ndvi_files)


Cong_land_cover_water <- crop(land_cover_water, con_vect)
Cong_land_cover_water <- mask(Cong_land_cover_water, con_vect)

Cong_land_cover_built <- crop(land_cover_builts, con_vect)
Cong_land_cover_built <- mask(Cong_land_cover_built, con_vect)

Cong_land_cover_cropland <- crop(land_cover_cropland, con_vect)
Cong_land_cover_cropland <- mask(Cong_land_cover_cropland, con_vect)

Cong_land_cover_trees <- crop(land_cover_trees, con_vect)
Cong_land_cover_trees <- mask(Cong_land_cover_trees, con_vect)

Cong_land_cover_wetland <- crop(land_cover_wetland, con_vect)
Cong_land_cover_wetland <- mask(Cong_land_cover_wetland, con_vect)


Cong_rainfall_2021 <- crop(rainfall_2021, con_vect)
Cong_rainfall_2021 <- mask(Cong_rainfall_2021, con_vect)

Cong_temperature_2021 <- crop(annual_temperature_climate , con_vect)
Cong_temperature_2021 <- mask(Cong_temperature_2021, con_vect)


Cong_ITN_intervention_2021  <- crop(ITN_intervention_2021,  con_vect)
Cong_ITN_intervention_2021  <- mask(ITN_intervention_2021,  con_vect)

Cong_ITN_intervention_2021_df <- terra::extract(Cong_ITN_intervention_2021,  spatial_data_cong)
Cong_pop_density_2020_df <- terra::extract(Cong_pop,  spatial_data_cong)
Cong_rainfall_2021_df    <- terra::extract(Cong_rainfall_2021, spatial_data_cong)
Cong_temperature_2021_df <- terra::extract(Cong_temperature_2021, spatial_data_cong)
ndvi_vals_Cong_monthly <- terra::extract(Cong_ndvi_stack, spatial_data_cong)
ndvi_mean_Cong <- apply(ndvi_vals_Cong_monthly[,-1], 1, mean, na.rm = TRUE)
elevation_Cong_df <- terra::extract(elevation_Cong, spatial_data_cong)
Cong_land_cover_water_df <- terra::extract(Cong_land_cover_water, spatial_data_cong)
Cong_land_cover_built_df <- terra::extract(Cong_land_cover_built, spatial_data_cong)
Cong_land_cover_cropland_df <- terra::extract(Cong_land_cover_cropland, spatial_data_cong)
Cong_land_cover_trees_df <- terra::extract(Cong_land_cover_trees, spatial_data_cong)
Cong_land_cover_wetland_df <- terra::extract(Cong_land_cover_wetland, spatial_data_cong)

shp_cong_country <- shp_con %>% 
  group_by(adm0_name) %>% 
  summarise(geometry = st_union(geometry))


bbox <- st_bbox(shp_con)


rivers <- st_read("data/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp",
                  wkt_filter = st_as_text(st_as_sfc(bbox)))

lakes  <- st_read("data/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp",
                  wkt_filter = st_as_text(st_as_sfc(bbox)))

rivers_valid <- st_make_valid(rivers)
lakes_valid  <- st_make_valid(lakes)

rivers_valid <- ms_simplify(rivers_valid, keep = 0.3, keep_shapes = TRUE)
lakes_valid  <- ms_simplify(lakes_valid,  keep = 0.3, keep_shapes = TRUE)

rivers_valid <- st_transform(rivers_valid, st_crs(shp_con))
lakes_valid  <- st_transform(lakes_valid,  st_crs(shp_con))


# nearest river segment for each point
nn_riv <- st_nn(spatial_data_cong, rivers_valid, k = 1, returnDist = TRUE)
dist_riv_min <- sapply(nn_riv$dist, `[`, 1)

# nearest lake polygon for each point
nn_lake <- st_nn(spatial_data_cong, lakes_valid, k = 1, returnDist = TRUE)
dist_lake_min <- sapply(nn_lake$dist, `[`, 1)

spatial_data_cong$dist_to_water <- pmin(dist_riv_min, dist_lake_min)

spatial_data_cong <- spatial_data_cong %>% 
  mutate(pop_density_2020 = Cong_pop_density_2020_df[,2],
         annual_rainfall_2021 = Cong_rainfall_2021_df[,2],
         annual_temperature_2021 = Cong_temperature_2021_df[,2],
         vegetation_index = ndvi_mean_Cong,
         elevation = elevation_Cong_df[,2],
         lc_water = Cong_land_cover_water_df[,2],
         lc_built = Cong_land_cover_built_df[,2],
         lc_trees = Cong_land_cover_trees_df[,2],
         lc_wetland = Cong_land_cover_wetland_df[,2],
         lc_cropland = Cong_land_cover_cropland_df[,2],
         annual_ITN_usage_2021 = Cong_ITN_intervention_2021_df[,2])

spatial_data_cong_df <- spatial_data_cong %>% 
  st_drop_geometry() 

st_write(spatial_data_cong, "data/cong_data/spatial_data_cong_sf.gpkg", append = FALSE)


############################ Gabon case ########################################

ndvi_path <- "data/gab_data/NDVI2021/"
# Load all 12 monthly NDVI tif files
ndvi_files <- list.files(
  ndvi_path,
  pattern    = "_1_km_monthly_NDVI.*\\.tif$",
  full.names = TRUE
)


ITN_intervention_2021 <- rast("data/ITN_2000/2025_GBD2024_Africa_ITN_2021.tif")

gab_pop <- rast("data/gab_data/gab_pd_2020_1km_UNadj.tif")

elevation_gab <- rast("data/gab_data/elevation/gab_elv_msk.tif")
gab_ndvi_stack<- rast(ndvi_files)

gab_land_cover_water <- crop(land_cover_water, gab_vect)
gab_land_cover_water <- mask(gab_land_cover_water, gab_vect)

gab_land_cover_built <- crop(land_cover_builts, gab_vect)
gab_land_cover_built <- mask(gab_land_cover_built, gab_vect)

gab_land_cover_cropland <- crop(land_cover_cropland, gab_vect)
gab_land_cover_cropland <- mask(gab_land_cover_cropland, gab_vect)

gab_land_cover_trees <- crop(land_cover_trees, gab_vect)
gab_land_cover_trees <- mask(gab_land_cover_trees, gab_vect)

gab_land_cover_wetland <- crop(land_cover_wetland, gab_vect)
gab_land_cover_wetland <- mask(gab_land_cover_wetland, gab_vect)

gab_rainfall_2021 <- crop(rainfall_2021, gab_vect)
gab_rainfall_2021 <- mask(gab_rainfall_2021, gab_vect)

gab_temperature_2021 <- crop(annual_temperature_climate , gab_vect)
gab_temperature_2021 <- mask(gab_temperature_2021, gab_vect)


gab_ITN_intervention_2021  <- crop(ITN_intervention_2021, gab_vect)
gab_ITN_intervention_2021  <- mask(gab_ITN_intervention_2021, gab_vect)


gab_ITN_intervention_2021_df <- terra::extract(gab_ITN_intervention_2021,  GAB_df_1_sf)
gab_pop_density_2020_df <- terra::extract(gab_pop,  GAB_df_1_sf)
gab_rainfall_2021_df    <- terra::extract(gab_rainfall_2021, GAB_df_1_sf)
gab_temperature_2021_df <- terra::extract(gab_temperature_2021, GAB_df_1_sf)
ndvi_vals_gab_monthly <- terra::extract(gab_ndvi_stack, GAB_df_1_sf)
ndvi_mean_gab <- apply(ndvi_vals_gab_monthly[,-1], 1, mean, na.rm = TRUE)
elevation_gab_df <- terra::extract(elevation_gab, GAB_df_1_sf)
gab_land_cover_water_df <- terra::extract(gab_land_cover_water, GAB_df_1_sf)
gab_land_cover_built_df <- terra::extract(gab_land_cover_built, GAB_df_1_sf)
gab_land_cover_cropland_df <- terra::extract(gab_land_cover_cropland, GAB_df_1_sf)
gab_land_cover_trees_df <- terra::extract(gab_land_cover_trees, GAB_df_1_sf)
gab_land_cover_wetland_df <- terra::extract(gab_land_cover_wetland, GAB_df_1_sf)

GAB_df_1_sf <- GAB_df_1_sf %>% 
  mutate(pop_density_2020 = gab_pop_density_2020_df[,2],
         annual_rainfall_2021 = gab_rainfall_2021_df[,2],
         annual_temperature_2021 = gab_temperature_2021_df[,2],
         vegetation_index = ndvi_mean_gab,
         elevation = elevation_gab_df[,2],
         lc_water = gab_land_cover_water_df[,2],
         lc_built = gab_land_cover_built_df[,2],
         lc_trees = gab_land_cover_trees_df[,2],
         lc_wetland = gab_land_cover_wetland_df[,2],
         lc_cropland = gab_land_cover_cropland_df[,2],
         annual_ITN_usage_2021 = gab_ITN_intervention_2021_df[,2])


# distance to water

bbox <- st_bbox(shp_gab)


rivers <- st_read("data/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp",
                  wkt_filter = st_as_text(st_as_sfc(bbox)))

lakes  <- st_read("data/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp",
                  wkt_filter = st_as_text(st_as_sfc(bbox)))

rivers_valid <- st_make_valid(rivers)
lakes_valid  <- st_make_valid(lakes)

rivers_valid <- ms_simplify(rivers_valid, keep = 0.3, keep_shapes = TRUE)
lakes_valid  <- ms_simplify(lakes_valid,  keep = 0.3, keep_shapes = TRUE)

rivers_valid <- st_transform(rivers_valid, st_crs(shp_gab))
lakes_valid  <- st_transform(lakes_valid,  st_crs(shp_gab))

# nearest river segment for each point
nn_riv <- st_nn(GAB_df_1_sf_sf, rivers_valid, k = 1, returnDist = TRUE)
dist_riv_min <- sapply(nn_riv$dist, `[`, 1)

# nearest lake polygon for each point
nn_lake <- st_nn(GAB_df_1_sf_sf, lakes_valid, k = 1, returnDist = TRUE)
dist_lake_min <- sapply(nn_lake$dist, `[`, 1)

GAB_df_1_sf$dist_to_water <- pmin(dist_riv_min, dist_lake_min)

st_write(GAB_df_1_sf, "data/gab_data/GAB_df_1_sf.gpkg", append = FALSE)


############################## Cameroon case ###################################

cmr_pop <- rast("data/cmr_data/cmr_pd_2011_1km_UNadj.tif")

elevation_cmr <- rast("data/cmr_data/elevation/CMR_elv_msk.tif")
ndvi_path <- "data/cmr_data/NDVI2021/"
# Load all 12 monthly NDVI tif files
ndvi_files <- list.files(
  ndvi_path,
  pattern    = "_1_km_monthly_NDVI.*\\.tif$",
  full.names = TRUE
)

cmr_ndvi_stack <- rast(ndvi_files)

cmr_land_cover_water <- crop(land_cover_water, cmr_vect)
cmr_land_cover_water <- mask(cmr_land_cover_water, cmr_vect)

cmr_land_cover_built <- crop(land_cover_builts, cmr_vect)
cmr_land_cover_built <- mask(cmr_land_cover_built, cmr_vect)

cmr_land_cover_cropland <- crop(land_cover_cropland, cmr_vect)
cmr_land_cover_cropland <- mask(cmr_land_cover_cropland, cmr_vect)

cmr_land_cover_trees <- crop(land_cover_trees, cmr_vect)
cmr_land_cover_trees <- mask(cmr_land_cover_trees, cmr_vect)

cmr_land_cover_wetland <- crop(land_cover_wetland, cmr_vect)
cmr_land_cover_wetland <- mask(cmr_land_cover_wetland, cmr_vect)

cmr_rainfall_2021 <- crop(rainfall_2021, cmr_vect)
cmr_rainfall_2021 <- mask(cmr_rainfall_2021, cmr_vect)

cmr_temperature_2021 <- crop(annual_temperature_climate , cmr_vect)
cmr_temperature_2021 <- mask(cmr_temperature_2021, cmr_vect)

cmr_ITN_intervention_2021  <- crop(ITN_intervention_2021,  cmr_vect)
cmr_ITN_intervention_2021  <- mask(ITN_intervention_2021,  cmr_vect)

cmr_ITN_intervention_2021_df <- terra::extract(cmr_ITN_intervention_2021,  CMR_df_1_sf)
cmr_pop_density_2020_df <- terra::extract(cmr_pop,  CMR_df_1_sf)
cmr_rainfall_2021_df    <- terra::extract(cmr_rainfall_2021, CMR_df_1_sf)
cmr_temperature_2021_df <- terra::extract(cmr_temperature_2021, CMR_df_1_sf)
ndvi_vals_cmr_monthly <- terra::extract(cmr_ndvi_stack, CMR_df_1_sf)
ndvi_mean_cmr <- apply(ndvi_vals_cmr_monthly[,-1], 1, mean, na.rm = TRUE)
elevation_cmr_df <- terra::extract(elevation_cmr, CMR_df_1_sf)
cmr_land_cover_water_df <- terra::extract(cmr_land_cover_water, CMR_df_1_sf)
cmr_land_cover_built_df <- terra::extract(cmr_land_cover_built, CMR_df_1_sf)
cmr_land_cover_cropland_df <- terra::extract(cmr_land_cover_cropland, CMR_df_1_sf)
cmr_land_cover_trees_df <- terra::extract(cmr_land_cover_trees, CMR_df_1_sf)
cmr_land_cover_wetland_df <- terra::extract(cmr_land_cover_wetland, CMR_df_1_sf)

shp_cmr_country <- shp_cmr %>% 
  mutate(country = "Cameroon")%>% 
  group_by(country) %>% 
  summarise(geometry = st_union(geometry))


bbox <- st_bbox(shp_cmr)


rivers <- st_read("data/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp",
                  wkt_filter = st_as_text(st_as_sfc(bbox)))

lakes  <- st_read("data/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp",
                  wkt_filter = st_as_text(st_as_sfc(bbox)))

rivers_valid <- st_make_valid(rivers)
lakes_valid  <- st_make_valid(lakes)

rivers_valid <- ms_simplify(rivers_valid, keep = 0.3, keep_shapes = TRUE)
lakes_valid  <- ms_simplify(lakes_valid,  keep = 0.3, keep_shapes = TRUE)


rivers_valid <- st_transform(rivers_valid, st_crs(shp_cmr))
lakes_valid  <- st_transform(lakes_valid,  st_crs(shp_cmr))


# nearest river segment for each point
nn_riv <- st_nn(CMR_df_1_sf, rivers_valid, k = 1, returnDist = TRUE)
dist_riv_min <- sapply(nn_riv$dist, `[`, 1)

# nearest lake polygon for each point
nn_lake <- st_nn(CMR_df_1_sf, lakes_valid, k = 1, returnDist = TRUE)
dist_lake_min <- sapply(nn_lake$dist, `[`, 1)

CMR_df_1_sf$dist_to_water <- pmin(dist_riv_min, dist_lake_min)

CMR_df_1_sf <- CMR_df_1_sf %>% 
  mutate(pop_density_2020 = cmr_pop_density_2020_df[,2],
         annual_rainfall_2021 = cmr_rainfall_2021_df[,2],
         annual_temperature_2021 = cmr_temperature_2021_df[,2],
         vegetation_index = ndvi_mean_cmr,
         elevation = elevation_cmr_df[,2],
         lc_water = cmr_land_cover_water_df[,2],
         lc_built = cmr_land_cover_built_df[,2],
         lc_trees = cmr_land_cover_trees_df[,2],
         lc_wetland = cmr_land_cover_wetland_df[,2],
         lc_cropland = cmr_land_cover_cropland_df[,2],
         annual_ITN_usage_2021 = cmr_ITN_intervention_2021_df[,2])

st_write(CMR_df_1_sf, "data/cmr_data/CMR_df_1_sf.gpkg", append = FALSE)


# Cameroon case with clusters of households

cmr_cluster_ITN_intervention_2021_df <- terra::extract(cmr_ITN_intervention_2021,
                                                       cluster_data_cmr_sf)
cmr_cluster_pop_density_2020_df <- terra::extract(cmr_pop,  cluster_data_cmr_sf)
cmr_cluster_rainfall_2021_df    <- terra::extract(cmr_rainfall_2021, cluster_data_cmr_sf)
cmr_cluster_temperature_2021_df <- terra::extract(cmr_temperature_2021, cluster_data_cmr_sf)
ndvi_vals_cluster_cmr_monthly <- terra::extract(cmr_ndvi_stack, cluster_data_cmr_sf)
ndvi_mean_cluster_cmr <- apply(ndvi_vals_cluster_cmr_monthly[,-1], 1, mean, na.rm = TRUE)
elevation_cluster_cmr_df <- terra::extract(elevation_cmr, cluster_data_cmr_sf)
cmr_cluster_land_cover_water_df <- terra::extract(cmr_land_cover_water, cluster_data_cmr_sf)
cmr_cluster_land_cover_built_df <- terra::extract(cmr_land_cover_built, cluster_data_cmr_sf)
cmr_cluster_land_cover_cropland_df <- terra::extract(cmr_land_cover_cropland, cluster_data_cmr_sf)
cmr_cluster_land_cover_trees_df <- terra::extract(cmr_land_cover_trees, cluster_data_cmr_sf)
cmr_cluster_land_cover_wetland_df <- terra::extract(cmr_land_cover_wetland, cluster_data_cmr_sf)


# nearest river segment for each point
nn_riv <- st_nn(cluster_data_cmr_sf, rivers_valid, k = 1, returnDist = TRUE)
dist_riv_min <- sapply(nn_riv$dist, `[`, 1)

# nearest lake polygon for each point
nn_lake <- st_nn(cluster_data_cmr_sf, lakes_valid, k = 1, returnDist = TRUE)
dist_lake_min <- sapply(nn_lake$dist, `[`, 1)

cluster_data_cmr_sf$dist_to_water <- pmin(dist_riv_min, dist_lake_min)

cluster_data_cmr_sf <- cluster_data_cmr_sf %>% 
  mutate(pop_density_2020 = cmr_cluster_pop_density_2020_df[,2],
         annual_rainfall_2021 = cmr_cluster_rainfall_2021_df[,2],
         annual_temperature_2021 = cmr_cluster_temperature_2021_df[,2],
         vegetation_index = ndvi_mean_cluster_cmr,
         elevation = elevation_cluster_cmr_df[,2],
         lc_water = cmr_cluster_land_cover_water_df[,2],
         lc_built = cmr_cluster_land_cover_built_df[,2],
         lc_trees = cmr_cluster_land_cover_trees_df[,2],
         lc_wetland = cmr_cluster_land_cover_wetland_df[,2],
         lc_cropland = cmr_cluster_land_cover_cropland_df[,2],
         annual_ITN_usage_2021 = cmr_cluster_ITN_intervention_2021_df[,2]) %>% 
  dplyr::select(-c(grid_id, lat, lng))


st_write(cluster_data_cmr_sf, "data/cmr_data/cmr_cluster_df_1_sf.gpkg", append = FALSE)

################################################################################
#                 case of grid cells at country level                          #
################################################################################
countries1_sf <- rbind(shp_ben %>% dplyr::select(COUNTRY, NAME_3, geometry) %>% 
                         rename(country = COUNTRY,
                                district = NAME_3), shp_con %>% dplyr::select(adm0_name,
                                                                              adm2_name,
                                                                              geometry) %>% 
                         rename(country = adm0_name,
                                district = adm2_name)) 


countries_sf <- rbind(countries1_sf, shp_cmr %>% dplyr::select(Admin2, geometry) %>% 
                        mutate(country = "Cameroon") %>% 
                        rename(district = Admin2))

countries_sf1 <- rbind(countries1_sf, shp_cmr %>% dplyr::select(Admin2, geometry) %>% 
                         mutate(country = "Cameroon") %>% 
                         rename(district = Admin2))

countries_sf1$country_id <- as.integer(factor(countries_sf1$country))


# 1. Reproject countries to a metric CRS
countries_m <- st_transform(countries_sf1, 3857)

# 2. Create bbox in the same CRS
bbox_m <- st_bbox(countries_m) %>% st_as_sfc()

# 3. Build 5 km grid (now 5000 means 5000 meters)
grid_5km <- st_make_grid(
  bbox_m,
  cellsize = 5000,
  what = "centers"
)

grid_sf <- st_sf(geometry = grid_5km)

# 4. Intersect
grid_sf <- st_intersection(grid_sf, countries_m)

grid_vect <- vect(countries_sf %>% st_transform(4326))

grid_pop <- raster::merge(ben_pop, Cong_pop, cmr_pop)
grid_elevation <- raster::merge(elevation_ben, elevation_cmr, elevation_Cong)
grid_NDVI <- raster::merge(ben_ndvi_stack, cmr_ndvi_stack, Cong_ndvi_stack)

grid_land_cover_water <- crop(land_cover_water, grid_vect)
grid_land_cover_water <- mask(grid_land_cover_water, grid_vect)

grid_land_cover_built <- crop(land_cover_builts, grid_vect)
grid_land_cover_built <- mask(grid_land_cover_built, grid_vect)

grid_land_cover_cropland <- crop(land_cover_cropland, grid_vect)
grid_land_cover_cropland <- mask(grid_land_cover_cropland, grid_vect)

grid_land_cover_trees <- crop(land_cover_trees, grid_vect)
grid_land_cover_trees <- mask(grid_land_cover_trees, grid_vect)

grid_land_cover_wetland <- crop(land_cover_wetland, grid_vect)
grid_land_cover_wetland <- mask(grid_land_cover_wetland, grid_vect)

grid_rainfall_2021 <- crop(rainfall_2021, grid_vect)
grid_rainfall_2021 <- mask(grid_rainfall_2021, grid_vect)

grid_temperature_2021 <- crop(annual_temperature_climate, grid_vect)
grid_temperature_2021 <- mask(grid_temperature_2021, grid_vect)


grid_ITN_usage_2021 <- crop(ITN_intervention_2021, grid_vect)
grid_ITN_usage_2021 <- mask(grid_ITN_usage_2021, grid_vect)

grid_ITN_usage_2021_df <- terra::extract(grid_ITN_usage_2021,  grid_sf%>% st_transform(4326))
grid_pop_density_2020_df <- terra::extract(grid_pop,  grid_sf%>% st_transform(4326))
grid_rainfall_2021_df    <- terra::extract(grid_rainfall_2021, grid_sf%>% st_transform(4326))
grid_temperature_2021_df <- terra::extract(grid_temperature_2021, grid_sf%>% st_transform(4326))
ndvi_vals_grid_monthly <- terra::extract(grid_NDVI, grid_sf%>% st_transform(4326))
ndvi_mean_grid <- apply(ndvi_vals_grid_monthly[,-1], 1, mean, na.rm = TRUE)
elevation_grid_df <- terra::extract(grid_elevation, grid_sf%>% st_transform(4326))
grid_land_cover_water_df <- terra::extract(grid_land_cover_water, grid_sf%>% st_transform(4326))
grid_land_cover_built_df <- terra::extract(grid_land_cover_built, grid_sf%>% st_transform(4326))
grid_land_cover_cropland_df <- terra::extract(grid_land_cover_cropland, grid_sf%>% st_transform(4326))
grid_land_cover_trees_df <- terra::extract(grid_land_cover_trees, grid_sf%>% st_transform(4326))
grid_land_cover_wetland_df <- terra::extract(grid_land_cover_wetland, grid_sf%>% st_transform(4326))

# Distance to water

bbox <- st_bbox(countries_sf)


rivers <- st_read("data/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp",
                  wkt_filter = st_as_text(st_as_sfc(bbox)))

lakes  <- st_read("data/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp",
                  wkt_filter = st_as_text(st_as_sfc(bbox)))

rivers_valid <- st_make_valid(rivers)
lakes_valid  <- st_make_valid(lakes)

rivers_valid <- ms_simplify(rivers_valid, keep = 0.3, keep_shapes = TRUE)
lakes_valid  <- ms_simplify(lakes_valid,  keep = 0.3, keep_shapes = TRUE)



# nearest river segment for each point
nn_riv <- st_nn(grid_sf %>% st_transform(4326), rivers_valid, k = 1, returnDist = TRUE)
dist_riv_min <- sapply(nn_riv$dist, `[`, 1)

# nearest lake polygon for each point
nn_lake <- st_nn(grid_sf%>% st_transform(4326), lakes_valid, k = 1, returnDist = TRUE)
dist_lake_min <- sapply(nn_lake$dist, `[`, 1)

grid_sf$dist_to_water <- pmin(dist_riv_min, dist_lake_min)

grid_sf <- grid_sf %>% 
  mutate(annual_rainfall_2021 = grid_rainfall_2021_df[,2],
         annual_temperature_2021 = grid_temperature_2021_df[,2],
         elevation               = elevation_grid_df[,2],
         vegetation_index = ndvi_mean_grid,
         pop_density_2020 = grid_pop_density_2020_df[,2],
         lc_cropland = grid_land_cover_cropland_df[,2],
         lc_water = grid_land_cover_water_df[,2],
         lc_wetland = grid_land_cover_wetland_df[,2],
         lc_built = grid_land_cover_built_df[,2],
         lc_trees = grid_land_cover_trees_df[,2],
         annual_ITN_usage_2021 = grid_ITN_usage_2021_df[,2]) 



numeric_cols <- names(grid_sf)[sapply(grid_sf, is.numeric)]

for (col in numeric_cols) {
  grid_sf <- fill_by_nearest(grid_sf, col)
}


st_write(grid_sf, "data/grid_data_for_prediction_sf.gpkg", append = FALSE)

################################################################################
#                 case of grid cells at continental level                      #
################################################################################

shp_africa <- st_read("data/Africa_data/shp_africa/Africa_simplified.shp") %>% st_transform(3857)

# Dissolve all polygons into one
shp_africa <- st_union(shp_africa)
shp_africa <- st_sf(geometry = shp_africa)

bbox <- st_bbox(shp_africa)

grid_5km_africa <- st_make_grid(
  st_as_sfc(bbox),
  cellsize = 5000,
  what = "centers"
)

grid_sf_africa <- st_sf(geometry = grid_5km_africa)
grid_sf_africa <- st_filter(grid_sf_africa, shp_africa)

ndvi_path <- "data/Africa_data/NDVI_africa/"
# Load all 12 monthly NDVI tif files
ndvi_files <- list.files(
  ndvi_path,
  pattern    = "_1_km_monthly_NDVI.*\\.tif$",
  full.names = TRUE
)

path1 <- "data/Africa_data/climateEngine_download.precipitation.tif"
path2 <- "data/Africa_data/ERA5_LAND_DAILY_temperature_2m_2021-01-01_2021-12-31_-39p5213to41p4698N_-22p2966to54p0315E.temperature_2m.tif"


ITN_intervention_2021 <- rast("data/ITN_2000/2025_GBD2024_Africa_ITN_2021.tif")
africa_rainfall_2021 <- rast(path1)
africa_temperature_2021 <- rast(path2)
africa_pop <- rast("data/Africa_data/pop_density/gpw_v4_population_density_rev11_2020_30_sec.tif")
africa_elevation <- rast("data/Africa_data/elevation/wc2.1_30s_elev.tif")
africa_ndvi_stack<- rast(ndvi_files)



africa_land_cover_water <- crop(land_cover_water, vect(shp_africa %>% st_transform(4326)))
africa_land_cover_water <- mask(africa_land_cover_water, vect(shp_africa %>% st_transform(4326)))

africa_land_cover_built <- crop(land_cover_builts, vect(shp_africa %>% st_transform(4326)))
africa_land_cover_built <- mask(africa_land_cover_built, vect(shp_africa %>% st_transform(4326)))

africa_land_cover_cropland <- crop(land_cover_cropland, vect(shp_africa %>% st_transform(4326)))
africa_land_cover_cropland <- mask(africa_land_cover_cropland, vect(shp_africa %>% st_transform(4326)))

africa_land_cover_trees <- crop(land_cover_trees, vect(shp_africa %>% st_transform(4326)))
africa_land_cover_trees <- mask(africa_land_cover_trees, vect(shp_africa %>% st_transform(4326)))

africa_land_cover_wetland <- crop(land_cover_wetland, vect(shp_africa %>% st_transform(4326)))
africa_land_cover_wetland <- mask(africa_land_cover_wetland, vect(shp_africa %>% st_transform(4326)))

africa_ITN_intervention_2021 <- crop(ITN_intervention_2021, vect(shp_africa %>% st_transform(4326)))
africa_ITN_intervention_2021 <- mask(africa_ITN_intervention_2021, vect(shp_africa %>% st_transform(4326)))

africa_pop_2020 <- crop(africa_pop, vect(shp_africa %>% st_transform(4326)))
africa_pop_2020 <- mask(africa_pop_2020, vect(shp_africa %>% st_transform(4326)))

africa_rainfall_2021 <- crop(africa_rainfall_2021, vect(shp_africa %>% st_transform(4326)))
africa_rainfall_2021 <- mask(africa_rainfall_2021, vect(shp_africa %>% st_transform(4326)))

africa_temperature_2021 <- crop(africa_temperature_2021, vect(shp_africa %>% st_transform(4326)))
africa_temperature_2021 <- mask(africa_temperature_2021, vect(shp_africa %>% st_transform(4326)))

africa_ndvi_stack <- crop(africa_ndvi_stack, vect(shp_africa %>% st_transform(4326)))
africa_ndvi_stack <- mask(africa_ndvi_stack, vect(shp_africa %>% st_transform(4326)))

africa_elevation <- crop(africa_elevation, vect(shp_africa %>% st_transform(4326)))
africa_elevation <- mask(africa_elevation, vect(shp_africa %>% st_transform(4326)))



# africa_ITN_intervention_2021  <- project(ITN_intervention_2021,  st_crs(grid_sf_africa)$wkt)
# africa_pop  <- project(africa_pop,  st_crs(grid_sf_africa)$wkt)
# africa_rainfall_2021 <- project(africa_rainfall_2021, st_crs(grid_sf_africa)$wkt)
# africa_temperature_2021 <- project(africa_temperature_2021, st_crs(grid_sf_africa)$wkt)
# africa_ndvi_stack <- project(africa_ndvi_stack, st_crs(grid_sf_africa)$wkt)
# africa_elevation   <- project(africa_elevation,   st_crs(grid_sf_africa)$wkt)
# africa_land_cover_water   <- project(africa_land_cover_water,   st_crs(grid_sf_africa)$wkt)
# africa_land_cover_built   <- project(africa_land_cover_built,   st_crs(grid_sf_africa)$wkt)
# africa_land_cover_cropland   <- project(africa_land_cover_cropland,   st_crs(grid_sf_africa)$wkt)
# africa_land_cover_trees   <- project(africa_land_cover_trees,   st_crs(grid_sf_africa)$wkt)
# africa_land_cover_wetland   <- project(africa_land_cover_wetland,   st_crs(grid_sf_africa)$wkt)

africa_ITN_intervention_2021_df <- terra::extract(africa_ITN_intervention_2021, 
                                                  grid_sf_africa %>% st_transform(4326))
africa_pop_density_2020_df <- terra::extract(africa_pop_2020,  
                                             grid_sf_africa %>% st_transform(4326))
africa_rainfall_2021_df    <- terra::extract(africa_rainfall_2021, 
                                             grid_sf_africa %>% st_transform(4326))
africa_temperature_2021_df <- terra::extract(africa_temperature_2021, 
                                             grid_sf_africa %>% st_transform(4326))

ndvi_vals_africa_monthly <- terra::extract(africa_ndvi_stack, 
                                           grid_sf_africa %>% st_transform(4326))
ndvi_mean_africa <- apply(ndvi_vals_africa_monthly[,-1], 1, mean, na.rm = TRUE)

africa_elevation_df <- terra::extract(africa_elevation,
                                      grid_sf_africa %>% st_transform(4326))
africa_land_cover_water_df <- terra::extract(africa_land_cover_water, 
                                             grid_sf_africa %>% st_transform(4326))
africa_land_cover_built_df <- terra::extract(africa_land_cover_built,
                                             grid_sf_africa %>% st_transform(4326))
africa_land_cover_cropland_df <- terra::extract(africa_land_cover_cropland,
                                                grid_sf_africa %>% st_transform(4326))
africa_land_cover_trees_df <- terra::extract(africa_land_cover_trees,
                                             grid_sf_africa %>% st_transform(4326))
africa_land_cover_wetland_df <- terra::extract(africa_land_cover_wetland,
                                               grid_sf_africa %>% st_transform(4326))

shp_africa_summarise <- shp_africa %>% 
  group_by(continent) %>% 
  summarise(geometry = st_union(geometry))

# distance to water

# bbox <- st_bbox(grid_sf_africa)

rivers <- vect("data/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp")
lakes <- vect("data/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp")


africa4326 <- st_transform(shp_africa, 4326)
africa4326 <- vect(africa4326)

lakes_africa <- crop(lakes, africa4326)

r_template <- rast(
  ext(vect(st_as_sf(shp_africa))),
  resolution = 5000,
  crs = "EPSG:3857")

rivers <- project(rivers, "EPSG:3857")
lakes_africa <- project(lakes_africa, "EPSG:3857")

rivers2 <- rivers[rivers$ORD_STRA >= 3, ]

riv_r <- rasterize(
  rivers2,
  r_template,
  field = 1,
  background = NA)


lake_r <- rasterize(
  lakes_africa,
  r_template,
  field = 1,
  background = NA)

water <- cover(riv_r, lake_r)

dist_water <- distance(
  water,
  filename = "data/dist_to_water_5km.tif",
  overwrite = TRUE)

dist_water <- rast("data/dist_to_water_5km.tif")

grid_sf_africa$dist_to_water <- terra::extract(
  dist_water,
  vect(grid_sf_africa))[,2]


grid_sf_africa <- grid_sf_africa %>% 
  mutate(pop_density_2020 = africa_pop_density_2020_df[,2],
         annual_rainfall_2021 = africa_rainfall_2021_df[,2],
         annual_temperature_2021 = africa_temperature_2021_df[,2],
         vegetation_index = ndvi_mean_africa,
         elevation = africa_elevation_df[,2],
         lc_water = africa_land_cover_water_df[,2],
         lc_built = africa_land_cover_built_df[,2],
         lc_trees = africa_land_cover_trees_df[,2],
         lc_wetland = africa_land_cover_wetland_df[,2],
         lc_cropland = africa_land_cover_cropland_df[,2],
         annual_ITN_usage_2021 = africa_ITN_intervention_2021_df[,2])

# filling missing values
numeric_cols <- names(grid_sf_africa)[sapply(grid_sf_africa, is.numeric)]

for (col in numeric_cols) {
  grid_sf_africa <- fill_by_nearest(grid_sf_africa, col)
}



st_write(grid_sf_africa, "data/Africa_data/grid_data_for_prediction_sf.gpkg", append = FALSE)

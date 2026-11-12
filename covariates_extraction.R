################################################################################
#        Covariates extraction at sampled and prediction locations             #
################################################################################

# Load necessary libraries
library(tidyverse)
library(sf)
library(terra)
library(lwgeom)
library(exactextractr)
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


# distance to water

# ------------------------------------------------------------------
# Idea: densify rivers and lake shorelines into points every
# `spacing_m` meters, then do a nearest-neighbor search with a KD-tree.
# Error <= spacing_m / 2 (250 m for 500 m spacing), memory is small.
# Points are converted to 3D coordinates on the globe, so distances
# are valid everywhere in Africa (no projection distortion).
# ------------------------------------------------------------------

spacing_m <- 500
min_order <- 1      # Strahler filter. 1 = all rivers (modeling choice!)
R <- 6371008.8      # mean Earth radius (m)

# put every location in a coordinate system where plain straight-line distance is correct everywhere
to_xyz <- function(lonlat) {
  lon <- lonlat[, 1] * pi / 180
  lat <- lonlat[, 2] * pi / 180
  cbind(R * cos(lat) * cos(lon), R * cos(lat) * sin(lon), R * sin(lat))
}

# correct for the fact that the KD-tree measures a straight 
# line through the Earth, not the path along its surface

chord_to_arc <- function(d) 2 * R * asin(pmin(d / (2 * R), 1))

# Densify geometries in chunks and return xyz points
# turn rivers and shorelines into many closely spaced 3D points, so that 
# "distance to a line" becomes "distance to the nearest point"

densify_xyz <- function(geom, spacing_m, chunk = 50000) {
  idx <- split(seq_along(geom), ceiling(seq_along(geom) / chunk))
  out <- lapply(idx, function(i) {
    g <- st_segmentize(geom[i], units::set_units(spacing_m, "m"))  # adds vertices so that no two are more than 500 m apart
    to_xyz(st_coordinates(g)[, 1:2, drop = FALSE])                 # extracts vertices
  })
  do.call(rbind, out)
}

# Nearest distance (m, great-circle) from query points to a point cloud
# for every query point (a centroid), find the distance to the closest point of the water cloud

nn_dist <- function(xyz_water, xyz_query, chunk = 100000) {
  idx <- split(seq_len(nrow(xyz_query)), ceiling(seq_len(nrow(xyz_query)) / chunk))
  d <- unlist(lapply(idx, function(i)
    nabor::knn(data = xyz_water,
               query = xyz_query[i, , drop = FALSE],
               k = 1)$nn.dists[, 1]))
  chord_to_arc(d)
}

zone_dist <- function(centres, crs_m, half = 2000, n = 5) {
  # centres: sf POINTS (centres of the 5x5 km zones), in the metric CRS crs_m
  xy   <- st_coordinates(centres)[, 1:2, drop = FALSE]
  off  <- seq(-half, half, length.out = n)          # -2000, -1000, 0, 1000, 2000 m
  offs <- expand.grid(dx = off, dy = off)           # 25 offsets
  k    <- nrow(offs)
  
  idx <- rep(seq_len(nrow(xy)), each = k)           # zone id of each sample point
  px  <- xy[idx, 1] + rep(offs$dx, times = nrow(xy))
  py  <- xy[idx, 2] + rep(offs$dy, times = nrow(xy))
  
  pp <- st_as_sf(data.frame(x = px, y = py), coords = c("x", "y"), crs = crs_m) |>
    st_transform(4326)
  xyz <- to_xyz(st_coordinates(pp)[, 1:2])
  
  d_river <- nn_dist(xyz_riv,  xyz)
  d_lake  <- nn_dist(xyz_lake, xyz)
  d       <- pmin(d_river, d_lake)
  
  # points inside a lake: distance 0
  d[lengths(st_intersects(pp, lakes_geom)) > 0] <- 0
  
  as.numeric(tapply(d, idx, mean))                  # mean over the 25 points of each zone
}

################################################################################
#                                Loading the data                              #
################################################################################

cluster_data_cmr_sf  <- st_read("data/cmr_data/spatial_data_sample_cmr_cluster.gpkg")
CMR_df_1_sf          <- st_read("data/cmr_data/spatial_data_sample_cmr.gpkg")

GAB_df_1_sf          <- st_read("data/gab_data/spatial_data_sample_gab.gpkg")

spatial_data_ben     <- st_read("data/ben_data/spatial_data_sample_ben.gpkg")

spatial_data_cong    <- st_read("data/cong_data/spatial_data_sample_cong.gpkg")

cmr_raw_pr_mix       <- read.csv("data/cmr_data/cmr_MAP_full_raw_prevalence_mix.csv")

cmr_raw_pr_mix_sf    <- st_as_sf(cmr_raw_pr_mix, 
                                 coords = c("longitude", "latitude"),
                                 crs = 4326)


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

crs_m <- 32631    # a projected CRS in metres for Gabon (UTM 32S); an equal-area CRS also works

# 1. Zones: 5x5 km squares for survey sites and for grid cells
pts_m   <- st_transform(spatial_data_ben, crs_m)
zones_s <- st_buffer(pts_m, dist = 2500, endCapStyle = "SQUARE")   # square centred on each site
zones_s <- zones_s |> st_transform(4326)
zones_v <- terra::vect(zones_s)

ben_ITN_intervention_2021_df <- terra::extract(ben_ITN_intervention_2021,  
                                               zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
ben_pop_density_2020_df <- terra::extract(ben_pop,  zones_v, fun = "mean", 
                                          na.rm = TRUE, weights = TRUE)
ben_rainfall_2021_df    <- terra::extract(ben_rainfall_2021, zones_v, fun = "mean", 
                                          na.rm = TRUE, weights = TRUE)
ben_temperature_2021_df <- terra::extract(ben_temperature_2021, zones_v, fun = "mean", 
                                          na.rm = TRUE, weights = TRUE)
ndvi_vals_ben_monthly <- terra::extract(ben_ndvi_stack, zones_v, fun = "mean", 
                                        na.rm = TRUE, weights = TRUE)
ndvi_mean_ben <- apply(ndvi_vals_ben_monthly[,-1], 1, mean, , na.rm = TRUE, weights = TRUE)
elevation_ben_df <- terra::extract(elevation_ben, zones_v, fun = "mean")
ben_land_cover_water_df <- terra::extract(ben_land_cover_water, zones_v, fun = "mean", 
                                          na.rm = TRUE, weights = TRUE)
ben_land_cover_built_df <- terra::extract(ben_land_cover_built, zones_v, fun = "mean", 
                                          na.rm = TRUE, weights = TRUE)
ben_land_cover_cropland_df <- terra::extract(ben_land_cover_cropland, zones_v, fun = "mean", 
                                             na.rm = TRUE, weights = TRUE)
ben_land_cover_trees_df <- terra::extract(ben_land_cover_trees, zones_v, fun = "mean", 
                                          na.rm = TRUE, weights = TRUE)
ben_land_cover_wetland_df <- terra::extract(ben_land_cover_wetland, zones_v, fun = "mean", 
                                            na.rm = TRUE, weights = TRUE)



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

# ------------------------------------------------------------------
# Data (all in EPSG:4326)
# ------------------------------------------------------------------
wkt <- st_as_text(st_as_sfc(st_bbox(shp_ben)))

rivers <- st_read("data/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp",
                  wkt_filter = wkt)
rivers <- rivers[rivers$ORD_STRA >= min_order, ] |> st_geometry()

lakes <- st_read("data/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp",
                 wkt_filter = wkt)
lakes_geom <- st_make_valid(st_geometry(lakes))
shores <- st_boundary(lakes_geom)             # turns each polygon into its outline (a line)
shores <- shores[!st_is_empty(shores)]
shores <- st_cast(shores, "MULTILINESTRING")  # makes all the outlines one geometry type

xyz_riv  <- densify_xyz(rivers, spacing_m)
xyz_lake <- densify_xyz(shores, spacing_m)


spatial_data_ben_sf$dist_to_water <- zone_dist(st_centroid(pts_m),  crs_m)


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

crs_m <- 6933    # a projected CRS in metres for Gabon (UTM 32S); an equal-area CRS also works

# 1. Zones: 5x5 km squares for survey sites and for grid cells
pts_m   <- st_transform(spatial_data_cong, crs_m)
zones_s <- st_buffer(pts_m, dist = 2500, endCapStyle = "SQUARE")   # square centred on each site
zones_s <- zones_s |> st_transform(4326)
zones_v <- terra::vect(zones_s)

Cong_ITN_intervention_2021_df <- terra::extract(Cong_ITN_intervention_2021,  zones_v, fun = "mean", 
                                                na.rm = TRUE, weights = TRUE)
Cong_pop_density_2020_df <- terra::extract(Cong_pop,  zones_v, fun = "mean", 
                                           na.rm = TRUE, weights = TRUE)
Cong_rainfall_2021_df    <- terra::extract(Cong_rainfall_2021, zones_v, fun = "mean", 
                                           na.rm = TRUE, weights = TRUE)
Cong_temperature_2021_df <- terra::extract(Cong_temperature_2021, zones_v, fun = "mean", 
                                           na.rm = TRUE, weights = TRUE)
ndvi_vals_Cong_monthly <- terra::extract(Cong_ndvi_stack, zones_v, fun = "mean", 
                                         na.rm = TRUE, weights = TRUE)
ndvi_mean_Cong <- apply(ndvi_vals_Cong_monthly[,-1], 1, mean, fun = "mean", 
                        na.rm = TRUE, weights = TRUE)
elevation_Cong_df <- terra::extract(elevation_Cong, zones_v, fun = "mean", 
                                    na.rm = TRUE, weights = TRUE)
Cong_land_cover_water_df <- terra::extract(Cong_land_cover_water, zones_v, fun = "mean", 
                                           na.rm = TRUE, weights = TRUE)
Cong_land_cover_built_df <- terra::extract(Cong_land_cover_built, zones_v, fun = "mean", 
                                           na.rm = TRUE, weights = TRUE)
Cong_land_cover_cropland_df <- terra::extract(Cong_land_cover_cropland, zones_v, fun = "mean", 
                                              na.rm = TRUE, weights = TRUE)
Cong_land_cover_trees_df <- terra::extract(Cong_land_cover_trees, zones_v, fun = "mean", 
                                           na.rm = TRUE, weights = TRUE)
Cong_land_cover_wetland_df <- terra::extract(Cong_land_cover_wetland, zones_v, fun = "mean", 
                                             na.rm = TRUE, weights = TRUE)

# distance to water

# ------------------------------------------------------------------
# Data (all in EPSG:4326)
# ------------------------------------------------------------------
wkt <- st_as_text(st_as_sfc(st_bbox(shp_con)))

rivers <- st_read("data/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp",
                  wkt_filter = wkt)
rivers <- rivers[rivers$ORD_STRA >= min_order, ] |> st_geometry()

lakes <- st_read("data/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp",
                 wkt_filter = wkt)
lakes_geom <- st_make_valid(st_geometry(lakes))
shores <- st_boundary(lakes_geom)             # turns each polygon into its outline (a line)
shores <- shores[!st_is_empty(shores)]
shores <- st_cast(shores, "MULTILINESTRING")  # makes all the outlines one geometry type

xyz_riv  <- densify_xyz(rivers, spacing_m)
xyz_lake <- densify_xyz(shores, spacing_m)


spatial_data_cong$dist_to_water <- zone_dist(st_centroid(pts_m),  crs_m)


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

numeric_cols <- names(spatial_data_cong)[sapply(spatial_data_cong, is.numeric)]

for (col in numeric_cols) {
  spatial_data_cong <- fill_by_nearest(spatial_data_cong, col)
}


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

crs_m <- 6933    # a projected CRS in metres for Gabon (UTM 32S); an equal-area CRS also works

# 1. Zones: 5x5 km squares for survey sites and for grid cells
pts_m   <- st_transform(GAB_df_1_sf, crs_m)
zones_s <- st_buffer(pts_m, dist = 2500, endCapStyle = "SQUARE")   # square centred on each site
zones_s <- zones_s |> st_transform(4326)
zones_v <- terra::vect(zones_s)

gab_ITN_intervention_2021_df <- terra::extract(gab_ITN_intervention_2021,  zones_v, fun = "mean", 
                                               na.rm = TRUE, weights = TRUE)
gab_pop_density_2020_df <- terra::extract(gab_pop,  zones_v, fun = "mean", 
                                          na.rm = TRUE, weights = TRUE)
gab_rainfall_2021_df    <- terra::extract(gab_rainfall_2021, zones_v, fun = "mean", 
                                          na.rm = TRUE, weights = TRUE)
gab_temperature_2021_df <- terra::extract(gab_temperature_2021, zones_v, fun = "mean", 
                                          na.rm = TRUE, weights = TRUE)
ndvi_vals_gab_monthly <- terra::extract(gab_ndvi_stack, zones_v, fun = "mean", 
                                        na.rm = TRUE, weights = TRUE)
ndvi_mean_gab <- apply(ndvi_vals_gab_monthly[,-1], 1, mean, 
                       na.rm = TRUE, weights = TRUE)
elevation_gab_df <- terra::extract(elevation_gab, zones_v, fun = "mean", 
                                   na.rm = TRUE, weights = TRUE)
gab_land_cover_water_df <- terra::extract(gab_land_cover_water, zones_v, fun = "mean", 
                                          na.rm = TRUE, weights = TRUE)
gab_land_cover_built_df <- terra::extract(gab_land_cover_built, zones_v, fun = "mean", 
                                          na.rm = TRUE, weights = TRUE)
gab_land_cover_cropland_df <- terra::extract(gab_land_cover_cropland, zones_v, fun = "mean", 
                                             na.rm = TRUE, weights = TRUE)
gab_land_cover_trees_df <- terra::extract(gab_land_cover_trees, zones_v, fun = "mean", 
                                          na.rm = TRUE, weights = TRUE)
gab_land_cover_wetland_df <- terra::extract(gab_land_cover_wetland, zones_v, fun = "mean", 
                                            na.rm = TRUE, weights = TRUE)

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

# ------------------------------------------------------------------
# Data (all in EPSG:4326)
# ------------------------------------------------------------------
wkt <- st_as_text(st_as_sfc(st_bbox(shp_gab)))

rivers <- st_read("data/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp",
                  wkt_filter = wkt)
rivers <- rivers[rivers$ORD_STRA >= min_order, ] |> st_geometry()

lakes <- st_read("data/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp",
                 wkt_filter = wkt)
lakes_geom <- st_make_valid(st_geometry(lakes))
shores <- st_boundary(lakes_geom)             # turns each polygon into its outline (a line)
shores <- shores[!st_is_empty(shores)]
shores <- st_cast(shores, "MULTILINESTRING")  # makes all the outlines one geometry type

xyz_riv  <- densify_xyz(rivers, spacing_m)
xyz_lake <- densify_xyz(shores, spacing_m)


GAB_df_1_sf$dist_to_water <- zone_dist(st_centroid(pts_m),  crs_m)

numeric_cols <- names(GAB_df_1_sf)[sapply(GAB_df_1_sf, is.numeric)]

for (col in numeric_cols) {
  GAB_df_1_sf <- fill_by_nearest(GAB_df_1_sf, col)
}

st_write(GAB_df_1_sf, "data/gab_data/GAB_df_1_sf.gpkg", append = FALSE)

######################## Prediction grid for Gabon #############################

# resolution in meters (5000 m = 5 km)
res <- 5000  

box_gab <- st_bbox(shp_gab %>% st_transform(3857)) %>% st_as_sfc()
# Create grid covering bounding box of Gabon
grid <- st_make_grid(
  box_gab,
  cellsize = res,
  what = "centers"
)

grid <- st_sf(geometry = grid)

# Keep only grid points inside Gabon
grid <- st_intersection(grid, shp_gab%>% st_transform(3857))

gab_vect <- vect(grid %>% st_transform(4326))

ITN_intervention_2021 <- rast("data/ITN_2000/2025_GBD2024_Africa_ITN_2021.tif")

elevation_grid <- rast("data/gab_data/elevation/gab_elv_msk.tif")
grid_ndvi_stack<- rast(ndvi_files)
# grid_ITN_usage_2021  <- project(ITN_intervention_2021,  st_crs(grid)$wkt)
grid_pop  <- crop(gab_pop,  gab_vect)
grid_pop  <- mask(gab_pop,  gab_vect)

grid_land_cover_water <- crop(land_cover_water, gab_vect)
grid_land_cover_water <- mask(grid_land_cover_water, gab_vect)

grid_land_cover_built <- crop(land_cover_builts, gab_vect)
grid_land_cover_built <- mask(grid_land_cover_built, gab_vect)

grid_land_cover_cropland <- crop(land_cover_cropland, gab_vect)
grid_land_cover_cropland <- mask(grid_land_cover_cropland, gab_vect)

grid_land_cover_trees <- crop(land_cover_trees, gab_vect)
grid_land_cover_trees <- mask(grid_land_cover_trees, gab_vect)

grid_land_cover_wetland <- crop(land_cover_wetland, gab_vect)
grid_land_cover_wetland <- mask(grid_land_cover_wetland, gab_vect)

grid_rainfall_2021 <- crop(rainfall_2021, gab_vect)
grid_rainfall_2021 <- mask(grid_rainfall_2021, gab_vect)

grid_temperature_2021 <- crop(annual_temperature_climate , gab_vect)
grid_temperature_2021 <- mask(grid_temperature_2021, gab_vect)


grid_ITN_intervention_2021  <- crop(ITN_intervention_2021, gab_vect)
grid_ITN_intervention_2021  <- mask(grid_ITN_intervention_2021, gab_vect)

itn_est <- grid_ITN_intervention_2021[["Data"]]

crs_m <- 6933    # a projected CRS in metres for Gabon (UTM 32S); an equal-area CRS also works

# 1. Zones: 5x5 km squares for survey sites and for grid cells
pts_m   <- st_transform(grid, crs_m)
zones_s <- st_buffer(pts_m, dist = 2500, endCapStyle = "SQUARE")   # square centred on each site
zones_s <- zones_s |> st_transform(4326)


# grid_ITN_usage_2021_df <- terra::extract(grid_ITN_usage_2021,  grid)
grid$pop_density_2020 <- exact_extract(grid_pop, zones_s, fun = "mean", progress = TRUE)
grid$annual_rainfall_2021    <- exact_extract(grid_rainfall_2021, zones_s, fun = "mean", progress = TRUE)
grid$annual_temperature_2021 <- exact_extract(grid_temperature_2021, zones_s, fun = "mean", progress = TRUE)
ndvi_vals_grid_monthly <- exact_extract(grid_ndvi_stack, zones_s, fun = "mean", progress = TRUE)
ndvi_mean_grid <- apply(ndvi_vals_grid_monthly[,-1], 1, mean, na.rm = TRUE, weights = TRUE)
grid$vegetation_index <- ndvi_mean_grid
grid$elevation <- exact_extract(elevation_grid, zones_s, fun = "mean", progress = TRUE)
grid$lc_water <- exact_extract(grid_land_cover_water, zones_s, fun = "mean", progress = TRUE)
grid$lc_built <- exact_extract(grid_land_cover_built, zones_s, fun = "mean", progress = TRUE)
grid$lc_cropland <- exact_extract(grid_land_cover_cropland, zones_s, fun = "mean", progress = TRUE)
grid$lc_trees <- exact_extract(grid_land_cover_trees, zones_s, fun = "mean", progress = TRUE)
grid$lc_wetland <- exact_extract(grid_land_cover_wetland, zones_s, fun = "mean", progress = TRUE)
grid$annual_ITN_usage_2021 <- exact_extract(itn_est, zones_s, fun = "mean", progress = TRUE)


# distance to water

# ------------------------------------------------------------------
# Data (all in EPSG:4326)
# ------------------------------------------------------------------
wkt <- st_as_text(st_as_sfc(st_bbox(shp_gab)))

rivers <- st_read("data/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp",
                  wkt_filter = wkt)
rivers <- rivers[rivers$ORD_STRA >= min_order, ] |> st_geometry()

lakes <- st_read("data/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp",
                 wkt_filter = wkt)
lakes_geom <- st_make_valid(st_geometry(lakes))
shores <- st_boundary(lakes_geom)             # turns each polygon into its outline (a line)
shores <- shores[!st_is_empty(shores)]
shores <- st_cast(shores, "MULTILINESTRING")  # makes all the outlines one geometry type

xyz_riv  <- densify_xyz(rivers, spacing_m)
xyz_lake <- densify_xyz(shores, spacing_m)


grid$dist_to_water <- zone_dist(st_centroid(pts_m),  crs_m)

numeric_cols <- names(grid)[sapply(grid, is.numeric)]

for (col in numeric_cols) {
  grid <- fill_by_nearest(grid, col)
}

st_write(grid, "data/gab_data/gab_map_grid.gpkg", append = FALSE)

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

crs_m <- 6933    # a projected CRS in metres for Gabon (UTM 32S); an equal-area CRS also works

# 1. Zones: 5x5 km squares for survey sites and for grid cells
pts_m   <- st_transform(CMR_df_1_sf, crs_m)
zones_s <- st_buffer(pts_m, dist = 2500, endCapStyle = "SQUARE")   # square centred on each site
zones_s <- zones_s |> st_transform(4326)
zones_v <- terra::vect(zones_s)

cmr_ITN_intervention_2021_df <- terra::extract(cmr_ITN_intervention_2021,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_pop_density_2020_df <- terra::extract(cmr_pop,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_rainfall_2021_df    <- terra::extract(cmr_rainfall_2021, zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_temperature_2021_df <- terra::extract(cmr_temperature_2021, zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
ndvi_vals_cmr_monthly <- terra::extract(cmr_ndvi_stack, zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
ndvi_mean_cmr <- apply(ndvi_vals_cmr_monthly[,-1], 1, mean, na.rm = TRUE, weights = TRUE)
elevation_cmr_df <- terra::extract(elevation_cmr, zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_land_cover_water_df <- terra::extract(cmr_land_cover_water, zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_land_cover_built_df <- terra::extract(cmr_land_cover_built, zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_land_cover_cropland_df <- terra::extract(cmr_land_cover_cropland, zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_land_cover_trees_df <- terra::extract(cmr_land_cover_trees, zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_land_cover_wetland_df <- terra::extract(cmr_land_cover_wetland, zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)

# ------------------------------------------------------------------
# Data (all in EPSG:4326)
# ------------------------------------------------------------------
wkt <- st_as_text(st_as_sfc(st_bbox(shp_cmr)))

rivers <- st_read("data/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp",
                  wkt_filter = wkt)
rivers <- rivers[rivers$ORD_STRA >= min_order, ] |> st_geometry()

lakes <- st_read("data/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp",
                 wkt_filter = wkt)
lakes_geom <- st_make_valid(st_geometry(lakes))
shores <- st_boundary(lakes_geom)             # turns each polygon into its outline (a line)
shores <- shores[!st_is_empty(shores)]
shores <- st_cast(shores, "MULTILINESTRING")  # makes all the outlines one geometry type

xyz_riv  <- densify_xyz(rivers, spacing_m)
xyz_lake <- densify_xyz(shores, spacing_m)


CMR_df_1_sf$dist_to_water <- zone_dist(st_centroid(pts_m),  crs_m)



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

numeric_cols <- names(grid)[sapply(grid, is.numeric)]

for (col in numeric_cols) {
  grid <- fill_by_nearest(grid, col)
}

st_write(CMR_df_1_sf, "data/cmr_data/CMR_df_1_sf.gpkg", append = FALSE)


# Cameroon case with clusters of households

# 1. Zones: 5x5 km squares for survey sites and for grid cells
crs_m <- 6933
pts_m   <- st_transform(cluster_data_cmr_sf, crs_m)
zones_s <- st_buffer(pts_m, dist = 2500, endCapStyle = "SQUARE")   # square centred on each site
zones_s <- zones_s |> st_transform(4326)
zones_v <- terra::vect(zones_s)

cmr_cluster_ITN_intervention_2021_df <- terra::extract(cmr_ITN_intervention_2021,
                                                       zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_cluster_pop_density_2020_df <- terra::extract(cmr_pop,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_cluster_rainfall_2021_df    <- terra::extract(cmr_rainfall_2021, zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_cluster_temperature_2021_df <- terra::extract(cmr_temperature_2021, zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
ndvi_vals_cluster_cmr_monthly <- terra::extract(cmr_ndvi_stack, zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
ndvi_mean_cluster_cmr <- apply(ndvi_vals_cluster_cmr_monthly[,-1], 1, mean, na.rm = TRUE, weights = TRUE)
elevation_cluster_cmr_df <- terra::extract(elevation_cmr, zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_cluster_land_cover_water_df <- terra::extract(cmr_land_cover_water, zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_cluster_land_cover_built_df <- terra::extract(cmr_land_cover_built, zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_cluster_land_cover_cropland_df <- terra::extract(cmr_land_cover_cropland, zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_cluster_land_cover_trees_df <- terra::extract(cmr_land_cover_trees, zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_cluster_land_cover_wetland_df <- terra::extract(cmr_land_cover_wetland, zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)


# ------------------------------------------------------------------
# Data (all in EPSG:4326)
# ------------------------------------------------------------------
wkt <- st_as_text(st_as_sfc(st_bbox(shp_cmr)))

rivers <- st_read("data/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp",
                  wkt_filter = wkt)
rivers <- rivers[rivers$ORD_STRA >= min_order, ] |> st_geometry()

lakes <- st_read("data/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp",
                 wkt_filter = wkt)
lakes_geom <- st_make_valid(st_geometry(lakes))
shores <- st_boundary(lakes_geom)             # turns each polygon into its outline (a line)
shores <- shores[!st_is_empty(shores)]
shores <- st_cast(shores, "MULTILINESTRING")  # makes all the outlines one geometry type

xyz_riv  <- densify_xyz(rivers, spacing_m)
xyz_lake <- densify_xyz(shores, spacing_m)


cluster_data_cmr_sf$dist_to_water <- zone_dist(st_centroid(pts_m),  crs_m)

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

numeric_cols <- names(cluster_data_cmr_sf)[sapply(cluster_data_cmr_sf, is.numeric)]

for (col in numeric_cols) {
  cluster_data_cmr_sf <- fill_by_nearest(cluster_data_cmr_sf, col)
}


st_write(cluster_data_cmr_sf, "data/cmr_data/cmr_cluster_df_1_sf.gpkg", append = FALSE)


####################### Cameroon case benchmark ################################

map_2011_cmr <- st_read("data/cmr_data/cmr_MAP_full_raw_prevalence_2011.csv")

map_2011_cmr_sf <- st_as_sf(map_2011_cmr, coords = c("longitude", "latitude"), crs = 4326)


path3 <- "data/cmr_data/chirps_2011/CHIRPS_PENTAD_precipitation_2011-01-01_2011-12-31_-5p8264to13p9448N_0p6930to19p1696E.precipitation.tif"
path4 <- "data/cmr_data/ERA5_temperature/ERA5_LAND_DAILY_temperature_2m_2011-01-01_2011-12-31_-5p4765to13p9448N_0p6930to19p3454E.temperature_2m.tif"
ndvi_path <- "data/cmr_data/NDVI2011/"
# Load all 12 monthly NDVI tif files
ndvi_files <- list.files(
  ndvi_path,
  pattern    = "_1_km_monthly_NDVI.*\\.tif$",
  full.names = TRUE
)

cmr_ndvi_stack<- rast(ndvi_files)


cmr_pop <- rast("data/cmr_data/cmr_pd_2011_1km_UNadj.tif")
cmr_rainfall_2011 <- rast(path3)
cmr_temperature_2011 <- rast(path4)
elevation_cmr <- rast("data/cmr_data/elevation/cmr_elv_msk.tif")
cmr_ITN_2011  <- rast("data/ITN_2000/2025_GBD2024_Africa_ITN_2011.tif")
cmr_ITN_2011  <- rast("data/ITN_2000_2/2026_GBD2025_Africa_ITN_Use_Rate_2011.tif")
cmr_IRS_2011  <- rast("data/IRS_2000/2025_GBD2024_Africa_IRS_2011.tif")
cmr_effective_treat_2011 <- rast("data/Effective_treat_2000_2025/2026_GBD2025_Global_Antimalarial_EFT_2011.tif")

cmr_pop  <- crop(cmr_pop,  cmr_vect)
cmr_pop  <- mask(cmr_pop,  cmr_vect)

cmr_rainfall_2011 <- crop(cmr_rainfall_2011, cmr_vect)
cmr_rainfall_2011 <- mask(cmr_rainfall_2011, cmr_vect)

cmr_temperature_2011 <- crop(cmr_temperature_2011, cmr_vect)
cmr_temperature_2011 <- mask(cmr_temperature_2011, cmr_vect)

cmr_ndvi_stack <- crop(cmr_ndvi_stack, cmr_vect)
cmr_ndvi_stack <- mask(cmr_ndvi_stack, cmr_vect)

elevation_cmr   <- crop(elevation_cmr, cmr_vect)
elevation_cmr   <- mask(elevation_cmr, cmr_vect)

cmr_land_cover_water     <- crop(land_cover_water, cmr_vect)
cmr_land_cover_water     <- mask(cmr_land_cover_water, cmr_vect)

cmr_land_cover_built     <- crop(land_cover_builts, cmr_vect)
cmr_land_cover_built     <- mask(cmr_land_cover_built, cmr_vect)

cmr_land_cover_cropland  <- crop(land_cover_cropland, cmr_vect)
cmr_land_cover_cropland  <- mask(cmr_land_cover_cropland, cmr_vect)

cmr_land_cover_trees     <- crop(land_cover_trees, cmr_vect)
cmr_land_cover_trees     <- mask(cmr_land_cover_trees, cmr_vect)

cmr_land_cover_wetland   <- crop(land_cover_wetland, cmr_vect)
cmr_land_cover_wetland   <- mask(cmr_land_cover_wetland, cmr_vect)

cmr_ITN_2011   <- crop(cmr_ITN_2011, cmr_vect)
cmr_ITN_2011   <- mask(cmr_ITN_2011, cmr_vect)

cmr_IRS_2011   <- crop(cmr_IRS_2011, cmr_vect)
cmr_IRS_2011   <- mask(cmr_IRS_2011, cmr_vect)

cmr_effective_treat_2011   <- crop(cmr_effective_treat_2011, cmr_vect)
cmr_effective_treat_2011   <- mask(cmr_effective_treat_2011, cmr_vect)

# 1. Zones: 5x5 km squares for survey sites and for grid cells
crs_m <- 6933
pts_m   <- st_transform(map_2011_cmr_sf, crs_m)
zones_s <- st_buffer(pts_m, dist = 2500, endCapStyle = "SQUARE")   # square centred on each site
zones_s <- zones_s |> st_transform(4326)
zones_v <- terra::vect(zones_s)

cmr_pop_density_2011_df <- terra::extract(cmr_pop,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_rainfall_2011_df    <- terra::extract(cmr_rainfall_2011,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_temperature_2011_df <- terra::extract(cmr_temperature_2011,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
ndvi_vals_cmr_monthly <- terra::extract(cmr_ndvi_stack,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
ndvi_mean_cmr <- apply(ndvi_vals_cmr_monthly[,-1], 1, mean, na.rm = TRUE, weights = TRUE)
elevation_cmr_df <- terra::extract(elevation_cmr,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_land_cover_water_df <- terra::extract(cmr_land_cover_water,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_land_cover_built_df <- terra::extract(cmr_land_cover_built,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_land_cover_cropland_df <- terra::extract(cmr_land_cover_cropland,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_land_cover_trees_df <- terra::extract(cmr_land_cover_trees,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_land_cover_wetland_df <- terra::extract(cmr_land_cover_wetland,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_ITN_2011_df   <- terra::extract(cmr_ITN_2011,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_IRS_2011_df   <- terra::extract(cmr_IRS_2011,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_effective_treat_2011_df   <- terra::extract(cmr_effective_treat_2011,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)


# distance to water

# ------------------------------------------------------------------
# Data (all in EPSG:4326)
# ------------------------------------------------------------------
wkt <- st_as_text(st_as_sfc(st_bbox(shp_cmr)))

rivers <- st_read("data/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp",
                  wkt_filter = wkt)
rivers <- rivers[rivers$ORD_STRA >= min_order, ] |> st_geometry()

lakes <- st_read("data/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp",
                 wkt_filter = wkt)
lakes_geom <- st_make_valid(st_geometry(lakes))
shores <- st_boundary(lakes_geom)             # turns each polygon into its outline (a line)
shores <- shores[!st_is_empty(shores)]
shores <- st_cast(shores, "MULTILINESTRING")  # makes all the outlines one geometry type

xyz_riv  <- densify_xyz(rivers, spacing_m)
xyz_lake <- densify_xyz(shores, spacing_m)


map_2011_cmr_sf$dist_to_water <- zone_dist(st_centroid(pts_m),  crs_m)

cmr_pf_data_2011 <- map_2011_cmr_sf %>% 
  dplyr::select(site_id, lower_age, upper_age, examined, positive, pr, species, dist_to_water) %>% 
  mutate(pop_density = cmr_pop_density_2011_df[,2],
         rainfall = cmr_rainfall_2011_df[,2],
         temperature = cmr_temperature_2011_df[,2],
         elevation = elevation_cmr_df[,2],
         vegetation_index = ndvi_mean_cmr,
         landcover_water = cmr_land_cover_water_df[,2],
         landcover_trees = cmr_land_cover_trees_df[,2],
         landcover_builts = cmr_land_cover_built_df[,2],
         landcover_cropland = cmr_land_cover_cropland_df[,2],
         landcover_wetland = cmr_land_cover_wetland_df[,2],
         ITN_2011         = cmr_ITN_2011_df[,2],
         effective_treat  = cmr_effective_treat_2011_df[,2]) 

numeric_cols <- names(cmr_pf_data_2011)[sapply(cmr_pf_data_2011, is.numeric)]

for (col in numeric_cols) {
  cmr_pf_data_2011 <- fill_by_nearest(cmr_pf_data_2011, col)
}

st_write(cmr_pf_data_2011, 
         "data/cmr_data/benchmarking_spatial_data_cmr_sf.gpkg", append = FALSE)


################# Cameroon benchmark mix #######################################

# 1. Zones: 5x5 km squares for survey sites and for grid cells
crs_m <- 6933
pts_m   <- st_transform(cmr_raw_pr_mix_sf, crs_m)
zones_s <- st_buffer(pts_m, dist = 2500, endCapStyle = "SQUARE")   # square centred on each site
zones_s <- zones_s |> st_transform(4326)
zones_v <- terra::vect(zones_s)


cmr_pop_density_2011_df <- terra::extract(cmr_pop,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_rainfall_2011_df    <- terra::extract(cmr_rainfall_2011,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_temperature_2011_df <- terra::extract(cmr_temperature_2011,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
ndvi_vals_cmr_monthly <- terra::extract(cmr_ndvi_stack,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
ndvi_mean_cmr <- apply(ndvi_vals_cmr_monthly[,-1], 1, mean, na.rm = TRUE, weights = TRUE)
elevation_cmr_df <- terra::extract(elevation_cmr,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_land_cover_water_df <- terra::extract(cmr_land_cover_water,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_land_cover_built_df <- terra::extract(cmr_land_cover_built,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_land_cover_cropland_df <- terra::extract(cmr_land_cover_cropland,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_land_cover_trees_df <- terra::extract(cmr_land_cover_trees,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_land_cover_wetland_df <- terra::extract(cmr_land_cover_wetland,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_ITN_2011_df   <- terra::extract(cmr_ITN_2011,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_IRS_2011_df   <- terra::extract(cmr_IRS_2011,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)
cmr_effective_treat_2011_df   <- terra::extract(cmr_effective_treat_2011,  zones_v, fun = "mean", na.rm = TRUE, weights = TRUE)


# distance to water

# ------------------------------------------------------------------
# Data (all in EPSG:4326)
# ------------------------------------------------------------------
wkt <- st_as_text(st_as_sfc(st_bbox(shp_cmr)))

rivers <- st_read("data/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp",
                  wkt_filter = wkt)
rivers <- rivers[rivers$ORD_STRA >= min_order, ] |> st_geometry()

lakes <- st_read("data/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp",
                 wkt_filter = wkt)
lakes_geom <- st_make_valid(st_geometry(lakes))
shores <- st_boundary(lakes_geom)             # turns each polygon into its outline (a line)
shores <- shores[!st_is_empty(shores)]
shores <- st_cast(shores, "MULTILINESTRING")  # makes all the outlines one geometry type

xyz_riv  <- densify_xyz(rivers, spacing_m)
xyz_lake <- densify_xyz(shores, spacing_m)


cmr_raw_pr_mix_sf$dist_to_water <- zone_dist(st_centroid(pts_m),  crs_m)

cmr_raw_pr_mix_sf_final <- cmr_raw_pr_mix_sf %>% 
  dplyr::select(site_id, lower_age, upper_age, examined, positive, pr, species, dist_to_water) %>% 
  mutate(pop_density        = cmr_pop_density_2011_df[,2],
         rainfall           = cmr_rainfall_2011_df[,2],
         temperature        = cmr_temperature_2011_df[,2],
         elevation          = elevation_cmr_df[,2],
         vegetation_index   = ndvi_mean_cmr,
         landcover_water    = cmr_land_cover_water_df[,2],
         landcover_trees    = cmr_land_cover_trees_df[,2],
         landcover_builts   = cmr_land_cover_built_df[,2],
         landcover_cropland = cmr_land_cover_cropland_df[,2],
         landcover_wetland  = cmr_land_cover_wetland_df[,2],
         ITN_2011           = cmr_ITN_2011_df[,2],
         IRS_2011           = cmr_IRS_2011_df[,2],
         effective_treat    = cmr_effective_treat_2011_df[,2]) 


numeric_cols <- names(cmr_raw_pr_mix_sf_final)[sapply(cmr_raw_pr_mix_sf_final, is.numeric)]

for (col in numeric_cols) {
  cmr_raw_pr_mix_sf_final <- fill_by_nearest(cmr_raw_pr_mix_sf_final, col)
}

st_write(cmr_raw_pr_mix_sf_final, 
         "data/cmr_data/benchmarking_spatial_data_cmr_sf_mix.gpkg", append = FALSE)



################ Prediction grid for benchmark exercise Cameroon ###############

# resolution in meters (5000 m = 5 km)
res <- 5000  

box_cmr <- st_bbox(shp_cmr %>% st_transform(3857)) %>% st_as_sfc()
# Create grid covering bounding box of Cameroon
grid <- st_make_grid(
  box_cmr,
  cellsize = res,
  what = "centers"
)

grid <- st_sf(geometry = grid)

# Keep only grid points inside Cameroon
grid <- st_intersection(grid, shp_cmr%>% st_transform(3857))

cmr_vect <- vect(grid %>% st_transform(4326))

# grid_ITN_usage_2021  <- project(ITN_intervention_2021,  st_crs(grid)$wkt)
grid_pop  <- crop(cmr_pop,  cmr_vect)
grid_pop  <- mask(cmr_pop,  cmr_vect)

grid_rainfall_2021 <- crop(cmr_rainfall_2011, cmr_vect)
grid_rainfall_2021 <- mask(cmr_rainfall_2011, cmr_vect)

grid_temperature_2021 <- crop(cmr_temperature_2011, cmr_vect)
grid_temperature_2021 <- mask(cmr_temperature_2011, cmr_vect)

grid_ndvi_stack <- crop(cmr_ndvi_stack, cmr_vect)
grid_ndvi_stack <- mask(cmr_ndvi_stack, cmr_vect)

elevation_grid   <- crop(elevation_cmr, cmr_vect)
elevation_grid   <- mask(elevation_cmr, cmr_vect)

grid_land_cover_water   <- crop(land_cover_water, cmr_vect)
grid_land_cover_water   <- mask(grid_land_cover_water, cmr_vect)

grid_land_cover_built   <- crop(land_cover_builts, cmr_vect)
grid_land_cover_built   <- mask(grid_land_cover_built, cmr_vect)

grid_land_cover_cropland   <- crop(land_cover_cropland, cmr_vect)
grid_land_cover_cropland   <- mask(grid_land_cover_cropland, cmr_vect)

grid_land_cover_trees   <- crop(land_cover_trees, cmr_vect)
grid_land_cover_trees   <- mask(grid_land_cover_trees, cmr_vect)

grid_land_cover_wetland   <- crop(land_cover_wetland, cmr_vect)
grid_land_cover_wetland   <- mask(grid_land_cover_wetland, cmr_vect)

grid_ITN   <- crop(cmr_ITN_2011, cmr_vect)
grid_ITN   <- mask(grid_ITN, cmr_vect)

grid_IRS   <- crop(cmr_IRS_2011, cmr_vect)
grid_IRS   <- mask(grid_IRS, cmr_vect)

itn_est <- grid_ITN[["Data"]]
irs_est <- grid_IRS[["Data"]]

grid_effective_treat   <- crop(cmr_effective_treat_2011, cmr_vect)
grid_effective_treat   <- mask(grid_effective_treat, cmr_vect)

# 1. Zones: 5x5 km squares for survey sites and for grid cells
crs_m <- 6933
pts_m   <- st_transform(grid, crs_m)
zones_s <- st_buffer(pts_m, dist = 2500, endCapStyle = "SQUARE")   # square centred on each site
zones_s <- zones_s |> st_transform(4326)
zones_v <- terra::vect(zones_s)


# grid_ITN_usage_2021_df <- terra::extract(grid_ITN_usage_2021,  grid)
grid_pop_density_2020_df <- exact_extract(grid_pop,  zones_s, fun = "mean", progress = TRUE)
grid_rainfall_2021_df    <- exact_extract(grid_rainfall_2021, zones_s, fun = "mean", progress = TRUE)
grid_temperature_2021_df <- exact_extract(grid_temperature_2021, zones_s, fun = "mean", progress = TRUE)
ndvi_vals_grid_monthly <- exact_extract(grid_ndvi_stack, zones_s, fun = "mean", progress = TRUE)
ndvi_mean_grid <- apply(ndvi_vals_grid_monthly[,-1], 1, mean, na.rm = TRUE, weights = TRUE)
elevation_grid_df <- exact_extract(elevation_grid, zones_s, fun = "mean", progress = TRUE)
grid_land_cover_water_df <- exact_extract(grid_land_cover_water, zones_s, fun = "mean", progress = TRUE)
grid_land_cover_built_df <- exact_extract(grid_land_cover_built, zones_s, fun = "mean", progress = TRUE)
grid_land_cover_cropland_df <- exact_extract(grid_land_cover_cropland, zones_s, fun = "mean", progress = TRUE)
grid_land_cover_trees_df <- exact_extract(grid_land_cover_trees, zones_s, fun = "mean", progress = TRUE)
grid_land_cover_wetland_df <- exact_extract(grid_land_cover_wetland, zones_s, fun = "mean", progress = TRUE)
grid_ITN_df                <- exact_extract(itn_est, zones_s, fun = "mean", progress = TRUE)
grid_effective_treat_df    <- exact_extract(grid_effective_treat, zones_s, fun = "mean", progress = TRUE)


# Distance to water

# Data (all in EPSG:4326)
# ------------------------------------------------------------------
wkt <- st_as_text(st_as_sfc(st_bbox(shp_cmr)))

rivers <- st_read("data/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp",
                  wkt_filter = wkt)
rivers <- rivers[rivers$ORD_STRA >= min_order, ] |> st_geometry()

lakes <- st_read("data/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp",
                 wkt_filter = wkt)
lakes_geom <- st_make_valid(st_geometry(lakes))
shores <- st_boundary(lakes_geom)             # turns each polygon into its outline (a line)
shores <- shores[!st_is_empty(shores)]
shores <- st_cast(shores, "MULTILINESTRING")  # makes all the outlines one geometry type

xyz_riv  <- densify_xyz(rivers, spacing_m)
xyz_lake <- densify_xyz(shores, spacing_m)


grid$dist_to_water <- zone_dist(st_centroid(pts_m),  crs_m)

grid <- grid %>% 
  mutate(rainfall            = grid_rainfall_2021_df,
         temperature         = grid_temperature_2021_df,
         pop_density         = grid_pop_density_2020_df,
         elevation           = elevation_grid_df,
         vegetation_index    = ndvi_mean_grid,
         landcover_cropland  = grid_land_cover_cropland_df,
         landcover_water     = grid_land_cover_water_df,
         landcover_wetland   = grid_land_cover_wetland_df,
         landcover_builts    = grid_land_cover_built_df,
         landcover_trees     = grid_land_cover_trees_df,
         ITN_2011            = grid_ITN_df,
         effective_treat     = grid_effective_treat_df) 

numeric_cols <- names(grid)[sapply(grid, is.numeric)]

for (col in numeric_cols) {
  grid <- fill_by_nearest(grid, col)
}

st_write(grid, "data/cmr_data/cmr_map_grid.gpkg", append = FALSE)


############ Prediction grid at the continental level benchmarking #############

africa_pf_PR_with_coord_join_2021 <- read.csv("data/Africa_data/africa_MAP_full_raw_prevalence_2021.csv")

africa_pf_PR_with_coord_join_2021_sf <- st_as_sf(africa_pf_PR_with_coord_join_2021,
                                                 coords = c("longitude", "latitude"),
                                                 crs = 4326)

shp_africa <- st_read("data/Africa_data/shp_africa/Africa_simplified.shp") 

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
IRS_intervention_2021 <- rast("data/IRS_2000/2025_GBD2024_Africa_IRS_2021.tif")
africa_rainfall_2021 <- rast(path1)
africa_temperature_2021 <- rast(path2)
africa_pop <- rast("data/Africa_data/pop_density/gpw_v4_population_density_rev11_2020_30_sec.tif")
africa_elevation <- rast("data/Africa_data/elevation/wc2.1_30s_elev.tif")
africa_ndvi_stack<- rast(ndvi_files)

vect_africa <- vect(shp_africa)

africa_land_cover_water <- crop(land_cover_water, vect_africa)
africa_land_cover_water <- mask(africa_land_cover_water, vect_africa)

africa_land_cover_built <- crop(land_cover_builts, vect_africa)
africa_land_cover_built <- mask(africa_land_cover_built, vect_africa)

africa_land_cover_cropland <- crop(land_cover_cropland, vect_africa)
africa_land_cover_cropland <- mask(africa_land_cover_cropland, vect_africa)

africa_land_cover_trees <- crop(land_cover_trees, vect_africa)
africa_land_cover_trees <- mask(africa_land_cover_trees, vect_africa)

africa_land_cover_wetland <- crop(land_cover_wetland, vect_africa)
africa_land_cover_wetland <- mask(africa_land_cover_wetland, vect_africa)

africa_ITN_intervention_2021 <- crop(ITN_intervention_2021, vect_africa)
africa_ITN_intervention_2021 <- mask(africa_ITN_intervention_2021, vect_africa)

itn_est <- africa_ITN_intervention_2021[["Data"]]

africa_IRS_intervention_2021 <- crop(IRS_intervention_2021, vect_africa)
africa_IRS_intervention_2021 <- mask(africa_IRS_intervention_2021, vect_africa)


africa_pop_2020 <- crop(africa_pop, vect_africa)
africa_pop_2020 <- mask(africa_pop_2020, vect_africa)

africa_rainfall_2021 <- crop(africa_rainfall_2021, vect_africa)
africa_rainfall_2021 <- mask(africa_rainfall_2021, vect_africa)

africa_temperature_2021 <- crop(africa_temperature_2021, vect_africa)
africa_temperature_2021 <- mask(africa_temperature_2021, vect_africa)

africa_ndvi_stack <- crop(africa_ndvi_stack, vect_africa)
africa_ndvi_stack <- mask(africa_ndvi_stack, vect_africa)

africa_elevation <- crop(africa_elevation, vect_africa)
africa_elevation <- mask(africa_elevation, vect_africa)

# 1. Zones: 5x5 km squares for survey sites and for grid cells
crs_m <- "+proj=laea +lat_0=2 +lon_0=20 +datum=WGS84 +units=m +no_defs"
pts_m   <- st_transform(africa_pf_PR_with_coord_join_2021_sf, crs_m)
zones_s <- st_buffer(pts_m, dist = 2500, endCapStyle = "SQUARE")   # square centred on each site
zones_s <- zones_s |> st_transform(4326)

africa_ITN_intervention_2021_df <- exact_extract(itn_est, 
                                                 zones_s, fun = "mean", progress = TRUE)
africa_IRS_intervention_2021_df <- exact_extract(africa_IRS_intervention_2021, 
                                                 zones_s, fun = "mean", progress = TRUE)
africa_pop_density_2020_df <- exact_extract(africa_pop_2020,  
                                            zones_s, fun = "mean", progress = TRUE)
africa_rainfall_2021_df    <- exact_extract(africa_rainfall_2021, 
                                            zones_s, fun = "mean", progress = TRUE)
africa_temperature_2021_df <- exact_extract(africa_temperature_2021, 
                                            zones_s, fun = "mean", progress = TRUE)

ndvi_vals_africa_monthly <- exact_extract(africa_ndvi_stack, 
                                          zones_s, fun = "mean", progress = TRUE)
ndvi_mean_africa <- apply(ndvi_vals_africa_monthly[,-1], 1, mean, na.rm = TRUE)

africa_elevation_df <- exact_extract(africa_elevation,
                                     zones_s, fun = "mean", progress = TRUE)
africa_land_cover_water_df <- exact_extract(africa_land_cover_water, 
                                            zones_s, fun = "mean", progress = TRUE)
africa_land_cover_built_df <- exact_extract(africa_land_cover_built,
                                            zones_s, fun = "mean", progress = TRUE)
africa_land_cover_cropland_df <- exact_extract(africa_land_cover_cropland,
                                               zones_s, fun = "mean", progress = TRUE)
africa_land_cover_trees_df <- exact_extract(africa_land_cover_trees,
                                            zones_s, fun = "mean", progress = TRUE)
africa_land_cover_wetland_df <- exact_extract(africa_land_cover_wetland,
                                              zones_s, fun = "mean", progress = TRUE)

# Distance to water

shp_africa <- st_make_valid(shp_africa)

# Data (all in EPSG:4326)
# ------------------------------------------------------------------
wkt <- st_as_text(st_as_sfc(st_bbox(shp_africa)))

rivers <- st_read("data/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp",
                  wkt_filter = wkt)
rivers <- rivers[rivers$ORD_STRA >= min_order, ] |> st_geometry()

lakes <- st_read("data/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp",
                 wkt_filter = wkt)
lakes_geom <- st_make_valid(st_geometry(lakes))
shores <- st_boundary(lakes_geom)             # turns each polygon into its outline (a line)
shores <- shores[!st_is_empty(shores)]
shores <- st_cast(shores, "MULTILINESTRING")  # makes all the outlines one geometry type

xyz_riv  <- densify_xyz(rivers, spacing_m)
xyz_lake <- densify_xyz(shores, spacing_m)


africa_pf_PR_with_coord_join_2021_sf$dist_to_water <- zone_dist(st_centroid(pts_m),  crs_m)

africa_pf_PR_with_coord_join_2021_sf <- africa_pf_PR_with_coord_join_2021_sf %>% 
  mutate(pop_density_2020 = africa_pop_density_2020_df,
         annual_rainfall_2021 = africa_rainfall_2021_df,
         annual_temperature_2021 = africa_temperature_2021_df,
         vegetation_index = ndvi_mean_africa,
         elevation = africa_elevation_df,
         lc_water = africa_land_cover_water_df,
         lc_built = africa_land_cover_built_df,
         lc_trees = africa_land_cover_trees_df,
         lc_wetland = africa_land_cover_wetland_df,
         lc_cropland = africa_land_cover_cropland_df,
         annual_ITN_usage_2021 = africa_ITN_intervention_2021_df,
         annual_IRS_usage_2021 = africa_IRS_intervention_2021_df)

# filling missing values
numeric_cols <- names(africa_pf_PR_with_coord_join_2021_sf)[sapply(africa_pf_PR_with_coord_join_2021_sf, is.numeric)]

for (col in numeric_cols) {
  africa_pf_PR_with_coord_join_2021_sf <- fill_by_nearest(africa_pf_PR_with_coord_join_2021_sf, col)
}



st_write(africa_pf_PR_with_coord_join_2021_sf, 
         "data/Africa_data/africa_pf_PR_with_coord_join_2021_sf.gpkg", append = FALSE)

################################################################################
#                 case of grid cells at country level                          #
################################################################################
countries1 <- rbind(shp_ben %>% dplyr::select(COUNTRY, NAME_3, geometry) %>% 
                         rename(country = COUNTRY,
                                district = NAME_3), shp_con %>% dplyr::select(adm0_name,
                                                                              adm2_name,
                                                                              geometry) %>% 
                         rename(country = adm0_name,
                                district = adm2_name)) 


countries2 <- rbind(countries1, shp_cmr %>% dplyr::select(Admin2, geometry) %>% 
                        mutate(country = "Cameroon") %>% 
                        rename(district = Admin2))

countries3 <- rbind(countries2, shp_gab %>% dplyr::select(COUNTRY, NAME_2, geometry) %>% 
                         rename(country = COUNTRY,
                                district = NAME_2))

countries3$country_id <- as.integer(factor(countries3$country))


# 1. Reproject countries to a metric CRS
countries_m <- st_transform(countries3, 3857)

# 2. Create bbox in the same CRS
bbox_m <- st_bbox(countries_m) %>% st_as_sfc()

# 3. Build 5 km grid 
grid_5km <- st_make_grid(
  bbox_m,
  cellsize = 5000,
  what = "centers"
)

grid_sf <- st_sf(geometry = grid_5km)

# 4. Intersect
grid_sf <- st_intersection(grid_sf, countries_m)

grid_vect <- vect(grid_sf %>% st_transform(4326))

grid_pop <- raster::merge(ben_pop, Cong_pop, cmr_pop, gab_pop)
grid_elevation <- raster::merge(elevation_ben, elevation_cmr, elevation_Cong, elevation_gab)
grid_NDVI <- raster::merge(ben_ndvi_stack, cmr_ndvi_stack, Cong_ndvi_stack, gab_ndvi_stack)

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

itn_est <- grid_ITN_usage_2021[["Data"]]

grid_IRS_usage_2021 <- crop(IRS_intervention_2021, grid_vect)
grid_IRS_usage_2021 <- mask(grid_IRS_usage_2021, grid_vect)

# 1. Zones: 5x5 km squares for survey sites and for grid cells
crs_m <- "+proj=laea +lat_0=2 +lon_0=20 +datum=WGS84 +units=m +no_defs"
pts_m   <- st_transform(grid_sf, crs_m)
zones_s <- st_buffer(pts_m, dist = 2500, endCapStyle = "SQUARE")   # square centred on each site
zones_s <- zones_s |> st_transform(4326)

grid_ITN_usage_2021_df <- exact_extract(itn_est,  zones_s, fun = "mean", progress = TRUE)
grid_IRS_usage_2021_df <- exact_extract(grid_IRS_usage_2021,  zones_s, fun = "mean", progress = TRUE)
grid_pop_density_2020_df <- exact_extract(grid_pop,  zones_s, fun = "mean", progress = TRUE)
grid_rainfall_2021_df    <- exact_extract(grid_rainfall_2021, zones_s, fun = "mean", progress = TRUE)
grid_temperature_2021_df <- exact_extract(grid_temperature_2021, zones_s, fun = "mean", progress = TRUE)
ndvi_vals_grid_monthly <- exact_extract(grid_NDVI, zones_s, fun = "mean", progress = TRUE)
ndvi_mean_grid <- apply(ndvi_vals_grid_monthly[,-1], 1, mean, na.rm = TRUE, weights = TRUE)
elevation_grid_df <- exact_extract(grid_elevation, zones_s, fun = "mean", progress = TRUE)
grid_land_cover_water_df <- exact_extract(grid_land_cover_water, zones_s, fun = "mean", progress = TRUE)
grid_land_cover_built_df <- exact_extract(grid_land_cover_built, zones_s, fun = "mean", progress = TRUE)
grid_land_cover_cropland_df <- exact_extract(grid_land_cover_cropland, zones_s, fun = "mean", progress = TRUE)
grid_land_cover_trees_df <- exact_extract(grid_land_cover_trees, zones_s, fun = "mean", progress = TRUE)
grid_land_cover_wetland_df <- exact_extract(grid_land_cover_wetland, zones_s, fun = "mean", progress = TRUE)

# Distance to water

# Data (all in EPSG:4326)
# ------------------------------------------------------------------
wkt <- st_as_text(st_as_sfc(st_bbox(countries3)))

rivers <- st_read("data/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp",
                  wkt_filter = wkt)
rivers <- rivers[rivers$ORD_STRA >= min_order, ] |> st_geometry()

lakes <- st_read("data/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp",
                 wkt_filter = wkt)
lakes_geom <- st_make_valid(st_geometry(lakes))
shores <- st_boundary(lakes_geom)             # turns each polygon into its outline (a line)
shores <- shores[!st_is_empty(shores)]
shores <- st_cast(shores, "MULTILINESTRING")  # makes all the outlines one geometry type

xyz_riv  <- densify_xyz(rivers, spacing_m)
xyz_lake <- densify_xyz(shores, spacing_m)


grid_sf$dist_to_water <- zone_dist(st_centroid(pts_m),  crs_m)

grid_sf <- grid_sf %>% 
  mutate(annual_rainfall_2021 = grid_rainfall_2021_df,
         annual_temperature_2021 = grid_temperature_2021_df,
         elevation               = elevation_grid_df,
         vegetation_index = ndvi_mean_grid,
         pop_density_2020 = grid_pop_density_2020_df,
         lc_cropland = grid_land_cover_cropland_df,
         lc_water = grid_land_cover_water_df,
         lc_wetland = grid_land_cover_wetland_df,
         lc_built = grid_land_cover_built_df,
         lc_trees = grid_land_cover_trees_df,
         annual_ITN_usage_2021 = grid_ITN_usage_2021_df,
         annual_IRS_usage_2021 = grid_IRS_usage_2021_df) 



numeric_cols <- names(grid_sf)[sapply(grid_sf, is.numeric)]

for (col in numeric_cols) {
  grid_sf <- fill_by_nearest(grid_sf, col)
}


st_write(grid_sf, "data/Africa_data/Country_level_data/grid_data_for_prediction_sf.gpkg",
         append = FALSE)

################################################################################
#                 case of grid cells at continental level                      #
################################################################################

shp_africa <- st_read("data/Africa_data/shp_africa/Africa_simplified.shp") |> st_transform(3857)


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
IRS_intervention_2021 <- rast("data/IRS_2000/2025_GBD2024_Africa_IRS_2021.tif")
africa_rainfall_2021 <- rast(path1)
africa_temperature_2021 <- rast(path2)
africa_pop <- rast("data/Africa_data/pop_density/gpw_v4_population_density_rev11_2020_30_sec.tif")
africa_elevation <- rast("data/Africa_data/elevation/wc2.1_30s_elev.tif")
africa_ndvi_stack<- rast(ndvi_files)

vect_grid <- vect(grid_sf_africa |> st_transform(4326))

africa_land_cover_water <- crop(land_cover_water, vect_grid)
africa_land_cover_water <- mask(africa_land_cover_water, vect_grid)

africa_land_cover_built <- crop(land_cover_builts, vect_grid)
africa_land_cover_built <- mask(africa_land_cover_built, vect_grid)

africa_land_cover_cropland <- crop(land_cover_cropland, vect_grid)
africa_land_cover_cropland <- mask(africa_land_cover_cropland, vect_grid)

africa_land_cover_trees <- crop(land_cover_trees, vect_grid)
africa_land_cover_trees <- mask(africa_land_cover_trees, vect_grid)

africa_land_cover_wetland <- crop(land_cover_wetland, vect_grid)
africa_land_cover_wetland <- mask(africa_land_cover_wetland, vect_grid)

africa_ITN_intervention_2021 <- crop(ITN_intervention_2021, vect_grid)
africa_ITN_intervention_2021 <- mask(africa_ITN_intervention_2021, vect_grid)

africa_IRS_intervention_2021 <- crop(IRS_intervention_2021, vect_grid)
africa_IRS_intervention_2021 <- mask(africa_IRS_intervention_2021, vect_grid)

africa_pop_2020 <- crop(africa_pop, vect_grid)
africa_pop_2020 <- mask(africa_pop_2020, vect_grid)

africa_rainfall_2021 <- crop(africa_rainfall_2021, vect_grid)
africa_rainfall_2021 <- mask(africa_rainfall_2021, vect_grid)

africa_temperature_2021 <- crop(africa_temperature_2021, vect_grid)
africa_temperature_2021 <- mask(africa_temperature_2021, vect_grid)

africa_ndvi_stack <- crop(africa_ndvi_stack, vect_grid)
africa_ndvi_stack <- mask(africa_ndvi_stack, vect_grid)

africa_elevation <- crop(africa_elevation, vect_grid)
africa_elevation <- mask(africa_elevation, vect_grid)

# 1. Zones: 5x5 km squares for survey sites and for grid cells

crs_m <- "+proj=laea +lat_0=2 +lon_0=20 +datum=WGS84 +units=m +no_defs"
pts_m   <- st_transform(grid_sf_africa, crs_m)
zones_s <- st_buffer(pts_m, dist = 2500, endCapStyle = "SQUARE")   # square centred on each site
zones_s <- zones_s |> st_transform(4326)

africa_ITN_intervention_2021_df <- exact_extract(africa_ITN_intervention_2021[[1]], 
                                                 zones_s, fun = "mean", progress = TRUE)
africa_IRS_intervention_2021_df <- exact_extract(africa_IRS_intervention_2021, 
                                                 zones_s, fun = "mean", progress = TRUE)
africa_pop_density_2020_df <- exact_extract(africa_pop_2020,  
                                            zones_s, fun = "mean", progress = TRUE)
africa_rainfall_2021_df    <- exact_extract(africa_rainfall_2021, 
                                            zones_s, fun = "mean", progress = TRUE)
africa_temperature_2021_df <- exact_extract(africa_temperature_2021, 
                                            zones_s, fun = "mean", progress = TRUE)

ndvi_vals_africa_monthly <- exact_extract(africa_ndvi_stack, 
                                          zones_s, fun = "mean", progress = TRUE)
ndvi_mean_africa <- apply(ndvi_vals_africa_monthly[,-1], 1, mean, na.rm = TRUE, weights = TRUE)

africa_elevation_df <- exact_extract(africa_elevation,
                                     zones_s, fun = "mean", progress = TRUE)
africa_land_cover_water_df <- exact_extract(africa_land_cover_water, 
                                            zones_s, fun = "mean", progress = TRUE)
africa_land_cover_built_df <- exact_extract(africa_land_cover_built,
                                             zones_s, fun = "mean", progress = TRUE)
africa_land_cover_cropland_df <- exact_extract(africa_land_cover_cropland,
                                               zones_s, fun = "mean", progress = TRUE)
africa_land_cover_trees_df <- exact_extract(africa_land_cover_trees,
                                            zones_s, fun = "mean", progress = TRUE)
africa_land_cover_wetland_df <- exact_extract(africa_land_cover_wetland,
                                              zones_s, fun = "mean", progress = TRUE)

# distance to water


# ------------------------------------------------------------------
# Data (all in EPSG:4326)
# ------------------------------------------------------------------
africa <- st_transform(shp_africa, 4326)
wkt <- st_as_text(st_as_sfc(st_bbox(africa)))

rivers <- st_read("data/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp",
                  wkt_filter = wkt)
rivers <- rivers[rivers$ORD_STRA >= min_order, ] |> st_geometry()

lakes <- st_read("data/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp",
                 wkt_filter = wkt)
lakes_geom <- st_make_valid(st_geometry(lakes))
shores <- st_boundary(lakes_geom)             # turns each polygon into its outline (a line)
shores <- shores[!st_is_empty(shores)]
shores <- st_cast(shores, "MULTILINESTRING")  # makes all the outlines one geometry type

xyz_riv  <- densify_xyz(rivers, spacing_m)
xyz_lake <- densify_xyz(shores, spacing_m)


grid_sf_africa$dist_to_water <- zone_dist(st_centroid(pts_m),  crs_m)


grid_sf_africa <- grid_sf_africa %>% 
  mutate(pop_density_2020 = africa_pop_density_2020_df,
         annual_rainfall_2021 = africa_rainfall_2021_df,
         annual_temperature_2021 = africa_temperature_2021_df,
         vegetation_index = ndvi_mean_africa,
         elevation = africa_elevation_df,
         lc_water = africa_land_cover_water_df,
         lc_built = africa_land_cover_built_df,
         lc_trees = africa_land_cover_trees_df,
         lc_wetland = africa_land_cover_wetland_df,
         lc_cropland = africa_land_cover_cropland_df,
         annual_ITN_usage_2021 = africa_ITN_intervention_2021_df,
         annual_IRS_usage_2021 = africa_IRS_intervention_2021_df)

# filling missing values
numeric_cols <- names(grid_sf_africa)[sapply(grid_sf_africa, is.numeric)]

for (col in numeric_cols) {
  grid_sf_africa <- fill_by_nearest(grid_sf_africa, col)
}

st_write(grid_sf_africa, "data/Africa_data/grid_data_for_prediction_sf.gpkg", append = FALSE)

################################################################################
#                   Visualization of the raster covariates                     #
################################################################################

################################################################################
#                       Load necessary libraries                               #
################################################################################

library(ggplot2)
library(terra)
library(sf)
library(tidyverse)
library(tidyterra)
library(RColorBrewer)

################################################################################
#                              loading the data                                #
################################################################################


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


precip_ratser_2021_data <- rast(precip_ratser_2021)
names(precip_ratser_2021_data) <- seq(as.Date("2022-01-01"),as.Date("2022-12-31"),
                                      by = "months") %>% 
  format("%b %Y")

precip_clean <- ifel(precip_ratser_2021_data < 0, NA, precip_ratser_2021_data)

precip_ratser_2021_data_annual <- sum(precip_clean, na.rm = TRUE)
names(precip_ratser_2021_data_annual) <- "annual_rainfall"

# global_spatial_limit_3857 <- project(global_spatial_limit, crs(shp_africa))

shp_africa <- read_sf("data/Africa_data/shp_africa/Africa_simplified.shp")

pf_limit_africa <- crop(global_spatial_limit, vect(shp_africa))
pf_limit_africa <- mask(pf_limit_africa, vect(shp_africa))

precip_ratser_2021_data_annual <- crop(precip_ratser_2021_data_annual, vect(shp_africa))
precip_ratser_2021_data_annual <- mask(precip_ratser_2021_data_annual, vect(shp_africa))

annual_temperature_climate  <- crop(annual_temperature_climate , vect(shp_africa))
annual_temperature_climate  <- mask(annual_temperature_climate , vect(shp_africa))
names(annual_temperature_climate) <- "temperature"

africa_water <- crop(land_cover_water, vect(shp_africa))
africa_water <- mask(africa_water, vect(shp_africa))

africa_trees <- crop(land_cover_trees, vect(shp_africa))
africa_trees <- mask(africa_trees, vect(shp_africa))

africa_builts <- crop(land_cover_builts, vect(shp_africa))
africa_builts <- mask(africa_builts, vect(shp_africa))

africa_cropland <- crop(land_cover_cropland, vect(shp_africa))
africa_cropland <- mask(africa_cropland, vect(shp_africa))

africa_wetland <- crop(land_cover_wetland, vect(shp_africa))
africa_wetland <- mask(africa_wetland, vect(shp_africa))

terra::plot(pf_limit_africa)
terra::plot(africa_wetland)
terra::plot(africa_cropland)
terra::plot(africa_builts)
terra::plot(africa_trees)
terra::plot(africa_water)



#asign a value for the class NA

pf_limit_africa2 <- classify(pf_limit_africa,
                             matrix(c(NA, NA, 0), ncol=3, byrow=TRUE))

lev <- data.frame(
  value = c(0, 1, 2),
  category = c("No transmission", "Unstable transmission", "Stable transmission")
)

levels(pf_limit_africa2) <- lev

pf_limit_africa2 <- crop(pf_limit_africa2, vect(shp_africa %>% st_transform(4326)))
pf_limit_africa2 <- mask(pf_limit_africa2, vect(shp_africa %>% st_transform(4326)))

terra::plot(pf_limit_africa2)

names(pf_limit_africa2) <- "Transmission level 2010 of Pf"

cats(pf_limit_africa2)

shp_africa_country <- shp_africa %>% 
  st_make_valid() %>% 
  group_by(country) %>% 
  summarise(geometry = st_union(geometry))

pf_limit_2010 <- ggplot(shp_africa_country)+
  geom_spatraster(data = pf_limit_africa2, mapping = aes(fill = `Transmission level 2010 of Pf`))+
  scale_fill_manual(values = c(
    "No transmission" = "grey80",
    "Unstable transmission" = "orange",
    "Stable transmission" = "red"), na.value =  "transparent", na.translate = FALSE)+
  geom_sf(fill = NA)+
  theme_void()+
  labs(title = "Limit of transmission of P. falciparum in Africa (2010)")

ggsave("Outputs/africa_output/limit_trans_pf.png", 
       plot = pf_limit_2010, width = 6, height = 4, dpi = 300)

annual_rainfall <- ggplot(shp_africa_country)+
  geom_spatraster(data = precip_ratser_2021_data_annual, mapping = aes(fill = annual_rainfall))+
  geom_sf(fill = NA)+
  scale_fill_stepsn(
    colours = brewer.pal(9, "Blues"),
    breaks = c(0, 250, 500, 1000, 1500, 2000, 3000, 5000),
    na.value = "transparent",
    name = "Annual rainfall (mm)")+
  theme_void()+
  labs(title = "annual rainfall in 2021")

ggsave("Outputs/africa_output/annual_rainfall_2021_chirps.png", 
       plot = annual_rainfall, width = 6, height = 4, dpi = 300)


annual_temperature_2021 <- ggplot(shp_africa_country)+
  geom_spatraster(data = annual_temperature_climate, mapping = aes(fill = temperature))+
  geom_sf(fill = NA)+
  scale_fill_stepsn(
    colours = rev(brewer.pal(9, "RdYlGn")),
    breaks = c(0, 10, 15, 20, 25, 30, 35),
    na.value = "transparent",
    name = "Annual mean temp (°C)")+
  theme_void()+
  labs(title = "annual temperature in 2021")

ggsave("Outputs/africa_output/annual_temperature_2021_climate.png", 
       plot = annual_temperature_2021, width = 6, height = 4, dpi = 300)


wetland_2021 <- ggplot(shp_africa_country)+
  geom_spatraster(data = africa_wetland, mapping = aes(fill = wetland))+
  geom_sf(fill = NA)+
  scale_fill_gradientn(
    colours = viridis::magma(10),
    trans = "log10",
    limits = c(1e-6, 1),
    oob = scales::squish,   # 0 to 1e-6
    na.value = "transparent",
    name = "Wetland fraction (log10)")+
  theme_void()+
  labs(title = "Wetland in 2021")

ggsave("Outputs/africa_output/wetland_2021.png", 
       plot = wetland_2021, width = 6, height = 4, dpi = 300)


africa_water_2021 <- ggplot(shp_africa_country)+
  geom_spatraster(data = africa_water , mapping = aes(fill = water))+
  geom_sf(fill = NA)+
  scale_fill_gradientn(
    colours = viridis::plasma(10),
    trans = "log10",
    limits = c(1e-6, 1),
    oob = scales::squish,
    na.value = "transparent",
    name = "Water fraction (log10)")+
  theme_void()+
  labs(title = "Water in 2021")

ggsave("Outputs/africa_output/africa_water_2021.png", 
       plot = africa_water_2021, width = 6, height = 4, dpi = 300)

africa_trees_2021 <- ggplot(shp_africa_country)+
  geom_spatraster(data = africa_trees , mapping = aes(fill = trees))+
  geom_sf(fill = NA)+
  scale_fill_gradientn(
    colours = viridis::cividis(10),
    trans = "log10",
    limits = c(1e-6, 1),
    oob = scales::squish,
    na.value = "transparent",
    name = "Tree fraction (log10)")+
  theme_void()+
  labs(title = "Trees in 2021")

ggsave("Outputs/africa_output/africa_trees_2021.png", 
       plot = africa_trees_2021, width = 6, height = 4, dpi = 300)

africa_builts_2021 <- ggplot(shp_africa_country)+
  geom_spatraster(data = africa_builts , mapping = aes(fill = built))+
  geom_sf(fill = NA)+
  scale_fill_gradientn(
    colours = viridis::rocket(10),
    trans = "log10",
    limits = c(1e-6, 1),
    oob = scales::squish,
    na.value = "transparent",
    name = "Built-up fraction (log10)")+
  theme_void()+
  labs(title = "Built-up in 2021")

ggsave("Outputs/africa_output/africa_builts_2021.png", 
       plot = africa_builts_2021, width = 6, height = 4, dpi = 300)

africa_cropland_2021 <- ggplot(shp_africa_country)+
  geom_spatraster(data = africa_cropland , mapping = aes(fill = cropland))+
  geom_sf(fill = NA)+
  scale_fill_gradientn(
    colours = viridis::inferno(10),
    trans = "log10",
    limits = c(1e-6, 1),
    oob = scales::squish,
    na.value = "transparent",
    name = "Cropland fraction (log10)")+
  theme_void()+
  labs(title = "Cropland in 2021")

ggsave("Outputs/africa_output/africa_cropland_2021.png", 
       plot = africa_cropland_2021, width = 6, height = 4, dpi = 300)

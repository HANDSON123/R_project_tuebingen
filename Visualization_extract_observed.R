################################################################################
#     Visualization of the extracted covariates at the sampled locations       # 
################################################################################

################################################################################
#                         Load necessary libraries                             #
################################################################################

library(sf)
library(reshape2)
library(ggcorrplot)
library(leaflet)
library(ggplot2)
library(RColorBrewer)
library(leafpop)
library(readr)
library(tidyverse)
library(MODISTools)
library(tidyterra)
library(mapview)
library(spData)
library(gganimate)
library(viridis)
library(tidytext)

################################################################################
#                            Loading the data                                  #  
################################################################################

spatial_data_ben_sf <- st_read("data/ben_data/spatial_data_ben_sf.gpkg")
spatial_data_cong   <- st_read("data/cong_data/spatial_data_cong_sf.gpkg")
spatial_data_cmr    <- st_read("data/cmr_data/cmr_cluster_df_1_sf.gpkg")
spatial_data_gab    <- st_read("data/gab_data/GAB_df_1_sf.gpkg") 

#loading the shape files
shp_ben <- read_sf("data/ben_data/BEN_shape_file/gadm41_BEN_3.shp")
shp_con <- read_sf("data/cong_data/CON_shape_file/cog_admin2.shp")
shp_cmr <- read_sf("data/cmr_data/CMR_shape_file/cameroun_districts_2024.shp")
shp_gab <- read_sf("data/gab_data/gab_shape/gadm41_GAB_2.shp")


#renaming the districts

shp_ben$NAME_3[shp_ben$NAME_3=="Gakpè"]="Gakpe"
shp_ben$NAME_3[shp_ben$NAME_3=="Tokpa Domé"]="Tokpa-dome"
shp_ben$NAME_3[shp_ben$NAME_3=="Ouakpé-Daho"]="Houakpe-daho"
shp_ben$NAME_3[shp_ben$NAME_3=="Aganmalomé"]="Aganmanlome"
shp_ben$NAME_3[shp_ben$NAME_3=="Agonkanmè"]="Agonkamey"

sampled_site_cmr <- shp_cmr %>% 
  filter(Admin2 == "Bankim")

sampled_site_gab <- shp_gab %>% 
  filter(NAME_2 %in% c("Komo", "Ogooué et des Lacs"))

sampled_site_ben <- shp_ben %>% 
  filter(NAME_3 %in% c("Tokpa-dome", "Aganmanlome", "Savi", "Kpomassè", "Gakpe", 
                       "Arrondissement III", "Arrondissement I", "Houakpe-daho",
                       "Agonkamey", "Arrondissement II") & NAME_2 %in% c("Kpomassè", "Ouidah"))

sampled_site_cong <- shp_con %>% 
  filter(adm2_name %in% c("Brazzaville", "Goma-Tsetse"))


################################## Benin #######################################

dist_water <- ggplot() +
  geom_sf(data = sampled_site_ben, fill = "white", color = "black") +
  geom_sf(data = spatial_data_ben_sf, aes(color = dist_to_water/1000)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Distance to water (km) in sampled site (Benin)",
       color = "Distance to water (km)") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/ben_output/dist_water_viz_ben.pdf", plot = dist_water, width = 6, height = 4, dpi = 300)

pop_density <- ggplot() +
  geom_sf(data = sampled_site_ben, fill = "white", color = "black") +
  geom_sf(data = spatial_data_ben_sf, aes(color = pop_density_2020)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "population density in sampled site (Benin)",
       color = "population density") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/ben_output/pop_density_viz_ben.pdf", plot = pop_density, width = 6, height = 4, dpi = 300)

annual_rainfall <- ggplot() +
  geom_sf(data = sampled_site_ben, fill = "white", color = "black") +
  geom_sf(data = spatial_data_ben_sf, aes(color = annual_rainfall_2021)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Annual rainfall in sampled site (Benin)",
       color = "Annual rainfall") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/ben_output/annual_rainfall_viz_ben.pdf", plot = annual_rainfall, width = 6, height = 4, dpi = 300)

annual_temperature <- ggplot() +
  geom_sf(data = sampled_site_ben, fill = "white", color = "black") +
  geom_sf(data = spatial_data_ben_sf, aes(color = annual_temperature_2021)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Annual temperature in sampled site (Benin)",
       color = "Annual temperature") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/ben_output/annual_temperature_viz_ben.pdf", plot = annual_temperature, width = 6, height = 4, dpi = 300)


Vegetation_index <- ggplot() +
  geom_sf(data = sampled_site_ben, fill = "white", color = "black") +
  geom_sf(data = spatial_data_ben_sf, aes(color = vegetaion_index)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Vegetation index in sampled site (Benin)",
       color = "Vegetation index") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/ben_output/Vegetation_index_viz_ben.pdf", plot = Vegetation_index, width = 6, height = 4, dpi = 300)

elevation_ben <- ggplot() +
  geom_sf(data = sampled_site_ben, fill = "white", color = "black") +
  geom_sf(data = spatial_data_ben_sf, aes(color = elevation)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "elevation in sampled site (Benin)",
       color = "elevation") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/ben_output/elevation_viz_ben.pdf", plot = elevation_ben, width = 6, height = 4, dpi = 300)


################################## Congo #######################################

dist_water <- ggplot() +
  geom_sf(data = sampled_site_cong, fill = "white", color = "black") +
  geom_sf(data = spatial_data_cong, aes(color = dist_to_water/1000)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Distance to water (km) in sampled site (Congo)",
       color = "Distance to water (km)") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/cong_output/dist_water_viz_cong.pdf", plot = dist_water, width = 6, height = 4, dpi = 300)

pop_density <- ggplot() +
  geom_sf(data = sampled_site_cong, fill = "white", color = "black") +
  geom_sf(data = spatial_data_cong, aes(color = pop_density_2020)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "population density in sampled site (Congo)",
       color = "population density") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/cong_output/pop_density_viz_cong.pdf", plot = pop_density, width = 6, height = 4, dpi = 300)

annual_rainfall <- ggplot() +
  geom_sf(data = sampled_site_cong, fill = "white", color = "black") +
  geom_sf(data = spatial_data_cong, aes(color = annual_rainfall_2021)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Annual rainfall in sampled site (Congo)",
       color = "Annual rainfall") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/cong_output/annual_rainfall_viz_cong.pdf", plot = annual_rainfall, width = 6, height = 4, dpi = 300)

annual_temperature <- ggplot() +
  geom_sf(data = sampled_site_cong, fill = "white", color = "black") +
  geom_sf(data = spatial_data_cong, aes(color = annual_temperature_2021)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Annual temperature in sampled site (Congo)",
       color = "Annual temperature") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/cong_output/annual_temperature_vizcong.pdf", plot = annual_temperature, width = 6, height = 4, dpi = 300)


Vegetation_index <- ggplot() +
  geom_sf(data = sampled_site_cong, fill = "white", color = "black") +
  geom_sf(data = spatial_data_cong, aes(color = vegetaion_index)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Vegetation index in sampled site (Congo)",
       color = "Vegetation index") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/cong_output/Vegetation_index_viz_cong.pdf", plot = Vegetation_index, width = 6, height = 4, dpi = 300)


elevation_cong <- ggplot() +
  geom_sf(data = sampled_site_cong, fill = "white", color = "black") +
  geom_sf(data = spatial_data_cong, aes(color = elevation)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "elevation in sampled site (Congo)",
       color = "elevation") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/cong_output/elevation_viz_cong.pdf",
       plot = elevation_cong, width = 6, height = 4, dpi = 300)

############################## Cameroon ########################################

dist_water <- ggplot() +
  geom_sf(data = sampled_site_cmr, fill = "white", color = "black") +
  geom_sf(data = spatial_data_cmr, aes(color = dist_to_water/1000)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Distance to water (km) in sampled site (Cameroon)",
       color = "Distance to water (km)") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/cmr_output/dist_water_viz_cmr.pdf",
       plot = dist_water, width = 6, height = 4, dpi = 300)

pop_density <- ggplot() +
  geom_sf(data = sampled_site_cmr, fill = "white", color = "black") +
  geom_sf(data = spatial_data_cmr, aes(color = pop_density_2020)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "population density in sampled site (Cameroon)",
       color = "population density") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/cmr_output/pop_density_viz_cmr.pdf",
       plot = pop_density, width = 6, height = 4, dpi = 300)

annual_rainfall <- ggplot() +
  geom_sf(data = sampled_site_cmr, fill = "white", color = "black") +
  geom_sf(data = spatial_data_cmr, aes(color = annual_rainfall_2021)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Annual rainfall in sampled site (Cameroon)",
       color = "Annual rainfall") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/cmr_output/annual_rainfall_viz_cmr.pdf", 
       plot = annual_rainfall, width = 6, height = 4, dpi = 300)

annual_temperature <- ggplot() +
  geom_sf(data = sampled_site_cmr, fill = "white", color = "black") +
  geom_sf(data = spatial_data_cmr, aes(color = annual_temperature_2021)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Annual temperature in sampled site (Cameroon)",
       color = "Annual temperature") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/cmr_output/annual_temperature_viz_cmr.pdf",
       plot = annual_temperature, width = 6, height = 4, dpi = 300)


Vegetation_index <- ggplot() +
  geom_sf(data = sampled_site_cmr, fill = "white", color = "black") +
  geom_sf(data = spatial_data_cmr, aes(color = vegetaion_index)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Vegetation index in sampled site (Cameroon)",
       color = "Vegetation index") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/cmr_output/Vegetation_index_viz_cmr.pdf",
       plot = Vegetation_index, width = 6, height = 4, dpi = 300)


elevation_cmr <- ggplot() +
  geom_sf(data = sampled_site_cmr, fill = "white", color = "black") +
  geom_sf(data = spatial_data_cmr, aes(color = elevation)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "elevation in sampled site (Cameroon)",
       color = "elevation") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/cmr_output/elevation_viz_cmr.pdf",
       plot = elevation_cmr, width = 6, height = 4, dpi = 300)

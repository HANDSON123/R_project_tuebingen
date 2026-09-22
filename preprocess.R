################################################################################
#                             pre-process script                               #
################################################################################

################################################################################ 
#                      Loading necessary libraries                             #
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
library(readxl)
library(corrplot)
library(MODISTools)
library(terra)
library(tidyterra)
library(mapview)
library(spData)
library(gganimate)
library(htmlwidgets)
library(webshot)
library(chirps)
library(viridis)
library(patchwork)
library(GGally)
library(car)
library(gstat)
library(nngeo)
library(tidygeocoder)
library(malariaAtlas)
library(rmapshaper)
library(tidytext)

################################################################################
#                               loading data                                   #
################################################################################


BEN_df <- read_excel("data/COMAL data Benin&Congo.xlsx", sheet = 1)

CON_df <- read_excel("data/COMAL data Benin&Congo.xlsx", sheet = 2)

CMR_df <- read_excel("data/cmr_data/Cameroon2.xlsx")

GAB_df <- read_excel("data/gab_data/Gabon2_GPS.xlsx")

shp_ben <- read_sf("data/ben_data/BEN_shape_file/gadm41_BEN_3.shp")
shp_con <- read_sf("data/cong_data/CON_shape_file/cog_admin2.shp")
shp_cmr <- read_sf("data/cmr_data/CMR_shape_file/cameroun_districts_2024.shp")
shp_gab <- read_sf("data/gab_data/gab_shape/gadm41_GAB_2.shp")


############################ renaming the districts ############################

shp_ben$NAME_3[shp_ben$NAME_3=="Gakpè"]="Gakpe"
shp_ben$NAME_3[shp_ben$NAME_3=="Tokpa Domé"]="Tokpa-dome"
shp_ben$NAME_3[shp_ben$NAME_3=="Ouakpé-Daho"]="Houakpe-daho"
shp_ben$NAME_3[shp_ben$NAME_3=="Aganmalomé"]="Aganmanlome"
shp_ben$NAME_3[shp_ben$NAME_3=="Agonkanmè"]="Agonkamey"

ben_vect <- vect(shp_ben)
cmr_vect <- vect(shp_cmr)
con_vect <- vect(shp_con)
gab_vect <- vect(shp_gab)

################################################################################
#                            data cleaning                                     #
################################################################################


################################ Cameroon ######################################

## function to extract the GPS location for Cameroon data

convert_custom_dms <- function(x) {
  hemi <- substr(x, 1, 1)
  x <- substring(x, 2)
  
  apos <- regexpr("'", x)[1]
  
  if (hemi == "N") {
    deg <- as.numeric(substr(x, 1, apos - 4))
  } else {
    deg <- as.numeric(substr(x, 1, apos - 3))
  }
  
  min <- as.numeric(substr(x, apos - 2, apos - 1))
  sec <- as.numeric(gsub("\"", "", substr(x, apos + 1, nchar(x))))
  
  dec <- deg + min/60 + sec/3600
  
  if (hemi %in% c("S", "W")) dec <- -dec
  return(dec)
}

CMR_df$lat <- sapply(strsplit(CMR_df$gps, " "), function(x) convert_custom_dms(x[1]))
CMR_df$lon <- sapply(strsplit(CMR_df$gps, " "), function(x) convert_custom_dms(x[2]))


CMR_df_1 <- CMR_df %>% 
  dplyr::select(ID, lon, lat, gender, age, microscopy, pcr_pm, pcr_pf) %>% 
  mutate(
    lng = lon,
    tot_participants = 1,
    tot_infected_pm  = if_else(pcr_pm == "pos", 1, 0),
    tot_infected_pf  = if_else(pcr_pf == "pos", 1, 0),
    microscopy_pm    = case_when(microscopy == "negative" ~ "neg",
                                 microscopy %in% c("Pm monoinfection", 
                                                   "Pf-Pm coinfection") ~ "pos"),
    microscopy_pf    = case_when(microscopy == "negative" ~ "neg",
                                 microscopy %in% c("Pf monoinfection", 
                                                   "Pf-Pm coinfection",
                                                   "Pf-Po coinfection") ~ "pos"),
    tot_infected_microscopy_pm = if_else(microscopy_pm == "pos", 1, 0),
    tot_infected_microscopy_pf = if_else(microscopy_pf == "pos", 1, 0)
  )%>% 
  group_by(lng, lat) %>% 
  summarise(age_lower                            = min(age),
            age_upper                            = max(age),
            Number_of_participant_tot            = sum(tot_participants, na.rm = TRUE),
            Number_of_infected_tot_pm            = sum(tot_infected_pm, na.rm = TRUE),
            Number_of_infected_microscopy_tot_pm = sum(tot_infected_microscopy_pm,
                                                       na.rm = TRUE),
            PR_Pm_tot                            = Number_of_infected_tot_pm  / Number_of_participant_tot ,
            PR_Pm_tot_microscopy                 = Number_of_infected_microscopy_tot_pm/Number_of_participant_tot ,
            Number_of_infected_tot_pf            = sum(tot_infected_pf, na.rm = TRUE),
            Number_of_infected_microscopy_tot_pf = sum(tot_infected_microscopy_pf,
                                                       na.rm = TRUE),
            PR_Pf_tot                            = Number_of_infected_tot_pf  / Number_of_participant_tot ,
            PR_Pf_tot_microscopy                 = Number_of_infected_microscopy_tot_pf/Number_of_participant_tot ,
            .groups = "drop") %>% 
  mutate(country = "Cameroon")

# Group households into clusters

households <- st_as_sf(
  CMR_df_1,
  coords = c("lng","lat"),
  crs = 4326
)

households <- st_transform(households, 32632)

grid_cmr_cluster <- st_make_grid(
  households,
  cellsize = 1000,
  square = TRUE
)

grid_cmr_cluster <- st_sf(
  grid_id = seq_along(grid_cmr_cluster),
  geometry = grid_cmr_cluster
)

households_grid <- st_join(
  households,
  grid_cmr_cluster,
  join = st_intersects
)

cluster_data_cmr <-
  households_grid %>%
  st_drop_geometry() %>% 
  group_by(grid_id) %>%
  summarise(age_lower = min(age_lower),
            age_upper = max(age_upper),
            Number_of_participant_tot            = sum(Number_of_participant_tot, na.rm = TRUE),
            Number_of_infected_tot_pm            = sum(Number_of_infected_tot_pm, na.rm = TRUE),
            Number_of_infected_microscopy_tot_pm = sum(Number_of_infected_microscopy_tot_pm,
                                                       na.rm = TRUE),
            PR_Pm_tot                            = Number_of_infected_tot_pm  / Number_of_participant_tot ,
            PR_Pm_tot_microscopy                 = Number_of_infected_microscopy_tot_pm/Number_of_participant_tot ,
            Number_of_infected_tot_pf            = sum(Number_of_infected_tot_pf, na.rm = TRUE),
            Number_of_infected_microscopy_tot_pf = sum(Number_of_infected_microscopy_tot_pf,
                                                       na.rm = TRUE),
            PR_Pf_tot                            = Number_of_infected_tot_pf  / Number_of_participant_tot ,
            PR_Pf_tot_microscopy                 = Number_of_infected_microscopy_tot_pf/Number_of_participant_tot ,
            .groups = "drop") %>% 
  mutate(country = "Cameroon")


cluster_data_cmr_sf <-
  grid_cmr_cluster  %>% 
  inner_join(cluster_data_cmr, by = "grid_id")

cluster_data_cmr_sf <- st_transform(cluster_data_cmr_sf, 4326)
cluster_data_cmr_sf <- st_centroid(cluster_data_cmr_sf)

coords_cmr <- st_coordinates(cluster_data_cmr_sf)

cluster_data_cmr_sf$lng <- coords_cmr[,1]
cluster_data_cmr_sf$lat <- coords_cmr[,2]



CMR_df_1_sf <- st_as_sf(CMR_df_1, coords = c("lng", "lat"), crs = 4326)
st_write(CMR_df_1_sf, "data/cmr_data/spatial_data_sample_cmr.gpkg", append = FALSE)

cluster_data_cmr_sf <- st_as_sf(cluster_data_cmr_sf, coords = c("lng", "lat"), crs = 4326)
st_write(cluster_data_cmr_sf, "data/cmr_data/spatial_data_sample_cmr_cluster.gpkg", append = FALSE)



raw_pr_cmr <- ggplot() +
  geom_sf(data = shp_cmr, fill = "white", color = "black") +
  geom_point(data = cluster_data_cmr_sf, aes(x = lng, y = lat, 
                                             size = Number_of_participant_tot,
                                             color = PR_Pm_tot)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Raw P. malariae prevalence",
       color = "Prevalence", size = "N tested") +
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/cmr_output/raw_obs_PR_cmr.pdf", plot = raw_pr_cmr, width = 6, height = 4, dpi = 300)

sampled_site_cmr <- shp_cmr %>% 
  filter(Admin2 == "Bankim")

sample_district <- ggplot() +
  geom_sf(data = sampled_site_cmr, fill = "white", color = "black") +
  geom_point(data = cluster_data_cmr_sf, aes(x = lng, y = lat, 
                                             size = Number_of_participant_tot, 
                                             color = PR_Pm_tot)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Raw P. malariae prevalence (Bankim)",
       color = "Prevalence", size = "N tested") +
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/cmr_output/raw_obs_PR_sampled_location_cmr.pdf", 
       plot = sample_district, width = 6, height = 4, dpi = 300)


mapview(shp_cmr,
        col.regions = "white",
        alpha.regions = 0.2,
        color = "black") +
  mapview(CMR_df_1_sf,
          zcol = "PR_Pm_tot",
          cex = "Number_tot_participants",
          col.regions = colorRampPalette(c("blue", "orange")),
          at = NULL,
          layer.name = "P. malariae prevalence")


################################## Gabon #######################################

GAB_df_1 <- GAB_df %>% 
  dplyr::select(ID, long, lat, gender, age, pcr_pm, pcr_pf, microscopy) %>% 
  mutate(
    lng = long,
    tot_participants = 1,
    tot_infected_pm  = if_else(pcr_pm == "pos", 1, 0),
    tot_infected_pf  = if_else(pcr_pf == "pos", 1, 0),
    microscopy_pm    = case_when(microscopy == "negative" ~ "neg",
                                 microscopy %in% c("Pm monoinfection", 
                                                   "Pf-Pm coinfection") ~ "pos"),
    microscopy_pf    = case_when(microscopy == "negative" ~ "neg",
                                 microscopy %in% c("Pf monoinfection", 
                                                   "Pf-Pm coinfection",
                                                   "Pf-Po coinfection") ~ "pos"),
    tot_infected_microscopy_pm = if_else(microscopy_pm == "pos", 1, 0),
    tot_infected_microscopy_pf = if_else(microscopy_pf == "pos", 1, 0)
  )%>% 
  group_by(lng, lat) %>% 
  summarise(age_lower                            = min(age),
            age_upper                            = max(age),
            Number_of_participant_tot            = sum(tot_participants, na.rm = TRUE),
            Number_of_infected_tot_pm            = sum(tot_infected_pm, na.rm = TRUE),
            Number_of_infected_microscopy_tot_pm = sum(tot_infected_microscopy_pm,
                                                       na.rm = TRUE),
            PR_Pm_tot                            = Number_of_infected_tot_pm  / Number_of_participant_tot ,
            PR_Pm_tot_microscopy                 = Number_of_infected_microscopy_tot_pm/Number_of_participant_tot ,
            Number_of_infected_tot_pf            = sum(tot_infected_pf, na.rm = TRUE),
            Number_of_infected_microscopy_tot_pf = sum(tot_infected_microscopy_pf,
                                                       na.rm = TRUE),
            PR_Pf_tot                            = Number_of_infected_tot_pf  / Number_of_participant_tot ,
            PR_Pf_tot_microscopy                 = Number_of_infected_microscopy_tot_pf/Number_of_participant_tot ,
            .groups = "drop") %>% 
  mutate(country = "Gabon")



GAB_df_1 <- st_as_sf(GAB_df_1, coords = c("lng", "lat"), crs = 4326)
GAB_df_1_sf <- st_transform(GAB_df_1, st_crs(shp_gab))
st_write(GAB_df_1_sf, "data/gab_data/spatial_data_sample_gab.gpkg", append = FALSE)

raw_pr_gab <- ggplot() +
  geom_sf(data = shp_gab, fill = "white", color = "black") +
  geom_sf(data = GAB_df_1,
          aes(size = Number_of_participant_tot, color = PR_Pm_tot)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Raw P. malariae prevalence",
       color = "Prevalence", size = "N tested") +
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/gab_output/raw_obs_PR_gab.pdf", plot = raw_pr_gab, width = 6, height = 4, dpi = 300)

sampled_site_gab <- shp_gab %>% 
  filter(NAME_2 %in% c("Komo", "Ogooué et des Lacs"))

sample_district_gab <- ggplot() +
  geom_sf(data = sampled_site_gab, fill = "white", color = "black") +
  geom_sf(data = GAB_df_1, aes(size = Number_of_participant_tot, color = PR_Pm_tot)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Raw P. malariae prevalence ",
       color = "Prevalence", size = "N tested") +
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/gab_output/raw_obs_PR_sampled_location_gab.pdf",
       plot = sample_district_gab, width = 6, height = 4, dpi = 300)


mapview(shp_gab,
        col.regions = "white",
        alpha.regions = 0.2,
        color = "black") +
  mapview(GAB_df_1_sf,
          zcol = "PR_Pm_tot",
          cex = "Number_of_participant_tot",
          col.regions = colorRampPalette(c("blue", "orange")),
          at = NULL)

################################ Benin #########################################

village_coords_ben <- tribble(
  ~Municipality, ~District,     ~Village,            ~Match_confidence,   ~lat,     ~lng,
  "Ouidah",      "Gakpe",       "Lokohoue",           "fuzzy_match",      6.443353,  2.132558,#
  "Ouidah",      "Gakpe",       "Kindjitokpa",        "fuzzy_match",      6.461650,  2.148980,#
  "Kpomasse",    "Tokpa-dome",  "Lokogbo-Gnonwa",     "exact_match",      6.471386,  1.984651,#
  "Ouidah",      "Houakpe-daho","Gbezoume",           "exact_match",      6.344216,  1.962808,#
  "Ouidah",      "Gakpe",       "Sehlinhoue",         "fuzzy_match",      6.416572,  2.173980,#
  "Kpomasse",    "Aganmanlome", "Hessa",              "fuzzy_match",      6.495617,  2.053485,#
  "Kpomasse",    "Aganmanlome", "KOUZOUME",           "fuzzy_match",      6.462216,  2.017743,#
  "Ouidah",      "Gakpe",       "Agossouhouihoué",    "fuzzy_match",      6.421688,  2.144019,#
  "Ouidah",      "Gakpe",       "Ahossihoue",         "fuzzy_match",      6.452980,  2.167627,#
  "Ouidah",      "Houakpe-daho","Gbehonou",           "exact_match",      6.339606,  2.005831,#
  "Kpomasse",    "Aganmanlome", "KOUGBEDJI",          "fuzzy_match",      6.473206,  2.077346,#
  "Ouidah",      "Houakpe-daho","Toligbé",            "exact_match",      6.327061,  2.027671,#
  "Ouidah",      "Houakpe-daho","Houakpe-daho centre","exact_match",      6.329397,  2.047797,#
  "Kpomasse",    "Aganmanlome", "Aganmanlome Centre", "fuzzy_match",      6.477114,  2.049987,#
  "Ouidah",      "Gakpe",       "TOKOLI",             "fuzzy_match",      6.431493,  2.156038,#
  "Kpomasse",    "Aganmanlome", "Lokossa",            "fuzzy_match",      6.478105,  2.027037,#
  "Kpomasse",    "Agonkamey",   "Assogbenou Kpevi",   "exact_match",      6.366429,  2.027881,#
  "Ouidah",      "Houakpe-daho","Seyigbe",            "exact_match",      6.324844,  2.041879,#
  "Kpomasse",    "Agonkamey",   "Assogbenou daho",    "exact_match",      6.355088,  2.001053 #
)

BEN_1 <- BEN_df %>% 
  mutate(Village = case_when(Village == "HESSA" ~ "Hessa",
                             TRUE ~ Village)) %>% 
  dplyr::select(Village, Sex, `Age (year)`, `qRT-qPCR Pm positive`, `qRT-qPCR Pf positive`,
                `Microscopic Malaria infection status (P.f, P.m, P.o)`) %>%
  rename(age        = `Age (year)`,
         microscopy = `Microscopic Malaria infection status (P.f, P.m, P.o)`) %>% 
  mutate(Pm_positive = case_when( `qRT-qPCR Pm positive`=="x"~"yes",
                                  TRUE~"no"),
         Pf_positive = case_when( `qRT-qPCR Pf positive`=="x"~"yes",
                                  TRUE~"no"),
         participants = 1,
         infected_pm  = case_when(Pm_positive=="yes"~1,
                                  TRUE~0),
         infected_pf  = case_when(Pf_positive=="yes"~1,
                                  TRUE~0),
         infected_microscopy_pm  = case_when(microscopy %in% c("Pf/Pm", "Pm") ~ 1,
                                             TRUE~0),
         infected_microscopy_pf  = case_when(microscopy %in% c("Pf", "Pf/Pm", "Pf/Po")~1,
                                             TRUE~0)) %>% 
  group_by(Village) %>% 
  summarise(age_lower                         = min(age, na.rm = TRUE),
            age_upper                            = max(age, na.rm = TRUE),
            Number_of_participant_tot            = sum(participants, na.rm = TRUE),
            Number_of_infected_tot_pm            = sum(infected_pm, na.rm = TRUE),
            Number_of_infected_microscopy_tot_pm = sum(infected_microscopy_pm,
                                                       na.rm = TRUE),
            PR_Pm_tot                            = Number_of_infected_tot_pm/Number_of_participant_tot,
            PR_Pm_tot_microscopy                 = Number_of_infected_microscopy_tot_pm/Number_of_participant_tot,
            Number_of_infected_tot_pf            = sum(infected_pf, na.rm = TRUE),
            Number_of_infected_microscopy_tot_pf = sum(infected_microscopy_pf,
                                                       na.rm = TRUE),
            PR_Pf_tot                            = Number_of_infected_tot_pf/Number_of_participant_tot,
            PR_Pf_tot_microscopy                 = Number_of_infected_microscopy_tot_pf/Number_of_participant_tot,
            .groups = "drop") %>%  
  distinct()



BEN_1_GPS <- BEN_1 %>% 
  left_join(village_coords_ben %>% dplyr::select(Village, lat, lng, Match_confidence))

spatial_data_ben <- st_as_sf(BEN_1_GPS, coords = c("lng", "lat"), crs = 4326)
st_write(spatial_data_ben, "data/ben_data/spatial_data_sample_ben.gpkg", append = FALSE)

PR_obs_ben <- ggplot() +
  geom_sf(data = shp_ben, fill = "white", color = "black") +
  geom_point(data = BEN_1_GPS, aes(x = lng, y = lat, size = Number_of_participant_tot, color = PR_Pm_tot)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Raw P. malariae prevalence",
       color = "Prevalence", size = "N tested") +
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/ben_output/PR_p_malariae_observed.pdf", plot = PR_obs_ben, width = 6, height = 4, dpi = 300)


sampled_site_ben <- shp_ben %>% 
  filter(NAME_3 %in% c("Tokpa-dome", "Aganmanlome", "Savi", "Kpomassè", "Gakpe", 
                       "Arrondissement III", "Arrondissement I", "Houakpe-daho",
                       "Agonkamey", "Arrondissement II") & NAME_2 %in% c("Kpomassè", "Ouidah"))

Raw_PR_sampled_location_ben <- ggplot() +
  geom_sf(data = sampled_site_ben, fill = "white", color = "black") +
  geom_point(data = BEN_1_GPS, aes(x = lng, y = lat, size = Number_of_participant_tot, color = PR_Pm_tot)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Raw P. malariae prevalence",
       color = "Prevalence", size = "N tested") +
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/ben_output/Raw_PR_sampled_location.pdf", plot = Raw_PR_sampled_location_ben,
       width = 6, height = 4, dpi = 300)


mapview(shp_ben,
        col.regions = "white",
        alpha.regions = 0.2,
        color = "black") +
  mapview(spatial_data_ben,
          zcol = "PR_Pm_tot",
          cex = "Number_of_participant_tot",
          col.regions = colorRampPalette(c("blue", "orange")),
          at = NULL,
          layer.name = "P. malariae prevalence")


################################## Congo #######################################

village_coords_cong <- tribble(
  ~District,     ~Match_confidence,     ~lat,                  ~lng,
  "Ntoula",      "fuzzy_match",      -4.36,                 15.15,#
  "Djoumouna",   "exact_match",      -4.37601066518807,     15.15963007611599,#
  "Mayanga",     "exact_match",      -4.284712200560734,    15.186363375868131,#
)

CONG_1 <- CON_df %>% 
  dplyr::select(District, Sex, `Age (year)`, `qRT-qPCR Pm positive`, `qRT-qPCR Pf positive`,
                `Microscopic Malaria infection status (P.f, P.m, P.o)`) %>% 
  rename(age        = `Age (year)`,
         microscopy = `Microscopic Malaria infection status (P.f, P.m, P.o)`) %>% 
  mutate(Pm_positive = case_when( `qRT-qPCR Pm positive`=="x"~"yes",
                                  TRUE~"no"),
         Pf_positive = case_when( `qRT-qPCR Pf positive`=="x"~"yes",
                                  TRUE~"no"),
         participants = 1,
         infected_pm = case_when(Pm_positive=="yes"~1,
                                 TRUE~0),
         infected_pf = case_when(Pf_positive=="yes"~1,
                                 TRUE ~ 0),
         infected_microscopy_pm = case_when(microscopy %in% c("Pf+Pm", "Pm", 
                                                              "Pf +Pm") ~ 1,
                                            TRUE ~ 0),
         infected_microscopy_pf = case_when(microscopy %in% c("Pf", "Pf+Pm",
                                                              "Pf +Pm", "PF", "pf") ~ 1,
                                            TRUE ~ 0)) %>% 
  group_by(District) %>% 
  summarise(age_lower                         = min(age, na.rm = TRUE),
            age_upper                            = max(age, na.rm = TRUE),
            Number_of_participant_tot            = sum(participants, na.rm = TRUE),
            Number_of_infected_tot_pm            = sum(infected_pm, na.rm = TRUE),
            Number_of_infected_microscopy_tot_pm = sum(infected_microscopy_pm,
                                                       na.rm = TRUE),
            PR_Pm_tot                            = Number_of_infected_tot_pm/Number_of_participant_tot,
            PR_Pm_tot_microscopy                 = Number_of_infected_microscopy_tot_pm/Number_of_participant_tot,
            Number_of_infected_tot_pf            = sum(infected_pf, na.rm = TRUE),
            Number_of_infected_microscopy_tot_pf = sum(infected_microscopy_pf,
                                                       na.rm = TRUE),
            PR_Pf_tot                            = Number_of_infected_tot_pf/Number_of_participant_tot,
            PR_Pf_tot_microscopy                 = Number_of_infected_microscopy_tot_pf/Number_of_participant_tot,
            .groups = "drop") %>% 
  distinct()



CONG_1_GPS <- CONG_1 %>% 
  left_join(village_coords_cong %>% dplyr::select(District, lat, lng, Match_confidence))

spatial_data_cong <- st_as_sf(CONG_1_GPS, coords = c("lng", "lat"), crs = 4326)
st_write(spatial_data_cong, "data/cong_data/spatial_data_sample_cong.gpkg", append = FALSE)

raw_pr_cong <- ggplot() +
  geom_sf(data = shp_con, fill = "white", color = "black") +
  geom_point(data = CONG_1_GPS, aes(x = lng, y = lat,
                                    size = Number_of_participant_tot, 
                                    color = PR_Pm_tot)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Raw P. malariae prevalence",
       color = "Prevalence", size = "N tested") +
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/cong_output/raw_obs_PR.pdf", 
       plot = raw_pr_cong, width = 6, height = 4, dpi = 300)

sampled_site_cong <- shp_con %>% 
  filter(adm2_name %in% c("Brazzaville", "Goma-Tsetse"))

sampled_location_pr_cong <- ggplot() +
  geom_sf(data = sampled_site_cong, fill = "white", color = "black") +
  geom_point(data = CONG_1_GPS, aes(x = lng, y = lat, 
                                    size = Number_of_participant_tot, 
                                    color = PR_Pm_tot)) +
  scale_color_viridis_c(option = "C") +
  labs(title = "Raw P. malariae prevalence",
       color = "Prevalence", size = "N tested") +
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave("Outputs/cong_output/raw_obs_PR_sampled_location.pdf", 
       plot = sampled_location_pr_cong, width = 6, height = 4, dpi = 300)

mapview(shp_con,
        col.regions = "white",
        alpha.regions = 0.2,
        color = "black") +
  mapview(spatial_data_cong,
          zcol = "PR_Pm_tot",
          cex = "Number_of_participant_tot",
          col.regions = colorRampPalette(c("blue", "red")),
          at = NULL,
          layer.name = "P. malariae prevalence")


################################################################################
#                 Simplification of the African map                            #
################################################################################

```{r}
shp_global <- read_sf("data/Africa_data/shape_files_1/geoBoundariesCGAZ_ADM2.shp")

custom_match <- c(
  XKX = "Europe",
  ESH = "Africa",
  `111` = "Africa",        # Abyei
  `112` = "Asia",          # Aksai Chin
  `113` = "Asia",          # China–India border
  `114` = "Asia",          # Demchok
  `115` = "Europe",        # Dragonja
  `116` = "Africa",        # Dramana-Shakatoe
  `117` = "South America", # Falkland Islands
  `118` = "Asia",          # Gaza Strip
  `119` = "Asia",          # Kalapani
  `120` = "South America", # Isla Brasilera
  `121` = "Asia",          # Siachen-Saltoro
  `122` = "Africa",        # Koualou
  `123` = "Asia",          # Liancourt Rocks
  `124` = "Asia",          # No Man's Land
  `125` = "Asia",          # Paracel Islands
  `126` = "Asia",          # Sanafir & Tiran Islands
  `127` = "Asia",          # Senkaku Islands
  `128` = "Asia",          # Spratly Islands
  `129` = "Asia"           # West Bank
)

shp_global$continent <- countrycode(
  shp_global$shapeGroup,
  origin = "iso3c",
  destination = "continent",
  custom_match = custom_match
)

shp_africa <- shp_global %>% 
  filter(continent=="Africa")

custom_match <- c(
  `111` = "Abyei",                   # Abyei
  `116` = "Dramana-Shakatoe",        # Dramana-Shakatoe
  `122` = "Koualou"                  # Koualou
)


shp_africa$country <- countrycode(
  shp_africa$shapeGroup,
  origin = "iso3c",
  destination = "country.name",
  custom_match = custom_match
)

shp_africa_simple <- ms_simplify(
  shp_africa,
  keep = 0.05,          # keep 5% of the vertices
  keep_shapes = TRUE
)

st_write(
  shp_africa_simple,
  "data/Africa_data/shp_africa/Africa_simplified.shp",
  delete_layer = TRUE,
  append = TRUE
)


ggplot()+
  geom_sf(data = shp_africa_simple, fill = "grey", color = "black")

```

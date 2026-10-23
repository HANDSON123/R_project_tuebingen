################################################################################
#          Variability and relationship between covariates and prevalence      #  
################################################################################

####################### Import necessary libraries #############################

library(tidyverse)
library(sf)
library(terra)
library(lwgeom)

###################### Loading the data ########################################

spatial_data_cmr     <- st_read("data/cmr_data/cmr_cluster_df_1_sf.gpkg")
spatial_data_gab     <- st_read("data/gab_data/GAB_df_1_sf.gpkg")
spatial_data_ben     <- st_read("data/ben_data/spatial_data_ben_sf.gpkg")
spatial_data_cong    <- st_read("data/cong_data/spatial_data_cong_sf.gpkg")

sampled_data_join_sf <- st_read("data/Africa_data/sampled_data_join_sf.gpkg")



# Names of your covariates
covariates <- c(
  "vegetation_index",
  "pop_density_2020",
  "annual_temperature_2021",
  "annual_rainfall_2021",
  "elevation",
  "annual_ITN_usage_2021",
  "dist_to_water",
  "lc_water",
  "lc_trees",
  "lc_built",
  "lc_cropland",
  "lc_wetland")

# Compute summary statistics

# compute_cov_stats <- function(data, country_col, covariates){
#   
#   # Remove geometry
#   data <- st_drop_geometry(data)
#   
#   # Loop over countries
#   results <- lapply(unique(data[[country_col]]), function(ctry){
#     
#     dat_country <- data[data[[country_col]] == ctry, ]
#     
#     means <- sapply(covariates,
#                     function(v) mean(dat_country[[v]], na.rm = TRUE))
#     
#     vars <- sapply(covariates,
#                    function(v) var(dat_country[[v]], na.rm = TRUE))
#     
#     sds <- sapply(covariates,
#                   function(v) sd(dat_country[[v]], na.rm = TRUE))
#     
#     ranges <- sapply(covariates, function(v){
#       x <- dat_country[[v]]
#       diff(range(x, na.rm = TRUE))
#     })
#     
#     std_sd <- sds / (ranges)
#     
#     data.frame(
#       Country = ctry,
#       Covariate = covariates,
#       Mean = means,
#       Standard_sd = std_sd
#     )
#   })
#   
#   results <- bind_rows(results)
#   
#   results
# }
# 
# # Run
# cov_summary_all <- compute_cov_stats(
#   data = sampled_data_join_sf,
#   country_col = "country",
#   covariates = covariates
# )
# 
# 
# cv_table <- cov_summary_all %>%
#   select(Country, Covariate, StandardizedVariance) %>%
#   pivot_wider(
#     names_from = Covariate,
#     values_from = StandardizedVariance)
# 
# 
# top5_covariates <- cov_summary_all %>%
#   group_by(Country) %>%
#   slice_max(StandardizedVariance, n = 5)
# 
# cv_table <- top5_covariates %>%
#   select(Country, Covariate, StandardizedVariance) %>%
#   pivot_wider(
#     names_from = Covariate,
#     values_from = StandardizedVariance)


plot_data <- spatial_data_gab %>%
  st_drop_geometry() %>%
  select(PR_Pm_tot, all_of(covariates)) %>%
  pivot_longer(
    cols = all_of(covariates),
    names_to = "Covariate",
    values_to = "Value"
  )

ggplot(plot_data,
       aes(x = PR_Pm_tot,
           y = Value)) +
  geom_point(size = 2.5, alpha = 0.7) +
  facet_wrap(~Covariate,
             scales = "free_x",
             ncol = 3) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    linewidth = 1.1)+
  theme_classic()

ggplot(spatial_data_gab,
       aes(x = PR_Pm_tot,
           y = annual_rainfall_2021)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    linewidth = 1.1)+
  theme_classic()



scaled_data <- spatial_data_gab %>% 
  mutate(annual_rainfall_2021    = as.numeric(scale(annual_rainfall_2021)),
         annual_temperature_2021 = as.numeric(scale(annual_temperature_2021)),
         pop_density_2020        = as.numeric(scale(pop_density_2020)),
         dist_to_water           = as.numeric(scale(dist_to_water)),
         elevation               = as.numeric(scale(elevation)),
         vegetation_index        = as.numeric(scale(vegetation_index)),
         lc_water                = as.numeric(scale(lc_water)),
         lc_wetland              = as.numeric(scale(lc_wetland)),
         lc_cropland             = as.numeric(scale(lc_cropland)),
         lc_built                = as.numeric(scale(lc_built)),
         lc_trees                = as.numeric(scale(lc_trees)))


scaled_data1 <- spatial_data_gab %>% 
  mutate(annual_rainfall_2021    = as.numeric(scale(annual_rainfall_2021)),
         annual_temperature_2021 = as.numeric(scale(annual_temperature_2021)),
         pop_density_2020        = as.numeric(scale(pop_density_2020)),
         dist_to_water           = as.numeric(scale(dist_to_water)))

plot_data <- scaled_data1 %>%
  pivot_longer(
    cols = all_of(covariates),
    names_to = "Covariate",
    values_to = "Value")

unbounded <- c("annual_rainfall_2021", "annual_temperature_2021",
               "pop_density_2020", "dist_to_water")


ggplot(plot_data %>% filter(Covariate %in% unbounded),
       aes(x = Covariate,
           y = Value)) +
  geom_boxplot(fill = "lightblue") +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkslategrey")+
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

bounded <- c(
  "vegetation_index",
  "lc_water",
  "lc_trees",
  "lc_built",
  "lc_cropland",
  "lc_wetland")

ggplot(plot_data %>% filter((Covariate %in% bounded)),
       aes(x = Covariate,
           y = Value)) +
  geom_boxplot(fill = "lightblue") +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkslategrey")+
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

selected_covariate_final <- c("annual_rainfall_2021", "annual_temperature_2021",
                              "vegetation_index", "dist_to_water", "lc_built", "PR_Pm_tot")

spatial_data_gab_1 <- spatial_data_gab[, selected_covariate_final]

plot_data <- spatial_data_gab_1 %>% 
  st_drop_geometry()%>%
  pivot_longer(
    cols = 1:4,
    names_to = "Covariate",
    values_to = "Value")


ggplot(plot_data,
       aes(x = PR_Pm_tot,
           y = Value)) +
  geom_point(size = 2.5, alpha = 0.7) +
  facet_wrap(~Covariate,
             scales = "free_x",
             ncol = 3) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    linewidth = 1.1)+
  theme_classic()

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

# sampled_data_join_sf <- spatial_data_cmr |> 
#   bind_rows(spatial_data_ben |> mutate(country = "Benin") |> select(-Village)) |> 
#   bind_rows(spatial_data_cong|> mutate(country = "Congo") |> select(-District)) |> 
#   bind_rows(spatial_data_gab)
# 
# st_write(sampled_data_join_sf, 
#          "data/Africa_data/sampled_data_join_sf.gpkg", append = FALSE)


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

################################################################################
#                              Gabon case                                      #
################################################################################

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
  theme_bw(base_size = 25) +
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
  theme_bw(base_size = 25) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

selected_covariate_final <- c("annual_rainfall_2021", "annual_temperature_2021",
                              "vegetation_index", "dist_to_water", "pop_density_2020",
                               "PR_Pm_tot")

spatial_data_gab_1 <- scaled_data[, selected_covariate_final]

plot_data <- spatial_data_gab_1 %>% 
  st_drop_geometry()%>%
  pivot_longer(
    cols = 1:5,
    names_to = "Covariate",
    values_to = "Value")


ggplot(plot_data,
       aes(x = Value,
           y = PR_Pm_tot)) +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "red",
              linewidth = 1) +
  geom_point(size = 2.5, alpha = 0.7) +
  facet_wrap(~Covariate,
             scales = "free_x",
             ncol = 3) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    linewidth = 1.1)+
  theme_classic(base_size = 25)+
  labs(x = "Covariates", y = "Prevalence")

# selected_covariate_final1 <- c(
#   "annual_rainfall_2021",
#   "annual_temperature_2021",
#   "vegetation_index",
#   "dist_to_water",
#   "pop_density_2020"
# )
# 
# outcome <- "PR_Pm_tot"
# 
# plots <- lapply(selected_covariate_final1, function(var) {
#   
#   data_plot <- spatial_data_gab_1[, c(outcome, var)]
#   data_plot <- na.omit(data_plot)
#   
#   cor_test <- cor.test(
#     data_plot[[outcome]],
#     data_plot[[var]],
#     method = "pearson"
#   )
#   
#   r <- cor_test$estimate
#   p <- cor_test$p.value
#   
#   ggplot(data_plot, aes(x = .data[[var]], y = .data[[outcome]])) +
#     geom_point(alpha = 0.5, size = 1.5) +
#     geom_smooth(method = "lm", se = TRUE, color = "red") +
#     annotate(
#       "text",
#       x = Inf, y = Inf,
#       label = sprintf("r = %.3f\np = %.3g", r, p),
#       hjust = 1.2, vjust = 1.5
#     ) +
#     labs(x = var, y = outcome) +
#     theme_bw()
# })
# 
# wrap_plots(plots, ncol = 2)

################################################################################
#                           Cameroon case                                      #
################################################################################

scaled_data <- spatial_data_cmr %>% 
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


scaled_data1 <- spatial_data_cmr %>% 
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
  theme_bw(base_size = 25) +
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
  theme_bw(base_size = 25) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

selected_covariate_final <- c("annual_rainfall_2021", "annual_temperature_2021",
                              "lc_trees", "lc_water", "dist_to_water", "pop_density_2020",
                              "PR_Pm_tot")

spatial_data_cmr_1 <- scaled_data[, selected_covariate_final]

plot_data <- spatial_data_cmr_1 %>% 
  st_drop_geometry()%>%
  pivot_longer(
    cols = 1:6,
    names_to = "Covariate",
    values_to = "Value")


ggplot(plot_data,
       aes(x = Value,
           y = PR_Pm_tot)) +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "red",
              linewidth = 1) +
  geom_point(size = 2.5, alpha = 0.7) +
  facet_wrap(~Covariate,
             scales = "free_x",
             ncol = 3) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    linewidth = 1.1)+
  theme_classic(base_size = 25)+
  labs(x = "Covariates", y = "Prevalence")


################################################################################
#                              Benin case                                      #
################################################################################

scaled_data <- spatial_data_ben %>% 
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


scaled_data1 <- spatial_data_ben %>% 
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
  theme_bw(base_size = 25) +
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
  theme_bw(base_size = 25) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

selected_covariate_final <- c("annual_rainfall_2021", "annual_temperature_2021",
                              "lc_trees", "lc_water", "dist_to_water", 
                              "pop_density_2020", "vegetation_index", "lc_wetland",
                              "PR_Pm_tot")

spatial_data_ben_1 <- scaled_data[, selected_covariate_final]

plot_data <- spatial_data_ben_1 %>% 
  st_drop_geometry()%>%
  pivot_longer(
    cols = 1:8,
    names_to = "Covariate",
    values_to = "Value")


ggplot(plot_data,
       aes(x = Value,
           y = PR_Pm_tot)) +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "red",
              linewidth = 1) +
  geom_point(size = 2.5, alpha = 0.7) +
  facet_wrap(~Covariate,
             scales = "free_x",
             ncol = 3) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    linewidth = 1.1)+
  theme_classic(base_size = 25)+
  labs(x = "Covariates", y = "Prevalence")


################################################################################
#                              Congo case                                      #
################################################################################

scaled_data <- spatial_data_cong %>% 
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


scaled_data1 <- spatial_data_cong %>% 
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
  theme_bw(base_size = 25) +
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
  theme_bw(base_size = 25) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

selected_covariate_final <- c("annual_rainfall_2021", "annual_temperature_2021",
                              "lc_trees", "lc_water", "dist_to_water", "lc_built",
                              "pop_density_2020", "PR_Pm_tot")

spatial_data_cong_1 <- scaled_data[, selected_covariate_final]

plot_data <- spatial_data_cong_1 %>% 
  st_drop_geometry()%>%
  pivot_longer(
    cols = 1:7,
    names_to = "Covariate",
    values_to = "Value")


ggplot(plot_data,
       aes(x = Value,
           y = PR_Pm_tot)) +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "red",
              linewidth = 1) +
  geom_point(size = 2.5, alpha = 0.7) +
  facet_wrap(~Covariate,
             scales = "free_x",
             ncol = 3) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    linewidth = 1.1)+
  theme_classic(base_size = 25)+
  labs(x = "Covariates", y = "Prevalence")


################################################################################
#                        Continental case                                      #
################################################################################

scaled_data <- sampled_data_join_sf %>% 
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


scaled_data1 <- sampled_data_join_sf %>% 
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
  theme_bw(base_size = 25) +
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
  theme_bw(base_size = 25) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

selected_covariate_final <- c("annual_rainfall_2021", "annual_temperature_2021",
                              "lc_trees", "lc_water", "dist_to_water", "lc_built",
                              "pop_density_2020", "PR_Pm_tot")

spatial_data_cong_1 <- scaled_data[, selected_covariate_final]

plot_data <- spatial_data_cong_1 %>% 
  st_drop_geometry()%>%
  pivot_longer(
    cols = 1:7,
    names_to = "Covariate",
    values_to = "Value")


ggplot(plot_data,
       aes(x = Value,
           y = PR_Pm_tot)) +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "red",
              linewidth = 1) +
  geom_point(size = 2.5, alpha = 0.7) +
  facet_wrap(~Covariate,
             scales = "free_x",
             ncol = 3) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    linewidth = 1.1)+
  theme_classic(base_size = 25)+
  labs(x = "Covariates", y = "Prevalence")

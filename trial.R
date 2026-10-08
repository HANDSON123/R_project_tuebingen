countries1_sf <- rbind(shp_gab %>% dplyr::select(COUNTRY, NAME_2, geometry) %>%
                         rename(country = COUNTRY,
                                district = NAME_2), shp_con %>% dplyr::select(adm0_name,
                                                                              adm2_name,
                                                                              geometry) %>%
                         rename(country = adm0_name,
                                district = adm2_name))


countries_sf <- rbind(countries1_sf, shp_cmr %>% dplyr::select(Admin2, geometry) %>%
                        mutate(country = "Cameroon") %>%
                        rename(district = Admin2))


sampled_data_join_sf$country_id <- as.integer(factor(sampled_data_join_sf$country))  # 1, 2, 3, 4

sf_data <- st_as_sf(sampled_data_join_sf, coords = c("lon", "lat"), crs = 4326)
sf_data <- st_transform(sf_data, 3857)
sf_data <- sf_data %>%
  mutate(country_id = case_when(country_id %in% c(2, 3, 4) ~ 2,
                                TRUE ~ country_id))

coords <- st_coordinates(sf_data)

coords_A <- coords[sf_data$country_id == 1, ]
coords_B <- coords[sf_data$country_id == 2, ]


country_A <- st_transform(shp_ben, 3857)
country_B <- st_transform(countries_sf, 3857)

country_A_union <- st_union(country_A)
country_B_union <- st_union(country_B)

country_A_simplified <- st_simplify(country_A_union, dTolerance = 10000)
country_B_simplified <- st_simplify(country_B_union, dTolerance = 10000)


# Extract geometry
geom_A <- st_geometry(country_A_simplified)[[1]]

# If MULTIPOLYGON : take the first polygon
if (inherits(geom_A, "MULTIPOLYGON")) {
  geom_A <- geom_A[[1]]
}

# Extract the OUTER ring only (first ring)
outer_ring_A <- geom_A[[1]]

# Convert to matrix with X,Y
boundary_A <- as.matrix(outer_ring_A)[, 1:2, drop = FALSE]



# country_B is an sf polygon or multipolygon
geom_B <- st_geometry(country_B_simplified)[[1]]

# If MULTIPOLYGON : take the first polygon
if (inherits(geom_B, "MULTIPOLYGON")) {
  geom_B <- geom_B[[1]]
}

# Extract the OUTER ring only (first ring)
outer_ring_B <- geom_B[[1]]

# Convert to matrix with X,Y
boundary_B <- as.matrix(outer_ring_B)[, 1:2, drop = FALSE]



# mesh_A <- inla.mesh.2d(loc = coords_A, max.edge = c(5000, 20000))
# mesh_B <- inla.mesh.2d(loc = coords_B, max.edge = c(5000, 20000))
# mesh_C <- inla.mesh.2d(loc = coords_C, max.edge = c(5000, 20000))

mesh_A <- inla.mesh.2d(
  loc      = coords_A,
  boundary = list(boundary_A),
  max.edge = c(50000, 300000),   # inner ~50 km, outer ~1000 km
  cutoff   = 30000,
  offset   = c(100000, 200000)
)

mesh_B <- inla.mesh.2d(
  loc = coords_B,
  boundary   = list(boundary_B) ,
  max.edge = c(50000, 300000),   # inner ~30 km, outer ~80 km
  cutoff   = 30000,
  offset   = c(100000, 200000)
)

# Extract triangle vertex indices
tri <- mesh_B$graph$tv


# Build edges from triangles
edges_list <- list()

for (i in 1:nrow(tri)) {
  idx <- tri[i, ]
  # Each triangle has 3 edges
  edges_list[[length(edges_list) + 1]] <- st_linestring(mesh_B$loc[idx[c(1,2)], 1:2])
  edges_list[[length(edges_list) + 1]] <- st_linestring(mesh_B$loc[idx[c(2,3)], 1:2])
  edges_list[[length(edges_list) + 1]] <- st_linestring(mesh_B$loc[idx[c(3,1)], 1:2])
}

mesh_edges_sf <- st_as_sf(st_sfc(edges_list, crs = 3857))

coords_df_B <- as.data.frame(coords_B)
colnames(coords_df_B) <- c("X", "Y")

ggplot() +
  geom_sf(data = country_B_simplified, fill = "white", color = "black", size = 0.8) +
  geom_sf(data = mesh_edges_sf, color = "grey60", size = 0.3) +
  geom_point(data = coords_df_B, aes(X, Y), color = "red", size = 1.5) +
  theme_minimal() +
  ggtitle("Mesh B with Country Boundary")



# Build Matérn SPDE model of the spatial field on mesh

spde_A <- inla.spde2.pcmatern(
  mesh = mesh_A,
  prior.range = c(50000, 0.5),
  prior.sigma = c(0.1, 0.5)
)

spde_B <- inla.spde2.pcmatern(
  mesh = mesh_B,
  prior.range = c(50000, 0.5),
  prior.sigma = c(0.1, 0.5)   # 0.1 came from the squared root of the sill because it
)                              # represent the standard deviance while the sill represent the variance.

# Computes the interpolation weights that map your SPDE spatial field (defined on mesh nodes) to your actual malaria data locations

# 1. Create indicator weights (1 if the row belongs to the country, 0 if not)
# 'coords' must be the full 36-row matrix of ALL coordinates
w_A <- as.numeric(sf_data$country_id == 1)
w_B <- as.numeric(sf_data$country_id == 2)

# 2. Build the full matrices automatically.
# Passing the full 'coords' matrix ensures every matrix has exactly 36 rows.
# The 'weights' argument forces rows from other countries to be completely zero.
A_A_full <- inla.spde.make.A(mesh = mesh_A, loc = coords, weights = w_A)
A_B_full <- inla.spde.make.A(mesh = mesh_B, loc = coords, weights = w_B)

A_fixed <- Diagonal(n = nrow(sf_data))
A_spatial <- list(A_A_full, A_B_full)


A_A_full <- as(A_A_full, "dgCMatrix")
A_B_full <- as(A_B_full, "dgCMatrix")


#1. Variables to be scale
vars_to_scale <- c("annual_rainfall_2021", "annual_temperature_2021",
                   "pop_density_2020", "dist_to_water", "vegetation_index",
                   "elevation")

# 2. Create a named list to store the training parameters (Means and SDs)
scaling_stats <- list()

# 3. Scale the observed data and extract the parameters
for (var in vars_to_scale) {
  # Calculate mean and standard deviation from the 36 observed points
  mean_val <- mean(sf_data[[var]], na.rm = TRUE)
  sd_val   <- sd(sf_data[[var]], na.rm = TRUE)

  # Store them for later use on the grid
  scaling_stats[[var]] <- c(mean = mean_val, sd = sd_val)

  # Apply the scaling to the observed data
  sf_data[[var]] <- (sf_data[[var]] - mean_val) / sd_val
}

################################################################################
#        checking for over dispersion in our model
################################################################################

# 1. Calculate observed prevalence for each of the 36 villages
village_p <- sf_data$Number_of_infected_tot_pf / sf_data$Number_of_participant_tot

# 2. Calculate the global overall prevalence
global_p <- sum(sf_data$Number_of_infected_tot_pf) / sum(sf_data$Number_of_participant_tot)

# 3. Get the average number of participants per village
mean_N <- mean(sf_data$Number_of_participant_tot)

# 4. Calculate the Observed vs. Theoretical Variance
observed_variance    <- var(village_p)
theoretical_variance <- (global_p * (1 - global_p)) / mean_N

# 5. Compute the Raw Dispersion Ratio
raw_dispersion_ratio <- observed_variance / theoretical_variance
print(raw_dispersion_ratio)

# sf_data <- sf_data %>%
#   mutate(pop_density_2020 = as.numeric(scale(pop_density_2020)),
#          annual_rainfall_2021 = as.numeric(scale(annual_rainfall_2021)),
#          annual_temperature_2021 = as.numeric(scale(annual_temperature_2021)),
#          vegetaion_index = as.numeric(scale(vegetaion_index)),
#          dist_to_water  = as.numeric(scale(dist_to_water)),
#          elevation = as.numeric(scale(elevation)))

# Model building

effects_fixed_rainfall <- data.frame(
  intercept             = 1,
  annual_rainfall_2021  = sf_data$annual_rainfall_2021,
  vegetation_index      = sf_data$vegetation_index,
  pop_density_2020      = sf_data$pop_density_2020,
  dist_to_water         = sf_data$dist_to_water,
  lc_water              = sf_data$lc_water,
  lc_trees              = sf_data$lc_trees,
  lc_built              = sf_data$lc_built,
  lc_cropland           = sf_data$lc_cropland,
  lc_wetland            = sf_data$lc_wetland
  # annual_ITN_usage_2021 = sf_data$annual_ITN_usage_2021
)

effects_fixed_temperature <- data.frame(
  intercept             = 1,
  annual_temperature_2021  = sf_data$annual_temperature_2021,
  vegetation_index      = sf_data$vegetation_index,
  pop_density_2020      = sf_data$pop_density_2020,
  dist_to_water         = sf_data$dist_to_water,
  lc_water              = sf_data$lc_water,
  lc_trees              = sf_data$lc_trees,
  lc_built              = sf_data$lc_built,
  lc_cropland           = sf_data$lc_cropland,
  lc_wetland            = sf_data$lc_wetland
  # annual_ITN_usage_2021 = sf_data$annual_ITN_usage_2021
)


effects_elevation      <- data.frame(
  intercept             = 1,
  elevation             = sf_data$elevation,
  vegetation_index      = sf_data$vegetation_index,
  pop_density_2020      = sf_data$pop_density_2020,
  dist_to_water         = sf_data$dist_to_water,
  lc_water              = sf_data$lc_water,
  lc_trees              = sf_data$lc_trees,
  lc_built              = sf_data$lc_built,
  lc_cropland           = sf_data$lc_cropland,
  lc_wetland            = sf_data$lc_wetland
  # annual_ITN_usage_2021 = sf_data$annual_ITN_usage_2021
)

effects_fixed_baseline <- data.frame(
  intercept             = 1,
  vegetation_index      = sf_data$vegetation_index,
  pop_density_2020      = sf_data$pop_density_2020,
  dist_to_water         = sf_data$dist_to_water,
  lc_water              = sf_data$lc_water,
  lc_trees              = sf_data$lc_trees,
  lc_built              = sf_data$lc_built,
  lc_cropland           = sf_data$lc_cropland,
  lc_wetland            = sf_data$lc_wetland
  # annual_ITN_usage_2021 = sf_data$annual_ITN_usage_2021
)


effects_country <- list(
  country_id = sf_data$country_id
)

effects_A <- list(field_A = 1:spde_A$n.spde)
effects_B <- list(field_B = 1:spde_B$n.spde)

A_fixed   <- Diagonal(nrow(sf_data))
A_country <- Diagonal(nrow(sf_data))


stack_test_rainfall <- inla.stack(
  data = list(
    y       = sf_data$Number_of_infected_tot_pf,
    Ntrials = sf_data$Number_of_participant_tot
  ),
  A = list(
    A_fixed,
    A_country,
    A_A_full,
    A_B_full),
  effects = list(
    effects_fixed_rainfall,
    effects_country,
    effects_A,
    effects_B),
  tag = "est"
)

stack_test_temperature <- inla.stack(
  data = list(
    y       = sf_data$Number_of_infected_tot_pf,
    Ntrials = sf_data$Number_of_participant_tot
  ),
  A = list(
    A_fixed,
    A_country,
    A_A_full,
    A_B_full),
  effects = list(
    effects_fixed_temperature,
    effects_country,
    effects_A,
    effects_B),
  tag = "est"
)


stack_test_elevation <- inla.stack(
  data = list(
    y       = sf_data$Number_of_infected_tot_pf,
    Ntrials = sf_data$Number_of_participant_tot
  ),
  A = list(
    A_fixed,
    A_country,
    A_A_full,
    A_B_full),
  effects = list(
    effects_elevation,
    effects_country,
    effects_A,
    effects_B),
  tag = "est"
)


stack_test_baseline <- inla.stack(
  data = list(
    y       = sf_data$Number_of_infected_tot_pf,
    Ntrials = sf_data$Number_of_participant_tot
  ),
  A = list(
    A_fixed,
    A_country,
    A_A_full,
    A_B_full),
  effects = list(
    effects_fixed_baseline,
    effects_country,
    effects_A,
    effects_B),
  tag = "est"
)


# Penalize complexity which control how much countries differ from each other through
# the country level random effect
hyper_country <- list(
  prec = list(prior = "pc.prec", param = c(0.5, 0.5))
)

formula_rainfall <- y ~
  -1 +
  intercept + annual_rainfall_2021 + pop_density_2020 + vegetation_index +
  dist_to_water + lc_water + lc_trees + lc_built + lc_cropland + lc_wetland +
  f(country_id, model="iid", hyper=hyper_country) +
  f(field_A, model=spde_A) +
  f(field_B, model=spde_B)

formula_temperature <- y ~
  -1 +
  intercept + annual_temperature_2021 + pop_density_2020 + vegetation_index +
  dist_to_water + lc_water + lc_trees + lc_built + lc_cropland + lc_wetland +
  f(country_id, model="iid", hyper=hyper_country) +
  f(field_A, model=spde_A) +
  f(field_B, model=spde_B)

formula_elevation <- y ~
  -1 +
  intercept + elevation + pop_density_2020 + dist_to_water + vegetation_index +
  lc_water + lc_trees + lc_built + lc_cropland + lc_wetland +
  f(country_id, model="iid", hyper=hyper_country) +
  f(field_A, model=spde_A) +
  f(field_B, model=spde_B)

formula_baseline <- y ~
  -1 +
  intercept + pop_density_2020 + dist_to_water + vegetation_index +
  lc_water + lc_trees + lc_built + lc_cropland + lc_wetland +
  f(country_id, model="iid", hyper=hyper_country) +
  f(field_A, model=spde_A) +
  f(field_B, model=spde_B)



res_rainfall <- inla(
  formula_rainfall,
  family = "binomial",
  data = inla.stack.data(stack_test_rainfall),
  Ntrials = inla.stack.data(stack_test_rainfall)$Ntrials,
  control.predictor = list(
    A = inla.stack.A(stack_test_rainfall),
    compute = TRUE
  ),
  control.fixed = list(
    mean = 0,
    prec = 1
  ),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)


res_temperature <- inla(
  formula_temperature,
  family = "binomial",
  data = inla.stack.data(stack_test_temperature),
  Ntrials = inla.stack.data(stack_test_temperature)$Ntrials,
  control.predictor = list(
    A = inla.stack.A(stack_test_temperature),
    compute = TRUE
  ),
  control.fixed = list(
    mean = 0,
    prec = 1
  ),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)


res_elevation <- inla(
  formula_elevation,
  family = "binomial",
  data = inla.stack.data(stack_test_elevation),
  Ntrials = inla.stack.data(stack_test_elevation)$Ntrials,
  control.predictor = list(
    A = inla.stack.A(stack_test_elevation),
    compute = TRUE
  ),
  control.fixed = list(
    mean = 0,
    prec = 1
  ),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

res_baseline <- inla(
  formula_baseline,
  family = "binomial",
  data = inla.stack.data(stack_test_baseline),
  Ntrials = inla.stack.data(stack_test_baseline)$Ntrials,
  control.predictor = list(
    A = inla.stack.A(stack_test_baseline),
    compute = TRUE
  ),
  control.fixed = list(
    mean = 0,
    prec = 1
  ),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

# ## diasgnostic on the distance to water
#
# res_simple <- inla(
#   Number_of_infected_tot_pm ~ dist_to_water,
#   family = "binomial",
#   data = sf_data,
#   Ntrials = sf_data$Number_of_participant_tot
# )
#
# res_simple$summary.fixed

# The model’s hyperparameters show that country-level differences are
# almost negligible, meaning the covariates already explain most variation
# between Benin and Congo. The spatial field in Benin is clearly
# identifiable: it has a short correlation range (\~7 km) and a moderate
# spatial variance, so spatial structure genuinely improves the fit there.
# In contrast, Congo’s field is extremely uncertain—with only three
# locations, the model cannot reliably estimate its spatial range or
# variance, so any spatial effect in Congo should be interpreted
# cautiously. The practical consequence is that predictions and spatial
# maps will be trustworthy in Benin, while Congo spatial surface will be
# weak, noisy, and driven mostly by priors rather than data.

effects_fixed_rainfall <- data.frame(
  intercept             = 1,
  annual_rainfall_2021  = sf_data$annual_rainfall_2021,
  vegetation_index      = sf_data$vegetation_index,
  pop_density_2020      = sf_data$pop_density_2020,
  dist_to_water         = sf_data$dist_to_water,
  lc_water              = sf_data$lc_water,
  lc_trees              = sf_data$lc_trees,
  lc_built              = sf_data$lc_built,
  lc_cropland           = sf_data$lc_cropland,
  lc_wetland            = sf_data$lc_wetland
  # annual_ITN_usage_2021 = sf_data$annual_ITN_usage_2021
)

stack_test_rainfall <- inla.stack(
  data = list(
    y       = sf_data$Number_of_infected_tot_pm,
    Ntrials = sf_data$Number_of_participant_tot
  ),
  A = list(
    A_fixed,
    A_country,
    A_A_full,
    A_B_full),
  effects = list(
    effects_fixed_rainfall,
    effects_country,
    effects_A,
    effects_B),
  tag = "est"
)

formula_rainfall <- y ~
  -1 +
  intercept + annual_rainfall_2021 + vegetation_index + pop_density_2020 +
  dist_to_water + lc_water + lc_trees + lc_built + lc_cropland + lc_wetland +
  f(country_id, model="iid", hyper=hyper_country) +
  f(field_A, model=spde_A) +
  f(field_B, model=spde_B)

res_rainfall1 <- inla(
  formula_rainfall,
  family = "betabinomial",
  data = inla.stack.data(stack_test_rainfall),
  Ntrials = inla.stack.data(stack_test_rainfall)$Ntrials,
  control.predictor = list(
    A = inla.stack.A(stack_test_rainfall),
    compute = TRUE
  ),
  control.fixed = list(
    mean = 0,
    prec = 1
  ),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

hist(res_rainfall1$cpo$pit,
     breaks = 20,
     main = "PIT Histogram",
     xlab = "PIT")

res_rainfall1$dic$dic
res_rainfall1$summary.fixed
res_rainfall1$summary.fitted.values



# Comparing the different models


# Model 1:  Rainfall

# estimated values for the different covariates

res_rainfall$summary.fixed

# Hyper parameters examination

res_rainfall$summary.hyperpar

# DIC

res_rainfall$dic$dic

# WAIC

res_rainfall$waic$waic

# CPO

sum(log(res_rainfall$cpo$cpo), na.rm = TRUE)

summary(res_rainfall$cpo$cpo)

# PIT(probability integral transform) histogram

hist(res_rainfall$cpo$pit,
     breaks = 20,
     main = "PIT Histogram",
     xlab = "PIT")

# observed vs predicted prevalence

obs <- inla.stack.data(stack_test_rainfall)$y /
       inla.stack.data(stack_test_rainfall)$Ntrials

idx.est <- inla.stack.index(stack_test_rainfall, "est")$data

pred.est <- res_rainfall$summary.fitted.values$mean[idx.est]

predicted_observed <- data.frame(obs = obs, pred = pred.est)

ggplot(data.frame(obs = obs, pred = pred.est),
       aes(x = obs, y = pred)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1) +
  labs(title = "Observed vs predicted, rainfall model",
    x = "Observed prevalence",
    y = "Predicted prevalence"
  ) +
  theme_minimal()

# residual plot

residuals <- obs - pred.est

ggplot(data.frame(pred = pred.est, resid = residuals),
       aes(x = pred, y = resid)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", size = 1) +
  labs(title = "residual rainfall model",
    x = "Predicted",
    y = "Residuals"
  ) +
  theme_minimal()


# Model 2:  temperature

# estimated values for the different covariates

res_temperature$summary.fixed

# Hyper parameters examination

res_temperature$summary.hyperpar

# DIC

res_temperature$dic$dic

# WAIC

res_temperature$waic$waic

# CPO

sum(log(res_temperature$cpo$cpo), na.rm = TRUE)

summary(res_temperature$cpo$cpo)

# PIT histogram

hist(res_temperature$cpo$pit,
     breaks = 20,
     main = "PIT Histogram",
     xlab = "PIT")

# observed vs predicted prevalence

idx.est <- inla.stack.index(stack_test_temperature, "est")$data

pred.est <- res_temperature$summary.fitted.values$mean[idx.est]

predicted_observed <- data.frame(obs = obs, pred = pred.est)

ggplot(data.frame(obs = obs, pred = pred.est),
       aes(x = obs, y = pred)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1) +
  labs(title =  "Observed vs predicted, temperature model",
    x = "Observed prevalence",
    y = "Predicted prevalence"
  ) +
  theme_minimal()

# residual plot

residuals <- obs - pred.est

ggplot(data.frame(pred = pred.est, resid = residuals),
       aes(x = pred, y = resid)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", size = 1) +
  labs(title = "residual temperature model",
    x = "Predicted",
    y = "Residuals"
  ) +
  theme_minimal()


# Model 4:  elevation

# estimated values for the different covariates

res_elevation$summary.fixed

# Hyper parameters examination

res_elevation$summary.hyperpar

# DIC

res_elevation$dic$dic

# WAIC

res_elevation$waic$waic

# CPO

sum(log(res_elevation$cpo$cpo), na.rm = TRUE)

summary(res_elevation$cpo$cpo)

# PIT histogram

hist(res_elevation$cpo$pit,
     breaks = 20,
     main = "PIT Histogram",
     xlab = "PIT")

# observed vs predicted prevalence

idx.est <- inla.stack.index(stack_test_elevation, "est")$data

pred.est <- res_elevation$summary.fitted.values$mean[idx.est]

predicted_observed <- data.frame(obs = obs, pred = pred.est)

ggplot(data.frame(obs = obs, pred = pred.est),
       aes(x = obs, y = pred)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1) +
  labs(title = "Observed vs predicted, elevation model",
    x = "Observed prevalence",
    y = "Predicted prevalence"
  ) +
  theme_minimal()

# residual plot

residuals <- obs - pred.est

ggplot(data.frame(pred = pred.est, resid = residuals),
       aes(x = pred, y = resid)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", size = 1) +
  labs(title = "residual elevation model",
    x = "Predicted",
    y = "Residuals"
  ) +
  theme_minimal()


# Model 5:  baseline model

# estimated values for the different covariates

res_baseline$summary.fixed

# Hyper parameters examination

res_baseline$summary.hyperpar

# DIC

res_baseline$dic$dic

# WAIC

res_baseline$waic$waic

# CPO

sum(log(res_baseline$cpo$cpo), na.rm = TRUE)

summary(res_baseline$cpo$cpo)

# PIT histogram

hist(res_baseline$cpo$pit,
     breaks = 20,
     main = "PIT Histogram",
     xlab = "PIT")

# observed vs predicted prevalence

idx.est <- inla.stack.index(stack_test_baseline, "est")$data

pred.est <- res_baseline$summary.fitted.values$mean[idx.est]

predicted_observed <- data.frame(obs = obs, pred = pred.est)

ggplot(data.frame(obs = obs, pred = pred.est),
       aes(x = obs, y = pred)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1) +
  labs(title = "Observed vs predicted, baseline model",
    x = "Observed prevalence",
    y = "Predicted prevalence"
  ) +
  theme_minimal()

# residual plot

residuals <- obs - pred.est

ggplot(data.frame(pred = pred.est, resid = residuals),
       aes(x = pred, y = resid)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", size = 1) +
  labs(title = "residual baseline model",
    x = "Predicted",
    y = "Residuals"
  ) +
  theme_minimal()


# ## Conclusion on models comparison (P. malaria model)
#
# The rainfall model is the right choice because it consistently shows the
# strongest performance across all diagnostics. It has the lowest DIC and
# WAIC, the best CPO, the tightest observed‑vs‑predicted fit, and it
# reduces the spatial field variance more than any other covariate,
# meaning rainfall explains a large part of the spatial structure instead
# of leaving it to the latent field. In short, rainfall improves
# prediction, reduces unexplained spatial variation, and aligns with
# malaria biology making it the most reliable and scientifically justified
# model to continue with.


# ## *P. f* model examination
#
# -   The observed vs fitted values plot shows that the smooth, shrunk
# estimates of the true prevalence, which are generally consistent
# with the observed but less extreme. the number of points above and
# below the red line are roughly the same which means the model is not
# over estimating nor under estimating.
#
# -   The residual plot shows a rough horizontal cloud indicating that the
# model is stable and most residuals are between -0.1 and 0.1 meaning
# that the model fit well. we also have some moderate residuals around
# 0.3 and -0.3 meaning the model smooth appropriately. Finally, the
# number of point above and below the horizontal line are almost the
# same, which means the model is unbiased.
#
# -   The country_id precision is approximately \~15 which means there is
# some meaningful difference between countries after accounting for
# covariates + spatial fields.
#
# -   For all the three field, the range is \~7km to \~9km with a very
# wide credible interval, meaning areas withing these distances are
# spatially correlated or influence each other but the model is unsure
# to which extend because the uncertainty is high. This means that the
# prediction will be smooth but not overly flat.
#
# -   All the three spatial field also have a very small standard
# deviation 0.021 to 0.027, which means spatial variation exist but
# covariates explain most of the variation and spatial field refine
# the prediction instead of dominating it.
#

grid_sf <- read_sf("data/grid_data_for_prediction_sf.gpkg")

# Scale the prediction grid using the EXACT same parameters as in the fitting model
# (This ensures a grid cell is measured relative to our 36 observed points)
for (var in vars_to_scale) {
  mean_val <- scaling_stats[[var]]["mean"]
  sd_val   <- scaling_stats[[var]]["sd"]

  grid_sf[[var]] <- (grid_sf[[var]] - mean_val) / sd_val
}


# coords_pred <- st_coordinates(grid_sf)
#
# coords_pred_A <- coords_pred[grid_sf$country_id == 1, ]
# coords_pred_B <- coords_pred[grid_sf$country_id == 2, ]
# coords_pred_C <- coords_pred[grid_sf$country_id == 3, ]
#
# A_pred_A <- inla.spde.make.A(mesh = mesh_A, loc = coords_pred_A)
# A_pred_B <- inla.spde.make.A(mesh = mesh_B, loc = coords_pred_B)
# A_pred_C <- inla.spde.make.A(mesh = mesh_C, loc = coords_pred_C)
#
# A_pred_A_full <- matrix(0, nrow = nrow(grid_sf), ncol = spde_A$n.spde)
# A_pred_B_full <- matrix(0, nrow = nrow(grid_sf), ncol = spde_B$n.spde)
# A_pred_C_full <- matrix(0, nrow = nrow(grid_sf), ncol = spde_C$n.spde)
#
# A_pred_A_full[grid_sf$country_id == 1, ] <- as.matrix(A_pred_A)
# A_pred_B_full[grid_sf$country_id == 2, ] <- as.matrix(A_pred_B)
# A_pred_C_full[grid_sf$country_id == 3, ] <- as.matrix(A_pred_C)


# 1. Get the full coordinates matrix for all 37,515 grid cells
coords_pred <- st_coordinates(grid_sf)

# 2. Create 0/1 indicator weights across ALL 37,515 rows
w_pred_A <- as.numeric(grid_sf$country_id == 1)
w_pred_B <- as.numeric(grid_sf$country_id == 2)

# 3. Build native, optimized Sparse Matrices (dgCMatrix) automatically.
# Passing 'coords_pred' means every matrix instantly has exactly 37,515 rows.
# 'weights' forces cells outside the country to stay zero efficiently.
A_pred_A_full <- inla.spde.make.A(mesh = mesh_A, loc = coords_pred, weights = w_pred_A)
A_pred_B_full <- inla.spde.make.A(mesh = mesh_B, loc = coords_pred, weights = w_pred_B)

# Verify they are sparse structures
class(A_pred_A_full) # Should output "dgCMatrix"

effects_fixed_pred <- data.frame(
  intercept               = 1,
  # annual_rainfall_2021    = grid_sf$annual_rainfall_2021,
  annual_temperature_2021 = grid_sf$annual_temperature_2021,
  pop_density_2020        = grid_sf$pop_density_2020,
  vegetation_index        = sf_data$vegetation_index,
  # annual_ITN_usage_2021   = grid_sf$annual_ITN_usage_2021,
  dist_to_water           = grid_sf$dist_to_water,
  lc_water                = grid_sf$lc_water,
  lc_trees                = grid_sf$lc_trees,
  lc_built                = grid_sf$lc_built,
  lc_cropland             = grid_sf$lc_cropland,
  lc_wetland              = grid_sf$lc_wetland)

effects_country_pred <- list(country_id = grid_sf$country_id)

effects_A_pred <- list(field_A = 1:spde_A$n.spde)
effects_B_pred <- list(field_B = 1:spde_B$n.spde)

stack_pred <- inla.stack(
  data = list(y = NA,
    Ntrials = NA),
  A = list(
    Diagonal(nrow(grid_sf)),      # fixed effects
    Diagonal(nrow(grid_sf)),      # country iid
    A_pred_A_full,                # spatial field A
    A_pred_B_full),
  effects = list(
    effects_fixed_pred,
    effects_country_pred,
    effects_A_pred,
    effects_B_pred),
  tag = "pred"
)

stack_full <- inla.stack(stack_test_rainfall, stack_pred)

res_pred <- inla(
  formula_temperature,
  family = "betabinomial",
  data = inla.stack.data(stack_full),
  Ntrials = inla.stack.data(stack_full)$Ntrials,
  control.predictor = list(
    A = inla.stack.A(stack_full),
    compute = TRUE,
    link = 1
  ),
  control.fixed = list(mean = 0, prec = 1),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,
    config = TRUE)
)

#idx_pred <- inla.stack.index(stack_full, "pred")$data


set.seed(123)
samples <- inla.posterior.sample(100, res_pred)

idx_pred <- inla.stack.index(stack_full, tag = "pred")$data

pred_mat <- sapply(samples, function(s) {
  lp <- s$latent[idx_pred]                   # linear predictor for grid cells
  plogis(lp)                                 # convert logit  prevalence
})


grid_sf$prev_mean   <- rowMeans(pred_mat)
grid_sf$prev_sd     <- apply(pred_mat, 1, sd)
grid_sf$prev_q025   <- apply(pred_mat, 1, quantile, 0.025)
grid_sf$prev_q975   <- apply(pred_mat, 1, quantile, 0.975)



# v1 <- res_pred$summary.fixed
# v2 <- res_pred$summary.fixed


## P. malaria model interpretation

# ### clustered model
# 
# Across the three‑country study area, *Plasmodium malariae* prevalence is
# strongly influenced by environmental conditions. Rainfall and wetland
# show clear positive associations with infection risk, indicating that
# humid and water‑rich ecosystems provide favorable breeding conditions
# for malaria vectors. Built‑up areas display a strong negative effect,
# reflecting reduced transmission in urbanized landscapes where mosquito
# habitats are disrupted. Distance to water shows a small but
# statistically meaningful positive effect. In contrast, cropland,
# population density, forest cover, and land‑cover water all have credible
# intervals that include zero, meaning the data do not provide sufficient
# evidence to classify their influence as either protective or
# risk‑enhancing. Overall, rainfall and wetlands emerge as key drivers of
# *P. malariae* transmission, with urbanization acting as an important
# suppressor, while several other covariates show uncertain or negligible
# effects.
# 
# ### non clustered model
# 
# The updated model produces results that are broadly consistent with the
# previous analysis, identifying rainfall and wetland environments as
# strong positive drivers of *Plasmodium malariae* prevalence, and
# built‑up areas as a strong negative predictor. However, the new model
# shows slightly reduced effect sizes for rainfall and increased
# uncertainty for several covariates. Notably, the effect of distance to
# water, previously statistically positive, now has a credible interval
# that includes zero, indicating insufficient evidence for a clear
# association. Cropland, population density, forest cover, and land‑cover
# water continue to show wide credible intervals that cross zero, meaning
# their effects remain statistically uncertain. Overall, the main
# difference between the two models lies in the reduced certainty around
# distance to water and slightly adjusted effect sizes, while the key
# environmental drivers rainfall, wetlands, and urbanization remain
# robust.


################################################################################
#          Leave-One-Cluster-Out-Cross validation                              #   
################################################################################

#-----------------------------------------------------------
# PREPARE DATA
#-----------------------------------------------------------
sampled_data_join_sf <- st_read("data/Africa_data/sampled_data_join_sf.gpkg")
shp_africa <- st_read("data/Africa_data/shp_africa/Africa_simplified.shp") %>% st_transform(3857)

Africa_union <- st_union(shp_africa )

sampled_data_join_sf <- sampled_data_join_sf %>% st_transform(3857)

sampled_data_join_sf_sp     <- as(sampled_data_join_sf, "Spatial")
Africa_union_sp <- as(Africa_union, "Spatial")

coords_Africa <- coordinates(sampled_data_join_sf_sp)
sampled_data_join_sf$lon <- coords_Africa[,1]
sampled_data_join_sf$lat <- coords_Africa[,2]
#-----------------------------------------------------------
# BUILD ONE GLOBAL MESH
#-----------------------------------------------------------

mesh <- inla.mesh.2d(
  loc      = coords_Africa,
  boundary = inla.sp2segment(Africa_union_sp),
  max.edge = c(25e3, 2000e3),   
  cutoff   = 20e3,
  offset   = c(100e3, 2000e3))


#-----------------------------------------------------------
# DEFINE SPDE MODEL
#-----------------------------------------------------------

spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(40000, 0.5),
  prior.sigma = c(0.2, 0.5))

s.index <- inla.spde.make.index(
  name = "spatial",
  n.spde = spde$n.spde
)

#-----------------------------------------------------------
# DEFINE FOLDS
#-----------------------------------------------------------

folds <- unique(sampled_data_join_sf$country)

cv_results <- list()

all_predictions <- data.frame()

#-----------------------------------------------------------
# LEAVE-ONE-COUNTRY-OUT CROSS VALIDATION
#-----------------------------------------------------------

for(i in seq_along(folds)) {
  
  country_out <- folds[i]
  
  cat("\n=====================================\n")
  cat("LEAVING OUT:", country_out, "\n")
  cat("=====================================\n")
  
  train <- sampled_data_join_sf %>%
    filter(country != country_out)
  
  test <- sampled_data_join_sf %>%
    filter(country == country_out)
  
  train <- sf::st_drop_geometry(train)
  test <- sf::st_drop_geometry(test)
  #---------------------------------------------------------
  # PROJECTOR MATRICES
  #---------------------------------------------------------
  
  A_train <- inla.spde.make.A(
    mesh = mesh,
    loc = as.matrix(train[, c("lon","lat")])
  )
  
  A_test <- inla.spde.make.A(
    mesh = mesh,
    loc = as.matrix(test[, c("lon","lat")])
  )
  
  #---------------------------------------------------------
  # ESTIMATION STACK
  #---------------------------------------------------------
  
  stack.est <- inla.stack(
    data = list(
      y = train$Number_of_infected_tot_pm,
      Ntrials = train$Number_of_participant_tot
    ),
    
    A = list(
      A_train,
      1
    ),
    
    effects = list(
      spatial = s.index,
      
      data.frame(
        intercept            = 1,
        annual_rainfall_2021 = train$annual_rainfall_2021,
        vegetation_index     = train$vegetation_index,
        # pop_density_2020     = train$pop_density_2020,
        dist_to_water        = train$dist_to_water,
        lc_water             = train$lc_water,
        lc_trees             = train$lc_trees,
        lc_built             = train$lc_built,
        lc_cropland          = train$lc_cropland,
        lc_wetland           = train$lc_wetland
      )
    ),
    
    tag = "est"
  )
  
  #---------------------------------------------------------
  # PREDICTION STACK
  #---------------------------------------------------------
  
  stack.pred <- inla.stack(
    data = list(
      y = NA,
      Ntrials = test$Number_of_participant_tot
    ),
    
    A = list(
      A_test,
      1
    ),
    
    effects = list(
      spatial = s.index,
      
      data.frame(
        intercept            = 1,
        annual_rainfall_2021 = test$annual_rainfall_2021,
        vegetation_index     = test$vegetation_index,
        # pop_density_2020     = test$pop_density_2020,
        dist_to_water        = test$dist_to_water,
        lc_water             = test$lc_water,
        lc_trees             = test$lc_trees,
        lc_built             = test$lc_built,
        lc_cropland          = test$lc_cropland,
        lc_wetland           = test$lc_wetland
      )
    ),
    
    tag = "pred"
  )
  
  #---------------------------------------------------------
  # COMBINED STACK
  #---------------------------------------------------------
  
  stk <- inla.stack(
    stack.est,
    stack.pred
  )
  
  #---------------------------------------------------------
  # MODEL FORMULA
  #---------------------------------------------------------
  
  formula <-
    y ~ -1 + intercept +
    annual_rainfall_2021 +
    vegetation_index +
    dist_to_water +
    lc_water +
    lc_trees +
    lc_built +
    lc_cropland +
    lc_wetland +
    f(spatial, model = spde)
  
  #---------------------------------------------------------
  # FIT MODEL
  #---------------------------------------------------------
  
  fit <- inla(
    formula,
    family = "betabinomial",
    data = inla.stack.data(stk),
    
    Ntrials = inla.stack.data(stk)$Ntrials,
    
    control.predictor = list(
      A = inla.stack.A(stk),
      compute = TRUE
    ),
    control.fixed = list(
      mean.intercept = -3,
      prec.intercept = 0.8,
      mean = 0,
      prec = 0.2
    ),
    control.compute = list(
      dic = TRUE,
      waic = TRUE,
      cpo = TRUE,
      config = TRUE
    )
  )
  
  
  #=========================================================
  # POSTERIOR SAMPLES
  #=========================================================
  
  samples <- inla.posterior.sample(
    n = 200,
    result = fit
  )
  
  #=========================================================
  # PREDICTIONS
  #=========================================================
  
  idx.pred <- inla.stack.index(stk,"pred")$data
  
  
  eta_mat <- inla.posterior.sample.eval(
    function() APredictor[idx.pred],
    inla.posterior.sample(200, fit)
  )
  
  prob_mat <- plogis(eta_mat)
  
  pred_mean    <- rowMeans(prob_mat)
  pred_sd      <- apply(prob_mat, 1, sd)
  pred_lower   <- apply(prob_mat, 1, quantile, 0.025)
  pred_upper   <- apply(prob_mat, 1, quantile, 0.975)
  pred_median  <- apply(prob_mat, 1, median)
  
  
  obs_prev <- test$PR_Pm_tot 
  
  #---------------------------------------------------------
  # PERFORMANCE METRICS
  #---------------------------------------------------------
  
  rmse <- sqrt(
    mean(
      (obs_prev - pred_mean)^2,
      na.rm = TRUE
    )
  )
  
  mae <- mean(
    abs(obs_prev - pred_mean),
    na.rm = TRUE
  )
  
  bias <- mean(
    pred_mean - obs_prev,
    na.rm = TRUE
  )
  
  cor_val <- cor(
    obs_prev,
    pred_mean,
    use = "complete.obs"
  )
  
  coverage <- mean(
    obs_prev >= pred_lower &
      obs_prev <= pred_upper,
    na.rm = TRUE
  )
  
  cv_results[[i]] <- data.frame(
    Country = country_out,
    RMSE = rmse,
    MAE = mae,
    Bias = bias,
    Correlation = cor_val,
    Coverage95 = coverage
  )
  
  pred_df <- data.frame(
    Country = country_out,
    Observed = obs_prev,
    Predicted = pred_mean,
    Lower = pred_lower,
    Upper = pred_upper
  )
  
  all_predictions <- rbind(
    all_predictions,
    pred_df
  )
  
  cat(
    "Fold:", country_out,
    " n_test =", nrow(test),
    " predictions =", length(pred_mean),
    "\n"
  )
}

cv_summary <- do.call(
  rbind,
  cv_results
)

write.csv(
  cv_summary,
  "Outputs/africa_output/Spatial_CV_Summary_country.csv",
  row.names = FALSE
)

write.csv(
  all_predictions,
  "Outputs/africa_output/Spatial_CV_Predictions_country.csv",
  row.names = FALSE
)


################################################################################
#          CREATE SPATIAL FOLDS WITHIN COUNTRIES                               #   
################################################################################

set.seed(123)

sampled_data_join_sf1 <- sampled_data_join_sf %>%
  group_by(country) %>%
  group_modify(~{
    
    nloc <- nrow(.x)
    
    # Congo has only 3 locations
    if(nloc <= 3){
      
      .x$spatial_fold <- 1
      
    } else {
      
      k <- ifelse(nloc >= 18, 4, 3)
      
      km <- kmeans(
        as.matrix(
          cbind(.x$lon, .x$lat)
        ),
        centers = k,
        nstart = 100
      )
      
      .x$spatial_fold <- km$cluster
    }
    
    .x
    
  }) %>%
  ungroup()

#=========================================================
# CREATE FOLD IDENTIFIERS
#=========================================================

sampled_data_join_sf1$fold_id <- paste(
  sampled_data_join_sf1$country,
  sampled_data_join_sf1$spatial_fold,
  sep = "_"
)

#=========================================================
# EXCLUDE CONGO FOLDS FROM VALIDATION
#=========================================================

folds <- unique(
  sampled_data_join_sf1$fold_id
)

#=========================================================
# STORAGE OBJECTS
#=========================================================

all_predictions <- data.frame()

cv_results <- list()

#=========================================================
# SPATIAL BLOCK CROSS VALIDATION
#=========================================================

for(i in seq_along(folds)){
  
  fold_out <- folds[i]
  
  cat(
    "\n=====================================\n"
  )
  
  cat(
    "LEAVING OUT:",
    fold_out,
    "\n"
  )
  
  cat(
    "=====================================\n"
  )
  
  #-------------------------------------------------------
  # SPLIT DATA
  #-------------------------------------------------------
  
  train <- sampled_data_join_sf1 %>%
    filter(fold_id != fold_out)
  
  test <- sampled_data_join_sf1 %>%
    filter(fold_id == fold_out)
  
  train <- st_drop_geometry(train)
  test <- st_drop_geometry(test)
  
  #-------------------------------------------------------
  # PROJECTOR MATRICES
  #-------------------------------------------------------
  
  A_train <- inla.spde.make.A(
    mesh = mesh,
    loc = as.matrix(
      train[,c("lon","lat")]
    )
  )
  
  A_test <- inla.spde.make.A(
    mesh = mesh,
    loc = as.matrix(
      test[,c("lon","lat")]
    )
  )
  
  #-------------------------------------------------------
  # STACK ESTIMATION
  #-------------------------------------------------------
  
  stack.est <- inla.stack(
    
    data = list(
      y = train$Number_of_infected_tot_pm,
      Ntrials = train$Number_of_participant_tot
    ),
    
    A = list(
      A_train,
      1
    ),
    
    effects = list(
      
      spatial = s.index,
      
      data.frame(
        intercept            = 1,
        annual_rainfall_2021 = train$annual_rainfall_2021,
        vegetation_index     = train$vegetation_index,
        dist_to_water        = train$dist_to_water,
        lc_water             = train$lc_water,
        lc_trees             = train$lc_trees,
        lc_built             = train$lc_built,
        lc_cropland          = train$lc_cropland,
        lc_wetland           = train$lc_wetland
      )
    ),
    
    tag = "est"
  )
  
  #-------------------------------------------------------
  # STACK PREDICTION
  #-------------------------------------------------------
  
  stack.pred <- inla.stack(
    
    data = list(
      y = NA,
      Ntrials = test$Number_of_participant_tot
    ),
    
    A = list(
      A_test,
      1
    ),
    
    effects = list(
      
      spatial = s.index,
      
      data.frame(
        intercept            = 1,
        annual_rainfall_2021 = test$annual_rainfall_2021,
        vegetation_index     = test$vegetation_index,
        dist_to_water        = test$dist_to_water,
        lc_water             = test$lc_water,
        lc_trees             = test$lc_trees,
        lc_built             = test$lc_built,
        lc_cropland          = test$lc_cropland,
        lc_wetland           = test$lc_wetland
      )
    ),
    
    tag = "pred"
  )
  
  stk <- inla.stack(
    stack.est,
    stack.pred
  )
  
  #-------------------------------------------------------
  # MODEL
  #-------------------------------------------------------
  
  formula <-
    y ~ -1 +
    intercept +
    annual_rainfall_2021 +
    vegetation_index +
    dist_to_water +
    lc_water +
    lc_trees +
    lc_built +
    lc_cropland +
    lc_wetland +
    f(
      spatial,
      model = spde
    )
  
  fit <- inla(
    
    formula,
    
    family = "betabinomial",
    
    data = inla.stack.data(stk),
    
    Ntrials =
      inla.stack.data(stk)$Ntrials,
    
    control.predictor = list(
      A = inla.stack.A(stk),
      compute = TRUE
    ),
    
    control.fixed = list(
      mean.intercept = -3,
      prec.intercept = 0.8,
      mean = 0,
      prec = 0.2
    ),
    
    control.compute = list(
      dic = TRUE,
      waic = TRUE,
      cpo = TRUE,
      config = TRUE
    )
  )
  
  #=========================================================
  # POSTERIOR SAMPLES
  #=========================================================
  
  samples <- inla.posterior.sample(
    n = 200,
    result = fit
  )
  
  #=========================================================
  # PREDICTIONS
  #=========================================================
  
  idx.pred <- inla.stack.index(stk,"pred")$data
  
  
  eta_mat <- inla.posterior.sample.eval(
    function() APredictor[idx.pred],
    inla.posterior.sample(200, fit)
  )
  
  prob_mat <- plogis(eta_mat)
  
  pred_mean    <- rowMeans(prob_mat)
  pred_sd      <- apply(prob_mat, 1, sd)
  pred_lower   <- apply(prob_mat, 1, quantile, 0.025)
  pred_upper   <- apply(prob_mat, 1, quantile, 0.975)
  pred_median  <- apply(prob_mat, 1, median)
  obs_prev <- test$PR_Pm_tot
  
  #-------------------------------------------------------
  # METRICS
  #-------------------------------------------------------
  
  rmse <- sqrt(
    mean(
      (obs_prev - pred_mean)^2,
      na.rm = TRUE
    )
  )
  
  mae <- mean(
    abs(obs_prev - pred_mean),
    na.rm = TRUE
  )
  
  cor_val <- cor(
    obs_prev,
    pred_mean,
    use = "complete.obs"
  )
  
  coverage <- mean(
    obs_prev >= pred_lower &
      obs_prev <= pred_upper,
    na.rm = TRUE
  )
  
  cv_results[[i]] <- data.frame(
    Fold = fold_out,
    RMSE = rmse,
    MAE = mae,
    Correlation = cor_val,
    Coverage95 = coverage
  )
  
  pred_df <- data.frame(
    Fold = fold_out,
    Country = test$country,
    Observed = obs_prev,
    Predicted = pred_mean,
    Lower95 = pred_lower,
    Upper95 = pred_upper
  )
  
  all_predictions <- rbind(
    all_predictions,
    pred_df
  )
  
}

#=========================================================
# RESULTS
#=========================================================

cv_summary <- do.call(
  rbind,
  cv_results
)

print(cv_summary)

#=========================================================
# OVERALL PERFORMANCE
#=========================================================

overall_results <- data.frame(
  Mean_RMSE =
    mean(cv_summary$RMSE, na.rm = TRUE),
  
  Mean_MAE =
    mean(cv_summary$MAE, na.rm = TRUE),
  
  Mean_Correlation = mean(cv_summary$Correlation, na.rm = TRUE),  
  Mean_Coverage95 = mean(cv_summary$Coverage95, na.rm = TRUE))

print(overall_results)

#=========================================================
# CHECK PREDICTIONS
#=========================================================

cat("\n")
cat("Total predictions:", nrow(all_predictions), "\n")
cat("\n")

print(
  table(all_predictions$Country)
)

#=========================================================
# OBSERVED VS PREDICTED
#=========================================================

ggplot(
  all_predictions, aes(x = Observed, y = Predicted, colour = Country)) +
  geom_point(size = 3) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Within-country spatial cross-validation",
    x = "Observed prevalence",
    y = "Predicted prevalence")

#=========================================================
# SAVE RESULTS
#=========================================================

write.csv(
  cv_summary,
  "Outputs/africa_output/Spatial_CV_Summary_within_country.csv",
  row.names = FALSE
)

write.csv(
  all_predictions,
  "Outputs/africa_output/Spatial_CV_Predictions_within_country.csv",
  row.names = FALSE
)

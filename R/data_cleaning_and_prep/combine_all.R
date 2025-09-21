library(data.table)
library(tidyverse)
library(sf)
library(tidylog)
library(ggrepel)
library(MetBrewer)

# Flux data ----------------
dt_flux_t <- fread("data/processed_data/preliminary_data/prelim_fluxes.csv") %>% 
  filter(!is.na(flux_best)) %>% 
  dplyr::select(type, elevation, treatment, temperature, flux_best, tier, site, plot_id, date, year) %>% 
  unique() %>% 
  mutate(month = month(date))
summary(dt_flux_t)
table(dt_flux_t[dt_flux_t$treatment == "c", ]$tier)
table(dt_flux_t[dt_flux_t$treatment == "c" & tier == "Peru_2018", ]$month)

dt_flux <- fread("data/processed_data/preliminary_data/prelim_fluxes.csv") %>% 
  dplyr::select(type, elevation, treatment, temperature, flux_best, tier, site, plot_id, date, year, par) %>%
  unique() %>% 
  mutate(temperature_nee = ifelse(type == "nee", temperature, NA),
         temperature_reco = ifelse(type == "reco", temperature, NA), 
         month = month(date),
         country = case_when(
           grepl("China", tier) ~ "China", 
           grepl("Colorado", tier) ~ "USA", 
           grepl("Norway", tier) ~ "Norway", 
           grepl("Peru", tier) ~ "Peru", 
           grepl("South_Africa", tier) ~ "South Africa", 
           grepl("Svalbard", tier) ~ "Svalbard"
         ), 
         country = as.factor(country)) %>% 
  pivot_wider(names_from = type, values_from = flux_best, values_fn = mean) %>% 
  filter(!(tier == "Peru_2019" & month == 7)) %>% #remove the dry season
  filter(!(tier == "Svalbard_2018" & site == "ITEX")) %>% #remove itex experiment
  group_by(country, tier, elevation, treatment,
           site, plot_id) %>% 
  summarize(nee = mean(nee, na.rm = T), #in case there are multiple measurements per plot and year
            reco = mean(reco, na.rm = T), 
            temperature_nee = mean(temperature_nee, na.rm = T), 
            temperature_reco = mean(temperature_reco, na.rm = T), 
            par = mean(par, na.rm =T) ) %>% 
  mutate(plot_id = as.character(plot_id),
         nee = ifelse(tier %in% c("South_Africa_2023"), nee, nee*-1), 
         reco = ifelse(tier %in% c("South_Africa_2023"), reco, reco*-1)) %>% 
  filter(tier %in% c("China_2016", "Colorado_2016", "Svalbard_2018",
                     "Norway_2022", "Peru_2018", "South_Africa_2023")) %>%
  filter(nee > -20 & nee < 10 & reco > 0 & reco < 20) %>% #assuming all else are wrong 
  #filter(treatment == "c") %>% # peru, we also include naturally burned sites (cause fire is also part of the system e.g., in SA)
  ungroup()
table(dt_flux_t$tier)
table(dt_flux$tier)

#Colorado 2018 and Peru 2019 have the most control plots, which is why we will go with them for now 
#But for Peru 2018 we have leaf traits, so we go for this one instead. 
# And for Colorado 2016 we have lead P and only three plots less, so we go for it :)

#uniqueDates <- dt_flux %>% ungroup() %>%  dplyr::select(tier, plot_id, date) %>% unique()


dt_trait_raw <- fread("data/processed_data/preliminary_data/prelim_traitstrap.csv") %>% 
  rename(tier = Tier, 
         site = Site, 
         plot_id = PlotID)

glimpse(dt_trait_raw)
dt_trait_raw[dt_trait_raw$tier == "Norway_2022" & Trait == "cn_ratio", ]
dt_trait_mean <- dt_trait_raw %>% dplyr::select(tier, site, plot_id, Trait, mean) %>% 
  pivot_wider(values_from = mean, names_from = Trait) %>% unique() %>% 
  mutate(sla_cm2_g = ifelse(is.infinite(sla_cm2_g), NA, sla_cm2_g)) %>% 
  dplyr::select(tier, site, plot_id, plant_height_cm, dry_mass_g, ldmc, leaf_area_cm2, sla_cm2_g, wet_mass_g,
                n_percent, cn_ratio, p_percent, c_percent, cp_ratio, np_ratio) %>% 
  group_by(tier, site, plot_id) %>% 
  summarize(
    dry_mass_g = mean(dry_mass_g, na.rm = T), 
    plant_height_cm = mean(plant_height_cm, na.rm = T), 
    ldmc = mean(ldmc, na.rm = T), 
    leaf_area_cm2 = mean(leaf_area_cm2, na.rm = T), 
    sla_cm2_g = mean(sla_cm2_g, na.rm = T), 
    wet_mass_g = mean(wet_mass_g, na.rm = T), 
    n_percent = mean(n_percent, na.rm = T), 
    cn_ratio = mean(cn_ratio, na.rm = T), 
    cp_ratio = mean(cp_ratio, na.rm = T), 
    np_ratio = mean(np_ratio, na.rm = T), 
    p_percent = mean(p_percent, na.rm = T), 
    c_percent = mean(c_percent, na.rm = T)
  ) %>% filter(!grepl("NA", plot_id))

table(dt_trait_mean[!is.na(dt_trait_mean$n_percent),]$tier)


dt_trait_mean %>% as.data.table()

dt_trait_mean
glimpse(dt_trait_mean)

plots_with_traits <- unique(dt_trait_mean[dt_trait_mean$plot_id %in% c(unique(dt_flux_t$plot_id)), ]$plot_id)
plots_without_traits <- unique(dt_flux_t[!dt_flux_t$plot_id %in% plots_with_traits, ]$plot_id)
n_distinct(dt_flux_t$plot_id)

scale(as.numeric(dt_trait_mean$sla_cm2_g))

#spatial predictors -------------------
dt_sps <- fread("data/processed_data/preliminary_data/spatial_predictors.csv") 

#remotely sensed covariates ---------------------
dt_rs_not_sa <- fread("data/environmental_data/predictors.csv") %>% 
  mutate(plot_id = as.character(plot_id), 
         date = as_date(date)) %>% 
  dplyr::select(-"flux_id") %>% 
  filter(!grepl("SA", plot_id))

dt_flux_dates <- fread("data/processed_data/preliminary_data/prelim_fluxes.csv") 
sa_dates <- dt_flux_dates[grepl("SA", dt_flux_dates$plot_id), c("plot_id", "date")] %>% unique() %>% mutate(date = as_date(date))

dt_rs_sa <- fread("data/environmental_data/predictors.csv") %>% 
  mutate(plot_id = as.character(plot_id), 
         date = as_date(date)) %>% 
  dplyr::select(-"flux_id", -"date") %>% 
  filter(grepl("SA", plot_id)) %>% 
  left_join(sa_dates)

dt_rs <- rbind(dt_rs_not_sa, dt_rs_sa) %>% 
  dplyr::select(T_diff1, T_diff7, T_diff14,
                T_diff28, plot_id) %>% 
  group_by(plot_id) %>% 
  summarize(
    temp_diff_day = mean(T_diff1, na.rm = T), 
    temp_diff_week = mean(T_diff7, na.rm = T), 
    temp_diff_fortnight = mean(T_diff14, na.rm = T), 
    temp_iff_month = mean(T_diff28, na.rm = T), 
    
  )
hist(dt_rs$temp_diff_week)

#get lon and lat --------------------
dt_loc <- read_sf("data/processed_data/preliminary_data/prelim_flux_loc.gpkg")

coords <- st_coordinates(dt_loc)
dt_coords <- dt_loc %>%
  as.data.table() %>% 
  mutate(lon = coords[,1], 
         lat = coords[,2], 
         geom = NULL) %>% 
  dplyr::select(plot_id, lon, lat) %>% 
  unique() 

# Add spatial blocks -------------------

p <- dt_coords %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

sf_use_s2(TRUE)
pb <- p %>% 
  st_buffer(50)

sf_use_s2(FALSE)

sbb <- st_union(pb) %>% 
  st_cast("POLYGON") %>% 
  st_as_sf() %>% 
  rownames_to_column("spatial_block")

dt_coords$spatial_block <- st_intersects(p, sbb) %>% as.numeric()

# get climate -----------------------

dt_clim <- fread("data/environmental_data/climate.csv") %>%
  rename(MAP = annual_precipitation, 
         MAT = annual_mean_temperature) %>% 
  mutate(MMP = MAP/12,
         npp_temp = 3000 / (1 + exp(1.315 - 0.119 * MAT)),
         npp_prec = 3000 * (1 - exp(-0.000664 * (MMP * 12))),
         miami_npp = pmin(npp_temp, npp_prec)) %>%
  group_by(plot_id) %>%
  summarize(
    mmp = mean(MMP, na.rm = TRUE),
    mat = mean(MAT, na.rm = TRUE),
    miami_npp = mean(miami_npp, na.rm = TRUE)
  )

## get biomass proxy (cover sum x height) ---------------------------

dt_ch <- fread("data/processed_data/preliminary_data/prelim_coverXheight.csv") %>% 
  group_by(plot_id) %>%
  summarize(
    cover_sum = mean(coverSum, na.rm = TRUE),
    veg_height = mean(vegHeight, na.rm = TRUE),
    height_x_cover = mean(HeightXCover , na.rm = TRUE)
  )

## get soil moisture ------------------------

dt_moist <- fread("data/environmental_data/soil_moist_temp.csv") %>% 
  mutate(Year = year(date)) %>% 
  group_by(plot_id) %>% 
  summarize(soil_moisture = mean(soil_moisture, na.rm = T)) 

n_distinct(dt_moist$plot_id)

## Woodiness & grassiness ------------------------

dt_wg <- fread("data/processed_data/preliminary_data/prelim_woodiness.csv") %>% 
  group_by(plot_id) %>% 
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE))) %>% 
  filter(plot_id != "")

plot(dt_wg$grassiness, dt_wg$woodiness)

## Species Richness ------------------------

dt_sr <- fread("data/processed_data/preliminary_data/prelim_species_richness.csv") %>% 
  group_by(tier, plot_id) %>% 
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE))) %>% 
  filter(plot_id != "") %>% 
  rename(species_richness = SpeciesRichness, 
         veg_plot_size_m2 = VegPlotSizeM2)

# Functional Diversity 

dt_fd <- fread("data/processed_data/preliminary_data/plant_functional_diversity.csv")


### combine everything but pcas --------------------
dt_mod_1 <- dt_flux %>%
  left_join(dt_trait_mean) %>%
 # left_join(dt_trait_var) %>%
 # left_join(dt_pca) %>%
  left_join(dt_clim) %>% 
  left_join(dt_moist) %>% 
  left_join(dt_fd) %>% 
  left_join(dt_ch) %>% 
  left_join(dt_sps) %>% 
  left_join(dt_rs) %>% 
  left_join(dt_wg) %>% 
  left_join(dt_sr) %>% 
  left_join(dt_coords) %>% 
  filter(!is.na(sla_cm2_g)) %>% 
  as.data.table() %>% 
  unique() %>% 
  mutate(
    gpp = reco-nee, 
    gpp = ifelse(gpp < 0, 0, gpp), 
    gpp_reco_ratio = gpp/reco, 
    temperature_mean = case_when(
      is.na(temperature_nee) ~ temperature_reco, 
      is.na(temperature_reco) ~ temperature_nee, 
      .default = (temperature_nee + temperature_reco)/2
    ),
    np_ratio = n_percent/p_percent,
    unique_id = paste0("flux", 1:nrow(.)), 
    country = as.factor(country), 
    gpp_reco_ratio = ifelse(abs(gpp_reco_ratio) > 10, NA, gpp_reco_ratio),
    abs_lat = abs(lat), 
    reco_gpp_ratio = reco/gpp,
    reco_gpp_ratio = ifelse(abs(reco_gpp_ratio) > 10, NA, reco_gpp_ratio),
    map = mmp*12)  %>% 
  mutate(gradient = case_when(
    country == "South Africa" ~ "Drakensberg", 
    country == "Peru" ~ "Central Andes", 
    country == "China" ~ "Eastern Himalayas", 
    country == "USA" ~ "Rocky Mountains", 
    country == "Norway" ~ "Southern Scandes", 
    country == "Svalbard" ~ "Svalbard" 
  ))


## PCAs -----------------------
#dt_mod %>% group_by(site) %>% summarize(mean_h = mean(plant_height_cm, na.rm = T)) %>% print(n = 50)

#morphological traits 
dt_m_t <- dt_mod_1 %>%
  ungroup() %>% 
  select(lat, gradient, plot_id, sla_cm2_g, ldmc, leaf_area_cm2, dry_mass_g, plant_height_cm) %>%
  filter(complete.cases(.)) %>% 
  unique() %>% 
  filter(plot_id != "") %>% 
  mutate(across(where(is.numeric), ~ as.numeric(scale(.))))
pr_m_t <- princomp(dt_m_t %>% select(-plot_id, -gradient, -lat))

dt_m_t$morph_traits_pc1 <- pr_m_t$scores[,1]
dt_m_t$morph_traits_pc2 <- pr_m_t$scores[,2]

#viz 
loadings_m_t <- as_tibble(pr_m_t$loadings[, 1:2], rownames = "Variable")
scores_m_t <- as_tibble(pr_m_t$scores[, 1:2]) %>%
  mutate(gradient = dt_m_t$gradient, 
         lat = dt_m_t$lat) %>% 
  mutate(gradient = fct_reorder(gradient, lat))

#importance 
#https://stackoverflow.com/questions/62422277/how-to-obtain-principal-component-variance-explained-in-r-prcomp-and-prepro
summ_m_t <- summary(pr_m_t)
summ_m_t
#Alternative (same results)
eigenvalues_m_t <- pr_m_t$sdev^2
variance_explained_m_t <- eigenvalues_m_t / sum(eigenvalues_m_t)
(percent_var_m_t <- round(variance_explained_m_t[1:2] * 100, 1))


p_m_t <- ggplot() +
  geom_point(data = scores_m_t, aes(x = Comp.1, y = Comp.2, color = gradient),  size = 1, alpha = 0.75) +
  geom_segment(data = loadings_m_t, aes(x = 0, y = 0, xend = Comp.1*5, yend = Comp.2*5),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_text_repel(data = loadings_m_t, aes(x = Comp.1*5, y = Comp.2*5, label = Variable),
                  size = 4, color = "black") +
  scale_color_met_d(name = "Archambault") +
  scale_fill_met_d(name = "Archambault") +
  labs(title = "a)", 
       subtitle = "PCA Morphological Traits",
       x = paste0("PC1 (", percent_var_m_t[1], "% variance)"),
       y = paste0("PC2 (", percent_var_m_t[2], "% variance)")) +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0, size = 14, face = "bold"),
        panel.grid = element_line(color = "snow2"), 
        axis.text = element_text(size = 12), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        panel.background = element_rect(fill = "snow", color = NA))
p_m_t

#chemical traits 
dt_c_t <- dt_mod_1 %>%
  ungroup() %>% 
  select(plot_id, n_percent, cn_ratio, c_percent, p_percent, cp_ratio, np_ratio) %>%
  filter(complete.cases(.)) %>% 
  unique() %>% 
  filter(plot_id != "") %>% 
  mutate(across(where(is.numeric),~ as.numeric(scale(.))))
pr_c_t <- princomp(dt_c_t %>% select(-plot_id))

dt_c_t$chem_traits_pc1 <- pr_c_t$scores[,1]
dt_c_t$chem_traits_pc2 <- pr_c_t$scores[,2]

#all traits 
dt_a_t <- dt_mod_1 %>%
  ungroup() %>% 
  select(lat, gradient, plot_id, sla_cm2_g, ldmc, leaf_area_cm2, dry_mass_g, plant_height_cm,
         n_percent, cn_ratio, c_percent, p_percent, cp_ratio, np_ratio) %>%
  filter(complete.cases(.)) %>% 
  unique() %>%
  filter(plot_id != "") %>% 
  mutate(across(where(is.numeric), ~ as.numeric(scale(.))))
pr_a_t <- princomp(dt_a_t %>% select(-plot_id, -lat, -gradient))

dt_a_t$all_traits_pc1 <- pr_a_t$scores[,1]
dt_a_t$all_traits_pc2 <- pr_a_t$scores[,2]

loadings_a_t <- as_tibble(pr_a_t$loadings[, 1:2], rownames = "Variable")

s_a_alibi_point <- data.table(
  Comp.1 = 0, 
  Comp.2 = 0,
  gradient = "Drakensberg", 
  lat = -2
)

norway_alibi_point <- data.table(
  Comp.1 = 0, 
  Comp.2 = 0,
  gradient = "Southern Scandes", 
  lat = 1
)

dt_a_t[dt_a_t$gradient == "Southern Scandes", ]
dt_a_t[dt_a_t$gradient == "Svalbard", ]$lat
dt_a_t[dt_a_t$gradient == "Rocky Mountains", ]$lat
dt_a_t[dt_a_t$gradient == "Central Andes", ]$lat



scores_a_t <- as_tibble(pr_a_t$scores[, 1:2]) %>%
  mutate(gradient = dt_a_t$gradient, 
         lat = dt_a_t$lat) %>% 
  rbind(s_a_alibi_point) %>% 
  rbind(norway_alibi_point) %>% 
  mutate(gradient = fct_reorder(gradient, lat))

#importance 
#https://stackoverflow.com/questions/62422277/how-to-obtain-principal-component-variance-explained-in-r-prcomp-and-prepro
summ_a_t <- summary(pr_a_t)
summ_a_t
#Alternative (same results)
eigenvalues_a_t <- pr_a_t$sdev^2
variance_explained_a_t <- eigenvalues_a_t / sum(eigenvalues_a_t)
(percent_var_a_t <- round(variance_explained_a_t[1:2] * 100, 1))


p_a_t <- ggplot() +
  geom_point(data = scores_a_t, aes(x = Comp.1, y = Comp.2, color = gradient),  size = 1, alpha = 0.75) +
  geom_segment(data = loadings_a_t, aes(x = 0, y = 0, xend = Comp.1*5, yend = Comp.2*5),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = loadings_a_t, aes(x = Comp.1*5, y = Comp.2*5, label = Variable),
                  size = 4, color = "black") +
  scale_color_met_d(name = "Archambault") +
  scale_fill_met_d(name = "Archambault") +
  labs(title = "b)",
       subtitle = "PCA All Traits",
       x = paste0("PC1 (", percent_var_a_t[1], "% variance)"),
       y = paste0("PC2 (", percent_var_a_t[2], "% variance)")) +
  theme_minimal() +
  labs(color = "Gradient") +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = 0, size = 14, face = "bold"),
        panel.grid = element_line(color = "snow2"), 
        axis.text = element_text(size = 12), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        panel.background = element_rect(fill = "snow", color = NA))
p_a_t

library(patchwork)
p_pca <- p_m_t | p_a_t
ggsave(plot = p_pca, "builds/plots/supplement/pca_viz.png", dpi = 600, height = 5, width = 10)
#compile PCA dataframe 

dt_pca <- dt_m_t %>% 
  select(plot_id, morph_traits_pc1, morph_traits_pc2) %>% 
  unique() %>% 
  left_join(dt_a_t[, c("plot_id", "all_traits_pc1", "all_traits_pc2")]) %>% 
  unique() %>% 
  left_join(dt_c_t[, c("plot_id", "chem_traits_pc1", "chem_traits_pc2")]) %>% 
  distinct() 


table(dt_a_t$plot_id)
# Combine actually everything -------

dt_mod <- dt_mod_1 %>% 
  left_join(dt_pca) %>% 
 #get country level means 
  group_by(country) %>% 
  mutate(
    gpp_country_mean = mean(gpp, na.rm = TRUE), 
    reco_country_mean = mean(reco, na.rm = TRUE),
    nee_country_mean = mean(nee, na.rm = TRUE),
    mmp_country_mean = mean(mmp, na.rm = TRUE), 
    map_country_mean = mean(map, na.rm = TRUE), 
    mat_country_mean = mean(mat, na.rm = TRUE), 
    sla_country_mean = mean(sla_cm2_g, na.rm = TRUE), 
    leaf_area_country_mean = mean(leaf_area_cm2, na.rm = TRUE),
    plant_height_country_mean = mean(plant_height_cm, na.rm = TRUE),
    ldmc_country_mean = mean(ldmc, na.rm = TRUE), 
    # sla_sd_country_mean = mean(sla_cm2_g_sd, na.rm = TRUE), 
    # leaf_area_sd_country_mean = mean(leaf_area_cm2_sd, na.rm = TRUE),
    # ldmc_sd_country_mean = mean(ldmc_sd, na.rm = TRUE), 
    height_x_cover_country_mean = mean(height_x_cover, na.rm = TRUE), 
    temperature_nee_country_mean = mean(temperature_nee, na.rm = TRUE), 
    temperature_mean_country_mean = mean(temperature_mean, na.rm = TRUE), 
    temperature_reco_country_mean = mean(temperature_reco, na.rm = TRUE), 
    elevation_country_mean = mean(elevation, na.rm = TRUE),
    
    # PCA scores: traits
    all_traits_pc1_country_mean = mean(all_traits_pc1, na.rm = TRUE),
    all_traits_pc2_country_mean = mean(all_traits_pc2, na.rm = TRUE),
    chem_traits_pc1_country_mean = mean(chem_traits_pc1, na.rm = TRUE),
    chem_traits_pc2_country_mean = mean(chem_traits_pc2, na.rm = TRUE),
    morph_traits_pc1_country_mean = mean(morph_traits_pc1, na.rm = TRUE),
    morph_traits_pc2_country_mean = mean(morph_traits_pc2, na.rm = TRUE),
    
    # Alternative hypotheses
    n_percent_country_mean = mean(n_percent, na.rm = TRUE),
    p_percent_country_mean = mean(p_percent, na.rm = TRUE),
    cn_ratio_country_mean = mean(cn_ratio, na.rm = TRUE),
    cp_ratio_country_mean = mean(cp_ratio, na.rm = TRUE),
    np_ratio_country_mean = mean(np_ratio, na.rm = TRUE),
    c_percent_country_mean = mean(c_percent, na.rm = TRUE),
    par_country_mean = mean(par, na.rm = TRUE),
    soil_moisture_country_mean = mean(soil_moisture, na.rm = TRUE),
    woodiness_country_mean = mean(woodiness, na.rm = TRUE),
    grassiness_country_mean = mean(grassiness, na.rm = TRUE),
    species_richness_country_mean = mean(species_richness, na.rm = TRUE), 
    functional_diversity_q1_country_mean = mean(functional_diversity_q1, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  mutate(
    gpp_anomaly_country = gpp - gpp_country_mean,
    reco_anomaly_country = reco - reco_country_mean,
    nee_anomaly_country = nee - nee_country_mean,
    sla_anomaly_country = sla_cm2_g - sla_country_mean,
    leaf_area_anomaly_country = leaf_area_cm2 - leaf_area_country_mean,
    plant_height_anomaly_country = plant_height_cm - plant_height_country_mean,
    ldmc_country_anomaly_country = ldmc - ldmc_country_mean,
    # sla_sd_anomaly_country = sla_cm2_g_sd - sla_sd_country_mean,
    # leaf_area_sd_anomaly_country = leaf_area_cm2_sd - leaf_area_sd_country_mean,
    # ldmc_sd_anomaly_country = ldmc_sd - ldmc_sd_country_mean,
    height_x_cover_anomaly_country = height_x_cover - height_x_cover_country_mean,
    temperature_nee_anomaly_country = temperature_nee - temperature_nee_country_mean,
    temperature_mean_anomaly_country = temperature_mean - temperature_mean_country_mean,
    temperature_reco_anomaly_country = temperature_reco - temperature_reco_country_mean, 
    map_anomaly_country = map - map_country_mean, 
    mat_anomaly_country = mat - mat_country_mean,
    elevation_anomaly_country = elevation - elevation_country_mean, 

    # PCA scores: traits anomalies
    all_traits_pc1_anomaly_country = all_traits_pc1 - all_traits_pc1_country_mean,
    all_traits_pc2_anomaly_country = all_traits_pc2 - all_traits_pc2_country_mean,
    chem_traits_pc1_anomaly_country = chem_traits_pc1 - chem_traits_pc1_country_mean,
    chem_traits_pc2_anomaly_country = chem_traits_pc2 - chem_traits_pc2_country_mean,
    morph_traits_pc1_anomaly_country = morph_traits_pc1 - morph_traits_pc1_country_mean,
    morph_traits_pc2_anomaly_country = morph_traits_pc2 - morph_traits_pc2_country_mean, 
    
    #Alternative hypotheses 
    species_richness_anomaly_country = species_richness - species_richness_country_mean,
    functional_diversity_q1_anomaly_country = functional_diversity_q1 - functional_diversity_q1_country_mean,
    n_percent_anomaly_country = n_percent - n_percent_country_mean,
    p_percent_anomaly_country = p_percent - p_percent_country_mean,
    cn_ratio_anomaly_country = cn_ratio - cn_ratio_country_mean,
    cp_ratio_anomaly_country = cp_ratio - cp_ratio_country_mean,
    np_ratio_anomaly_country = np_ratio - np_ratio_country_mean,
    c_percent_anomaly_country = c_percent - c_percent_country_mean,
    par_anomaly_country = par - par_country_mean,
    soil_moisture_anomaly_country = soil_moisture - soil_moisture_country_mean,
    woodiness_anomaly_country = woodiness - woodiness_country_mean,
    grassiness_anomaly_country = grassiness - grassiness_country_mean
    

  ) %>% 
  rename(
    temperature_gpp = temperature_mean, 
    temperature_gpp_anomaly_country = temperature_mean_anomaly_country, 
    temperature_gpp_country_mean = temperature_mean_country_mean, 
  )
dt_mod[(is.na(dt_mod$morph_traits_pc1_anomaly_country)), ]$plot_id


fwrite(dt_mod, "data/processed_data/clean_data/global_fluxes_main_data.csv")


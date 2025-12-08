library(data.table)
library(tidyverse)
library(ggcorrplot)
library(car)
library(GGally)


dt_raw <- fread("data/processed_data/clean_data/global_fluxes_main_data.csv") %>% 
  dplyr::select(
    # identifiers 
    country, site, plot_id, treatment, gradient, 
    
    # fluxes
    nee, reco, gpp,
    
    # environmental 
    elevation, map, downscaled_temp,
    temperature_nee, temperature_reco, temperature_gpp,
    
    # trait means 
    sla_cm2_g, ldmc, leaf_area_cm2, dry_mass_g, plant_height_cm,
    n_percent, cn_ratio, cp_ratio, np_ratio, c_percent, p_percent,
    
    # others 
    height_x_cover, species_richness,
    functional_diversity_q1, lat, par_nee,
    
    #pca axis 
    all_traits_pc1, all_traits_pc2,
    chem_traits_pc1, chem_traits_pc2,
    morph_traits_pc1, morph_traits_pc2,

    # Anomalies
    par_nee_anomaly_country,
    gpp_anomaly_country,
    nee_anomaly_country,
    reco_anomaly_country,
    temperature_gpp_anomaly_country,
    temperature_nee_anomaly_country,
    temperature_reco_anomaly_country,
    downscaled_temp_anomaly_country,
    elevation_anomaly_country,
    
    leaf_area_anomaly_country, 
    sla_anomaly_country, 
    plant_height_anomaly_country,
    height_x_cover_anomaly_country, 
    all_traits_pc1_anomaly_country, all_traits_pc2_anomaly_country,
    chem_traits_pc1_anomaly_country, chem_traits_pc2_anomaly_country,
    morph_traits_pc1_anomaly_country, morph_traits_pc2_anomaly_country,
    
    functional_diversity_q1_anomaly_country,
    species_richness_anomaly_country,
    n_percent_anomaly_country,
    p_percent_anomaly_country,
    cn_ratio_anomaly_country, cp_ratio_anomaly_country, np_ratio_anomaly_country,
    c_percent_anomaly_country,
    soil_moisture_anomaly_country,
    woodiness_anomaly_country,
    grassiness_anomaly_country
  ) #%>% filter(complete.cases(.))

### PCA and traits ------------------------

morph_traits_data <- dt_raw %>%
  select(
    PC1 = morph_traits_pc1,
    PC2 = morph_traits_pc2,
    SLA = sla_cm2_g,
    LDMC = ldmc,
    `Leaf Area` = leaf_area_cm2,
    `Dry Mass` = dry_mass_g,
    Height = plant_height_cm
  )

p_morph_trait_pairs <- ggpairs(morph_traits_data,
        title = "Morphological Traits",
        upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("smooth", alpha = 0.3)),
        diag = list(continuous = "densityDiag")) +
  theme_minimal() +
  theme(legend.position = "none",
        legend.box="vertical",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid = element_line(color = "snow"), 
        #axis.title.x = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.text.x = element_text(size = 10, angle = 0, hjust = 1), 
        panel.background = element_rect(fill = "snow", color = NA), 
        strip.text.x = element_text(size = 12), 
        strip.text.y = element_text(size = 12), 
        strip.background = element_rect(fill = "linen", color = "linen") )

p_morph_trait_pairs
ggsave(plot = p_morph_trait_pairs, "builds/plots/supplement/trait_pca_and_morph_traits_pairs.png",
       dpi = 600, height = 8, width = 8)



# all traits
all_traits_data <- dt_raw %>%
  select(
    PC1         = all_traits_pc1,
    PC2         = all_traits_pc2,
    N          = n_percent,
    `C:N`    = cn_ratio,
    `C:P`    = cp_ratio,
    `N:P`    = np_ratio,
    C          = c_percent,
    P          = p_percent,
    SLA         = sla_cm2_g,
    LDMC        = ldmc,
    `Leaf Area` = leaf_area_cm2,
    `Dry Mass`  = dry_mass_g,
    Height      = plant_height_cm
  )

p_all_trait_pairs <- ggpairs(all_traits_data,
                              title = "All Traits",
                              upper = list(continuous = wrap("cor", size = 3)),
                              lower = list(continuous = wrap("smooth", alpha = 0.3)),
                              diag = list(continuous = "densityDiag")) +
  theme_minimal() +
  theme(legend.position = "none",
        legend.box="vertical",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid = element_line(color = "snow"), 
        #axis.title.x = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1), 
        panel.background = element_rect(fill = "snow", color = NA), 
        strip.text.x = element_text(size = 10), 
        strip.text.y = element_text(size = 10), 
        strip.background = element_rect(fill = "linen", color = "linen") )

p_all_trait_pairs
ggsave(plot = p_all_trait_pairs, "builds/plots/supplement/trait_pca_and_all_traits_pairs.png",
       dpi = 600, height = 10, width = 10)


## correlation of "regular" variables ------
dt_corr_reg <- dt_raw %>% 
  dplyr::select(
    # Fluxes
    NEE  = nee,
    Reco = reco,
    GPP  = gpp,
    
    # Environmental
    Elevation = elevation,
    MAP = map,
    `Growing Season Temp.` = downscaled_temp,
    `PAR`  = par_nee,
    `Temp NEE` = temperature_nee,
    `Temp Reco` = temperature_reco,
    `Temp GPP` = temperature_gpp,
    
    # Trait means
    SLA = sla_cm2_g,
    LDMC = ldmc,
    `Leaf Area` = leaf_area_cm2,
    `Dry Mass`= dry_mass_g,
    Height = plant_height_cm,
    
    # Others
    `Biomass` = height_x_cover
  ) %>% 
  filter(complete.cases(.))


corr_reg <- round(cor(dt_corr_reg), 1)
p_reg <- ggcorrplot(corr_reg, hc.order = F, type = "lower",
                       lab = TRUE) +
  labs(title = "Full Dataset")

p_reg

#all traits 
dt_corr_reg_at <- dt_raw %>% 
  dplyr::select(
    # Fluxes
    NEE  = nee,
    Reco = reco,
    GPP  = gpp,
    # Environmental
    Elevation   = elevation,
    `Downscaled Temp` = downscaled_temp,
    `Temp NEE` = temperature_nee,
    `Temp Reco` = temperature_reco,
    `Temp GPP` = temperature_gpp,
    `PAR`  = par_nee,
    
    # Trait means
    SLA = sla_cm2_g,
    LDMC = ldmc,
    `Leaf Area` = leaf_area_cm2,
    `Dry Mass`  = dry_mass_g,
    Height = plant_height_cm,
    N  = n_percent,
    `C:N` = cn_ratio,
    `C:P`= cp_ratio,
    `N:P` = np_ratio,
    C = c_percent,
    P = p_percent,
    # Others
    `Biomass` = height_x_cover
  ) %>% 
  filter(complete.cases(.))

corr_reg_at <- round(cor(dt_corr_reg_at), 1)
p_reg_at <- ggcorrplot(corr_reg_at, hc.order = F, type = "lower",
           lab = TRUE) +
  labs(title = "Subset With Leaf Chemical Traits ")

p_reg_at

library(patchwork)

p_reg_corr <- p_reg / p_reg_at

ggsave(plot = p_reg_corr, "builds/plots/supplement/correlation_regular_vars.png", dpi = 600, height = 15, width = 10)

# correlation of country level anomalies variables ------
dt_corr_country <- dt_raw %>% 
  dplyr::select(
    # Flux anomalies
    `GPP Anomaly` = gpp_anomaly_country,
    `NEE Anomaly` = nee_anomaly_country,
    `Reco Anomaly` = reco_anomaly_country,
    `Temp GPP Anomaly` = temperature_gpp_anomaly_country,
    `Temp NEE Anomaly` = temperature_nee_anomaly_country,
    `Temp Reco Anomaly` = temperature_reco_anomaly_country,
    `Growing Season Temp Anomaly`  = downscaled_temp_anomaly_country,
    `Elevation Anomaly`  = elevation_anomaly_country,
    `PAR Anomaly`  = par_nee_anomaly_country,
    
    
    # Trait anomalies
    `Leaf Area Anomaly` = leaf_area_anomaly_country,
    `SLA Anomaly`= sla_anomaly_country,
    `Biomass Anomaly`   = height_x_cover_anomaly_country,
  ) %>% 
  filter(complete.cases(.))

names(dt_corr_country) <- gsub("country", "gradient", names(dt_corr_country))

corr_country <- round(cor(dt_corr_country), 1)
p_country <- ggcorrplot(corr_country, hc.order = F, type = "lower",
                     lab = TRUE) +
  labs(title = "Full Dataset")
p_country

dt_corr_country_at <- dt_raw %>% 
  dplyr::select(
    # Flux anomalies
    `GPP Anomaly` = gpp_anomaly_country,
    `NEE Anomaly` = nee_anomaly_country,
    `Reco Anomaly` = reco_anomaly_country,
    `Temp GPP Anomaly` = temperature_gpp_anomaly_country,
    `Temp NEE Anomaly` = temperature_nee_anomaly_country,
    `Temp Reco Anomaly` = temperature_reco_anomaly_country,
    `Growing Season Temp Anomaly`  = downscaled_temp_anomaly_country,
    `Elevation Anomaly` = elevation_anomaly_country,
    `PAR Anomaly` = par_nee_anomaly_country,
    
    # Trait anomalies
    `Leaf Area Anomaly` = leaf_area_anomaly_country,
    `SLA Anomaly` = sla_anomaly_country,
    `Biomass Anomaly` = height_x_cover_anomaly_country,
    `N Anomaly` = n_percent_anomaly_country,
    `P Anomaly` = p_percent_anomaly_country,
    `Height Anomaly` = plant_height_anomaly_country,
    `C:N Anomaly`  = cn_ratio_anomaly_country,
    `C:P Anomaly` = cp_ratio_anomaly_country,
    `N:P Anomaly` = np_ratio_anomaly_country,
    `C Anomaly` = c_percent_anomaly_country
  ) %>% 
  filter(complete.cases(.))


names(dt_corr_country_at) <- gsub("country", "gradient", names(dt_corr_country_at))


corr_country_at <- round(cor(dt_corr_country_at), 1)
p_country_at <- ggcorrplot(corr_country_at, hc.order = F, type = "lower",
                        lab = TRUE) +
  labs(title = "Subset With Leaf Chemical Traits ")
p_country_at


p_country_corr <- p_country / p_country_at
ggsave(plot = p_country_corr, "builds/plots/supplement/correlation_country_level_vars.png", dpi = 600, height = 15, width = 10)


cor.test(dt_raw$height_x_cover_anomaly_country, dt_raw$plant_height_anomaly_country, na.rm = T)

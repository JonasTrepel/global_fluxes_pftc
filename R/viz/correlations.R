library(data.table)
library(tidyverse)
library(ggcorrplot)
library(car)
library(GGally)


dt_raw <- fread("data/processed_data/clean_data/global_fluxes_main_data.csv") %>% 
  dplyr::select(
    # identifiers 
    country, site, plot_id, treatment,
    
    # fluxes
    nee, reco, gpp,
    
    # environmental 
    elevation, map, mat,
    temperature_nee, temperature_reco, temperature_gpp,
    
    # trait means 
    sla_cm2_g, ldmc, leaf_area_cm2, dry_mass_g, plant_height_cm,
    n_percent, cn_ratio, cp_ratio, np_ratio, c_percent, p_percent,
    
    # others 
    height_x_cover, species_richness,
    functional_diversity_q1, lat,
    
    #pca axis 
    all_traits_pc1, all_traits_pc2,
    chem_traits_pc1, chem_traits_pc2,
    morph_traits_pc1, morph_traits_pc2,
    
    # Anomalies
    gpp_anomaly_country,
    nee_anomaly_country,
    reco_anomaly_country,
    temperature_gpp_anomaly_country,
    temperature_nee_anomaly_country,
    temperature_reco_anomaly_country,
    mat_anomaly_country,
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
    par_anomaly_country,
    soil_moisture_anomaly_country,
    woodiness_anomaly_country,
    grassiness_anomaly_country
  ) #%>% filter(complete.cases(.))

### PCA and traits ------------------------

morph_traits_data <- dt_raw %>%
  select(morph_traits_pc1, morph_traits_pc2, #traits_pc3, traits_pc4,
         sla_cm2_g, ldmc, leaf_area_cm2, dry_mass_g, plant_height_cm)

p_morph_trait_pairs <- ggpairs(morph_traits_data,
        title = "Morphological Traits",
        upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("smooth", alpha = 0.3)),
        diag = list(continuous = "densityDiag"))

p_morph_trait_pairs
ggsave(plot = p_morph_trait_pairs, "builds/plots/supplement/trait_pca_and_morph_traits_pairs.png",
       dpi = 600, height = 8, width = 8)

# Select variable pcas vs traits + height × cover
chem_traits_data <- dt_raw %>%
  select(chem_traits_pc1, chem_traits_pc2, 
         n_percent, cn_ratio, cp_ratio, np_ratio, c_percent, p_percent)

p_chem_trait_pairs <- ggpairs(chem_traits_data,
                               title = "Chemical Traits",
                               upper = list(continuous = wrap("cor", size = 3)),
                               lower = list(continuous = wrap("smooth", alpha = 0.3)),
                               diag = list(continuous = "densityDiag"))

p_chem_trait_pairs
ggsave(plot = p_chem_trait_pairs, "builds/plots/supplement/trait_pca_and_chem_traits_pairs.png",
       dpi = 600, height = 8, width = 8)

# all traits
all_traits_data <- dt_raw %>%
  select(all_traits_pc1, all_traits_pc2, 
         n_percent, cn_ratio, cp_ratio, np_ratio, c_percent, p_percent, 
         sla_cm2_g, ldmc, leaf_area_cm2, dry_mass_g, plant_height_cm)

p_all_trait_pairs <- ggpairs(all_traits_data,
                              title = "All Traits",
                              upper = list(continuous = wrap("cor", size = 3)),
                              lower = list(continuous = wrap("smooth", alpha = 0.3)),
                              diag = list(continuous = "densityDiag"))

p_all_trait_pairs
ggsave(plot = p_all_trait_pairs, "builds/plots/supplement/trait_pca_and_all_traits_pairs.png",
       dpi = 600, height = 10, width = 10)


## correlation of "regular" variables ------
dt_corr_reg <- dt_raw %>% 
  dplyr::select(
    # fluxes
    nee, reco, gpp,
    
    # environmental 
    elevation, map, mat,
    temperature_nee, temperature_reco, temperature_gpp,
    
    # trait means 
    sla_cm2_g, ldmc, leaf_area_cm2, dry_mass_g, plant_height_cm,
    n_percent, cn_ratio, cp_ratio, np_ratio, c_percent, p_percent,
    
    # others 
     height_x_cover, species_richness,
     functional_diversity_q1, 

    ) %>% filter(complete.cases(.))

corr_reg <- round(cor(dt_corr_reg), 1)
p_reg <- ggcorrplot(corr_reg, hc.order = F, type = "lower",
           lab = TRUE)

p_reg

ggsave(plot = p_reg, "builds/plots/supplement/correlation_regular_vars.png", dpi = 600, height = 10, width = 10)

# correlation of country level anomalies variables ------
dt_corr_country <- dt_raw %>% 
  dplyr::select(# fluxes
    gpp_anomaly_country,
    nee_anomaly_country,
    reco_anomaly_country,
    temperature_gpp_anomaly_country,
    temperature_nee_anomaly_country,
    temperature_reco_anomaly_country,
    mat_anomaly_country,
    elevation_anomaly_country,
    
    leaf_area_anomaly_country, 
    sla_anomaly_country, 
    height_x_cover_anomaly_country, 
    n_percent_anomaly_country,
    p_percent_anomaly_country,
    plant_height_anomaly_country,
    cn_ratio_anomaly_country, 
    cp_ratio_anomaly_country,
    np_ratio_anomaly_country,
    c_percent_anomaly_country) %>% 
  filter(complete.cases(.))

corr_country <- round(cor(dt_corr_country), 1)
p_country <- ggcorrplot(corr_country, hc.order = F, type = "lower",
                     lab = TRUE)
p_country
ggsave(plot = p_country, "builds/plots/supplement/correlation_country_level_vars.png", dpi = 600, height = 10, width = 10)



cor.test(dt_raw$height_x_cover_anomaly_country, dt_raw$plant_height_anomaly_country, na.rm = T)

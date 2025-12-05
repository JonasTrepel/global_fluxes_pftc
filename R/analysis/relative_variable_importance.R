library(glmm.hp)
library(tidyverse)
library(data.table)
library(glmmTMB)

dt_raw <- fread("data/processed_data/clean_data/global_fluxes_main_data.csv") %>% 
  dplyr::select(
    # identifiers 
    country, site, plot_id,# treatment,
    
    # fluxes
    nee, reco, gpp,
    
    # environmental 
    downscaled_temp, map, mat,
    temperature_nee, temperature_reco, temperature_gpp,
    
    # trait means 
    sla_cm2_g, ldmc, leaf_area_cm2, dry_mass_g, plant_height_cm,
    n_percent, cn_ratio, c_percent, p_percent,
    
    # others 
    height_x_cover, species_richness,
    functional_diversity_q1, lat,
    
    #pca axis 
    all_traits_pc1,  all_traits_pc2,
    chem_traits_pc1, chem_traits_pc2,
    morph_traits_pc1,  morph_traits_pc2,
    
    # Anomalies
    gpp_anomaly_country,
    nee_anomaly_country,
    reco_anomaly_country,
    temperature_gpp_anomaly_country,
    temperature_nee_anomaly_country,
    temperature_reco_anomaly_country,
    downscaled_temp_anomaly_country,
    downscaled_temp_anomaly_country,
    
    leaf_area_anomaly_country, 
    sla_anomaly_country, 
    height_x_cover_anomaly_country, 
    all_traits_pc1_anomaly_country,  all_traits_pc2_anomaly_country,
    chem_traits_pc1_anomaly_country, chem_traits_pc2_anomaly_country,
    morph_traits_pc1_anomaly_country,  morph_traits_pc2_anomaly_country,
    
    # functional_diversity_q1_anomaly_country,
    # species_richness_t1_anomaly_country,
    # # n_percent_anomaly_country,
    # # p_percent_anomaly_country,
    # # cn_ratio_anomaly_country,
    # # c_percent_anomaly_country,
     par_nee_anomaly_country,
    # soil_moisture_anomaly_country,
    # woodiness_t1_anomaly_country,
    # grassiness_t1_anomaly_country
  ) %>% filter(!is.na(morph_traits_pc1_anomaly_country))

m_dat <- dt_raw %>% 
  dplyr::select(-c(nee, reco, gpp, lat)) %>% 
  mutate(leaf_area_cm2 = log(leaf_area_cm2),
         species_richness = log(species_richness),
         functional_diversity_q1 = log(functional_diversity_q1),
         height_x_cover = log(height_x_cover)
         
  ) %>% 
  mutate(across(where(is.numeric), ~as.numeric(scale(.x)))) %>% 
  left_join(dt_raw[, c("plot_id", "nee", "reco", "gpp", "lat")]) %>% 
  distinct()





## Prepare models -------

fm_t4_nee_1 <- glmmTMB(nee_anomaly_country ~
                         temperature_nee_anomaly_country +
                         downscaled_temp_anomaly_country +
                         height_x_cover_anomaly_country +
                         morph_traits_pc1_anomaly_country +
                         morph_traits_pc2_anomaly_country +
                         par_nee_anomaly_country +
                         ( 1 | site), na.action = na.omit,
                       data = m_dat)
summary(fm_t4_nee_1)

fm_t4_nee_2 <- glmmTMB(nee_anomaly_country ~
                         temperature_nee_anomaly_country +
                         downscaled_temp_anomaly_country +
                         height_x_cover_anomaly_country +
                         leaf_area_anomaly_country +
                         sla_anomaly_country +
                         par_nee_anomaly_country +
                         ( 1 | site), na.action = na.omit,
                       data = m_dat)
summary(fm_t4_nee_2)

fm_t4_nee_3 <- glmmTMB(nee_anomaly_country ~
                         temperature_nee_anomaly_country +
                         downscaled_temp_anomaly_country +
                         height_x_cover_anomaly_country +
                         all_traits_pc1_anomaly_country +
                         all_traits_pc2_anomaly_country +
                         par_nee_anomaly_country +
                         ( 1 | site), na.action = na.omit,
                       data = m_dat)
summary(fm_t4_nee_3)



#Reco

fm_t4_reco_1 <- glmmTMB(reco_anomaly_country ~
                         temperature_reco_anomaly_country +
                         downscaled_temp_anomaly_country +
                         height_x_cover_anomaly_country +
                         morph_traits_pc1_anomaly_country +
                         morph_traits_pc2_anomaly_country +
                         ( 1 | site), na.action = na.omit,
                       data = m_dat)
summary(fm_t4_reco_1)

fm_t4_reco_2 <- glmmTMB(reco_anomaly_country ~
                         temperature_reco_anomaly_country +
                         downscaled_temp_anomaly_country +
                         height_x_cover_anomaly_country +
                         leaf_area_anomaly_country +
                         sla_anomaly_country +
                         ( 1 | site), na.action = na.omit,
                       data = m_dat)
summary(fm_t4_reco_2)

fm_t4_reco_3 <- glmmTMB(reco_anomaly_country ~
                         temperature_reco_anomaly_country +
                         downscaled_temp_anomaly_country +
                         height_x_cover_anomaly_country +
                         all_traits_pc1_anomaly_country +
                         all_traits_pc2_anomaly_country +
                         ( 1 | site), na.action = na.omit,
                       data = m_dat)
summary(fm_t4_reco_3)

#Gpp
fm_t4_gpp_1 <- glmmTMB(gpp_anomaly_country ~
                         temperature_gpp_anomaly_country +
                         downscaled_temp_anomaly_country +
                         height_x_cover_anomaly_country +
                         morph_traits_pc1_anomaly_country +
                         morph_traits_pc2_anomaly_country +
                         par_nee_anomaly_country + 
                         ( 1 | site), na.action = na.omit,
                       data = m_dat)
summary(fm_t4_gpp_1)

fm_t4_gpp_2 <- glmmTMB(gpp_anomaly_country ~
                         temperature_gpp_anomaly_country +
                         downscaled_temp_anomaly_country +
                         height_x_cover_anomaly_country +
                         leaf_area_anomaly_country +
                         sla_anomaly_country +
                         par_nee_anomaly_country + 
                         ( 1 | site), na.action = na.omit,
                       data = m_dat)
summary(fm_t4_gpp_2)

fm_t4_gpp_3 <- glmmTMB(gpp_anomaly_country ~
                         temperature_gpp_anomaly_country +
                         downscaled_temp_anomaly_country +
                         height_x_cover_anomaly_country +
                         all_traits_pc1_anomaly_country +
                         all_traits_pc2_anomaly_country +
                         par_nee_anomaly_country + 
                         ( 1 | site), na.action = na.omit,
                       data = m_dat)
summary(fm_t4_gpp_3)

#extract AIC etc 
model_list_fm <- list(
  
  fm_t4_nee_1 = fm_t4_nee_1,
  fm_t4_nee_2 = fm_t4_nee_2,
  fm_t4_nee_3 = fm_t4_nee_3,
  
  fm_t4_reco_1 = fm_t4_reco_1,
  fm_t4_reco_2 = fm_t4_reco_2,
  fm_t4_reco_3 = fm_t4_reco_3,
  
  fm_t4_gpp_1 = fm_t4_gpp_1,
  fm_t4_gpp_2 = fm_t4_gpp_2,
  fm_t4_gpp_3 = fm_t4_gpp_3
  
)

dt_res_fm <- data.frame()

for(i in 1:length(model_list_fm)){
  
  
  m <- model_list_fm[[i]]
  
  m_name <- names(model_list_fm)[i]
  
  v_par <- glmm.hp(m)
  
  tmp_dt <- v_par$hierarchical.partitioning %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "predictor") %>% 
    mutate( 
    model_name = m_name,
    rsq_m = as.numeric(r.squaredGLMM(m)[1]), 
    rsq_c = as.numeric(r.squaredGLMM(m)[2]))
  
  dt_res_fm <- rbind(dt_res_fm, tmp_dt)
  
  print(paste0(i, " done"))
}

dt_vp <- dt_res_fm %>% 
  mutate(
    predictor_tier = case_when(
      .default = "not_applicable",
      grepl("_1", model_name) ~ "morph_traits_pca",
      grepl("_2", model_name) ~ "traits_separately", 
      grepl("_3", model_name) ~ "all_traits_pca"), 
    response = case_when(
      grepl("reco", model_name) ~ "Reco",
      grepl("nee", model_name) ~ "NEE",
      grepl("gpp", model_name)~ "GPP"), 
  clean_term =  case_when(
    .default = predictor,
    grepl("temperature_", predictor) ~ "Inst. Temperature",
    predictor == "downscaled_temp_anomaly_country" ~ "Growing Season Temp.",
    predictor == "height_x_cover_anomaly_country"      ~ "'Biomass'",
    predictor == "par_nee_anomaly_country"      ~ "PAR",
    predictor == "morph_traits_pc1_anomaly_country"          ~ "Morph. Traits PC1",
    predictor == "morph_traits_pc2_anomaly_country"          ~ "Morph. Traits PC2",
    predictor == "all_traits_pc1_anomaly_country"          ~ "All Traits PC1",
    predictor == "all_traits_pc2_anomaly_country"          ~ "All Traits PC2",
    predictor == "sla_anomaly_country"                 ~ "SLA",
    predictor == "leaf_area_anomaly_country"           ~ "Leaf Area"), 
  label = paste0(response, "\n(R²m = ", round(rsq_m, 2), ")")) %>% 
  rename(var_imp = `I.perc(%)`) %>% 
  mutate(
    clean_term = factor(clean_term,
                        levels = c("'Biomass'",
                                   "Leaf Area",
                                   "SLA",
                                   "Chem. Traits PC2", "Chem. Traits PC1",
                                   "All Traits PC2", "All Traits PC1",
                                   "Morph. Traits PC2", "Morph. Traits PC1",
                                   "Inst. Temperature",
                                   "Growing Season Temp.",
                                   "PAR")))



p_vp_ts <- dt_vp %>% 
  filter(predictor_tier == "traits_separately") %>% 
  ggplot() +
  geom_col(aes(x = var_imp, y = clean_term)) +
  facet_wrap(~label) +
  labs(x = "Percentage of Marginal R²", y = "") +
  theme_minimal() +
  theme(legend.position = "none",
        legend.box="vertical",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid = element_line(color = "snow"), 
        #axis.title.x = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.text.x = element_text(size = 10, angle = 0, hjust = 1), 
        panel.background = element_rect(fill = "snow", color = NA), 
        strip.text.x = element_text(size = 12), 
        strip.text.y = element_text(size = 12, face = "bold"), 
        strip.background = element_rect(fill = "linen", color = "linen") )

p_vp_ts

ggsave(plot = p_vp_ts, "builds/plots/variable_importance_traits_separately.png", dpi = 600, height = 3, width = 10)

p_vp_mt <- dt_vp %>% 
  filter(predictor_tier == "morph_traits_pca") %>% 
  ggplot() +
  geom_col(aes(x = var_imp, y = clean_term)) +
  facet_wrap(~label) +
  labs(x = "Percentage of Marginal R²", y = "") +
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
        strip.text.y = element_text(size = 12, face = "bold"), 
        strip.background = element_rect(fill = "linen", color = "linen") )

p_vp_mt


p_vp_at <- dt_vp %>% 
  filter(predictor_tier == "all_traits_pca") %>% 
  ggplot() +
  geom_col(aes(x = var_imp, y = clean_term)) +
  facet_wrap(~label) +
  labs(x = "Percentage of Marginal R²", y = "") +
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
        strip.text.y = element_text(size = 12, face = "bold"), 
        strip.background = element_rect(fill = "linen", color = "linen") )

p_vp_at
ggsave(plot = p_vp_at, "builds/plots/supplement/variable_importance_all_traits_pca.png", dpi = 600, height = 3, width = 7)

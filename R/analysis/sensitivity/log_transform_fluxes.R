library(broom.mixed)
library(MuMIn)
library(DHARMa)
library(tidyverse)
library(data.table)
library(glmmTMB)


### We may want to check how estimates change if we use -1/kT and ln(flux) instead 

dt_raw <- fread("data/processed_data/clean_data/global_fluxes_main_data.csv") %>% 
  dplyr::select(
    # identifiers 
    country, site, plot_id, gradient,
    
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

    par_nee_anomaly_country,
  ) %>% filter(!is.na(morph_traits_pc1_anomaly_country))

hist(dt$nee)
hist(dt$reco)
hist(dt$gpp)


dt_raw %>% 
  mutate() %>% 
  ggplot() +
  geom_point(aes(x = log(reco), y = temperature_reco))

dt_mod <- dt_raw %>% 
  # Create absolute temperature (K) and transformed temperature for NEE, RECO, GPP
  mutate(
    temp_k_nee  = temperature_nee  + 273.15,
    temp_k_reco = temperature_reco + 273.15,
    temp_k_gpp  = temperature_gpp  + 273.15,
    
    temp_k_nee_trans  = -1 / ((1.38e-23) * temp_k_nee),
    temp_k_reco_trans = -1 / ((1.38e-23) * temp_k_reco),
    temp_k_gpp_trans  = -1 / ((1.38e-23) * temp_k_gpp)
  ) %>%
  group_by(gradient) %>% 
  mutate(
    temp_k_nee_trans_mean  = mean(temp_k_nee_trans,  na.rm = TRUE),
    temp_k_reco_trans_mean = mean(temp_k_reco_trans, na.rm = TRUE),
    temp_k_gpp_trans_mean  = mean(temp_k_gpp_trans,  na.rm = TRUE),
    
    temp_k_nee_trans_anomaly  = temp_k_nee_trans  - temp_k_nee_trans_mean,
    temp_k_reco_trans_anomaly = temp_k_reco_trans - temp_k_reco_trans_mean,
    temp_k_gpp_trans_anomaly  = temp_k_gpp_trans  - temp_k_gpp_trans_mean,
    
    positive_nee_anomaly_country  = nee_anomaly_country  + abs(min(nee_anomaly_country,  na.rm = TRUE)) + 1,
    positive_reco_anomaly_country = reco_anomaly_country + abs(min(reco_anomaly_country, na.rm = TRUE)) + 1,
    positive_gpp_anomaly_country  = gpp_anomaly_country  + abs(min(gpp_anomaly_country,  na.rm = TRUE)) + 1
  ) %>% 
  ungroup()

library(cowplot)

(p_nee_flux <- dt_mod %>% 
  ggplot() +
  facet_wrap(~gradient, ncol = 6) +
  geom_point(aes(x = positive_nee_anomaly_country,y = nee_anomaly_country), alpha = 0.7) +
  labs(title = "NEE", x = "Positive NEE anomaly", y = "NEE anomaly"))

(p_reco_flux <- dt_mod %>% 
  ggplot() +
  facet_wrap(~gradient, ncol = 6) +
  geom_point(aes(x = positive_reco_anomaly_country, y = reco_anomaly_country), alpha = 0.7) +
  labs(title = "Reco", x = "Positive Reco anomaly", y = "Reco anomaly"))

(p_gpp_flux <- dt_mod %>% 
  ggplot() +
  facet_wrap(~gradient, ncol = 6) +
  geom_point(aes(x = positive_gpp_anomaly_country, y = gpp_anomaly_country), alpha = 0.7) +
  labs(title = "GPP",x = "Positive GPP anomaly",y = "GPP anomaly"))

flux_panel <- plot_grid(
  p_nee_flux, p_reco_flux, p_gpp_flux,
  ncol = 1,
)
flux_panel

(p_nee_temp <- dt_mod %>% 
    ggplot() +
    facet_wrap(~gradient, ncol = 6) +
    geom_point(aes(x = temp_k_nee_trans_anomaly,y = temperature_nee_anomaly_country),alpha = 0.7) +
    labs(title = "NEE", x = "Transformed temperature anomaly (NEE)",y = "Temperature anomaly (NEE)"))

(p_reco_temp <- dt_mod %>% 
    ggplot() +
    facet_wrap(~gradient, ncol = 6) +
    geom_point(aes(x = temp_k_reco_trans_anomaly,y = temperature_reco_anomaly_country),alpha = 0.7) +
    labs(title = "Reco", x = "Transformed temperature anomaly (Reco)",y = "Temperature anomaly (Reco)"))

(p_gpp_temp <- dt_mod %>% 
    ggplot() +
    facet_wrap(~gradient, ncol = 6) +
    geom_point(aes(x = temp_k_gpp_trans_anomaly, y = temperature_gpp_anomaly_country), alpha = 0.7) +
    labs(title = "GPP", x = "Transformed temperature anomaly (GPP)", y = "Temperature anomaly (GPP)"))

temp_panel <- plot_grid(
  p_nee_temp, p_reco_temp, p_gpp_temp,
  ncol = 1
)

temp_panel

m_dat <- dt_mod %>% 
  dplyr::select(-c(nee, reco, gpp, lat)) %>% 
  mutate(leaf_area_cm2 = log(leaf_area_cm2),
         species_richness = log(species_richness),
         functional_diversity_q1 = log(functional_diversity_q1),
         height_x_cover = log(height_x_cover), 
         log_positive_nee_anomaly_country = log(positive_nee_anomaly_country),
         log_positive_gpp_anomaly_country = log(positive_gpp_anomaly_country),
         log_positive_reco_anomaly_country = log(positive_reco_anomaly_country)
  ) %>% 
  mutate(across(where(is.numeric), ~as.numeric(scale(.x)))) %>% 
  left_join(dt_raw[, c("plot_id", "nee", "reco", "gpp", "lat")]) %>% 
  distinct()

m_dat %>% 
  ggplot(aes(x = log_positive_nee_anomaly_country, y = positive_nee_anomaly_country)) +
  geom_point() +
  facet_wrap(~gradient)
## Prepare models -------
# NEEE
m_nee_t1_1 <- glmmTMB(positive_nee_anomaly_country ~
                         temperature_nee_anomaly_country +
                         par_nee_anomaly_country +
                         ( 1 | site), na.action = na.omit,
                       data = m_dat)
summary(m_nee_t1_1)

m_nee_t1_2 <- glmmTMB(log_positive_nee_anomaly_country ~
                        temp_k_nee_trans_anomaly +
                        par_nee_anomaly_country +
                        ( 1 | site), na.action = na.omit,
                      data = m_dat)
summary(m_nee_t1_2)

m_nee_t4_1 <- glmmTMB(positive_nee_anomaly_country ~
                     temperature_nee_anomaly_country +
                     downscaled_temp_anomaly_country +
                     height_x_cover_anomaly_country +
                     leaf_area_anomaly_country +
                     sla_anomaly_country +
                     par_nee_anomaly_country +
                     ( 1 | site), na.action = na.omit,
                   data = m_dat)
summary(m_nee_t4_1)

m_nee_t4_2 <- glmmTMB(log_positive_nee_anomaly_country ~
                        temp_k_nee_trans_anomaly +
                        downscaled_temp_anomaly_country +
                        height_x_cover_anomaly_country +
                        leaf_area_anomaly_country +
                        sla_anomaly_country +
                        par_nee_anomaly_country +
                        ( 1 | site), na.action = na.omit,
                      data = m_dat)
summary(m_nee_t4_2)


# GPP
m_gpp_t1_1 <- glmmTMB(positive_gpp_anomaly_country ~
                        temperature_gpp_anomaly_country +
                        par_nee_anomaly_country +
                        ( 1 | site), na.action = na.omit,
                      data = m_dat)
summary(m_gpp_t1_1)

m_gpp_t1_2 <- glmmTMB(log_positive_gpp_anomaly_country ~
                        temp_k_gpp_trans_anomaly +
                        par_nee_anomaly_country +
                        ( 1 | site), na.action = na.omit,
                      data = m_dat)
summary(m_gpp_t1_2)

m_gpp_t4_1 <- glmmTMB(positive_gpp_anomaly_country ~
                        temperature_gpp_anomaly_country +
                        downscaled_temp_anomaly_country +
                        height_x_cover_anomaly_country +
                        leaf_area_anomaly_country +
                        sla_anomaly_country +
                        par_nee_anomaly_country +
                        ( 1 | site), na.action = na.omit,
                      data = m_dat)
summary(m_gpp_t4_1)

m_gpp_t4_2 <- glmmTMB(log_positive_gpp_anomaly_country ~
                        temp_k_gpp_trans_anomaly +
                        downscaled_temp_anomaly_country +
                        height_x_cover_anomaly_country +
                        leaf_area_anomaly_country +
                        sla_anomaly_country +
                        par_nee_anomaly_country +
                        ( 1 | site), na.action = na.omit,
                      data = m_dat)
summary(m_gpp_t4_2)

# Reco
m_reco_t1_1 <- glmmTMB(positive_reco_anomaly_country ~
                        temperature_reco_anomaly_country +
                        ( 1 | site), na.action = na.omit,
                      data = m_dat)
summary(m_reco_t1_1)

m_reco_t1_2 <- glmmTMB(log_positive_reco_anomaly_country ~
                        temp_k_reco_trans_anomaly +
                        ( 1 | site), na.action = na.omit,
                      data = m_dat)
summary(m_reco_t1_2)

m_reco_t4_1 <- glmmTMB(positive_reco_anomaly_country ~
                        temperature_reco_anomaly_country +
                        downscaled_temp_anomaly_country +
                        height_x_cover_anomaly_country +
                        leaf_area_anomaly_country +
                        sla_anomaly_country +
                        ( 1 | site), na.action = na.omit,
                      data = m_dat)
summary(m_reco_t4_1)

m_reco_t4_2 <- glmmTMB(log_positive_reco_anomaly_country ~
                        temp_k_reco_trans_anomaly +
                        downscaled_temp_anomaly_country +
                        height_x_cover_anomaly_country +
                        leaf_area_anomaly_country +
                        sla_anomaly_country +
                        ( 1 | site), na.action = na.omit,
                      data = m_dat)
summary(m_reco_t4_2)


#extract estimates and Rsq etc 
model_list_fm <- list(
  
  m_nee_t1_1 = m_nee_t1_1,
  m_nee_t1_2 = m_nee_t1_2,
  m_nee_t4_1 = m_nee_t4_1,
  m_nee_t4_2 = m_nee_t4_2,
  
  m_gpp_t1_1 = m_gpp_t1_1,
  m_gpp_t1_2 = m_gpp_t1_2,
  m_gpp_t4_1 = m_gpp_t4_1,
  m_gpp_t4_2 = m_gpp_t4_2,
  
  m_reco_t1_1 = m_reco_t1_1,
  m_reco_t1_2 = m_reco_t1_2,
  m_reco_t4_1 = m_reco_t4_1,
  m_reco_t4_2 = m_reco_t4_2
  
)

dt_res_fm <- data.frame()

for(i in 1:length(model_list_fm)){
  
  
  m <- model_list_fm[[i]]
  
  m_name <- names(model_list_fm)[i]
  
  tmp_tidy <- tidy(m) %>% 
    filter(effect == "fixed")
  
  tmp_dt <- tmp_tidy %>% 
    mutate( 
      model_name = m_name,
      rsq_m = as.numeric(r.squaredGLMM(m)[1,1]), 
      rsq_c = as.numeric(r.squaredGLMM(m)[1,2]), 
      ci_lb = estimate-1.96*std.error, 
      ci_ub = estimate+1.96*std.error)
  
  dt_res_fm <- rbind(dt_res_fm, tmp_dt)
  
  print(paste0(i, " done"))
}

dt_p <- dt_res_fm %>% 
  mutate(
    predictor_tier = case_when(
      .default = "not_applicable",
      grepl("t1_1", model_name) ~ "simple_og",
      grepl("t1_2", model_name) ~ "simple_transformed",
      grepl("t4_1", model_name) ~ "full_og",
      grepl("t4_2", model_name) ~ "full_transformed",
    ), 
    clean_tier = case_when(
      .default = "not_applicable",
      grepl("t1_1", model_name) ~ "Original",
      grepl("t1_2", model_name) ~ "Transformed",
      grepl("t4_1", model_name) ~ "Original",
      grepl("t4_2", model_name) ~ "Transformed",
    ),
    response = case_when(
      grepl("reco", model_name) ~ "Reco",
      grepl("nee", model_name) ~ "NEE",
      grepl("gpp", model_name)~ "GPP"), 
    clean_term =  case_when(
      .default = term,
      term == "downscaled_temp_anomaly_country" ~ "Growing Season Temp.",
      grepl("temp", term) ~ "Inst. Temperature",
      term == "height_x_cover_anomaly_country"      ~ "'Biomass'",
      term == "par_nee_anomaly_country"      ~ "PAR",
      term == "morph_traits_pc1_anomaly_country"          ~ "Morph. Traits PC1",
      term == "morph_traits_pc2_anomaly_country"          ~ "Morph. Traits PC2",
      term == "all_traits_pc1_anomaly_country"          ~ "All Traits PC1",
      term == "all_traits_pc2_anomaly_country"          ~ "All Traits PC2",
      term == "sla_anomaly_country"                 ~ "SLA",
      term == "leaf_area_anomaly_country"           ~ "Leaf Area"), 
    label = paste0(response, "\n(R²m = ", round(rsq_m, 2), ")")) %>% 
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

p_s = dt_p %>% 
  as.data.table() %>% 
  filter(grepl("simple", predictor_tier)) %>% 
  filter(!grepl("ntercept", term)) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_pointrange(aes(x = estimate, y = clean_term, xmin = ci_lb, xmax = ci_ub, color = clean_tier),
                  position = position_dodge(width = .5)) +
  labs(y = NULL, x = "Estimate", color = "Tier") +
  facet_wrap(~response)
p_s  


p_f = dt_p %>% 
  as.data.table() %>% 
  filter(!grepl("simple", predictor_tier)) %>% 
  filter(!grepl("ntercept", term)) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_pointrange(aes(x = estimate, y = clean_term, xmin = ci_lb, xmax = ci_ub, color = clean_tier),
                  position = position_dodge(width = .5)) +
  labs(y = NULL, x = "Estimate", color = "Tier") +
  facet_wrap(~response)
p_f 


library(patchwork)
(p_comb = p_s / p_f)
ggsave(plot = p_comb, "builds/plots/supplement/log_transformed_fluxes.png", dpi = 600, height = 6, width = 8)

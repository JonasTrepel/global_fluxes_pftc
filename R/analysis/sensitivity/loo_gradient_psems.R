library(scico)
library("piecewiseSEM")
library(data.table)
library(tidyverse)
library(lme4)
library(MuMIn)
library(glmmTMB)
library("semPlot")
library(performance)
library(gridExtra)
library(DHARMa)

set.seed(161)

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
    par_nee_anomaly_country,
    
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
    # par_nee_anomaly_country,
    # soil_moisture_anomaly_country,
    # woodiness_t1_anomaly_country,
    # grassiness_t1_anomaly_country
  ) %>% filter(!is.na(morph_traits_pc1_anomaly_country))

dt <- dt_raw %>% 
  dplyr::select(-c(nee, reco, gpp, lat, 
                   gpp_anomaly_country, nee_anomaly_country, reco_anomaly_country)) %>% 
  mutate(leaf_area_cm2 = log(leaf_area_cm2),
         species_richness = log(species_richness),
         functional_diversity_q1 = log(functional_diversity_q1),
         height_x_cover = log(height_x_cover)
         
  ) %>% 
  mutate(across(where(is.numeric), ~as.numeric(scale(.x)))) %>% 
  left_join(dt_raw[, c("plot_id", "nee", "reco", "gpp", "lat", 
                       "gpp_anomaly_country", "nee_anomaly_country", "reco_anomaly_country")]) %>% 
  distinct()


cor.test(dt$par_nee_anomaly_country, dt$temperature_nee_anomaly_country)

gradients = unique(dt$gradient)
dt_res <- data.frame()
for(grad in gradients){

dt_sub = dt %>% 
  filter(!gradient %in% c(grad))
  
  
# 1. NEE ------------------------------------------

m_t4_nee_1 <- psem(
  
  # model for temperature
  glmmTMB(temperature_nee_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_sub),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_sub),
  
  # model for sla / ldmc axis 
  glmmTMB(morph_traits_pc1_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_sub),
  
  # model for lead area / dry mass axis 
  glmmTMB(morph_traits_pc2_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_sub),
  
  #model for flux
  glmmTMB(nee_anomaly_country ~
            temperature_nee_anomaly_country +
            height_x_cover_anomaly_country +
            morph_traits_pc1_anomaly_country + 
            morph_traits_pc2_anomaly_country +
            downscaled_temp_anomaly_country +
            par_nee_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_sub),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% morph_traits_pc1_anomaly_country,
  height_x_cover_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  morph_traits_pc1_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  temperature_nee_anomaly_country %~~% height_x_cover_anomaly_country,
  temperature_nee_anomaly_country %~~% morph_traits_pc1_anomaly_country,
  temperature_nee_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  downscaled_temp_anomaly_country %~~% par_nee_anomaly_country,
  temperature_nee_anomaly_country %~~% par_nee_anomaly_country,
  height_x_cover_anomaly_country %~~% par_nee_anomaly_country,
  morph_traits_pc1_anomaly_country %~~% par_nee_anomaly_country,
  morph_traits_pc2_anomaly_country %~~% par_nee_anomaly_country,
  data = dt_sub
)
# ss_t4_nee_1 <- summary(m_t4_nee_1)
# ss_t4_nee_1
# dSep(m_t4_nee_1)
# LLchisq(m_t4_nee_1)
# AIC(m_t4_nee_1)
# plot(m_t4_nee_1)
# anova(m_t4_nee_1)
# 

## Traits separately --------------------
m_t4_nee_2 <- psem(
  
  # model for temperature
  glmmTMB(temperature_nee_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_sub),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_sub),
  
  # model for sla 
  glmmTMB(sla_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_sub),
  
  # model for lead area
  glmmTMB(leaf_area_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_sub),
  
  #model for flux
  glmmTMB(nee_anomaly_country ~
            temperature_nee_anomaly_country +
            height_x_cover_anomaly_country +
            sla_anomaly_country + 
            leaf_area_anomaly_country +
            downscaled_temp_anomaly_country +
            par_nee_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_sub),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% sla_anomaly_country,
  height_x_cover_anomaly_country %~~% leaf_area_anomaly_country,
  sla_anomaly_country %~~% leaf_area_anomaly_country,
  temperature_nee_anomaly_country %~~% height_x_cover_anomaly_country,
  temperature_nee_anomaly_country %~~% sla_anomaly_country,
  temperature_nee_anomaly_country %~~% leaf_area_anomaly_country,
  downscaled_temp_anomaly_country %~~% par_nee_anomaly_country,
  temperature_nee_anomaly_country %~~% par_nee_anomaly_country,
  height_x_cover_anomaly_country %~~% par_nee_anomaly_country,
  leaf_area_anomaly_country %~~% par_nee_anomaly_country,
  sla_anomaly_country %~~% par_nee_anomaly_country,
  
  data = dt_sub
)
# ss_t4_nee_2 <- summary(m_t4_nee_2)
# ss_t4_nee_2
# dSep(m_t4_nee_2)
# LLchisq(m_t4_nee_2)
# AIC(m_t4_nee_2)
# plot(m_t4_nee_2)
# anova(m_t4_nee_2)

# 2. Reco ------------------------------------------
## Moprhological Trait PCA --------------------
m_t4_reco_1 <- psem(
  
  # model for temperature
  glmmTMB(temperature_reco_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_sub),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_sub),
  
  # model for sla / ldmc axis 
  glmmTMB(morph_traits_pc1_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_sub),
  
  # model for lead area / dry mass axis 
  glmmTMB(morph_traits_pc2_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_sub),
  
  #model for flux
  glmmTMB(reco_anomaly_country ~
            temperature_reco_anomaly_country +
            height_x_cover_anomaly_country +
            morph_traits_pc1_anomaly_country + 
            morph_traits_pc2_anomaly_country +
            downscaled_temp_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_sub),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% morph_traits_pc1_anomaly_country,
  height_x_cover_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  morph_traits_pc1_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  temperature_reco_anomaly_country %~~% height_x_cover_anomaly_country,
  temperature_reco_anomaly_country %~~% morph_traits_pc1_anomaly_country,
  temperature_reco_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  
  data = dt_sub
)
# ss_t4_reco_1 <- summary(m_t4_reco_1)
# ss_t4_reco_1
# dSep(m_t4_reco_1)
# LLchisq(m_t4_reco_1)
# AIC(m_t4_reco_1)
# plot(m_t4_reco_1)
# anova(m_t4_reco_1)


## Traits separately --------------------
m_t4_reco_2 <- psem(
  
  # model for temperature
  glmmTMB(temperature_reco_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_sub),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_sub),
  
  # model for sla 
  glmmTMB(sla_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_sub),
  
  # model for lead area
  glmmTMB(leaf_area_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_sub),
  
  #model for flux
  glmmTMB(reco_anomaly_country ~
            temperature_reco_anomaly_country +
            height_x_cover_anomaly_country +
            sla_anomaly_country + 
            leaf_area_anomaly_country +
            downscaled_temp_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_sub),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% sla_anomaly_country,
  height_x_cover_anomaly_country %~~% leaf_area_anomaly_country,
  sla_anomaly_country %~~% leaf_area_anomaly_country,
  temperature_reco_anomaly_country %~~% height_x_cover_anomaly_country,
  temperature_reco_anomaly_country %~~% sla_anomaly_country,
  temperature_reco_anomaly_country %~~% leaf_area_anomaly_country,
  
  data = dt_sub
)
# ss_t4_reco_2 <- summary(m_t4_reco_2)
# ss_t4_reco_2
# dSep(m_t4_reco_2)
# LLchisq(m_t4_reco_2)
# AIC(m_t4_reco_2)
# plot(m_t4_reco_2)
# anova(m_t4_reco_2)



# 3. GPP ------------------------------------------
## Moprhological Trait PCA --------------------
m_t4_gpp_1 <- psem(
  
  # model for temperature
  glmmTMB(temperature_gpp_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_sub),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_sub),
  
  # model for sla / ldmc axis 
  glmmTMB(morph_traits_pc1_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_sub),
  
  # model for lead area / dry mass axis 
  glmmTMB(morph_traits_pc2_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_sub),
  
  #model for flux
  glmmTMB(gpp_anomaly_country ~
            temperature_gpp_anomaly_country +
            height_x_cover_anomaly_country +
            morph_traits_pc1_anomaly_country + 
            morph_traits_pc2_anomaly_country +
            downscaled_temp_anomaly_country +
            par_nee_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_sub),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% morph_traits_pc1_anomaly_country,
  height_x_cover_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  morph_traits_pc1_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  temperature_gpp_anomaly_country %~~% height_x_cover_anomaly_country,
  temperature_gpp_anomaly_country %~~% morph_traits_pc1_anomaly_country,
  temperature_gpp_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  downscaled_temp_anomaly_country %~~% par_nee_anomaly_country,
  temperature_gpp_anomaly_country %~~% par_nee_anomaly_country,
  height_x_cover_anomaly_country %~~% par_nee_anomaly_country,
  morph_traits_pc1_anomaly_country %~~% par_nee_anomaly_country,
  morph_traits_pc2_anomaly_country %~~% par_nee_anomaly_country,
  
  data = dt_sub
)
# ss_t4_gpp_1 <- summary(m_t4_gpp_1)
# ss_t4_gpp_1
# dSep(m_t4_gpp_1)
# LLchisq(m_t4_gpp_1)
# AIC(m_t4_gpp_1)
# plot(m_t4_gpp_1)
# anova(m_t4_gpp_1)


## Traits separately --------------------
m_t4_gpp_2 <- psem(
  
  # model for temperature
  glmmTMB(temperature_gpp_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_sub),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_sub),
  
  # model for sla 
  glmmTMB(sla_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_sub),
  
  # model for lead area
  glmmTMB(leaf_area_anomaly_country ~ 
            downscaled_temp_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_sub),
  
  #model for flux
  glmmTMB(gpp_anomaly_country ~
            temperature_gpp_anomaly_country +
            height_x_cover_anomaly_country +
            sla_anomaly_country + 
            leaf_area_anomaly_country +
            downscaled_temp_anomaly_country +
            par_nee_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_sub),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% sla_anomaly_country,
  height_x_cover_anomaly_country %~~% leaf_area_anomaly_country,
  sla_anomaly_country %~~% leaf_area_anomaly_country,
  temperature_gpp_anomaly_country %~~% height_x_cover_anomaly_country,
  temperature_gpp_anomaly_country %~~% sla_anomaly_country,
  temperature_gpp_anomaly_country %~~% leaf_area_anomaly_country,
  downscaled_temp_anomaly_country %~~% par_nee_anomaly_country,
  temperature_gpp_anomaly_country %~~% par_nee_anomaly_country,
  height_x_cover_anomaly_country %~~% par_nee_anomaly_country,
  sla_anomaly_country %~~% par_nee_anomaly_country,
  leaf_area_anomaly_country %~~% par_nee_anomaly_country,
  
  data = dt_sub
)
# ss_t4_gpp_2 <- summary(m_t4_gpp_2)
# ss_t4_gpp_2
# dSep(m_t4_gpp_2)
# LLchisq(m_t4_gpp_2)
# AIC(m_t4_gpp_2)
# plot(m_t4_gpp_2)
# anova(m_t4_gpp_2)

## Extract model estimates ---------------------------------------

model_list <- list(
  m_t4_nee_1 = m_t4_nee_1,
  m_t4_nee_2 = m_t4_nee_2,

  m_t4_reco_1 = m_t4_reco_1,
  m_t4_reco_2 = m_t4_reco_2,

  m_t4_gpp_1 = m_t4_gpp_1,
  m_t4_gpp_2 = m_t4_gpp_2
)

dt_sub_res <- data.frame()

for(i in 1:length(model_list)){
  
  m <- model_list[[i]]
  
  m_name <- names(model_list)[i]
  
  m_s <- summary(m)
  
  tmp_covs <- m_s$coefficients %>% 
    as.data.table()
  names(tmp_covs)[9] <- "sig_symbol"  
  
  tmp_r2 <- m_s$R2
  
  tmp_dt_sub <- tmp_covs %>% 
    left_join(tmp_r2) %>% 
    filter(!grepl("~~", Response)) %>% 
    mutate(model_name = m_name, 
           excluded_gradient = grad, 
           aic = as.numeric(AIC_psem(m)[1]),
           aicc = as.numeric(AIC_psem(m)[2]))
  
  dt_sub_res <- rbind(dt_sub_res, tmp_dt_sub)
  
}

dt_res <- rbind(dt_sub_res, dt_res)

}


dt_est <- dt_res %>% 
  rename(p_value = P.Value, 
         std_error = Std.Error, 
         rsq_c = Conditional, 
         rsq_m = Marginal, 
         predictor = Predictor, 
         response = Response, 
         estimate = Estimate) %>% 
  mutate(
    std_error = as.numeric(std_error),
    model_tier = case_when(
      grepl("t1", model_name) ~ "tier_1",
      grepl("t2", model_name) ~ "tier_2",
      grepl("t3", model_name) ~ "tier_3",
      grepl("t4", model_name) ~ "tier_4"),
    predictor_tier = case_when(
      grepl("_1", model_name) ~ "morph_traits_pca",
      grepl("_2", model_name) ~ "traits_separately",
      grepl("_3", model_name) ~ "all_traits_pca", 
      grepl("_4", model_name) ~ "leaf_traits_with_all_traits_data"),
    response_tier = case_when(
      .default = "veg",
      grepl("reco", response) & !grepl("temperature", response) ~ "flux",
      grepl("nee", response) & !grepl("temperature", response) ~ "flux",
      grepl("gpp", response) & !grepl("temperature", response) ~ "flux"),
    rsq_label = paste0("(R²m = ", round(rsq_m, 2), "; R²c = ", round(rsq_c,2), ")"), 
    # clean_response = case_when(
    #   response == "nee_anomaly_country" ~ paste0("NEE\n", rsq_label),
    #   response == "gpp_anomaly_country" ~ paste0("GPP\n", rsq_label),
    #   response == "reco_anomaly_country" ~ paste0("Reco\n", rsq_label),
    #   response == "height_x_cover_anomaly_country" ~ paste0("'Biomass'\n", rsq_label), 
    #   response == "morph_traits_pc1_anomaly_country" ~ paste0("Morph. Traits PCA 1\n", rsq_label),
    #   response == "morph_traits_pc2_anomaly_country" ~ paste0("Morph. Traits PCA 2\n", rsq_label),
    #   response == "all_traits_pc1_anomaly_country" ~ paste0("All Traits PCA 1\n", rsq_label),
    #   response == "all_traits_pc2_anomaly_country" ~ paste0("All Traits PCA 2\n", rsq_label),
    #   response == "chem_traits_pc1_anomaly_country" ~ paste0("Chem. Traits PCA 1\n", rsq_label),
    #   response == "chem_traits_pc2_anomaly_country" ~ paste0("Chem. Traits PCA 2\n", rsq_label),
    #   grepl("temperature", response) ~ paste0("Inst. Temp.\n", rsq_label),
    #   response == "sla_anomaly_country" ~ paste0("SLA\n", rsq_label),
    #   response == "leaf_area_anomaly_country" ~ paste0("Leaf Area\n", rsq_label),
    # ), 
    clean_response = case_when(
      response == "nee_anomaly_country" ~ "NEE",
      response == "gpp_anomaly_country" ~ "GPP",
      response == "reco_anomaly_country" ~ "Reco",
      response == "height_x_cover_anomaly_country" ~ "'Biomass'", 
      response == "morph_traits_pc1_anomaly_country" ~ "Morph. Traits PCA 1",
      response == "morph_traits_pc2_anomaly_country" ~ "Morph. Traits PCA 2",
      response == "all_traits_pc1_anomaly_country" ~ "All Traits PCA 1",
      response == "all_traits_pc2_anomaly_country" ~ "All Traits PCA 2", 
      response == "chem_traits_pc1_anomaly_country" ~ "Chem. Traits PCA 1", 
      response == "chem_traits_pc2_anomaly_country" ~ "Chem. Traits PCA 2", 
      grepl("temperature", response) ~ "Inst. Temp.", 
      response == "sla_anomaly_country" ~ "SLA", 
      response == "leaf_area_anomaly_country" ~ "Leaf Area",
    ), 
    ci_lb = estimate - 1.96*std_error,
    ci_ub = estimate + 1.96*std_error, 
    clean_term =  case_when(
      .default = predictor,
      grepl("temperature_", predictor) ~ "Inst. Temp.",
      predictor == "downscaled_temp_anomaly_country" ~ "Growing Season Temp.",
      predictor == "height_x_cover_anomaly_country"      ~ "'Biomass'",
      predictor == "morph_traits_pc1_anomaly_country"          ~ "Morph. Traits PC1",
      predictor == "morph_traits_pc2_anomaly_country"          ~ "Morph. Traits PC2",
      predictor == "all_traits_pc1_anomaly_country"          ~ "All Traits PC1",
      predictor == "all_traits_pc2_anomaly_country"          ~ "All Traits PC2",
      predictor == "chem_traits_pc1_anomaly_country"          ~ "Chem. Traits PC1",
      predictor == "chem_traits_pc2_anomaly_country"          ~ "Chem. Traits PC2",
      predictor == "sla_anomaly_country"                 ~ "SLA",
      predictor == "par_nee_anomaly_country"                 ~ "PAR",
      predictor == "leaf_area_anomaly_country"           ~ "Leaf Area"), 
    significance = case_when(
      .default = "Non significant", 
      ci_lb > 0 ~ "Significantly positive", 
      ci_ub < 0 ~ "Significantly negative", 
    )) %>% 
  mutate(
    clean_term = factor(clean_term,
                        levels = c("'Biomass'",
                                   "Leaf Area",
                                   "SLA",
                                   "Chem. Traits PC2", "Chem. Traits PC1",
                                   "All Traits PC2", "All Traits PC1",
                                   "Morph. Traits PC2", "Morph. Traits PC1",
                                   "Inst. Temp.",
                                   "Growing Season Temp.", 
                                   "PAR")),
    response = factor(response,
                      levels = c("gpp_anomaly_country",
                                 "nee_anomaly_country",
                                 "reco_anomaly_country",
                                 "temperature_gpp_anomaly_country",
                                 "temperature_reco_anomaly_country",
                                 "temperature_nee_anomaly_country",
                                 "morph_traits_pc1_anomaly_country",
                                 "morph_traits_pc2_anomaly_country",
                                 "all_traits_pc1_anomaly_country",
                                 "all_traits_pc2_anomaly_country",
                                 "chem_traits_pc1_anomaly_country",
                                 "chem_traits_pc2_anomaly_country",
                                 "sla_anomaly_country",
                                 "leaf_area_anomaly_country",
                                 "height_x_cover_anomaly_country")
    )
  )

resp_clean_comb <- dt_est %>% 
  select(response, clean_response) %>%
  unique() %>% 
  arrange(response) %>% 
  select(clean_response) %>%
  unique() 

dt_est$clean_response <- factor(dt_est$clean_response, levels = c(
  resp_clean_comb$clean_response
))
# Plots -------------
scico(9, palette = 'bam')
#"#001959" "#818231" "#F9CCF9"
c(MetBrewer::met.brewer(name = "Archambault", n = 6))
#"#88a0dc" "#7c4b73" "#ed968c" "#ab3329" "#e78429" "#f9d14a"
palette = c("Drakensberg" = "#88a0dc",
            "Central Andes" = "#7c4b73",
            "Eastern Himalayas" = "#ed968c", 
            "Rocky Mountains" = "#ab3329", 
            "Southern Scandes" = "#e78429", 
            "Svalbard" = "#f9d14a")

theme_est <-   theme(
                     legend.box="vertical",
                     plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                     panel.grid = element_line(color = "snow"), 
                     #axis.title.x = element_blank(), 
                     axis.text = element_text(size = 12), 
                     axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
                     panel.border = element_rect(color = NA), 
                     panel.background = element_rect(fill = "snow"), 
                     strip.text.x = element_text(size = 12), 
                     strip.text.y = element_text(size = 12, face = "bold"), 
                     strip.background = element_rect(fill = "linen", color = "linen") )



# Tier 4--------------------------

# traits separately -------
p_t4_flux_ts <- dt_est %>% 
  filter(!grepl("ntercept", predictor)) %>% 
  filter(model_tier == "tier_4" & predictor_tier == "traits_separately" & response_tier == "flux") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = excluded_gradient),
                  linewidth = 1.1, size = 0.9, alpha = 0.9, 
                  position = position_dodge(width = .5)) +
  scale_color_manual(values = palette) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Coefficient Estimate (±95 % CI)", y = "", color = "Excluded\nGradient") +
  theme_est +
  theme(legend.position = "bottom")
p_t4_flux_ts

p_t4_pred_ts <- dt_est %>% 
  filter(!grepl("ntercept", predictor) & grepl("gpp", model_name)) %>% 
  filter(model_tier == "tier_4" & predictor_tier == "traits_separately" & response_tier == "veg") %>% 
  ggplot() +
  geom_vline(xintercept = 0,linetype = "dashed", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = excluded_gradient),
                  linewidth = 1.1, size = 0.9, alpha = 0.9, 
                  position = position_dodge(width = .6)) +
  scale_color_manual(values = palette) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Coefficient Estimate (±95 % CI)", y = "",) +
  theme_est +
  theme(legend.position = "none")
p_t4_pred_ts

p_t4_r2_ts <- dt_est %>% 
  filter(!grepl("ntercept", predictor)) %>% 
  filter(!response %in% c("temperature_nee_anomaly_country", "temperature_reco_anomaly_country")) %>% 
  filter(model_tier == "tier_4" & predictor_tier == "traits_separately") %>% 
  ggplot() +
  geom_col(aes(y = excluded_gradient, x = rsq_m,
                  color = excluded_gradient, fill = excluded_gradient), size = 0.25) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 7)  +
  theme_bw() +
  labs(x = "Marginal R²", y = "Excluded\nGradient",) +
  theme_est +
  theme(legend.position = "none")
p_t4_r2_ts

p_t4_ts <- grid.arrange(p_t4_pred_ts, p_t4_flux_ts, p_t4_r2_ts, heights = c(0.9, 2, 1))

ggsave(plot = p_t4_ts,
       "builds/plots/supplement/loo_t4_estimates_traits_separately.png",
       dpi = 600, 
       height = 12, width = 12)

#morph. traits PCA -------
p_t4_flux_mpca <- dt_est %>% 
  filter(!grepl("ntercept", predictor)) %>% 
  filter(model_tier == "tier_4" & predictor_tier == "morph_traits_pca" & response_tier == "flux") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = excluded_gradient),
                  linewidth = 1.1, size = 0.9, alpha = 0.9, 
                  position = position_dodge(width = .6)) +
  scale_color_manual(values = palette) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Coefficient Estimate (±95 % CI)", y = "", color = "Excluded\nGradient") +
  theme_est +
  theme(legend.position = "bottom")
p_t4_flux_mpca

p_t4_pred_mpca <- dt_est %>% 
  filter(!grepl("ntercept", predictor) & grepl("gpp", model_name)) %>% 
  filter(model_tier == "tier_4" & predictor_tier == "morph_traits_pca" & response_tier == "veg") %>% 
  ggplot() +
  geom_vline(xintercept = 0,linetype = "dashed", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = excluded_gradient),
                  linewidth = 1.1, size = 0.9, alpha = 0.9, 
                  position = position_dodge(width = .5)) +
  scale_color_manual(values = palette) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Coefficient Estimate (±95 % CI)", y = "",) +
  theme_est +
  theme(legend.position = "none")
p_t4_pred_mpca

p_t4_r2_mpca <- dt_est %>% 
  filter(!grepl("ntercept", predictor)) %>% 
  filter(!response %in% c("temperature_nee_anomaly_country", "temperature_reco_anomaly_country")) %>% 
  filter(model_tier == "tier_4" & predictor_tier == "morph_traits_pca") %>% 
  ggplot() +
  geom_col(aes(y = excluded_gradient, x = rsq_m,
               color = excluded_gradient, fill = excluded_gradient), size = 0.25) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 7)  +
  theme_bw() +
  labs(x = "Marginal R²", y = "Excluded\nGradient",) +
  theme_est +
  theme(legend.position = "none")
p_t4_r2_mpca

p_t4_mpca <- grid.arrange(p_t4_pred_mpca, p_t4_flux_mpca, p_t4_r2_mpca, heights = c(0.9, 2, 1))

ggsave(plot = p_t4_mpca,
       "builds/plots/supplement/loo_t4_estimates_morph_traits_pca.png",
       dpi = 600, 
       height = 12, width = 12)

